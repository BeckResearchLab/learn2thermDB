"""Wraps up use of the NCBI blast software for use with this package.

We are not using a known database, we are comparing lists of sequences, so we need tools
to support this.

Note BioPython provides a python API to the command line
"""
import tempfile
import os
import shutil
import re
import time
import subprocess
from collections import OrderedDict

import numpy as np
import pandas as pd
import duckdb as ddb

from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastpCommandline
from Bio.Blast import NCBIXML

from codecarbon import OfflineEmissionsTracker
import dask_jobqueue
from distributed import Client
import distributed
import dask

import learn2therm.io

from typing import List, Tuple, Dict, Any, Callable

import logging
logger = logging.getLogger(__name__)

class BlastFiles:
    """Temporary files for use with BLAST CLI.
    
    Blast expects two input FASTA and produces an XML. The FASTA are redundant to CSV
    we already have. These files are created for the context and removed after completion.

    Parameters
    ----------
    query_iterator : iterator of (seq id, sequence)
        sequences to be used as query
    subject_iterator : iterator of (seq id, sequence)
        sequences to be used as the "database"

    Returns
    -------
    query_filename : str, name of fasta file with query sequences
    subject_filename : str, name of fasta file with subject sequences
    output_filename : str, name of output file for blast to save reuslts, will be deleted out of context
    """
    def __init__(self, query_iterator, subject_iterator, dbtype: str = 'nucl', dev_sample_num: int = None):
        # we have to create the temporary fasta files
        logger.debug("Creating temporary files to deposit blast inputs and outputs.")
        os.makedirs('./tmp/', exist_ok=True)
        query_temp = tempfile.NamedTemporaryFile('w', delete=False, dir='./tmp/')
        logger.debug(f"query file: {query_temp.name}")
        if dev_sample_num is not None:
            logger.debug(f"Using only max {dev_sample_num} sequences from query and subject")
        self.qt = query_temp.name
        n = 0
        for id_, seq in query_iterator:
            if seq == 'None' or seq is None:
                continue
            query_temp.write(f">{id_}\n{seq}\n")
            n +=1
            if dev_sample_num is not None:
                if n >= dev_sample_num:
                    break
        query_temp.close()
        logger.debug(f"added {n} sequences to query file")

        # folder for subject DB after we make a fasta
        subject_folder = tempfile.mkdtemp(dir='./tmp/')
        self.st = subject_folder
        subject_fasta_file = subject_folder+'/subs.fasta'
        self.subject_fasta_file = subject_fasta_file
        n = 0
        file = open(subject_fasta_file, 'w')
        for id_, seq in subject_iterator:
            if seq == 'None' or seq is None:
                continue
            file.write(f">{id_}\n{seq}\n")
            n +=1
            if dev_sample_num is not None:
                if n >= dev_sample_num:
                    break
        file.close()
        logger.debug(f"added {n} sequences to db file")

        # create db for it
        NcbimakeblastdbCommandline(dbtype=dbtype, input_file=subject_fasta_file, parse_seqids=True)()
        logger.debug(f"created database")
        # create the output xml file
        out_temp = tempfile.NamedTemporaryFile('w', delete=False, dir='./tmp')
        self.ot = out_temp.name

    def __enter__(self):
        return self.qt, self.subject_fasta_file, self.ot

    def __exit__(self, type, value, traceback):
        logger.debug("Removing temporary files used by blast")
        os.remove(self.qt)
        shutil.rmtree(self.st)
        os.remove(self.ot)
        

class BlastMetrics:
    """Handles computation of metrics for each alignment in a blast record.

    The HSP with the largest average sequence coverage is used for local metrics.
    
    Parameters
    ----------
    blast_record : the record containing all hits for a query
    """
    def __init__(self, blast_record):
        self.record = blast_record
        self.qid = self.record.query.split(' ')[0]
        logger.debug(f"Query {self.qid} with {len(self.record.alignments)} alignments.")

    def id_hsp_best_cov(self, alignment):
        """Determine HSP with the most average coverage of both sequences.
        
        Returns
        -------
        Index of HSP with max average seq coverage
        Max average coverage
        """
        scores = []
        for hsp in alignment.hsps:
            scores.append(
                ((hsp.query_end +1 - hsp.query_start)/self.record.query_length + (hsp.sbjct_end +1 - hsp.sbjct_start)/alignment.length)/2)
        return np.argmax(scores), max(scores)

    def compute_metrics(self, metric_names: List[str]):
        """Compute the specified metrics for each alignment"""
        # Check if all metric names are valid
        for metric_name in metric_names:
            if not hasattr(self, metric_name):
                raise ValueError(f"No metric found with name : {metric_name}")

        logger.debug(f"Computing metrics {metric_names} for all alignments in query {self.qid}")

        outputs = []
        for alignment in self.record.alignments:
            hsp_id, _ = self.id_hsp_best_cov(alignment)
            hsp = alignment.hsps[hsp_id]

            # Compute all the metrics for the current alignment
            metric_values = [getattr(self, metric_name)(alignment, hsp) for metric_name in metric_names]

            # Append the query ID, subject ID, and metric values to the outputs list
            if  '|' in alignment.hit_id:
                hit_id = alignment.hit_id.split('|')[-2]
            else:
                hit_id = alignment.hit_id
            outputs.append((self.qid, hit_id, *metric_values))

        # Create the DataFrame with columns for query ID, subject ID, and each metric
        columns = ['query_id', 'subject_id'] + metric_names
        return pd.DataFrame(data=outputs, columns=columns)

    @staticmethod
    def raw_gap_excluding_percent_id(n_matches, n_gaps, n_columns):
        """Percent matches in sequence, excluding gaps.
        
        Parameters
        ----------
        n_matches : int, number of matches in match columns
        n_gaps : number of gaps in match columns
        n_columns : total number of alignment match columns
        """
        return n_matches / (n_columns - n_gaps)

    @staticmethod
    def raw_gap_including_percent_id(n_matches, n_columns):
        """Percent matches in sequence, including gaps gaps.
        
        Parameters
        ----------
        n_matches : int, number of matches in match columns
        n_columns : total number of alignment match columns
        """
        return n_matches / (n_columns)

    @staticmethod
    def raw_gap_compressed_percent_id(n_matches, n_gaps, n_columns, n_compressed_gaps):
        """Percent matches in sequence, including but compressing gaps.
        
        Parameters
        ----------
        n_matches : int, number of matches in match columns
        n_gaps : number of gaps in match columns
        n_columns : total number of alignment match columns
        n_compressed_gaps : number of compressed gaps in match columns
        """
        return n_matches / (n_columns - n_gaps + n_compressed_gaps)

    def local_gap_compressed_percent_id(self, alignment, hsp):
        """Percent matches in match sequence, including but compressing gaps.
        
        The largest local HSP score is used
        """
        n_matches = hsp.identities
        n_gaps = hsp.gaps
        n_columns = len(hsp.query)
        n_compressed_gaps = len(re.findall('-+', hsp.query))+len(re.findall('-+', hsp.sbjct))
        return self.raw_gap_compressed_percent_id(n_matches, n_gaps, n_columns, n_compressed_gaps)

    def scaled_local_query_percent_id(self, alignment, hsp):
        """Percent matches in query sequence based on best HSP."""
        return hsp.identities/self.record.query_length

    def scaled_local_symmetric_percent_id(self, alignment, hsp):
        """Percent matches compared to average seq length of query and subject based on best HSP"""
        return 2*hsp.identities/(self.record.query_length + alignment.length)

    def local_E_value(self, alignment, hsp):
        """E value of HSP with most identities."""
        return hsp.expect

    def query_align_start(self, alignment, hsp):
        """Start index of alignment in query."""
        return hsp.query_start

    def query_align_end(self, alignment, hsp):
        """End index of alignment in query."""
        return hsp.query_end
    
    def subject_align_end(self, alignment, hsp):
        """End index of alignment in subject."""
        return hsp.sbjct_end

    def subject_align_start(self, alignment, hsp):
        """Start index of alignment in subject."""
        return hsp.sbjct_start

    def query_align_len(self, alignment, hsp):
        """Length of AA on query string taken up by alignment"""
        return int(hsp.query_end +1 - hsp.query_start)

    def query_align_cov(self, alignment, hsp):
        """Fraction of AA on query string taken up by alignment"""
        return (hsp.query_end +1 - hsp.query_start)/self.record.query_length
    
    def subject_align_len(self, alignment, hsp):
        """Length of AA on query string taken up by alignment"""
        return int(hsp.sbjct_end +1 - hsp.sbjct_start)

    def subject_align_cov(self, alignment, hsp):
        """Fraction of AA on query string taken up by alignment"""
        return (hsp.sbjct_end +1 - hsp.sbjct_start)/alignment.length
    
    def bit_score(self, alignment, hsp):
        return hsp.score

class AlignmentHandler:
    """Abstract parent. Handles preparing, alignment, and post evaluation for two sets of protein sequences.
    
    High level steps:
    - prepare input and output files for the desired protein sets, uses BlastFiles
    - run alignment
    - apply metrics, uses BlastMetrics
    - give back results as a DataFrame

    Children must define `_call_alignment`

    Parameters
    ----------
    seqs_A : DataFrame
        DataFrame containing protein sequences for the first set
        Must have columns 'pid' and 'sequence'
    seqs_B : DataFrame
        DataFrame containing protein sequences for the second set
        Must have columns 'pid' and 'sequence'
    alignment_params : dict
        parameters to pass to alignment method
    metrics: dict
        metric names to compute
    """
    def __init__(
        self,
        seqs_A: pd.DataFrame,
        seqs_B: pd.DataFrame,
        alignment_params: dict = {},
        metrics: list = ['scaled_local_symmetric_percent_id'],
    ):
        self.seqs_A = seqs_A
        self.seqs_B = seqs_B
        
        if type(alignment_params) == dict:
            self.alignment_params = alignment_params
        else:
            raise ValueError(f"`alignment_params` must be dict, kwargs for align method")
        
        if type(metrics) == list:
            for mname in metrics:
                if not hasattr(BlastMetrics, mname):
                    raise ValueError(f"Specified metric {mname} not available.")
            self.metrics = metrics
        else:
            raise ValueError(f"`metrics` should be list")
        return

    def _call_alignment(self, query_file: str, subject_db: str, output_file: str):
        """Parse the alignment parameters and call the correct CL command
        
        Parameters
        ----------
        query_file : str
            path to temporary fasta file containing queries
        subject_db : str
            path of temporary db containing subjects
        output_file: str
            path to temporary file to deposit alignment outputs, XML format
        """
        raise NotImplemented()

    def _compute_metrics(self, blast_record):
        """computes each of the metrics asked for for a single record.
        
        Returns
        -------
        DataFrame of protein pair indexes and list of metrics computed for a single query sequence
        """
        metric_handler = BlastMetrics(blast_record)
        df_outputs = metric_handler.compute_metrics(self.metrics)
        return df_outputs

    def run(self):
        """Run the whole alignment workflow for this pair of sequence sets."""

        time0 = time.time()

        # check the total number of pairwise seq we are doing
        time1 = time.time()
        count_A = len(self.seqs_A)
        count_B = len(self.seqs_B)

        pairwise_space = count_A * count_B
        logger.info(
            f"Running alignment considering seqs totalling ({count_A},{count_B}). Total pairwise space {pairwise_space}")
        time2 = time.time()
        logger.debug(f"Parsed alignment, took {(time2-time1)/60}m")

        # create the iterators over the actual protein sequences
        def seq_id_iter(df):
            for _, row in df.iterrows():
                yield row['pid'], row['sequence']
        iter_A = seq_id_iter(self.seqs_A)
        iter_B = seq_id_iter(self.seqs_B)

        # create the temporary files that a blastlike algorithm expects
        logger.info(f"Using alignment parameters {self.alignment_params}")
        time1 = time.time()
        with BlastFiles(iter_A, iter_B, dbtype='prot') as (query_tmp, subject_tmp, out_tmp):
            time2 = time.time()
            logger.debug(f"Temporary file preparation complete at {(query_tmp, subject_tmp, out_tmp)}, took {(time2-time1)/60}m")

            # run the alignment itself
            time1 = time.time()
            self._call_alignment(query_tmp, subject_tmp, out_tmp)
            time2 = time.time()
            logger.debug(f"Alignment complete, took {(time2-time1)/60}m")

            # apply the metrics to it
            # blast record io
            xml_f = open(out_tmp, 'r')
            blast_result_records = NCBIXML.parse(xml_f)

            time1 = time.time()
            hits = 0 # tracks the number of successful alignments
            dataframes = []
            for record in blast_result_records:
                df = self._compute_metrics(record)
                dataframes.append(df)
                hits += len(df)
            if hits == 0:
                df_results = pd.DataFrame()
            else:
                df_results = pd.concat(dataframes, ignore_index=True)
            time2 = time.time()
            logger.debug(f"Computed metrics, saved to file, head of results DF: \n{df_results.head()} took {(time2-time1)/60}m")  
            out_metadata = {'pw_space': pairwise_space, 'hits': hits, 'execution_time': (time2-time0)/60}
            logger.debug(f"Metadata for this run: {out_metadata}")
        return df_results, out_metadata

class BlastAlignmentHandler(AlignmentHandler):
    """Alignment using blastp
    
    Biopython's NcbiBlasp API is used, but at command line that is actually called
    is `blastp`
    """
    def _call_alignment(self, query_file: str, subject_db: str, output_file: str):
        """Parse the alignment parameters and call the correct CL command
        
        Parameters
        ----------
        query_file : str
            path to temporary fast file containing queries
        subject_db : str
            path of temporary db containing subjects
        output_file: str
            path to temporary file to deposit alignment outputs, XML format
        """
        # parse parameters 
        NcbiblastpCommandline(
            query=query_file,
            db=subject_db,
            outfmt=5,
            out=output_file,
            max_target_seqs=10000000, # very large so we do not throw out any pairs. will have to increase if there is more than this num of mesos
            max_hsps=100,
            **self.alignment_params
        )()

class DiamondAlignmentHandler(AlignmentHandler):
    """Alignment using Diamond
    
    https://github.com/bbuchfink/diamond
    
    Note the database file first has to be updated. Currently we create a brand
    new diamond database directly from the fasta because the diamond DB conversion
    from blast is not yet conda installable. This is minimal extra overhead.
    
    Calls `diamond`.
    """
    def _call_alignment(self, query_file: str, subject_db: str, output_file: str):
        """Parse the alignment parameters and call the correct CL command
        
        Parameters
        ----------
        query_file : str
            path to temporary fast file containing queries
        subject_db : str
            path of temporary db containing subjects
        output_file: str
            path to temporary file to deposit alignment outputs, XML format
        """
        time0 = time.time()
        # diamond requires a special DB type
        command = f"diamond makedb --in {subject_db} -d {subject_db}.diamond"
        logger.debug(f"Updating DB for diamond")
        process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # wait for it to be done
        while process.poll() is None:
            time.sleep(5)
        stdout, stderr = process.communicate()
        logger.debug(stdout)
        if process.poll() != 0:
            logger.error(stderr)
            raise ValueError(f"Conversion to diamond db failed, check logs.")
        time1 = time.time()
        logger.debug(f"Updated DB for diamond took {(time1-time0)/60}m")
        
        
        command = f"diamond blastp -d {subject_db}.diamond -q {query_file} -o {output_file} --outfmt 5 --max-target-seqs 1000000 --max-hsps 100 --id 0 --masking 0"
        command = command + f" --{self.alignment_params['sensitivity']}"
        if self.alignment_params['iterate']:
            command = command + " --iterate"
        command = command + f" --matrix {self.alignment_params['matrix']}"
        command = command + f" --gapopen {self.alignment_params['gapopen']}"
        command = command + f" --gapextend {self.alignment_params['gapextend']}"
        command = command + f" --threads {self.alignment_params['num_threads']}"
        command = command + f" --evalue {self.alignment_params['evalue']}"
        command = command + f" --query-cover {self.alignment_params['hsp_cov']}"
        command = command + f" --subject-cover {self.alignment_params['hsp_cov']}"
        if self.alignment_params['global_ranking']:
            command = command + f" --global_ranking {self.alignment_params['global_ranking']}"
        
        logger.debug(f"Diamond command: {command}")
        
        process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # wait for it to be done
        while process.poll() is None:
            time.sleep(5)
        stdout, stderr = process.communicate()
        logger.debug(stdout)
        if process.poll() != 0:
            logger.error(stderr)
            raise ValueError(f"Diamond alignment failed, check logs.")
        return
            
class TaxaAlignmentWorker:
    """Runs alignment for a taxa pair.

    Produces an output file in the output directory with the following format:
        align_taxa_{meso_index}_{thermo_index}.parquet

    Note that the file is always created, even if the alignment fails.
    Additionally, if the file already exists, it will not be overwritten.
    
    Parameters
    ----------
    meso_index : int
        index of meso taxa
    thermo_index : int  
        index of thermo taxa
    protein_database : str
        path to protein database
    output_dir : str
        path to output directory
    max_seq_len : int
        maximum sequence length
    alignment_params : Dict[str, Any]
        parameters for alignment
    alignment_handler : AlignmentHandler
        class to handle alignment
    metrics : list
        list of metrics to compute
    restart : bool
        whether to restart if output file already exists
    """
    def __init__(
        self,
        meso_index: int,
        thermo_index: int,
        protein_database: str,
        output_dir: str,
        max_seq_len: int,
        alignment_params: Dict[str, Any],
        aligner_class: AlignmentHandler = BlastAlignmentHandler,
        metrics: list = ['scaled_local_symmetric_percent_id'],
        retry_counter: int = 3,
    ):
        self.meso_index = meso_index
        self.thermo_index = thermo_index
        self.protein_database = protein_database
        if not output_dir.endswith("/"):
            output_dir = output_dir + "/"
        self.output_dir = output_dir
        self.max_seq_len = max_seq_len
        self.alignment_params = alignment_params
        self.aligner_class = aligner_class
        self.metrics = metrics
        self.retry_counter = retry_counter

    @property
    def pair_indexes(self):
        return (self.thermo_index, self.meso_index)

    def run(self):
        # define output file
        output_filename = f"align_taxa_{self.thermo_index}-{self.meso_index}.parquet"

        # start carbon tracking
        tracker = OfflineEmissionsTracker( 
            project_name=output_filename,
            output_dir=self.output_dir,
            country_iso_code='USA',
            region='Washington'
        )
        tracker.start()
        # if we have already done this part, get the results so that we can compute global
        # metrics and escape
        if os.path.exists(self.output_dir+output_filename):
            emissions = tracker.stop()
            run_metadata =  {'pw_space': None, 'hits': None, 'execution_time': None, 'emissions': emissions, 'target': output_filename}
            logger.info(f"Worker already started for pair {self.pair_indexes}, skipping.")
        else:
            # save empty file to indicate that we are working on this pair
            out_columns = ['thermo_pid', 'meso_pid'] + self.metrics + ['thermo_taxid', 'meso_taxid']
            pd.DataFrame(columns=out_columns).astype(pd.StringDtype()).to_parquet(self.output_dir+output_filename)

            # connect to DB
            conn = ddb.connect(self.protein_database, read_only=True)
            # get sequences from DB, lower length format
            meso_sequences = conn.execute(f"SELECT pid, sequence FROM proteins WHERE taxid={self.meso_index}").df()
            thermo_sequences = conn.execute(f"SELECT pid, sequence FROM proteins WHERE taxid={self.thermo_index}").df()
            logger.info(f"Found {len(meso_sequences)} meso sequences and {len(thermo_sequences)} thermo sequences for pair {self.pair_indexes}")
            # filter by length
            meso_sequences = meso_sequences[meso_sequences['sequence'].str.len() <= self.max_seq_len]
            thermo_sequences = thermo_sequences[thermo_sequences['sequence'].str.len() <= self.max_seq_len]
            logger.info(f"After filtering by length, found {len(meso_sequences)} meso sequences and {len(thermo_sequences)} thermo sequences for pair {self.pair_indexes}")

            conn.close()
            # give sequences to alignment handler
            handler = self.aligner_class(
                seqs_A=thermo_sequences,
                seqs_B=meso_sequences,
                alignment_params=self.alignment_params,
                metrics=self.metrics)
            # run handler, get output
            results, run_metadata = handler.run()

            # if we got none, skip this part
            if results.empty:
                logger.info(f"No results for pair {self.pair_indexes}")
            else:
                # add come data and do some renaming
                results['thermo_taxid'] = self.thermo_index
                results['meso_taxid'] = self.meso_index
                results = results.rename(columns={'query_id': 'thermo_pid', 'subject_id': 'meso_pid'})
                results.to_parquet(self.output_dir+output_filename)
            # return metadata
            emissions = tracker.stop()
            run_metadata['emissions'] = emissions
            run_metadata['target'] = output_filename
        output_metadata = OrderedDict()
        output_metadata['pw_space'] = run_metadata['pw_space']
        output_metadata['hits'] = run_metadata['hits']
        output_metadata['execution_time'] = run_metadata['execution_time']
        output_metadata['emissions'] = run_metadata['emissions']
        output_metadata['target'] = run_metadata['target']
        return output_metadata


class TaxaAlignmentClusterState:
    """Prepares and tracks the state of execution for a set of taxa pair alignments.
    
    This class exists to ensure that the correct work is skipped in the case that 
    we want to allow task skipping, and that no excess compute is wasted to job times.
    
    The TaxaAlignmentWorker will skip work if the output files already exists. Unfortunately,
    they would skip taxa pairs that started but did not finish due to the cluster going down.
    We want to ensure that they only skip pairs that have either completed or could not be
    completed in worker time.
    
    Additionally, compute is wasted if a taxa pair starts on a worker that does not have enough
    time left to complete the job. This class can ensure each task gets its own worker, which shuts
    down after completion.
    
    Two primary modes should be used:
    killer_workers = False
        - use when tasks are fast. In this case, workers are not closed after each task
        - some percentage of tasks will get cutoff even if they are quick, and subsequent
          workers will skip them
    killer_workers = True
        - always use when tasks are long
        - when tasks are short, this adds non negligable overhead, so instead should be used
          as final sweep
    
    The state is tracked as jobs are run in a csv file located in the alignment deposit called `completion_state.metadat`

    Parameters
    ----------
    pairs : list of tuple of str
        eg [(thermo_index, meso_index), ...]
    client : dask client
    aligner_class : AlignmentHandler
    max_seq_len : int
        maximum length of proteins
    protein_database : str
        duckdb database of proteins
    output_dir : str
        path to directory to deposit completed metrics, files will be created of the form `align_taxa_X-Y.parquet`
        X is thermo taxa index, Y is meso
    alignment_params : dict
        parameters to pass to alignment method
    metrics: dict
        metric names to compute
    restart: bool, default True
        whether to completely restart all work regardless of execution state
    killer_workers: bool, default True
        Whether or not to kill worker after each task
    worker_function: Callable, default lambda aligner: aligner.run()
        Function taking one argument, the aligner, that runs the job
    """
    def __init__(
        self,
        pairs: List[Tuple[str]],
        client,
        protein_database: str,
        output_dir: str,
        max_seq_len: int = 100,
        aligner_class: AlignmentHandler = BlastAlignmentHandler,
        alignment_params: dict = {},
        metrics: list = ['scaled_local_symmetric_percent_id'],
        restart: bool = True,
        killer_workers: bool = True,
        aggregate_outputs: bool = True,
        worker_function: Callable = lambda aligner: aligner.run()
    ):
        self.killer_workers = killer_workers
        self.worker_function = worker_function
        self.aggregate_outputs = aggregate_outputs

        if not output_dir.endswith('/'):
            output_dir = output_dir + '/'
        self.output_dir = output_dir
        
        # prepare and empty deposit if either we asked for restart
        # or cannot find existing metadata
        if restart or not os.path.exists(output_dir+'completion_state.metadat'):
            logger.info(f"Starting execution of {len(pairs)} from scratch")
            shutil.rmtree(output_dir, ignore_errors=True)
            os.makedirs(output_dir)
            self.metadata=None
        
        # otherwise load the metadata
        else:
            self.metadata=pd.read_csv(output_dir+'completion_state.metadat')
            # get all pairs what we successfully got a hit for
            completed = self.metadata[~self.metadata['hits'].isna()]['target'].values
            logger.info(f"Found {len(self.metadata)} worker completions, and {len(completed)} finished jobs")
            logger.info(f"Completed pairs: {completed}")

            # check existing files and clean the ones that were not completed
            # or timed out, which will be stored in metadata
            existing_files = os.listdir(output_dir)
            existing_files = [f for f in existing_files if f.startswith('align_taxa')]
            cleanup_counter = 0
            
            for target in existing_files:
                if target in completed:
                    logger.debug(f"{target} already completed")
                else:
                    logger.info(f"{target} erroneously ended, cleaning up file")
                    os.remove(output_dir+target)
                    cleanup_counter += 1

            logger.info(f"Found {len(completed)} pairs already complete. Cleaned up {cleanup_counter} erroneous files.")
            
            # now cleanup the incoming pairs to only the ones that are not done
            new_pairs = []
            for pair in pairs:
                if f"align_taxa_{pair[0]}-{pair[1]}.parquet" in completed:
                    pass
                else:
                    new_pairs.append(pair)

            logger.info(f"Trying again for {len(new_pairs)} pairs.")
            pairs = new_pairs
            
        # create aligners and send out the job 
        aligners = [TaxaAlignmentWorker(
            meso_index=mi,
            thermo_index=ti,
            max_seq_len=max_seq_len,
            protein_database=protein_database,
            output_dir=output_dir,
            aligner_class=aligner_class,
            alignment_params=alignment_params,
            metrics=metrics,
        ) for (ti, mi) in pairs]
        self.aligners = aligners
        self.client = client
        
        # these futures will record metadata as they complete, as well as terminated
        # workers that just completed it if that option is specified
        def futures_modified(futures):
            for future in futures:
                result = future.result()
                logger.info(f"Completed pair: {result}")
                if self.killer_workers:
                    who_has = client.who_has(future)
                    closing = list(list(who_has.values())[0])
                    client.retire_workers(closing)
                    logger.info(f"Retiring worker at {closing} that completed a task.")
                
                # record that we completed one
                if self.metadata is None:
                    self.metadata = pd.DataFrame(data=[list(result.values())], columns=list(result.keys()))
                    self.metadata.to_csv(self.output_dir+'completion_state.metadat', index=False)
                else:
                    file = open(self.output_dir+'completion_state.metadat', 'a')
                    output_values = [str(v) for v in result.values()]
                    file.write(','.join(output_values)+'\n')
                    file.close()
                yield result

        self._futures = distributed.as_completed(client.map(worker_function, aligners))

        if len(aligners) > 0:
            self.futures_modified = futures_modified(
                self._futures
            )
        else:
            self.futures_modified = None

    def _aggregate_output_files(self):
        """Aggregate files of distinct taxa pair files into larger chunks.
        
        We lose the file viewable information about which taxa re contained, but
        too many distinct files is cumbersome for DVC and duckdb, so we aggregate
        into files of size 500000 pairs.
        """
        unaggregated_files = os.listdir(self.output_dir)
        unaggregated_files = [f for f in unaggregated_files if f.startswith('align_taxa')]
        aggregate_files = os.listdir(self.output_dir)
        aggregate_files = [f for f in aggregate_files if f.startswith('agg_chunk')]
        dfs = []
        file_size = 0
        i = len(aggregate_files)
        for f in unaggregated_files:
            df = pd.read_parquet(self.output_dir+f)
            dfs.append(df)
            file_size += len(df)
            if file_size >= 500000:
                agg = pd.concat(dfs, axis=0, ignore_index=True)
                agg.to_parquet(self.output_dir+f'agg_chunk_{i}.parquet')
                logger.info(f"Saved chuck composed of {len(dfs)} taxa pairs.")
                i += 1
                file_size = 0
                dfs = []
            os.remove(self.output_dir+f)
        logger.info("Aggregated protein pair files into chunks.")

    def _close(self):
        self.client.cancel(self._futures, force=True)
        time.sleep(15)
        self.client.close(timeout=15)
        logger.info(f"Canceled futures.")
        if self.aggregate_outputs:
            self._aggregate_output_files()

    def __enter__(self):
        return self.futures_modified

    def __exit__(self, type, value, traceback):
        self._close()
            