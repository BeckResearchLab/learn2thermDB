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

import numpy as np
import pandas as pd

from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastpCommandline
from Bio.Blast import NCBIXML

import dask_jobqueue
from distributed import Client
import distributed
import dask

import learn2therm.io

from typing import List, Tuple

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

    Notes
    -----
    metric names:
        local_X - metric uses only information from a single HSP
        scaled_local_X - metric uses single HSP information normalized by global sequence length. 
            This is an attempt to filter out small but accurate chunks without needing the global alignment/qcoverage.
            Eg if match is |, mismatch is X and gap is -
                |||||||||||||||X|||||||||-------------------------------------------------------------------------------------------------
            Could score very well using local alignment score even though there is a huge gap tail.
        global_X - metric uses global alignment, currently cannot do without HSP tiling
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
                ((hsp.query_end - hsp.query_start)/self.record.query_length + (hsp.sbjct_end - hsp.sbjct_start)/alignment.length)/2)
        return np.argmax(scores), max(scores)

    def compute_metric(self, metric_name: str):
        """Compute the metric with specified name for each alignment"""
        if not hasattr(self, metric_name):
            raise ValueError(f"No metric found with name : {metric_name}")
        else:
            metric = getattr(self, metric_name)
        
        logger.debug(f"Computing metric `{metric_name}` for all alignments in query {self.qid}")

        outputs = []
        for alignment in self.record.alignments:
            outputs.append((self.qid, alignment.hit_id.split('|')[-1], metric(alignment)))
        return pd.DataFrame(data=outputs, columns=['query_id', 'subject_id', metric_name])

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

    def local_gap_compressed_percent_id(self, alignment):
        """Percent matches in match sequence, including but compressing gaps.
        
        The largest local HSP score is used
        """
        best_hsp_idx, _ = self.id_hsp_best_cov(alignment)
        hsp = alignment.hsps[best_hsp_idx]

        n_matches = hsp.identities
        n_gaps = hsp.gaps
        n_columns = len(hsp.query)
        n_compressed_gaps = len(re.findall('-+', hsp.query))+len(re.findall('-+', hsp.sbjct))
        return self.raw_gap_compressed_percent_id(n_matches, n_gaps, n_columns, n_compressed_gaps)

    def scaled_local_query_percent_id(self, alignment):
        """Percent matches in query sequence based on best HSP."""
        best_hsp_idx, _ = self.id_hsp_best_cov(alignment)
        hsp = alignment.hsps[best_hsp_idx]
        return hsp.identities/self.record.query_length

    def scaled_local_symmetric_percent_id(self, alignment):
        """Percent matches compared to average seq length of query and subject based on best HSP"""
        best_hsp_idx, _ = self.id_hsp_best_cov(alignment)
        hsp = alignment.hsps[best_hsp_idx]
        return 2*hsp.identities/(self.record.query_length + alignment.length)

    def local_E_value(self, alignment):
        """E value of HSP with most identities."""
        best_hsp_idx, _ = self.id_hsp_best_cov(alignment)
        hsp = alignment.hsps[best_hsp_idx]
        return hsp.expect

    def local_average_coverage(self, alignment):
        """The coverage of the HSP averaged for query and subject"""
        return self.id_hsp_best_cov(alignment)[1]


class AlignmentHandler:
    """Handles preparing, parallel processing, and post evaluation of alignments.
    
    High level steps:
    - find the proteins from meso and thermo taxa
    - prepare temporary input and output files of the desired proteins, uses BlastFiles
    - run alignment
    - apply metrics, uses BlastMetrics
    - deposit results to file

    This class expects files containing proteins of a very specific format, any deviation
    will cause issues:
    - files containing proteins for each taxa must all be in the same directory
    - file name of the form `taxa_index_X.csv`
    - file contains columns, in order: seq_id, protein_seq, protein_desc, protein_len
    - column seperator is semicolon ;
    
    Notes
    -----
    Class will create an empty file for the output. Alternatively, if the file already exists,
    execution is skipped. This is to ensure that if a pair takes too long for a job, the pair
    is not repeatedly ran and failed over and over. This leaves the vulnerability of the cluster
    ending mid computation and future jobs skipping the pair. AlignmentClusterFutures handles this
    issue.
    
    Parameters
    ----------
    meso_index : str
        index of mesophile
    thermo_index : str
        index of thermophile
    max_seq_len : int
        maximum length og protein sequence to consider
    protein_deposit : str
        path to directory containing protein files. Each file is of the form `taxa_index_X.csv`
        X is taxa index
    alignment_score_deposit : str
        path to directory to deposit completed metrics, files will be created of the form `taxa_pair_X-Y.csv`
        X is thermo taxa index, Y is meso
    alignment_params : dict
        parameters to pass to alignment method
    metrics: dict
        metric names to compute
    """
    def __init__(
        self,
        meso_index: str,
        thermo_index: str,
        max_seq_len: int, 
        protein_deposit: str,
        alignment_score_deposit: str,
        alignment_params: dict = {},
        metrics: list = ['scaled_local_symmetric_percent_id'],
    ):
        self.meso = meso_index
        self.thermo = thermo_index
        self.max_seq_len = max_seq_len
        if type(protein_deposit) == str:
            self.protein_deposit = protein_deposit
        else:
            raise ValueError(f"`protein_deposit` must be a string of a filepath")
        
        if type(alignment_score_deposit) == str:
            self.alignment_score_deposit = alignment_score_deposit
        else:
            raise ValueError(f"`alignment_score_deposit` must be a string of a filepath")
        
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

        # prepare instance variables
        self.meso_input_path = None
        self.thermo_input_path = None
        self.output_path = None
        self.output_exists = None
        return
    
    @property
    def pair_indexes(self):
        return f"{self.thermo}-{self.meso}"

    def _prepare_input_output_state(self):
        """Check and setup the input and output file paths.
        
        We are inside the worker if we have gotten to this stage.
        """
        # inputs
        if not os.path.exists(self.protein_deposit):
            raise ValueError(f"{self.protein_deposit} does not exist, cannot look for proteins")
        else:
            if not self.protein_deposit.endswith('/'):
                self.protein_deposit = self.protein_deposit+'/'
            self.meso_input_path = self.protein_deposit+f"taxa_index_{self.meso}.csv"
            self.thermo_input_path = self.protein_deposit+f"taxa_index_{self.thermo}.csv"

            if not os.path.exists(self.meso_input_path):
                raise ValueError(f"Could not find protein file for taxa {self.meso}")
            if not os.path.exists(self.thermo_input_path):
                raise ValueError(f"Could not find protein file for taxa {self.thermo}")
        
        # output
        if not os.path.exists(self.alignment_score_deposit):
            raise ValueError(f"Could not find deposit for data {self.alignment_score_deposit}")
        else:
            if not self.alignment_score_deposit.endswith('/'):
                self.alignment_score_deposit = self.alignment_score_deposit+'/'
            self.output_path = self.alignment_score_deposit+f"taxa_pair_{self.thermo}-{self.meso}.csv"
            # check if this output has been created
            if os.path.exists(self.output_path):
                self.output_exists = True
            else:
                self.output_exists = False
                # create an empty file
                file = open(self.output_path, 'w')
                file.close()
        return

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
        raise NotImplemented()
    
    def _compute_metrics(self, blast_record):
        """computes each of the metrics asked for for a single record.
        
        Returns
        -------
        DataFrame of protein pair indexes and list of metrics computed for a single query sequence
        """
        metric_handler = BlastMetrics(blast_record)
        df_outputs = [metric_handler.compute_metric(m) for m in self.metrics]
        # each df has query and subject id and one metric
        # join them all
        joined_df = df_outputs[0].set_index(['query_id', 'subject_id'])
        for df in df_outputs[1:]:
            joined_df = joined_df.join(df.set_index(['query_id', 'subject_id']))
        out = joined_df.reset_index()
        # ensure ordering
        out = out[['query_id', 'subject_id']+self.metrics]
        return out

    def run(self):
        """Run the whole alignment workflow for this taxa pair."""
        time0 = time.time()
        # check and prepare the state
        logger.debug(f"Checking and preparing file state for {self.pair_indexes}")
        self._prepare_input_output_state()

        # check the total number of pairwise seq we are doing
        time1 = time.time()
        meso_lengths = pd.read_csv(self.meso_input_path, sep=';', usecols=[3])['protein_len']
        thermo_lengths = pd.read_csv(self.thermo_input_path, sep=';', usecols=[3])['protein_len']
        logger.info(
            f"Running alignment for {self.pair_indexes} with total proteins counts ({len(thermo_lengths)},{len(meso_lengths)})")
        meso_count = (meso_lengths<=self.max_seq_len).sum()
        thermo_count = (thermo_lengths<=self.max_seq_len).sum()
        pairwise_space = meso_count * thermo_count
        logger.info(
            f"Alignment for {self.pair_indexes}, considering only seqs <= {self.max_seq_len} long, totalling ({thermo_count},{meso_count}). Total pairwise space {pairwise_space}")
        time2 = time.time()
        logger.debug(f"Parsed alignment space for {self.pair_indexes}, took {(time2-time1)/60}m")

        # if we have already done this part, get the results so that we can compute global
        # metrics and escape
        if self.output_exists:
            try:
                hits = len(pd.read_csv(self.output_path, usecols=[0]))
            except pd.errors.EmptyDataError:
                hits = 0
            return {'pair': self.pair_indexes, 'pw_space': pairwise_space, 'hits':hits}

        # create the iterators over the actual protein sequences
        # these only exist so that we can load a small chunk of sequences into memory
        # at a time
        meso_iter = learn2therm.io.csv_id_seq_iterator(
            self.meso_input_path,
            max_seq_length=self.max_seq_len,
            seq_col="protein_seq",
            sep=';',
            index_col='seq_id',
            chunksize=100
        )
        thermo_iter = learn2therm.io.csv_id_seq_iterator(
            self.thermo_input_path,
            max_seq_length=self.max_seq_len,
            seq_col="protein_seq",
            sep=';',
            index_col='seq_id',
            chunksize=100
        )

        # create the temporary files that a blastlike algorithm expects
        logger.info(f"Pair {self.pair_indexes}, using alignment parameters {self.alignment_params}")
        time1 = time.time()
        with BlastFiles(thermo_iter, meso_iter, dbtype='prot') as (query_tmp, subject_tmp, out_tmp):
            time2 = time.time()
            logger.debug(f"Pair {self.pair_indexes} temporary file preparation complete at {(query_tmp, subject_tmp, out_tmp)}, took {(time2-time1)/60}m")

            # run the alignment itself
            time1 = time.time()
            self._call_alignment(query_tmp, subject_tmp, out_tmp)
            time2 = time.time()
            logger.debug(f"Pair {self.pair_indexes} alignment complete, took {(time2-time1)/60}m")

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
            to_deposit = pd.concat(dataframes, ignore_index=True)
            logger.debug(f"Head of results DF: \n{to_deposit.head()}")
            logger.debug(f"Saving results to file {self.output_path}")
            to_deposit.rename(columns={'query_id':'thermo_protein_id', "subject_id": "meso_protein_id"}, inplace=True)
            to_deposit.to_csv(self.output_path)
            time2 = time.time()
            logger.debug(f"Computing metrics for {self.pair_indexes} took {(time2-time1)/60}m")  
      
        return {'pair': self.pair_indexes, 'pw_space': pairwise_space, 'hits': hits, 'execution_time': (time2-time0)/60}

class BlastAlignmentHandler(AlignmentHandler):

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
            evalue=10000000, # very large so we do not lose any hits
            # the rest are tunable params
            qcov_hsp_perc=50,
            max_hsps=100,
            **self.alignment_params
        )()

class DiamondAlignmentHandler(AlignmentHandler):

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
        logger.debug(f"Updating DB for diamond for pair {self.pair_indexes}")
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
        logger.debug(f"Updated DB for diamond for pair {self.pair_indexes}, took {(time1-time0)/60}m")
        
        
        command = f"diamond blastp -d {subject_db}.diamond -q {query_file} -o {output_file} --outfmt 5 --max-target-seqs 1000000 --max-hsps 100 --evalue 10000000 --query-cover 50 --subject-cover 0 --id 0 --masking 0"
        command = command + f" --{self.alignment_params['sensitivity']}"
        if self.alignment_params['iterate']:
            command = command + " --iterate"
        command = command + f" --matrix {self.alignment_params['matrix']}"
        command = command + f" --gapopen {self.alignment_params['gapopen']}"
        command = command + f" --gapextend {self.alignment_params['gapextend']}"
        command = command + f" --threads {self.alignment_params['num_threads']}"
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

class AlignmentClusterFutures:
    """Prepares and tracks the state of execution for a set of taxa pair alignments.
    
    This class exists to ensure that the correct work is skipped in the case that 
    we want to allow task skipping, and that no excess compute is wasted to job times
    
    The aligner handlers will skip work if the output files already exists. Unfortunately,
    they would skip taxa pairs that started but did not finish due to the cluster going down.
    We want to ensure that they only skip pairs that have either completed or could not be
    completed in worker time.
    
    Additionally, compute is wasted if a taxa pair starts on a worker that does not have enough
    time left to complete the job. This class ensures each task gets its own worker.
    
    Parameters
    ----------
    pairs : list of tuple of str
        eg [(thermo_index, meso_index), ...]
    client : dask client
    aligner_class : AlignmentHandler
    max_seq_len : int
        maximum length of proteins
    protein_deposit : str
        directory containing protein sequence files with appropriate format
    alignment_score_deposit : str
        path to directory to deposit completed metrics, files will be created of the form `taxa_pair_X-Y.csv`
        X is thermo taxa index, Y is meso
    worker_function : callable
        Must take a single argument, an AlignmentHandler instance, and return a dict containing the key "pair"
    alignment_params : dict
        parameters to pass to alignment method
    metrics: dict
        metric names to compute
    restart: bool
        whether to completely restart all work regardless of execution state
    """
    def __init__(
        self,
        pairs: List[Tuple[str]],
        client,
        protein_deposit: str,
        alignment_score_deposit: str,
        max_seq_len: int = 100,
        aligner_class: AlignmentHandler = BlastAlignmentHandler,
        worker_function: callable = lambda handler: handler.run(),
        alignment_params: dict = {},
        metrics: list = ['scaled_local_symmetric_percent_id'],
        restart: bool = True
    ):
        if not alignment_score_deposit.endswith('/'):
            alignment_score_deposit = alignment_score_deposit + '/'
        self.alignment_score_deposit = alignment_score_deposit
        
        if restart or not os.path.exists(alignment_score_deposit+'completion_state.metadat'):
            logger.info(f"Starting execution of {len(pairs)} from scratch")
            shutil.rmtree(alignment_score_deposit, ignore_errors=True)
            os.makedirs(alignment_score_deposit)
            self.metadata=None

        else:
            self.metadata=pd.read_csv(alignment_score_deposit+'completion_state.metadat', index_col=0)
            completed = self.metadata['pair'].values
            logger.info(f"Completed pairs: {completed}")

            # check existing files
            existing_files = os.listdir(alignment_score_deposit)
            existing_files = [f for f in existing_files if f.startswith('taxa')]
            cleanup_counter = 0
            for filename in existing_files:
                pair = filename.split('_')[-1].split('.')[0]
                if pair in completed:
                    pass
                    logger.info(f"Pair {pair} already completed")
                else:
                    logger.info(f"Pair {pair} erroneously ended, cleaning up file")
                    os.remove(alignment_score_deposit+filename)
                    cleanup_counter += 1

            logger.info(f"Found {len(completed)} pairs already complete. Cleaned up {cleanup_counter} erroneous files.")
            
        # create aligners and send out the job 
        aligners = [aligner_class(
            meso_index=mi,
            thermo_index=ti,
            max_seq_len=max_seq_len,
            protein_deposit=protein_deposit,
            alignment_score_deposit=alignment_score_deposit,
            metrics=metrics,
            alignment_params=alignment_params,
        ) for (ti, mi) in pairs]
        self.aligners = aligners
        self.client = client
        
        # these futures will record metadata as they complete, as well as terminated
        # workers that just completed it
        def futures_modified(futures):
            for future in futures:
                result = future.result()
                logger.info(f"Completed pair: {result}")
                who_has = client.who_has(future)
                closing = list(list(who_has.values())[0])
                client.retire_workers(closing)
                logger.debug(f"Retiring worker at {closing} that completed a task.")
                
                # record that we completed one
                if self.metadata is None:
                    self.metadata = pd.DataFrame(data=[list(result.values())], columns=list(result.keys()))
                else:
                    self.metadata = pd.concat([self.metadata, pd.Series(result).to_frame().T], axis=0, ignore_index=True)
                self.metadata.to_csv(self.alignment_score_deposit+'completion_state.metadat')
                yield result

        self.futures_modified = futures_modified(
            distributed.as_completed(client.map(worker_function, aligners)))

    def __enter__(self):
        return self.futures_modified

    def __exit__(self, type, value, traceback):
        pass
            