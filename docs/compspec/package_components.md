# Component specification for importable classes and code

## `database.Learn2ThermDB`
Wrapper of sqlite database containing data. Allows for creation and splitting of datasets.

__Params__:
- `db_path`: str, location of sqlite database on disk
- `read_only`: bool, whether to open the database read only

### _Methods_:
#### `from_files`
Class method to create DB from files.
__Params__:
- `files_path`: str, location of data files produced by pipeling
- `db_path`: str, location on disk where db file will be stored
- `read_only`: bool, whether to open the db read only

***

## `io`

### `seq_io_gnuzipped`
Opens biopython files using biopython that are gnuzipped.

__Params__:
- `filepath`: str, location of zipped file
- `filetype`: str, type of data file, passed directly to biopython

__Returns__:
- `records`: list of SeqRecord

### `csv_id_seq_iterator`
Returns a one by one iterator of seq ids and sequences from a large csv to avoid OOM.

__Params__:
- `csv_filepath` : str
  - path to file containing data
- `seq_col` : str
  - name of column containing sequences
- `index_col`: str, default None
  - which column name is associated with the index, otherwise the 0th column will be used
- `id_filter` : Collection, Optional
  - If given, only return sequences with the provided indexes
- `chunksize` : int, default 512
  - Number of sequences that will be stored in memory at once.
- `max_seq_length` : int, default None
  - Maximum length of sequence to return
**kwargs passed to pandas read csv 

__Returns__:
- iterator of id, sequence (eg amino acids or nucleotides)

***

## `utils`

### `WorkerLogger`
Custom logger for logging one worker, adds dask worker information to log message

### `start_logger_if_necessary`
Start a logger if necessary on a dask worker

__Params__:
- `logger_name` : str
  - name of logger to start or retrieve
- `log_file` : str
  - path to file to log to
- `log_level`
  - log level to respect
- `worker`: str
  - name of worker using this logger
- `filemode` : str
  - mode to apply to log file eg "a" for append

***

## `bacdive.BacdiveClient`
Subclass of bacdive client to get bacdive recorc by ncbi id.

### _Methods_:
#### `getIDByNCBITaxID`
Get the bacdive id by NCBI id

__Params__:
- `tid`: int, ncbi id

***

## `blast`

### `BlastFiles`
Context manager to handle creating and cleaning up of temporary files for local alignment

__Params__:
- `query_iterator` : iterator of (seq id, sequence)
  - sequences to be used as query
- `subject_iterator` : iterator of (seq id, sequence)
  - sequences to be used as the "database"

__Returns__
- `query_filename` : str, name of fasta file with query sequences
- `subject_filename` : str, name of fasta file with subject sequences
- `output_filename` : str, name of output file for blast to save reuslts, will be deleted out of context

### `BlastMetrics`
Class that handles computing alignmnment metrics for the output of blast

__Params__:
- `blast_record`: str, record file containing all hits for a particular query

#### _Methods_:

#### `id_hsp_best_cov`
Finds the hsp for the alignment with the largest symmetric coverage
__Params__:
- `alignment`: biopython alignment object associated with a hit

__Returns__:
- `int`, positions of best HSP
- `float`, value of symmetric coverage for alignment

#### `compute_metric`
Compute a metric for all hits in the record considering the HSP for each hit with the highest average coverage.

__Params__:
- `metric_name`: str, must match a class method associated with a metric, options below

__Returns__:
- `dataframe`, columns are metric names and rows are hits for the query

#### Metric options
All can be computed for any partular HSP
- `local_gap_compressed_percent_id`: percent id relative to the total alignment length, considering gaps as -1 regardless of size
- `scaled_local_query_percent_id`: percent id of the query
- `scaled_local_symmetric_percent_id`: average percent id of query and subject
- `bit_score`: bit score of HSP according to the scoring matrix
- `local_E_value`: evalue of the HSP
- `query_align_start`: position of start of alignment in query
- `query_align_end`: position of end of alignment in query
- `query_align_len`: length of alignment over the query sequence
- `query_align_cov`: frac of alignment over the total query length
- `subject_align_start`: position of start of alignment in subject
- `subject_align_end`: position of end of alignment in subject
- `subject_align_len`: length of alignment over the subject sequence
- `subject_align_cov`: frac of alignment over the total subject length

### `AlignmentHandler`
Abstract parent. Handles preparing, alignment, and post evaluation for a taxa pair
    
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

Children must define `_call_alignment`

__Notes__:

Class will create an empty file for the output. Alternatively, if the file already exists,
execution is skipped. This is to ensure that if a pair takes too long for a job, the pair
is not repeatedly ran and failed over and over. This leaves the vulnerability of the cluster
ending mid computation and future jobs skipping the pair. AlignmentClusterState handles this
issue.

__Params__
- `meso_index` : str
  - index of mesophile
- `thermo_index` : str
  - index of thermophile
- `max_seq_len` : int
  - maximum length og protein sequence to consider
- `protein_deposit` : str
  - path to directory containing protein files. Each file is of the form `taxa_index_X.csv`
  - X is taxa index
- `alignment_score_deposit` : str
  - path to directory to deposit completed metrics, files will be created of the form `taxa_pair_X-Y.csv`
  - X is thermo taxa index, Y is meso
- `alignment_params` : dict
  - parameters to pass to alignment method
- `metrics`: dict
  - metric names to compute

#### _Methods_:
#### `run`
Executes the protocol. 

Once the alignment for a single taxa pair is complete, a single new csv file will be popualted with pairwise alignments found and associated metrics computed for them.

### `BlastAlignmentHandler(AlignmentHandler)`
Align proteins in a taxa pair with blastp

### `DiamondAlignmentHandler(AlignmentHandler)`
Align proteins in a taxa pair with Diamond

### `AlignmentClusterState`
Prepares and tracks the state of execution for a set of taxa pair alignments.

This class exists to ensure that the correct work is skipped in the case that 
we want to allow task skipping, and that no excess compute is wasted to job times.

The aligner handlers will skip work if the output files already exists. Unfortunately,
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
        handlers will skip them
killer_workers = True
    - always use when tasks are long
    - when tasks are short, this adds non negligable overhead, so instead should be used
        as final sweep

The state is tracked as jobs are run in a csv file located in the alignment deposit called `completion_state.metadat`

__Params__:
- `pairs` : list of tuple of str
  - eg [(thermo_index, meso_index), ...]
- `client` : dask client
- `aligner_class` : AlignmentHandler
- `max_seq_len` : int
  - maximum length of proteins
- `protein_deposit` : str
  - directory containing protein sequence files with appropriate format
- `alignment_score_deposit` : str
  - path to directory to deposit completed metrics, files will be created of the form `taxa_pair_X-Y.csv`
  - X is thermo taxa index, Y is meso
- `worker_function` : callable
  - Must take a single argument, an AlignmentHandler instance, and return a dict containing the key "pair"
- `alignment_params` : dict
  - parameters to pass to alignment method
- `metrics`: dict
  - metric names to compute
- `restart`: bool, default True
  - whether to completely restart all work regardless of execution state
- `killer_workers`: bool, default True
  - Whether or not to kill worker after each task


***

## `hmmer`
A set of importable functions for running hmmer via the package pyhmmer, parsing output, and saving results

### `hmmpress_hmms`
Presses the HMMs in the given HMM database and stores the resulting files in a specified directory.
Note to use this function, the user must run `s2.5_get_HMM_profiles.py` first to obtain the HMMs.

__Params__:
- `hmmdb_path` : str
  - Path to the HMM database.
- `pfam_data_folder` : str, optional
  - Path to the directory where the HMMs should be stored.

__Returns__:
- `None`

### `prefetch_targets`
Prefetch HMM profiles from a given HMM database.
This function is meant to be a mode of the `run_pyhmmer` function.

__Params__:
- pressed_path : str
  - Path to the pressed HMM database.

__Returns__:
- targets : pyhmmer.plan7.OptimizedProfileBlock
  - The HMM profiles loaded from the database.


### `save_to_digital_sequences`
Save protein sequences from a DataFrame to a digital sequence block.

__Params__:
- `dataframe` : pd.DataFrame
  - DataFrame containing PIDs (Protein IDs) and sequences.
    columns: 'pid', 'protein_seq'

__Returns__:
- DigitalSequenceBlock
  - A digital sequence block containing the converted sequences.

### `run_pyhmmer`
Run HMMER's hmmscan or hmmsearch program on a set of input sequences using with HMMs from a database.

__Params__:
- `seqs` : pyhmmer.easel.DigitalSequenceBlock
  - Path to the input sequence file.
- `hmms_path` : str
  - Path to the HMM database.
- `prefetch` : bool, optional
  - Specifies how the HMM are stored in meomry.
    Also, can be a pyhmmer.plan7.OptimizedProfileBlock object.
- `output_file` : str, optional
  - Path to the output file if the users wants to write the file.
- `cpu` : int, optional
  - The number of CPUs to use. Default is 4.
- `scan`: bool, optional
  - Whether to run hmmscan or hmmsearch. Default is True (hmmscan).
- `eval_con` : float, optional
  - E-value threshold for domain reporting. Default is 1e-10.

__Returns__:
- `all_hits` : pyhmmer.plan7.TopHits or domtblout file
  - If the output_file has a name, it will be written to a domtblout file.
    Otherwise, the user will get a list of pyhmmeer TopHits objects.

__Notes__:

This function runs HMMER's hmmscan/hmmsearch program on a set of input sequences
using HMMs from a given database.
The function supports two modes: normal mode and prefetching mode.
In normal mode, the HMMs are pressed and stored in a directory before execution.
In prefetching mode, the HMMs are kept in memory for faster search.

### `parse_pyhmmer`
Parses the TopHit pyhmmer object getting the query and accession IDs and saves to a DataFrame

__Params__:
- `all_hits` : list
  - A list of TopHit objects from pyhmmer.
- `chunk_query_ids` : list
  - A list of query IDs from the chunk.

__Returns__:
-  `pandas.DataFrame`
  - A dataframe containing the query and accession IDs.

### `TODO` parsing functions documentation