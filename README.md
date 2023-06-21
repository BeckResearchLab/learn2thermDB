# Pipeline to produce dataset of pairs low and high temperature proteins

## Getting started
### Environment
Create and activate the environment specified in `environment.yml`

```
conda env create --file environment.yml
conda activate learn2thermDB
pip install .
```

Ensure that the following are set in order access full functionality
- `ENV_EMAIL` - The email will be used to access NCBI FTP in `s0.0` and Bacdive in `s0.1`
- `NCBI_API_KEY` - API key from NCBI. Needed to get 16s sequences in `s0.0`. See [here](https://support.nlm.nih.gov/knowledgebase/article/KA-05317/en-us)
- `LOGLEVEL` (optional) - Specified logging level to run the package. eg 'INFO' or 'DEBUG'
- `FATCAT_EXEC` - Necessary for struvtura alignment steps of validation, `s2.11` and `s2.12`

### Config

#### DVC
Configure DVC as desired. This is not required, however here are a few recommendations:
- Add a remote data repository
- Use hardlinking instead of copying to the DVC cache to increase speed.
See [here](https://dvc.org/doc/command-reference/config)

### Parameters

`params.yaml` contains all of the tunable parameters used to run the pipeline. Modifying these parameters will product a different result than in the presented paper. Note that DVC tracks input and output states, thus by changing a parameter, DVC whill know which stages need to be run. The parameters are discussed in the context of their pipeline stages in `./docs/compspec/pipeline_components.md`. The params file itself has comments indicating what the parameter does.

#### Dask
Dask configuration is required for `s1.3`. This script stands up a Dask cluster of workers to conduct protein alignment. Each worker runs alignment between a one taxa pair at a time. In order to conduct this step, `.config/dask/jobqueue.yaml` must be updated. The pipeline was initially ran using slurm, and the config file is provided, however names and accounts will necessarly need to be changed for your cluster. If using a distributed scheduler other than SLURM, it must be supported by dask and the appropriate configurations made. See [here](https://jobqueue.dask.org/en/latest/api.html).


### Execution
Data Version Control (DVC) is used to track data, parameters, metrics, and execution pipelines.

To use a DVC remote, see the the [documentation](https://dvc.org/doc/command-reference/remote).

DVC tracked data, metrics, and models are found in `./data` while scripts and parameters can be found in `./pipeline`. To execute pipeline steps, run `dvc repro <stage-name>` where stages are listed below, and details on stagess can be found in `./docs/compspec/pipeline_components.md`:

- get_raw_data_taxa
- get_raw_data_proteins
- get_proteome_mdata
- parse_proteins
- label_taxa
- get_16s_blast_scores
- label_all_pairs
- protein_alignment
- make_database
- get_hait_pairs
- compare_to_Tm
- run_hait_alignment
- compare_hait_alignment
- get_HMM_profiles
- hmmer_hait
- run_hmmer
- parse_hmmer_result
- compare_hait_hmmer
- sample_data_for_structure
- structure_hait
- structure_l2t

Note that script execution is expected to occur with the top level as the current working directory, and paths are specified with respect to the repo top level.

The entire pipeline can be run in a single command using `dvc exp run`. Note however that many steps are resource intensive. At least 8 cores and 60GB of RAM is recommended. For `s1.3`, Dask must be configured. This stage requires a distributed computing cluster.

## Python package
Installable, importable code used in the pipeline is found in `learn2therm` and should be installed given the above steps in the __Environment__ section.


