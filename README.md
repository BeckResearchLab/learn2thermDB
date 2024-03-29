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
Dask configuration is required for `s1.3`. This script stands up a Dask cluster of workers to conduct protein alignment. Each worker runs alignment between a one taxa pair at a time. In order to conduct this step, `.config/dask/jobqueue.yaml` must be updated. The pipeline was initially ran using slurm, and the config file is provided, however names and accounts will necessarily need to be changed for your cluster. If using a distributed scheduler other than SLURM, it must be supported by dask and the appropriate configurations made. See [here](https://jobqueue.dask.org/en/latest/api.html).

If the pipeline is erroring out at this stage, it is likely an issue with the cluster configuration. Common issues experienced are workers not being able to find executables.
Ensure that the workers have appropriate environment setups. If using SLURM, the existing config file is a working example: note the sourceing of `bashrc`, environment activation, and environment variable exports in the job preludes/directives.


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

Dependancies between stages is shown below:
```mermaid
flowchart TD
        node1["../../analysis_pipeline/dvc.yaml:chosen_protein_search_space"]
        node2["../../analysis_pipeline/dvc.yaml:filtered_protein_search_space"]
        node3["../../analysis_pipeline/dvc.yaml:full_protein_search_space"]
        node4["../../dvc.yaml:compare_hait_alignment"]
        node5["../../dvc.yaml:compare_hait_hmmer"]
        node6["../../dvc.yaml:compare_to_Tm"]
        node7["../../dvc.yaml:get_16s_blast_scores"]
        node8["../../dvc.yaml:get_HMM_profiles"]
        node9["../../dvc.yaml:get_hait_pairs"]
        node10["../../dvc.yaml:get_proteome_mdata"]
        node11["../../dvc.yaml:get_raw_data_proteins"]
        node12["../../dvc.yaml:get_raw_data_taxa"]
        node13["../../dvc.yaml:hmmer_hait"]
        node14["../../dvc.yaml:label_all_pairs"]
        node15["../../dvc.yaml:label_taxa"]
        node16["../../dvc.yaml:make_database"]
        node17["../../dvc.yaml:parse_hmmer_result"]
        node18["../../dvc.yaml:parse_proteins"]
        node19["../../dvc.yaml:protein_alignment"]
        node20["../../dvc.yaml:run_hait_alignment"]
        node21["../../dvc.yaml:run_hmmer"]
        node22["../../dvc.yaml:sample_data_for_structure"]
        node23["../../dvc.yaml:structure_hait"]
        node24["../../dvc.yaml:structure_l2t"]
        node7-->node2
        node7-->node14
        node7-->node16
        node7-->node19
        node8-->node13
        node8-->node21
        node9-->node13
        node9-->node20
        node9-->node23
        node10-->node18
        node11-->node18
        node12-->node3
        node12-->node7
        node12-->node15
        node12-->node16
        node12-->node18
        node14-->node2
        node14-->node16
        node14-->node19
        node15-->node1
        node15-->node7
        node15-->node16
        node16-->node4
        node16-->node5
        node16-->node6
        node16-->node21
        node16-->node22
        node17-->node5
        node18-->node1
        node18-->node2
        node18-->node3
        node18-->node7
        node18-->node16
        node18-->node19
        node19-->node16
        node19-->node17
        node20-->node4
        node21-->node17
        node22-->node24
```

Note that script execution is expected to occur with the top level as the current working directory, and paths are specified with respect to the repo top level.

The entire pipeline can be run in a single command using `dvc exp run`. Note however that many steps are resource intensive. At least 8 cores and 60GB of RAM is recommended. For `s1.3`, Dask must be configured. This stage requires a distributed computing cluster.

## Python package
Installable, importable code used in the pipeline is found in `learn2therm` and should be installed given the above steps in the __Environment__ section.


