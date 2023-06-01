# Pipeline to produce dataset of pairs low and high temperature proteins

## Getting started
### Environment
Create and activate the environment specified in `environment.yml`

```
conda env create --file environment.yml
conda activate learn2thermDB
pip install .
```

Ensure that the following environmental variables are set for pipeline exacution:  
- `ENV_EMAIL` - The email will be used to access NCBI FTP in `s0.0`
- `NCBI_API_KEY` - NCBI api key associeted with the Entrex tool. Needed for `s0.0`
- `LOGLEVEL` (optional) - Specified logging level to run the package. eg 'INFO' or 'DEBUG'

### Config

#### DVC
Configure DVC as desired. This is not required, however here are a few recommendations:
- Add a remote data repository
- Use hardlinking instead of copying to the DVC cache to increase speed.
See [here](https://dvc.org/doc/command-reference/config)

#### Dask
Dask configuration is required for `s1.3`. This script stands up a Dask cluster of workers to conduct protein alignment. Each worker runs alignment between a one taxa pair at a time. In order to conduct this step, `.config/dask/jobqueue.yaml` must be updated. The pipeline was initially ran using slurm, and the config file is provided, however names and accounts will necessarly need to be changed for your cluster. If using a distributed scheduler other than SLURM, it must be supported by dask and the appropriate configurations made. See [here](https://jobqueue.dask.org/en/latest/api.html).


### Execution
Data Version Control (DVC) is used to track data, parameters, metrics, and execution pipelines.

To use a DVC remote, see the the [documentation](https://dvc.org/doc/command-reference/remote).

DVC tracked data, metrics, and models are found in `./data` while scripts and parameters can be found in `./pipeline`. To execute pipeline steps, run `dvc repro <stage-name>` where stages are listed below, and a detailed list can be found in `./docs/compspec/pipeline_components.md`:

- TODO

Note that script execution is expected to occur with the top level as the current working directory, and paths are specified with respect to the repo top level.

The entire pipeline can be run in a single command using `dvc exp run`. Note however that many steps are resource intensive. At least 8 cores and 60GB of RAM is recommended. For `s1.3`, Dask must be configured. This stage requires a distributed computing cluster.

### Python package
Installable, importable code used in the pipeline is found in `learn2therm` and should be installed given the above steps in the __Environment__ section.

## Directory
```
-data/                                      # Contains DVC tracked data, metrics
-learn2therm/                               # Contains git tracked importable code
-pipeline/                                  # DVC tracked executable pipeline steps. Steps to reproduce final result.
-analysis_pipeline/                         # DVC tracked executable pipeline steps. Optional, branches main pipeline
-notebooks/                                 # notebooks for testing and decision making
-environment.yml                            # Conda dependancies
-docs/                                      # repository documentation
```
