# Machine translation of mesophic protein to same-function thermophilic protein

## Getting started
### Environment
Create and activate the environment specified in `environment.yml`

```
conda env create --file environment.yml
conda activate learn2therm
pip install .
```

Set environmental variables to ensure everything runs smoothly:
- `ENV_EMAIL` : email used durring HTTP and FTP execution

### Execution
Data Version Control (DVC) is used to track data, parameters, metrics, and execution pipelines.

To use a DVC remote, see the the [documentation](https://dvc.org/doc/command-reference/remote).

DVC tracked data, metrics, and models are found in `./data` while scripts and parameters can be found in `./pipeline`. To execute pipeline steps, run `dvc repro <stage-name>` where stages are listed below:

- TODO

### Python package
Installable, importable code is found in `src` and should be installed given the above steps in the __Environemnt__ section.

## Directory
```
-data/                                      # Contains DVC tracked data, models, and metrics
-src/                                       # Contains git tracked importable code
-pipeline/                                  # Contains DVC tracked executable pipeline steps and parameters
-notebooks/                                 # notebooks for testing and decisionmaking
-environment.yml                            # Conda dependancies