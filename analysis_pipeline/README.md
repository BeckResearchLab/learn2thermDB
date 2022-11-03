# Analysis pipeline

This folder contains a dvc tracked pipeline. This pipeline is not part of the necessary steps to produce the final dataset, for those steps see the `pipeline` folder. Steps in this pipeline are tangential to the primary pipeline and are present to extract information or help make decisions.

Script labeling is of the format `t<Primary Script Number Necessary>_<Script Name>`

## Pipeline Directory

| Script | Function |
| ------ | -------- |
| t1.4_protein_alignment_resource_test.py | Run a small subset of taxa pairs through protein blasting in order to estimate cost |
