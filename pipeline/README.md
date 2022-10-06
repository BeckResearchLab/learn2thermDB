# Executables within data version controled pipeline

Script labeling is of the format `s<Sub Pipeline Number>.<Script Number>_<Script Name>`

## Pipeline Directory

| Script | Function |
| ------ | -------- |
| s0.0_get_ncbi_refseq.sh | Pull most recent NCBI refseq sequences for archaea and bacteria |
| s0.0_get_bacdive_ogt.py | Cross reference TAXID from NCBI to bacdive records, and attempt to extract information on OGT |

Parameters are grouped in terms of sub-pipelines:
| Sub-pipeline | Params |
| ------ | -------- |
| data extraction | s0_data_extraction_params.yaml |
| data processing | s1_data_processing_params.yaml |

See the script [compspec]('./docs/compspec/pipeline_components.md') for more details.

## ENV variables
Ensure that the following are set in order access full functionality
- `ENV_EMAIL` - The email will be used to access NCBI FTP in `s0.0` and Bacdive in `s0.1`
- `BACDIVE_PASSWORD` - Password associated with Bacdive account fo use in `s0.1`
- `LOGLEVEL` (optional) - Specified logging level to run the package. eg 'INFO' or 'DEBUG'