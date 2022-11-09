# Pipeline Documentation

This documentation will offer a comprehensive exposition on each script in the pipeline. For a high-level overview instead, see the script [compspecs]('learn2therm/docs/compspec/pipeline_components.md'). The purpose, inputs, steps, outputs, metrics, are all described for each script in this document.

1. `s0.0_get_ncbi_refseq.sh`
   _Objective_: extract refseq genome data for all bacteria and archaea.

   _Inputs_: there are no parameter inputs in this script as it is used purely for data extraction purposes.

   _Steps in the script_:
   a. The user will input his enviromental account credentials for NCBI
   b. The refseq genome data for both bateria and archaea will be downloaded from NCBI creating subdirectories for each species if necessary

   _Metrics_: there are no metrics being stored in this script.

   _Outputs_: the script will generate

2. `s0.1_get_bacdive_ogt.py`
