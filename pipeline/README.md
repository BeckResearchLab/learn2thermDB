# Executables within data version controled pipeline

Script labeling is of the format `s<Sub Pipeline Number>.<Script Number>_<Script Name>`

## Pipeline Directory

| Script | Function |
| ------ | -------- |
| s0.0_get_raw_data_taxa.py | Pull most recent NCBI 16s r RNA sequences, and OGT records from Enqvist |
| s0.1_get_raw_data_proteins.py | Download bacterial and archaeal UniProtKB proteins. |
| s0.2_get_proteome_mdata.py | Download Proteome metadata from Uniprot. A "best" proteome is selected for each species |
| s0.3_parse_proteins.py | Loop raw UniProt data and record protein seq and select DB identifiers if the protein belongs to a taxa with OGT. |
| s1.0_label_taxa.py | Assign thermophile labels to organisms |
| s1.1_get_16s_blast_scores.py | Run pairwise blastn on 16s sequences |
| s1.2_label_all_pairs.py | Use thresholds on 16s metrics to label organism pairs |

See the script [compspec]('./docs/compspec/pipeline_components.md') for more details.

## ENV variables
Ensure that the following are set in order access full functionality
- `ENV_EMAIL` - The email will be used to access NCBI FTP in `s0.0` and Bacdive in `s0.1`
- `NCBI_API_KEY` - API key from NCBI. Needed to get 16s sequences in `s0.0`. See [here](https://support.nlm.nih.gov/knowledgebase/article/KA-05317/en-us)
- `LOGLEVEL` (optional) - Specified logging level to run the package. eg 'INFO' or 'DEBUG'
- `FATCAT_EXEC` - Necessary for struvtura alignment steps of validation, `s2.11` and `s2.12`

## Notes
### s0.1

Extremely long runtime due to throttling by uniprot. Script does not actually DVC track the downloaded
files such that they can be skipped of the script is restarted. Thus, DVC cannot automatically determine
if this stage needs rerun - the files must be removed and the stage force ran to download newer version
of UniProt

## s0.3

Primary function of this script is twofold:
1. Vastly reduce the data size, keeping only primary sequence, database IDs eg alphafold and PDB, taxa ID, and proteome ID
2. Remove proteins that are unlabeled or redundant. If a protein belongs to an organism without OGT, it is skipped. If a protein is associated with a proteome but it is not the representative proteome selected for the organism, is it skipped.