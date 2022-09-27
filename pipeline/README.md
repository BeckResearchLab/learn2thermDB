# Executables within data version controled pipeline

Script labeling is of the format `s<Sub Pipeline Number>.<Script Number>_<Script Name>`

## Directory

| Script | Function |
| ------ | -------- |
| ------ | -------- |
| s0.0_get_ncbi_refseq.sh | Pull most recent NCBI refseq sequences for archaea and bacteria |

## Notes
- Ensure `ENV_EMAIL` is set as an environmental variable. The email will be used to access NCBI FTP in `s0.0`