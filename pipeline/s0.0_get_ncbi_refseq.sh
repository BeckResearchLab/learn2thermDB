#!/bin/bash
#get credentials
HOST=ftp.ncbi.nlm.nih.gov
USER=anonymous
PASSWORD=$ENV_EMAIL

# create subdirectory if necessary
mkdir -p ../data/refseq
mkdir -p ../data/refseq/archaea
mkdir -p ../data/refseq/bacteria

# download archaea data
cd ../data/refseq/archaea
lftp -c "open -u $USER,$PASSWORD $HOST; cd /genomes/refseq/archaea; mget -c -O './' /genomes/refseq/archaea/*/latest_assembly_versions/*/*.gbff.gz"
# unzip all
for FILE in *; do gzip -d $FILE; done
