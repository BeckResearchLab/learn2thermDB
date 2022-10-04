#!/bin/bash
#get credentials
HOST=ftp.ncbi.nlm.nih.gov
USER=anonymous
PASSWORD=$ENV_EMAIL

# create subdirectory if necessary
mkdir -p ./data/refseq
mkdir -p ./data/refseq/archaea
mkdir -p ./data/refseq/bacteria

# download archaea data
cd ./data/refseq/archaea
lftp -c "open -u $USER,$PASSWORD $HOST; cd /genomes/refseq/archaea; mget -P 15 -c -O './' /genomes/refseq/archaea/*/representative/*/*genomic.gbff.gz"

# download bacteria data
cd ../bacteria
lftp -c "open -u $USER,$PASSWORD $HOST; cd /genomes/refseq/bacteria; mget -P 15 -c -O './' /genomes/refseq/bacteria/*/representative/*/*genomic.gbff.gz"

# Notes
# The file format is .../refseq/organism_group(eg bacteria)/2 word genus group/<latest or all or reference>_assembly_version/<assembly accession>.<assembly version>_<assembbly name>
# Inside that folder is a series of files. I believe *_genomic.gbff.gz is the most complete info, we can extract proteins from it.
#    assembly_status.txt
#        A text file reporting the current status of the version of the assembly
#        for which data is provided. Any assembly anomalies are also reported.
#    *_assembly_report.txt file
#        Tab-delimited text file reporting the name, role and sequence 
#        accession.version for objects in the assembly. The file header contains 
#        meta-data for the assembly including: assembly name, assembly 
#        accession.version, scientific name of the organism and its taxonomy ID, 
#        assembly submitter, and sequence release date.
#    *_assembly_stats.txt file
#        Tab-delimited text file reporting statistics for the assembly including: 
#        total length, ungapped length, contig & scaffold counts, contig-N50, 
#        scaffold-L50, scaffold-N50, scaffold-N75, and scaffold-N90
#    *_assembly_regions.txt
#        Provided for assemblies that include alternate or patch assembly units. 
#        Tab-delimited text file reporting the location of genomic regions and 
#        listing the alt/patch scaffolds placed within those regions.
#    *_assembly_structure directory
#        This directory will only be present if the assembly has internal 
#        structure. When present, it will contain AGP files that define how 
#        component sequences are organized into scaffolds and/or chromosomes. 
#        Other files define how scaffolds and chromosomes are organized into 
#        non-nuclear and other assembly-units, and how any alternate or patch 
#        scaffolds are placed relative to the chromosomes. Refer to the README.txt
#        file in the assembly_structure directory for additional information.
#    *_cds_from_genomic.fna.gz
#        FASTA format of the nucleotide sequences corresponding to all CDS 
#        features annotated on the assembly, based on the genome sequence. See 
#        the "Description of files" section below for details of the file format.
#    *_feature_count.txt.gz
#        Tab-delimited text file reporting counts of gene, RNA, CDS, and similar
#        features, based on data reported in the *_feature_table.txt.gz file.
#        See the "Description of files" section below for details of the file 
#        format.
#    *_feature_table.txt.gz
#        Tab-delimited text file reporting locations and attributes for a subset 
#        of annotated features. Included feature types are: gene, CDS, RNA (all 
#        types), operon, C/V/N/S_region, and V/D/J_segment. Replaces the .ptt & 
#        .rnt format files that were provided in the old genomes FTP directories.
#        See the "Description of files" section below for details of the file 
#        format.
#    *_genomic.fna.gz file
#        FASTA format of the genomic sequence(s) in the assembly. Repetitive 
#        sequences in eukaryotes are masked to lower-case (see below).
#        The FASTA title is formatted as sequence accession.version plus 
#        description. The genomic.fna.gz file includes all top-level sequences in
#        the assembly (chromosomes, plasmids, organelles, unlocalized scaffolds,
#        unplaced scaffolds, and any alternate loci or patch scaffolds). Scaffolds
#        that are part of the chromosomes are not included because they are
#        redundant with the chromosome sequences; sequences for these placed 
#        scaffolds are provided under the assembly_structure directory.
#    *_genomic.gbff.gz file
#        GenBank flat file format of the genomic sequence(s) in the assembly. This
#        file includes both the genomic sequence and the CONTIG description (for 
#        CON records), hence, it replaces both the .gbk & .gbs format files that 
#        were provided in the old genomes FTP directories.
#    *_genomic.gff.gz file
#        Annotation of the genomic sequence(s) in Generic Feature Format Version 3
#        (GFF3). Sequence identifiers are provided as accession.version.
#        Additional information about NCBI's GFF files is available at 
#        ftp://ftp.ncbi.nlm.nih.gov/genomes/README_GFF3.txt.
#    *_genomic.gtf.gz file
#        Annotation of the genomic sequence(s) in Gene Transfer Format Version 2.2
#        (GTF2.2). Sequence identifiers are provided as accession.version.
#    *_genomic_gaps.txt.gz
#        Tab-delimited text file reporting the coordinates of all gaps in the 
#        top-level genomic sequences. The gaps reported include gaps specified in
#        the AGP files, gaps annotated on the component sequences, and any other 
#        run of 10 or more Ns in the sequences. See the "Description of files" 
#        section below for details of the file format.
#    *_protein.faa.gz file
#        FASTA format sequences of the accessioned protein products annotated on
#        the genome assembly. The FASTA title is formatted as sequence 
#        accession.version plus description.
#    *_protein.gpff.gz file
#        GenPept format of the accessioned protein products annotated on the 
#        genome assembly
#    *_rm.out.gz file
#        RepeatMasker output; 
#        Provided for Eukaryotes 
#    *_rm.run file
#        Documentation of the RepeatMasker version, parameters, and library; 
#        Provided for Eukaryotes 
#    *_rna.fna.gz file
#        FASTA format of accessioned RNA products annotated on the genome 
#        assembly; Provided for RefSeq assemblies as relevant (Note, RNA and mRNA 
#        products are not instantiated as a separate accessioned record in GenBank
#        but are provided for some RefSeq genomes, most notably the eukaryotes.)
#        The FASTA title is provided as sequence accession.version plus 
#        description.
#    *_rna.gbff.gz file
#        GenBank flat file format of RNA products annotated on the genome 
#        assembly; Provided for RefSeq assemblies as relevant
#    *_rna_from_genomic.fna.gz
#        FASTA format of the nucleotide sequences corresponding to all RNA 
#        features annotated on the assembly, based on the genome sequence. See 
#        the "Description of files" section below for details of the file format.
#    *_translated_cds.faa.gz
#        FASTA sequences of individual CDS features annotated on the genomic 
#        records, conceptually translated into protein sequence. The sequence 
#        corresponds to the translation of the nucleotide sequence provided in the
#        *_cds_from_genomic.fna.gz file. 
#    *_wgsmaster.gbff.gz
#        GenBank flat file format of the WGS master for the assembly (present only
#        if a WGS master record exists for the sequences in the assembly).
#    annotation_hashes.txt
#        Tab-delimited text file reporting hash values for different aspects
#        of the annotation data. See the "Description of files" section below 
#        for details of the file format.
#    md5checksums.txt file
#        file checksums are provided for all data files in the directory


# there are ~43k bacteria strains in here, and for many the "latest assembly versions" has multiple assemblies. 
# it would be nice to have only one assembly per species. Congruently "representative" contains an assembly labeled as representative
# of the strain, but not all species have it. There are only