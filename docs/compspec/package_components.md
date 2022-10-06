# Component specification for importable classes and code

### `Learn2ThermDB`
Wrapper of sqlite database containing data. Allows for creation and splitting of datasets.

__Params__:
- `db_loc`: str, location of sqlite database on disk
- `protein_blast_threshold`: float, blast score required to count a protein pair
- `temperature_window`: float, number of celcius on either side of the threshold to discard organisms
- `max_sequence_length`: int, maximum number of amino acids in a sequence 

#### _Methods_:
`from_files(db_loc: str, taxa_table: str, taxa_pair_table: str, protein_table: str, protein_pair_table: str)`
Class method to create DB from files.
__Params__:
- `db_loc`: str, location of sqlite database on disk to save
- `taxa_table`: str, location on disk of table of taxa with labels
- `taxa_pair_table`: str, location on disk of table of taxa pairs
- `protein_table`: str, location on disk of table of proteins
- `protein_pair_table`: str, location on disk of table of protein pairs


`tt_split(how: str, test_frac: float)`  
Determine indexes of protein pairs for each split.  
Options: random splitting, split by genus, split by 16s rRNA feature vector (?)

__Params__:
- `how` - str, determines which split type
- `test_frac` - float, approximately what fraction of sequence pairs will be for test

`get_HF_dataset_unlabeled()`  
Returns Hugging Face dataset object of unlabeled sequences, split.

`get_HF_dataset_pairs()`  
Returns Hugging Face dataset object of meso-thermo pairs, split.
