# check quantity of Tm data in database

Comparing to fireprotdb (s2.X_compare_to_Tm.py) yields only 6 proteins that we have OGT for and are in fireprotDB

The original paper that suggested the relationship between OGT and Tm as 24 + 0.9OGT (https://www.sciencedirect.com/science/article/pii/S0301462299001039?via%3Dihub) took data from here: 
(https://reader.elsevier.com/reader/sd/pii/S0022283697910421?token=67DBC7B2EBEE02DA1219994121E7F56164EB0695AB319BC7DA7A1C0082AC89FA19A6779CD49752531B99C4E523F4B7D3&originRegion=us-east-1&originCreation=20230517165058) has about 30. Let's see how many are in here

voigt_protein_pdbs.txt: contains pdb ids for proteins of prokaryotes from that old paper

look_for_voigt.py: check databse for those proteins and print out their labeled OGT
 - result: we only get 5 proteins back, why not more?
 - answer: we do, just some did not map to PDB. FOr example, if the same protein from 4 strains was found, and only one has PDB, but we took one of the other strains, we don't have the PDB id. 

## Conclusion:
- FireProtDB uses Uniprot ID. Not much we can do if they only have 6 proteins overlapping
