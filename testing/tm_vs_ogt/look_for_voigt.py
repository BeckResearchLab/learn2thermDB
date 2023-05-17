import duckdb as ddb

con = ddb.connect('../../tmp/database', read_only=True)

# load the PDB ids
with open('./voigt_protein_pdbs.txt', 'r') as f:
    pdbs = f.readlines()
pdbs = [s.upper().strip() for s in pdbs]
print(pdbs)

# check db
count = con.execute(f"SELECT COUNT(*) FROM proteins WHERE pdb_id IN {tuple(pdbs)}").df()
print("Count: ", count)
pdbs_found = con.execute(f"SELECT pdb_id FROM proteins WHERE pdb_id IN {tuple(pdbs)}").df()['pdb_id']
print("Found: ", pdbs_found)
pdbs_not_found = [pdb for pdb in pdbs if pdb not in pdbs_found]
print("Not found: ", pdbs_not_found)