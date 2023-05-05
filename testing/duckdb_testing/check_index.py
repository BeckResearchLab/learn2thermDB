import duckdb as ddb
import pandas as pd

conn = ddb.connect('./test.db')

conn.execute("CREATE OR REPLACE TABLE proteins AS SELECT * FROM read_parquet('../../data/proteins/uniprot_chunk_0.parquet')")
conn.commit()
conn.execute("CREATE INDEX taxid ON proteins(taxid)")
conn.commit()

assert len(conn.execute("SELECT * FROM duckdb_indexes()").df()) > 0

conn.close()
conn = ddb.connect('./test.db')

if len(conn.execute("SELECT * FROM duckdb_indexes()").df()) > 0:
    print('success')
else:
    print('fail')
