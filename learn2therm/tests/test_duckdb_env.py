"""Duckdb has trouble with index persistence. This test is to make sure that
the index is actually persisted with the current environment.
"""
import duckdb as ddb
import pandas as pd
import pytest
import tempfile

def test_duckdb_index_persistence():
    """Test that the index is persisted when the database is closed and reopened.
    """
    toy_data = pd.DataFrame({'a': [1,2,3], 'b': [4,5,6]})
    with tempfile.TemporaryDirectory(dir='./tmp/') as tmpdir:
        conn = ddb.connect(tmpdir+'/test.db')
        conn.execute("CREATE TABLE toy_data AS SELECT * FROM toy_data")
        conn.commit()
        conn.execute("CREATE INDEX a ON toy_data(a)")
        conn.commit()
        conn.close()
        conn = ddb.connect(tmpdir+'/test.db')
        assert len(conn.execute("SELECT * FROM duckdb_indexes()").df()) > 0
        
