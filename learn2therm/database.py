"""Wraps creation of a duckdb database for all of the data and quick reboot"""
import os
import sys
import time

import numpy as np
import pandas as pd
import duckdb as ddb

import logging
logger = logging.getLogger(__name__)

class L2TDatabase:
    
    """Python wrapper for database of all taxa, proteins, taxa pairs, and protein pairs.
    
    Parameters
    ----------
    db_path : str
        Path to duck database folder
    read_only : bool, default True
        Whether or not to leave the database as read only
    
    Attributes
    ----------
    history : List of Tuple
        (command string, execution time) of calls to the current connection
    read_only : bool
        Whether the current connection is read only
    transaction_open : bool
        Whether or not there is currently a transaction open
    """
    def __init__(self, db_path: str, read_only: bool=True):
        conn = ddb.connect(db_path, read_only=read_only)
        
        self.path = db_path
        self.conn = conn
        self.read_only = read_only
        
    @classmethod
    def from_files(
        cls,
        files_path: str,
        db_path: str,
        read_only: bool=True
    ):
        """Create the database files.
        
        Parameters
        ----------
        files_path : str
            Path to root dir of csv data files
        db_path : str
            Path where databse files will be created
        """
        conn = ddb.connect(db_path, read_only=False)
        
        if not files_path.endswith('/'):
            files_path.append('/')
        else:
            pass

        t1 = time.time()
        cls._create_taxa_table(conn, files_path)
        t2 = time.time()
        logger.info(f"Took {(t2-t1)/60}m to create taxa table")
        
        t1 = time.time()
        cls._create_proteins_table(conn, files_path)
        t2 = time.time()
        logger.info(f"Took {(t2-t1)/60}m to create protein table")

        t1 = time.time()
        cls._create_taxa_pairs_table(conn, files_path)
        t2 = time.time()
        logger.info(f"Took {(t2-t1)/60}m to create taxa pair table")

        t1 = time.time()
        cls._create_protein_pairs_table(conn, files_path)
        t2 = time.time()
        logger.info(f"Took {(t2-t1)/60}m to create protein pair table")

        # create indexes
        t1 = time.time()
        cls._create_indexes(conn)
        t2 = time.time()
        logger.info(f"Took {(t2-t1)/60}m to create speed indexes")

        conn.close()        
        return cls(db_path, read_only=read_only)

    @staticmethod
    def _create_taxa_table(conn, files_path):
        # schema for taxa
        conn.execute(f"""
            CREATE OR REPLACE TEMP TABLE taxa_raw AS SELECT 
                "taxid"::BIGINT AS taxid,
                "16s_seq"::VARCHAR AS "16s_seq",
                "16s_len"::BIGINT AS "16s_len",
                "temperature"::FLOAT AS temperature,
                "superkingdom"::DOUBLE AS superkingdom,
                "phylum"::DOUBLE AS phylum,
                "class"::DOUBLE AS class,
                "order"::DOUBLE AS order,
                "family"::DOUBLE AS family,
                "genus"::DOUBLE AS genus,
            FROM read_parquet('{files_path}taxa.parquet')""")
        conn.execute(f"""CREATE OR REPLACE TEMP TABLE taxa_labels AS SELECT 
            "taxid"::BIGINT AS taxid,
            "thermophile_label"::BOOLEAN AS thermophile_label,
        FROM read_parquet('{files_path}taxa_thermophile_labels.parquet')""")
        conn.execute("""CREATE OR REPLACE TABLE taxa AS 
            SELECT *
            FROM taxa_raw
            INNER JOIN taxa_labels ON (taxa_raw.taxid=taxa_labels.taxid)""")
        conn.execute("""
            ALTER TABLE taxa DROP COLUMN "taxid:1"
        """)
        conn.commit()
    
    @staticmethod
    def _create_proteins_table(conn, files_path):
        # create table
        conn.execute(f"""
            CREATE OR REPLACE TABLE proteins AS SELECT
                "pid"::VARCHAR AS pid,
                "taxid"::BIGINT AS taxid,
                "pdb_id"::VARCHAR AS pdb_id,
                "alphafold_id"::VARCHAR AS alphafold_id,
                "proteome"::VARCHAR AS proteome,
                "protein_seq"::VARCHAR AS protein_seq,
            FROM read_parquet('{files_path}proteins/*.parquet')
            """)

    @staticmethod
    def _create_taxa_pairs_table(conn, files_path):
        conn.execute(f"""
            CREATE OR REPLACE TEMP TABLE taxa_alignment AS SELECT *
            FROM read_parquet('{files_path}taxa_pairs/alignment/*.parquet')""")
        conn.execute(f"""CREATE OR REPLACE TEMP TABLE taxa_pair_labels AS 
            SELECT * FROM read_parquet('{files_path}taxa_pairs/pair_labels/*.parquet')""")
        conn.execute(f"""CREATE OR REPLACE TABLE taxa_pairs AS 
            SELECT * FROM taxa_alignment
            INNER JOIN taxa_pair_labels ON (taxa_alignment.__index_level_0__=taxa_pair_labels.__index_level_0__)""")
        conn.commit()
        conn.execute("ALTER TABLE taxa_pairs RENAME COLUMN query_id TO thermo_taxid")
        conn.execute("ALTER TABLE taxa_pairs RENAME COLUMN subject_id TO meso_taxid")
        conn.commit
        conn.execute("""
            ALTER TABLE taxa_pairs DROP COLUMN "__index_level_0__:1"
        """)
        conn.execute("ALTER TABLE taxa_pairs DROP COLUMN __index_level_0__")
        conn.execute("ALTER TABLE taxa_pairs ALTER COLUMN thermo_taxid TYPE BIGINT")
        conn.execute("ALTER TABLE taxa_pairs ALTER COLUMN meso_taxid TYPE BIGINT")
        conn.commit()
    
    @staticmethod
    def _create_protein_pairs_table(conn, files_path):
        conn.execute(f"""CREATE OR REPLACE TABLE pairs AS SELECT 
                meso_pid::VARCHAR AS meso_pid,
                thermo_pid::VARCHAR AS thermo_pid,
                meso_taxid::INT AS meso_taxid,
                thermo_taxid::INT AS thermo_taxid,
                local_gap_compressed_percent_id::FLOAT AS local_gap_compressed_percent_id,
                scaled_local_query_percent_id::FLOAT AS scaled_local_query_percent_id,
                scaled_local_symmetric_percent_id::FLOAT AS scaled_local_symmetric_percent_id,
                local_E_value::FLOAT AS local_E_value,
                query_align_start::INT AS query_align_start,
                query_align_end::INT AS query_align_end,
                subject_align_end::INT AS subject_align_end,
                subject_align_start::INT AS subject_align_start,
                query_align_len::INT AS query_align_len,
                query_align_cov::FLOAT AS query_align_cov,
                subject_align_len::INT AS subject_align_len,
                subject_align_cov::FLOAT AS subject_align_cov,
                bit_score::INT AS bit_score,  
            FROM read_parquet('{files_path}/protein_pairs/*.parquet')""")
        conn.commit()

    @staticmethod
    def _create_indexes(conn):

        # drop existing indexes
        conn.execute("DROP INDEX IF EXISTS taxa_primary")
        conn.execute("DROP INDEX IF EXISTS protein_primary")
        conn.execute("DROP INDEX IF EXISTS protein_to_taxa")
        conn.execute("DROP INDEX IF EXISTS taxa_pairs_to_taxa_meso")
        conn.execute("DROP INDEX IF EXISTS taxa_pairs_to_taxa_thermo")
        conn.execute("DROP INDEX IF EXISTS protein_pairs_to_proteins_meso")
        conn.execute("DROP INDEX IF EXISTS protein_pairs_to_proteins_thermo")
        conn.execute("DROP INDEX IF EXISTS protein_pairs_to_taxa_meso")
        conn.execute("DROP INDEX IF EXISTS protein_pairs_to_taxa_thermo")

        # primary indexes
        conn.execute("CREATE UNIQUE INDEX taxa_primary ON taxa (taxid)")
        logger.debug("Created taxa primary index")
        conn.execute("CREATE UNIQUE INDEX protein_primary ON proteins (pid)")
        logger.debug("Created protein primary index")

        # foreign indexes
        conn.execute("CREATE INDEX protein_to_taxa ON proteins (taxid)")
        logger.debug("Created protein to taxa index")
        conn.execute("CREATE INDEX taxa_pairs_to_taxa_meso ON taxa_pairs (meso_taxid)")
        logger.debug("Created taxa pairs to taxa meso index")
        conn.execute("CREATE INDEX taxa_pairs_to_taxa_thermo ON taxa_pairs (thermo_taxid)")
        logger.debug("Created taxa pairs to taxa thermo index")
        conn.execute("CREATE INDEX protein_pairs_to_proteins_meso ON pairs (meso_pid)")
        logger.debug("Created protein pairs to proteins meso index")
        conn.execute("CREATE INDEX protein_pairs_to_proteins_thermo ON pairs (thermo_pid)")
        logger.debug("Created protein pairs to proteins thermo index")
        conn.execute("CREATE INDEX protein_pairs_to_taxa_meso ON pairs (meso_taxid)")
        logger.debug("Created protein pairs to taxa pairs meso index")
        conn.execute("CREATE INDEX protein_pairs_to_taxa_thermo ON pairs (thermo_taxid)")
        logger.debug("Created protein pairs to taxa pairs thermo index")
        conn.commit()

    def execute(self, sql_statement: str):
        """Execute a sql statement on the database.
        
        Parameters
        ----------
        sql_statement : str
            statement to execute
        
        Returns
        -------
        DataFrame
        """
        return self.conn.execute(sql_statement).df()

    @property
    def table_schema(self):
        """View of tables in the schema"""
        return self.execute("SELECT * FROM information_schema.tables")

    @property
    def column_schema(self):
        """View of columns in the schema"""
        return self.execute("SELECT * FROM information_schema.columns")

    def print_metadata(self):
        output = f"Metadata for Learn2Therm database at: {self.path}"
        # get number of datapoints
        n_taxa = self.execute("SELECT COUNT(taxid) FROM taxa").iloc[0,0]
        n_taxa_with_OGT = self.execute("SELECT COUNT(taxid) FROM taxa WHERE temperature IS NOT NULL").iloc[0,0]
        n_proteins = self.execute("SELECT COUNT(proteins.pid) FROM proteins INNER JOIN taxa ON (proteins.taxid=taxa.taxid) WHERE taxa.temperature IS NOT NULL").iloc[0,0]
        n_taxa_pairs = self.execute("SELECT COUNT(*) FROM taxa_pairs WHERE is_pair").iloc[0,0]
        n_protein_pairs = self.execute("SELECT COUNT(*) FROM protein_pairs").iloc[0,0]

        # add to output string
        output = output+f"\nNumber of taxa: {n_taxa}"
        output = output+f"\nNumber of taxa with OGT labels: {n_taxa_with_OGT}"
        output = output+f"\nNumber proteins with labeled taxa: {n_proteins}"
        output = output+f"\nNumber of taxa pairs identified: {n_taxa_pairs}"
        output = output+f"\nNumber of protein pair hits meeting minumum requirements: {n_protein_pairs}"
        print(output)

