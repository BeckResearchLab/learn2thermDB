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

        t1 = time.time()
        cls._create_indexes(conn)
        t2 = time.time()
        logger.info(f"Took {(t2-t1)/60}m to create speed indexes")

        conn.close()        
        return cls(db_path, read_only=read_only)

    @staticmethod
    def _create_taxa_table(conn, files_path):
        # schema for taxa
        conn.execute("""
            CREATE OR REPLACE TABLE taxa(
                taxa_index INT PRIMARY KEY NOT NULL,
                ncbi_taxid INT NOT NULL,
                record_name STRING,
                filepath STRING,
                taxonomy STRING,
                organism STRING,
                bacdive_id INT,
                ogt_scraped_string STRING,
                seq_16srRNA STRING,
                len_16s INT,
                ogt FLOAT,
                thermophile_label BOOL,
            )""")
        conn.execute(f"""CREATE OR REPLACE TEMP TABLE taxa_tmp AS 
            SELECT
                "column0"::INT AS taxa_index,
                "taxid"::INT AS ncbi_taxid,
                "record_name"::STRING AS record_name,
                "filepath"::STRING AS filepath,
                "taxonomy"::STRING AS taxonomy,
                "organism"::STRING AS organism,
                "bacdive_id"::INT AS bacdive_id,
                "ogt_raw":: STRING AS ogt_scraped_string
            FROM read_csv_auto('{files_path}taxa/taxa_info_and_ogt.csv', header=True)""")
        conn.execute(f"""CREATE OR REPLACE TEMP TABLE taxa_16s_tmp AS 
            SELECT 
                "taxa_index"::INT AS taxa_index_1,
                "seq_16srRNA"::STRING AS seq_16srRNA,
                length(seq_16srRNA)::INT as len_16s
            FROM read_csv_auto('{files_path}taxa/16s_rRNA.csv', header=True, nullstr='None')""")
        conn.execute(f"""CREATE OR REPLACE TEMP TABLE taxa_labels_tmp AS 
            SELECT *
            FROM read_csv_auto('{files_path}taxa/labels.csv', header=True)""")
        conn.execute("""
            COPY (SELECT * EXCLUDE (column0, taxa_index_1) FROM taxa_tmp
                INNER JOIN taxa_16s_tmp ON (taxa_tmp.taxa_index=taxa_16s_tmp.taxa_index_1)
                INNER JOIN taxa_labels_tmp ON (taxa_tmp.taxa_index=taxa_labels_tmp.column0))
            TO './taxa_joined.csv' WITH (HEADER 1, DELIMITER '|')
        """)
        conn.execute("""
            COPY taxa FROM './taxa_joined.csv' ( HEADER, DELIMITER '|' )""")
        os.remove('./taxa_joined.csv')
    
    @staticmethod
    def _create_proteins_table(conn, files_path):
        # create table
        conn.execute(f"""
            CREATE OR REPLACE TABLE proteins AS SELECT 
                substr(seq_id, 0, strpos(seq_id, '.'))::INT AS taxa_index,
                "seq_id"::STRING AS protein_index,
                "protein_seq"::STRING AS protein_seq,
                "protein_desc"::STRING AS protein_desc,
                "protein_len"::INT AS protein_len
            FROM read_csv("{files_path}taxa/proteins/taxa_*.csv", auto_detect=False, header=True, sep=';', columns={{'seq_id': 'STRING', 'protein_seq': 'STRING', 'protein_desc': 'STRING', 'protein_len': 'INT'}})
            """)
        
        # at new integer index
        conn.execute("CREATE SEQUENCE protein_int_index_seq START 1")
        conn.execute("""
            ALTER TABLE proteins ADD COLUMN protein_int_index INT DEFAULT nextval('protein_int_index_seq')-1
        """)

    @staticmethod
    def _create_taxa_pairs_table(conn, files_path):
        conn.execute(f"""
            CREATE OR REPLACE TEMP TABLE taxa_pairs_tmp AS SELECT 
                *
            FROM read_csv_auto('{files_path}taxa_pairs/pairwise_16s_blast.csv', header=True)
        """)
        # now get labels, we might as well make this one table now
        conn.execute(f"""
            CREATE OR REPLACE TEMP TABLE taxa_pair_labels_tmp AS SELECT 
                column0:: INT AS taxa_pair_index,
                is_pair:: BOOL AS is_pair
            FROM read_csv_auto('{files_path}taxa_pairs/pair_labels.csv', header=True)
        """)
        # add new because it was only tracked in one of the files
        conn.execute("CREATE SEQUENCE taxa_pair_id_seq START 1")
        conn.execute("""
            ALTER TABLE taxa_pairs_tmp ADD COLUMN taxa_pair_index INT DEFAULT nextval('taxa_pair_id_seq')-1
        """)
        # join the two tables and create a true table
        conn.execute("""
            CREATE OR REPLACE TABLE taxa_pairs AS SELECT * FROM taxa_pairs_tmp
                INNER JOIN taxa_pair_labels_tmp ON (taxa_pairs_tmp.taxa_pair_index=taxa_pair_labels_tmp.taxa_pair_index)
        """)
        conn.execute("""
            ALTER TABLE taxa_pairs DROP "taxa_pair_index:1"
        """).df()
    
    @staticmethod
    def _create_protein_pairs_table(conn, files_path):
        conn.execute(f"""
            CREATE OR REPLACE TABLE protein_pairs AS SELECT  
                "thermo_protein_id"::STRING AS thermo_protein_index,
                "meso_protein_id"::STRING AS meso_protein_index,
                "local_gap_compressed_percent_id"::FLOAT AS local_gap_compressed_percent_id,
                "scaled_local_query_percent_id"::FLOAT AS scaled_local_query_percent_id,
                "scaled_local_symmetric_percent_id"::FLOAT AS scaled_local_symmetric_percent_id,
                "local_E_value"::FLOAT AS local_E_value,
                "query_align_start"::INT AS query_align_start,
                "query_align_end"::INT AS query_align_end,
                "subject_align_end"::INT AS subject_align_end,
                "subject_align_start"::INT AS subject_align_start,
                "query_align_len"::INT AS query_align_len,
                "query_align_cov"::FLOAT AS query_align_cov,
                "subject_align_len"::INT AS subject_align_len,
                "subject_align_cov"::FLOAT AS subject_align_cov,
                "bit_score"::INT AS bit_score,
                substr(thermo_protein_id, 0, strpos(thermo_protein_id, '.'))::INT AS thermo_index,
                substr(meso_protein_id, 0, strpos(meso_protein_id, '.'))::INT AS meso_index,
            FROM read_csv('{files_path}taxa_pairs/protein_alignment/taxa_pair*.csv', auto_detect=False, header=True, sep=',', columns={{
                'column0': 'INT',
                'thermo_protein_id': 'STRING',
                'meso_protein_id': 'STRING',
                'local_gap_compressed_percent_id': 'FLOAT',
                'scaled_local_query_percent_id': 'FLOAT',
                'scaled_local_symmetric_percent_id': 'FLOAT',
                'local_E_value': 'FLOAT',
                'query_align_start': 'INT',
                'query_align_end': 'INT',
                'subject_align_end': 'INT',
                'subject_align_start': 'INT',
                'query_align_len': 'INT',
                'query_align_cov': 'FLOAT',
                'subject_align_len': 'INT',
                'subject_align_cov': 'FLOAT',
                'bit_score': 'INT',
            }})""")
        conn.execute("CREATE SEQUENCE protein_pair_id_seq START 1")
        conn.execute("""
            ALTER TABLE protein_pairs ADD COLUMN prot_pair_index INT DEFAULT nextval('protein_pair_id_seq')-1
        """)
        # need to add the integer ids of meso and thermo proteins instead of just the slow string ones
        conn.execute("""
            ALTER TABLE protein_pairs ADD COLUMN meso_protein_int_index INT;
            ALTER TABLE protein_pairs ADD COLUMN thermo_protein_int_index INT;
            ALTER TABLE protein_pairs ADD COLUMN taxa_pair_index INT;
        """)
        conn.execute("""
            UPDATE protein_pairs SET meso_protein_int_index=(
                SELECT proteins.protein_int_index
                FROM proteins
                WHERE protein_pairs.meso_protein_index=proteins.protein_index
            )
        """)
        conn.execute("""
            UPDATE protein_pairs SET thermo_protein_int_index=(
                SELECT proteins.protein_int_index
                FROM proteins
                WHERE protein_pairs.thermo_protein_index=proteins.protein_index
            )
        """)
        conn.execute("""
            UPDATE protein_pairs SET taxa_pair_index=(
                SELECT taxa_pairs.taxa_pair_index
                FROM taxa_pairs
                WHERE protein_pairs.thermo_index=taxa_pairs.thermo_index AND protein_pairs.meso_index=taxa_pairs.meso_index
            )
        """)

    @staticmethod
    def _create_indexes(conn):
        # primary indexes
        conn.execute("CREATE UNIQUE INDEX taxa_primary ON taxa (taxa_index)")
        conn.execute("CREATE UNIQUE INDEX protein_primary ON proteins (protein_int_index)")
        conn.execute("CREATE UNIQUE INDEX taxa_pair_primary ON taxa_pairs (taxa_pair_index)")
        conn.execute("CREATE UNIQUE INDEX prot_pair_primary ON protein_pairs (prot_pair_index)")

        # foreign indexes
        conn.execute("CREATE INDEX protein_to_taxa ON proteins (taxa_index)")
        conn.execute("CREATE INDEX taxa_pair_to_meso ON taxa_pairs (meso_index)")
        conn.execute("CREATE INDEX taxa_pair_to_thermo ON taxa_pairs (thermo_index)")
        conn.execute("CREATE INDEX taxa_pair_both ON taxa_pairs (meso_index, thermo_index)")
        conn.execute("CREATE INDEX prot_pair_to_meso ON protein_pairs (meso_index)")
        conn.execute("CREATE INDEX prot_pair_to_thermo ON protein_pairs (thermo_index)")
        conn.execute("CREATE INDEX prot_pair_both ON protein_pairs (meso_index, thermo_index)")
        conn.execute("CREATE INDEX prot_pair_to_meso_prot ON protein_pairs (meso_protein_int_index)")
        conn.execute("CREATE INDEX prot_pair_to_thermo_prot ON protein_pairs (thermo_protein_int_index)")
        conn.execute("CREATE INDEX prot_pair_both_prot ON protein_pairs (meso_protein_int_index, thermo_protein_int_index)")
        conn.execute("CREATE INDEX prot_pair_to_taxa_pair ON protein_pairs (taxa_pair_index)")

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
        n_taxa = self.execute("SELECT COUNT(taxa_index) FROM taxa").iloc[0,0]
        n_taxa_with_OGT = self.execute("SELECT COUNT(taxa_index) FROM taxa WHERE ogt IS NOT NULL").iloc[0,0]
        n_proteins = self.execute("SELECT COUNT(proteins.protein_index) FROM proteins INNER JOIN taxa ON (proteins.taxa_index=taxa.taxa_index) WHERE taxa.ogt IS NOT NULL").iloc[0,0]
        n_taxa_pairs = self.execute("SELECT COUNT(taxa_pair_index) FROM taxa_pairs WHERE is_pair").iloc[0,0]
        n_protein_pairs = self.execute("SELECT COUNT(prot_pair_index) FROM protein_pairs").iloc[0,0]

        # add to output string
        output = output+f"\nNumber of taxa: {n_taxa}"
        output = output+f"\nNumber of taxa with OGT labels: {n_taxa_with_OGT}"
        output = output+f"\nNumber proteins with labeled taxa: {n_proteins}"
        output = output+f"\nNumber of taxa pairs identified: {n_taxa_pairs}"
        output = output+f"\nNumber of protein pair hits meeting minumum requirements: {n_protein_pairs}"
        print(output)

