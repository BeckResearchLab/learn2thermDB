import requests
import pandas as pd
import time
from Bio import SeqIO
from io import StringIO

from typing import List

def get_uniprot_sequences(uniprot_ids: List, chunksize=50) -> pd.DataFrame:
        """
        Retrieve uniprot sequences based on a list of uniprot sequence identifier.

        For large lists, they are chunked uo

        Parameters:
            uniprot_ids: List, list of uniprot identifier
            chunksize: int, number to query at once

        Returns:
            pd.DataFrame, pandas dataframe with uniprot id column and sequence
        """
        import requests

        def get_url(sublist):
            joined = ','.join(sublist)
            return f'https://rest.uniprot.org/uniprotkb/accessions?accessions={joined}&format=fasta'
        
        records = []
        print(len(uniprot_ids))
        for chunk in range(0, len(uniprot_ids), chunksize):
            sublist = uniprot_ids[chunk:chunk+chunksize]
            print(len(sublist))
            url = get_url(sublist)
            r = requests.get(url)
            r.raise_for_status()

            for record in SeqIO.parse(StringIO(r.text), 'fasta'):
                records.append(record)
            time.sleep(1)

        df = pd.DataFrame([{'id': record.id.split('|')[1], 'seq': str(record.seq)} for record in records])
        df = df.drop_duplicates(subset=['id'])
        return df