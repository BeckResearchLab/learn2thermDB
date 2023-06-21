# Code to create circos plot of learn2therm overlayed over ncbi taxonomy

Install requirements.txt in a new environment via `pip install -r requirements.txt`, python 3.9 was used here

Run the notebook to create the plot.

Note `snames_map.json` was created on 05.2023 by downloading the NCBI taxdump, here: 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
and parsing `names.dmp` at `filepath` with the following code:

```
import json
import os
from collections import defaultdict

col_delimiter = '\t|\t'
row_delimiter = '\t|\n'

other_names = defaultdict(set)
with open(filepath) as names_file:
    for line in names_file:
        line = line.rstrip(row_delimiter)
        values = line.split(col_delimiter)
        tax_id, name_txt, _, name_type = values[:4]
        if name_type == 'scientific name':
            scientific_names[tax_id] = name_txt
        else:
            other_names[tax_id].add(name_txt)
with open('snames_map.json'), 'w') as f:
    json.dump(scientific_names, f)
```