# get_ret_df_files.py
#
# This file maps locations of known genes from gene ontology databases
# Gene names are assigned from genes present at the https://genome.ucsc.edu
# using starting and ending base pair positions in hg38.
#
# Input file: /data/hg38length.txt
#
# Output files: <Chromosome>_ret_def.csv
#
# Later in the analysis we need to extract the names of important locations
# in the genome identified in the analysis. This file keeps track of those names.
#

import pandas as pd
import io
import numpy as np

def get_gene_name(start_p, end_p, chromosome):
    """
    Function to get genes associated with a given bin.

    """
    url = 'https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1355398003_PxotGFFkKyc9ASqRXr6WLvia6cuT'    
    session = requests.Session()
    res_str = str(chromosome + ':' + str(start_p) + '-' + str(end_p))
    params = {
        'hgsid': '1355398003_PxotGFFkKyc9ASqRXr6WLvia6cuT',
        'jsh_pageVertPos': '0',
        'clade': 'mammal',
        'org': 'Human',
        'db': 'hg38',
        'hgta_group': 'allTables',
        'hgta_table': 'refFlat',
        'hgta_regionType': 'range',
        'position': res_str,
        'hgta_outputType': 'primaryTable',
        'boolshad.sendToGalaxy': '0',
        'boolshad.sendToGreat': '0',
        'boolshad.sendToGenomeSpace': '0',
        'hgta_outFileName': '',
        'hgta_compressType': 'none',
        'hgta_doTopSubmit': 'get output'
    }
    response = session.post(url, data=params)
    gdf = pd.read_csv(io.StringIO(response.content.decode('utf-8')), sep='\t')
    return gdf

chr_length = pd.read_csv("./data/hg38length.txt", sep = '\t', header = None)
chr_length.columns = ['Chromosome', 'Total length (bp)', 'Nm', 'Nmm']
start_p = 0
divisor = 200
for i in range(1, 24):
    chr_i = i
    if chr_i == 23: # This is chromosome X
      chr_i = 'X'
    chromosome = "chr" + str(chr_i)

    end_p = (int(chr_length.loc[chr_length['Chromosome'] == chr_i, "Total length (bp)"])// divisor) * divisor
    gdf = get_gene_name(start_p, end_p, chromosome)
    start_df = pd.DataFrame()
    for t in np.unique(gdf['#geneName']):
        temp = gdf[gdf['#geneName'] == t]
        data = {'GeneName': [t], 'start': np.min(temp['txStart']), 'end': np.max(temp['txEnd'])}
        ret_df = pd.concat([start_df, pd.DataFrame(data)], axis = 0).reset_index(drop=True)
        start_df = ret_df

    ret_df.to_csv(chromosome + '_ret_df.csv', header=True, index=False, sep = "\t")
