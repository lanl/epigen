#for each chromosome get a list of known genes present at the https://genome.ucsc.edu and it is starting and ending bp position in the hg38 assembly
import pandas as pd
import os, io
import numpy as np

def get_gene_name(start_p, end_p, chromosome):
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
    if chr_i == 23: #this is chromosome X
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
