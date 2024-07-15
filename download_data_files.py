#download bigWig files for hg38 assembly from ENCODE
#input: genome_df38.csv file which has IDs to paste to ENCODE to download a file to the /hg38data directory
#how to run: takes sys.argv[1] and sys.argv[2] two arguments from the extrenal script to download ith-jth range of files from the genome_df38.csv

import requests
import json
import numpy as np
import six
import sys
import pandas as pd
from six.moves import urllib
df = pd.read_csv("./genome_df38.csv", delimiter=",")
for file_num in range(int(float(sys.argv[1])), int(float(sys.argv[2]))): #can be replaced to (0, df.shape[0])
    file_name = str(df['Accession'][file_num])
    print(file_num, file_name)
    requests.get("https://www.encodeproject.org")
    headers = {'accept': 'application/json'}
    url = str('https://www.encodeproject.org/annotations/') + str(file_name) + str('/')
    response = requests.get(url, headers=headers)
    out = response.json()
    url_file = out['files'][0]['cloud_metadata']['url']
    urllib.request.urlretrieve(url_file, str('./hg38data/') + str(file_name) + str('.bigWig'))


