# download_data_files.py: download bigWig files for hg38 assembly from ENCODE
#
# input: genome_df38.csv file which has IDs to paste to ENCODE
# to download a file to the /hg38data directory
#
# how to run: takes two arguments: i=sys.argv[1] and j=sys.argv[2]
# It then downloads ith-jth range of files from the genome_df38.csv
#
# Output: one .bigwig file per index, approx 0.5 GB per file.
# If storage is a concern, you can process these files in chunks
# using this script and get_values_for_chromosome.py
#

import requests
import sys
import pandas as pd

from six.moves import urllib

# load dataframe containing file information
df = pd.read_csv("./data/genome_df38.csv", delimiter=",")

# can be replaced to (0, df.shape[0])
for file_num in range(int(float(sys.argv[1])), int(float(sys.argv[2]))):

    file_name = str(df['Accession'][file_num]) # exact ID of file for encode
    print("Downloading: ", file_num, file_name)

    # ifnormation about where to get annoations from
    requests.get("https://www.encodeproject.org")
    headers = {'accept': 'application/json'}
    url = str('https://www.encodeproject.org/annotations/') + str(file_name) + str('/')

    # Find the actual URL for this accession number
    response = requests.get(url, headers=headers)
    out = response.json()
    url_file = out['files'][0]['cloud_metadata']['url']

    # Get and store the bigwig file.
    result_file = str('./hg38data/') + str(file_name) + str('.bigWig')
    urllib.request.urlretrieve(url_file, result_file)


