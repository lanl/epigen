# epigen
Code and scripts related to the paper. O# O4747

**General info**: 

**data** folder contains files necessary to execute the code, some large files needed to be downloaded as written in the scripts below.

This repository contains python scripts and bash scripts to obtain the paper results. 

All of these scripts related to the hg38-aligned dataset; it is easy to get files and results for the hg19-aligned dataset by changing 38 to 19 in the code and 1632 (#samples in the hg 38) to 1479 (#samples in the hg 19).

The data doesn't contain the chromosome Y.

**The order of execution**:

1. Download bigWig files from the ENCODE, run _run_download_data_files_script_ that calls _download_data_files.py_.
2. For each chromosome calculate values for each 200bp region, run _run_get_values_for_chromosome_script_ that calls _get_values_for_chromosome.py_.
3. For each chromosome obtain the dataframe that contains columns where each value (p-values here) is less than 0.05 in at least one sample (row): _create_pvdf.py_ that calls _run_create_pvdf_script_.
4. For each chromosome calculate the pairwise correlation matrix between samples, requires GPU, otherwise will take hours vs minutes per chromosome: correlation_matrix_calculation.py.
5. For each chromosome obtain a list of known genes and their bp position in the chromosome: get_ret_df_files.py.
6. Produce paper circular plot 2: get_newick_string.py and produce_circular_plot.py.
7. Produce paper umap plots 3, 4, and 5: produce_umap_plots.py.
8. Produce paper entropy plots 6 and 7: produce_entropy_plots.py.
9. Produce paper co-occurence plot 8: produce_cooccurence_plot.py.
10. Produce paper gene ontology resuls needed for plots 9 and 10: produce_go_results.py.
11. Produce paper datasets comparison plots 11 and 12: produce_comparison_plots_datasets19vs38.py.



© 2024. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

 
This program is Open-Source under the BSD-3 License.
 
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
