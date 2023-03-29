# Animal.SARS-CoV-2
***Acknowledgement:
The file titled "get_distance_to_focal_set_dsave.py" has been extracted from a block in Nextstrain (Hadfield et al. 2018) and used as a standalone scrirpt to generate
a proximity matrix. A few minor modifications have been added to add an additional input argument and save an intermideary file in the script. 

Nextstrain:
https://paperpile.com/c/t5H8a2/hIKC



cleanup_gisaid_animal_download.R: This script cleans up and filters the animal sequence data downloaded from GISAID. These filtering criteria include removing unwanted/ambiguous hosts, removing sequences with incomplete dates, removing sequences with too many n's.

prunne_iqtree_output.R: Thhis script prunes and cleans up the trees output by iqtree. The pruning is done using the IQR criterion.

transmission_permutation_test.R: This script does the permutation test on the transmission counts. First, it shuffles the tips of the tree and then finds the number of unfiltered transmission counts in the desired direction.
