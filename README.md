# TTsampler
This utility takes as input a phylogenetic tree and samples uniformly from the set of possible transmission trees, as decscribed [here](http://www.biorxiv.org/content/early/2017/07/08/160812).

It is primarily intended for those interested in exploring the structure of transmission tree space, and should not generally be used for phylodynamic inference.

Both an R package and a UNIX command-line script (TTSamplerCommandLine.R) are provided. Documentation for the former is [here](https://github.com/mdhall272/TTsampler/blob/master/TTsampler/TTsampler.pdf). The latter will not be updated further, but is retained and has equal functionality at present (October 2017).

TTSamplerCommandLine.R takes two compulsory arguments. The first should be the input tree in Newick format; it should be binary and since no manipulation is performed by the script, should be exactly the tree of interest. The second is simply a string identifying output files.

It should be run from the command line with Rscript. The following optional arguments are currently mutually exclusive. At present the script cannot handle all these variations at once, but this will be implemented eventually.

* -i Sample only transmission trees where hosts cease to be infected (and hence infectious) at sampling. 
* -I Used to specify a CSV file of minimum and maximum **heights** (not dates) for each host's infection. Assuming branch lengths are in calendar time, these should be given in backwards time with the zero point being the time of the youngest tip in the tree.
* -m Used to specify a CSV file that maps each tip to an identifier for the relevant host. This is for multiple sampling, and is unnecessary if only one or zero tips are present from any single host.
* -u Used to specify a number of unsampled hosts which had descendants within the sample.

The remaining optional arguments are:

* -s The number of samples to take. The default is 1.
* -d Use [ggtree](https://bioconductor.org/packages/release/bioc/html/ggtree.html) to draw each sample as an annotated tree in PDF format.
* -h Print the help message.

Standard output is two CSV files. The first, ending "\_tt.csv", lists every host (including, if -u is specified, unsampled hosts, which all begin "uh") in the first column, and the parent of each of these hosts for each sample replicate in the remaining columns. The second, ending "\_annotations.csv", if -u is not specified, gives node annotations (in the order used by ape when loading the tree file) for each sample; if -u is specified it also counts the number of unsampled individuals that occur in the transmission chain along the branch leading to each node (see p12 of the preprint).
