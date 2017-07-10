# TTsampler
This utility takes as input a phylogenetic tree and samples uniformly from the set of possible transmission trees, as decscribed [here](http://www.biorxiv.org/content/early/2017/07/08/160812).

It is primarily intended for those interested in exploring the structure of transmission tree space, and should not generally be used for phylodynamic inference.

The script takes two compulsory arguments. The first should be the input tree in Newick format; it should be binary and since no manipulation is performed by the script, should be exactly the tree of interest. The second is simply a string identifying output files.

The following optional arguments are currently mutually exclusive. At present the script cannot handle all these variations at once, but this will be implemented eventually.

* -i Sample only transmission trees where hosts cease to be infected (and hence infectious) at sampling. 
* -I Used to specify a CSV file of minimum and maximum **heights** (not dates) for each host's infection. At present these need to be given in units of branch lengths before the date of the last tip.
* -m Used to specify a CSV file that maps each tip to an identifier for the relevant host. This is for multiple sampling, and is unnecessary if only one or zero tips are present from any single host.
* -u Used to specify a number of unsampled hosts which had descendants within the sample.

The remaining optional arguments are:

* -s The number of samples to take.
* -d Use [ggtree](https://bioconductor.org/packages/release/bioc/html/ggtree.html) to draw each sample as an annotated tree in PDF format.
* -h Print the help message.

Standard output is two CSV files. The first, ending "\_tt.csv", lists every host (including, if -u is specified, unsampled hosts) in the first column, and the parent of each of these hosts for each sample replicate in the remaining columns. The second, ending "\_annotations.csv", if -u is not specified, gives node annotations (in the order used by ape when loading the tree file) for each sample; if -u is specified it also counts the number of unsampled individuals that occur in the transmission chain along the branch leading to that node.
