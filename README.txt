This program converts evolutionary similarity data into an evolutionary tree. The input data is in Phylip format, with each sequence having a similarity index comparing it to every other sequence.

The most similar groups are repeatedly merged using an agglomerative clustering algorithm. Each 'merge' represents a branch in the evolutionary tree.

The output data takes the form of a Newick tree. Each set of parentheses represents a 'merge,' and the numbers represent the evolutionary distance from the last branch. This machine-readable format is easily converted into a graphical evolutionary tree.


For visualizing a completed Newick tree: 
http://etetoolkit.org/treeview/

For more information on Phylip-formatted data:
https://en.wikipedia.org/wiki/PHYLIP

For more information on the Newick tree format: 
http://evolution.genetics.washington.edu/phylip/newicktree.html
