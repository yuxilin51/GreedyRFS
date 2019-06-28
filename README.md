# GreedyRFS
GreedyRFS is a heuristic method developed for the Robinson-Foulds Supertree problem. It takes a set of binary source trees which have overlapping leaf sets as input and outputs a binary supertree on the union of leaves. This heuristic relies on an exact algorithm for the Robinson-Foulds Supertree problem on two source trees. It repeatedly chooses two trees from the source trees (according to some greedy criterion) and merges them using the exact algorithm until only one supertree remains. (A paper detailing the algorithm and experimental results is submitted to Recomb-CG and will be available in the future.)

To use this method, first make sure that you have the following prerequisites for the code to compile and run: Python 3, two python packages -- [NetworkX](https://networkx.github.io/) and [Dendropy](https://dendropy.org/). Note that we require Python 3 as NetworkX only works with Python 3. 

You can either clone the git repository or download and extract the zip file. Then you can use the following command: 
```
Python GreedyRFS.py -t <input trees filename> -o <output filename>
```

The input trees have to be stored in one file: one tree per line in Newick format. The input trees have to be binary and have at least 4 overlapping leaves. 
