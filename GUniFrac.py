from optparse import OptionParser
import numpy as np
import pandas as pd
import preprocessing
from ete2 import Tree
import itertools
from scipy.spatial.distance import squareform

def importData(filename):
    table = pd.read_table(filename,index_col=0)
    relative_abundance = preprocessing.relative_abundance(table.values, axis=0)
    new_table = pd.DataFrame(relative_abundance, index=table.index, columns=table.columns)
    return new_table

def compute_GUniFrac(abundance,treefile, alpha=0.5, unweighted=False):
    n_samples = len(abundance.columns)
    n_distance = n_samples * (n_samples - 1) / 2
    d_array = np.zeros((n_distance))
    t = Tree(treefile,format=1)
    if set(t.get_leaf_names()) != set(abundance.index):
        print 'Error: OTU table contains unknown OTUs. All of OTU names in OTU table should be contained in tree file.'
        quit()
    
    for i,(sample1, sample2) in enumerate(itertools.combinations(abundance.columns, 2)):
        print 'calculating ',sample1,' vs. ',sample2,'...'
        denom = 0.0
        numer = 0.0
        for node in t.traverse():
            if node.is_root():
                continue
            else:
                p_a = 0.0
                p_b = 0.0
                for leaf in node.get_leaf_names():
                    if leaf in abundance.index:
                        p_a += abundance.loc[leaf,sample1]
                        p_b += abundance.loc[leaf,sample2]
            if p_a == 0.0 and p_b == 0.0:
                continue
            if unweighted:
                if p_a == 0.0 or p_b == 0.0:
                    numer += node.dist
                denom += node.dist
            else:
                denom += node.dist * (p_a + p_b) ** alpha
                numer += node.dist * (p_a + p_b) ** alpha * abs(p_a - p_b) / (p_a + p_b)
        d_array[i] = numer / denom
    return squareform(d_array)

if __name__ == '__main__':
    usage = 'usage: python %prog [options]'
    parser = OptionParser(usage)
    parser.add_option( "-f", "--file", action="store", dest="data_file", help="TAB-separated OTU table. rows are OTUs, columns are samples.")
    parser.add_option( "-t", "--tree", action="store", dest="tree_file", help="Rooted phylogenetic tree. NEWICK format.")
    parser.add_option( "-a", "--alpha", action="store", type="float", dest="alpha", default=0.5, help="alpha parameter of generalized UniFrac metric.")
    parser.add_option( "-u", "--unweighted", action="store_true", dest="unweighted", default=False, help="calculate unweighted UniFrac metric.")
    options, args = parser.parse_args()
    if options.data_file == None or options.tree_file == None:
        print "ERROR: requires options"
        parser.print_help()
        quit()
    datafile = options.data_file
    treefile = options.tree_file
    alpha = options.alpha
    unweighted = options.unweighted
    
    abundance = importData(datafile)
    distance_matrix = compute_GUniFrac(abundance, treefile, alpha=alpha, unweighted=unweighted)
    distance_table = pd.DataFrame(distance_matrix, index=abundance.columns, columns=abundance.columns)
    distance_table.to_csv('./result.dist',sep='\t')
