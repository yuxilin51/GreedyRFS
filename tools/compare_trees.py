
"""
Robust comparison of two trees using dendropy.

Written by EKM (molloy.erin.k@gmail.com) in October 2016.
"""
import sys
import os
import dendropy
from dendropy.calculate.treecompare \
    import false_positives_and_negatives


def compare_trees(tr1, tr2):
    from dendropy.calculate.treecompare \
        import false_positives_and_negatives

    lb1 = set([l.taxon.label for l in tr1.leaf_nodes()])
    lb2 = set([l.taxon.label for l in tr2.leaf_nodes()])

    com = lb1.intersection(lb2)
    if com != lb1 or com != lb2:
        com = list(com)
        tns = dendropy.TaxonNamespace(com)

        tr1.retain_taxa_with_labels(com)
        tr1.migrate_taxon_namespace(tns)

        tr2.retain_taxa_with_labels(com)
        tr2.migrate_taxon_namespace(tns)
    com = list(com)

    tr1.update_bipartitions()
    tr2.update_bipartitions()

    nl = len(com)
    ei1 = len(tr1.internal_edges(exclude_seed_edge=True))
    ei2 = len(tr2.internal_edges(exclude_seed_edge=True))

    [fp, fn] = false_positives_and_negatives(tr1, tr2)
    rf = float(fp + fn) / (ei1 + ei2)

    return(nl, ei1, ei2, fp, fn, rf)


if __name__ == "__main__":
    argc = len(sys.argv)
    assert (argc == 3), "Usage: assess_dataset.py tr1 tr2"

    tax = dendropy.TaxonNamespace()
    tr1 = dendropy.Tree.get(path=sys.argv[1],
                            schema='newick',
                            rooting='force-unrooted',
                            taxon_namespace=tax)

    tr2 = dendropy.Tree.get(path=sys.argv[2],
                            schema='newick',
                            rooting='force-unrooted',
                            taxon_namespace=tax)

    tr1.collapse_basal_bifurcation(set_as_unrooted_tree=True)
    tr2.collapse_basal_bifurcation(set_as_unrooted_tree=True)

    [nl, ei1, ei2, fp, fn, rf] = compare_trees(tr1, tr2)
    sys.stdout.write('%d %d %d %d %d %f' % (nl, ei1, ei2, fp, fn, rf))
    sys.stdout.flush()
    os._exit(0)  # CRITICAL ON BLUE WATERS LOGIN NODE
