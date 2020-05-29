import betatree.betatree as betatree
import numpy as np
from Bio import Phylo
from collections import Counter
import matplotlib.pyplot as plt
from collections import defaultdict
import scipy.stats as stats
from pathlib import Path
import datetime


def test_stat(function, alpha, length_of_selection=None):
    stat_pack = []
    ntrees = 1000
    for i in range(ntrees):
        my_stat = None
        while my_stat is None:
            myT = betatree.betatree(500, alpha = alpha, length_of_selection=length_of_selection)
            myT.coalesce()
            tree = myT.BioTree
            my_stat = function(tree)
        stat_pack.append(my_stat)
    stat_pack = np.array(stat_pack).flatten() # in case of u_stat, where I get lists of 2 values from each tree
    if length_of_selection is None:
        color = "blue"
    elif length_of_selection == 0:
        color ="red"
    elif length_of_selection == 1:
        color = "green"
    else:
        color = str(round(length_of_selection, 2))
    density = stats.gaussian_kde(stat_pack)
    n, x, _ = plt.hist(stat_pack, bins=50, color="white", alpha=0.2, label=length_of_selection, density=True, histtype="step")
    plt.plot(x, density(x), color=color)

def u_stat(tree):
    return [len(child.get_terminals())/len(tree.get_terminals()) for child in tree.root] if len(tree.root) == 2 else None

def sackin_index(tree):
    return np.mean([len(tree.get_path(leaf)) for leaf in tree.get_terminals()])

def number_of_cherries(tree):
    return sum([1 if clade.weight == 2 else 0 for clade in tree.find_clades()])

def colless_index(tree):
    return sum([abs(clade[0].weight-clade[1].weight) if len(clade) == 2 else 0 for clade in tree.find_clades()])

def allele_frequency_spectrum(alpha=2, color="red", length_of_selection = None):
    spectrums = []
    epochs = 1000
    n_leafs = 500
    for trial in range(epochs):
        myT = betatree.betatree(n_leafs, alpha=alpha, length_of_selection=length_of_selection)
        myT.coalesce()
        tree = myT.BioTree
        spectrum = ((len(clade.get_terminals()), clade.branch_length) for clade in tree.find_clades())
        spectrums.append(spectrum)
    afs = np.zeros(n_leafs+1)
    for spectrum in spectrums:
        for freq, prob in spectrum:
            afs[freq] += prob
    afs100 = np.zeros(101)
    for i in range(len(afs)):
        afs100[int(i//(n_leafs/100))] += afs[i]
    for_hist = []
    for i in range(len(afs100)):
        for_hist.extend([i]*int(afs100[i]))
    for_hist = np.array(for_hist)/n_leafs
    n, x, _ = plt.hist(for_hist, bins=50, color="white", alpha=0.2, label=length_of_selection, density=True, histtype="step")
    density = stats.gaussian_kde(for_hist)
    if length_of_selection == 0:
        color ="red"
    elif length_of_selection == 1:
        color = "green"
    plt.plot(x, density(x), color=color)
    plt.yscale("log")
    plt.xscale("log")

def afs_for_different_selection_periods():
    selections = np.linspace(0, 1, 11)
    # selections = [0]
    for selection in selections:
        print(round(selection, 2))
        color = str(selection)
        allele_frequency_spectrum(alpha=2, color=color, length_of_selection=selection)
    plt.yscale("log")
    plt.xscale("log")
    plt.show()

def stats_for_different_selection_periods():
    selections = np.linspace(0, 1, 11).round(2)
    statistics = [u_stat, sackin_index, number_of_cherries, colless_index]
    for function in statistics:
        print(function.__name__, " is running...")
        for selection in selections:
            print("With selection peroid of ", selection)
            test_stat(function, alpha=2, length_of_selection=selection)
        plt.show()
        plt.savefig(f'{function.__name__}.png')


afs_for_different_selection_periods()
stats_for_different_selection_periods()