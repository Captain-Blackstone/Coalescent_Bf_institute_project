#!/ebio/ag-neher/share/programs/EPD/bin/python
'''
author:     Taylor Kessinger & Richard Neher
date:       10/07/2014
content:    generate beta coalescent trees and calculate their SFS
'''
import os
import numpy as np
import random as rand
import scipy.special as sf
from Bio import Phylo
from betatree import *

def logit(x):
    return np.log(x/(1-x))

class SFS(betatree):
    '''
    class the generates many trees and accumulates an SFS.
    trees are generated by the inherited betatree class
    '''
    def __init__(self, sample_size, alpha=2, length_of_selection=None):
        betatree.__init__(self, sample_size, alpha)
        self.alleles = []
        self.sfs = None
        self.length_of_selection = length_of_selection

    def glob_trees(self, ntrees=10):
        '''
        generate many trees, accumulate the SFS
        parameters:
        ntrees -- number of trees to generate
        '''
        self.alleles=[]
        for ti in range(ntrees):
            self.coalesce()
            self.alleles.append([(clade.weight, clade.branch_length)
                                 for clade in self.BioTree.get_terminals()
                                 +self.BioTree.get_nonterminals()])

    def getSFS(self, ntrees = 10):
        '''
        calculate an SFS based on ntrees trees
        ntrees -- number of trees used to calculate the average SFS
        '''
        self.sfs = np.zeros(self.n+1)
        # loop over n trees and increment the sfs
        for ti in range(ntrees):
            self.coalesce()
            for clade in self.BioTree.get_terminals()+self.BioTree.get_nonterminals():
                self.sfs[clade.weight]+=clade.branch_length

        self.sfs/=ntrees


    def binSFS(self, mode = 'logit', bins=10):
        '''
        use the precalcutated SFS and bin it.
        mode -- one of linear, log, or logit. a binning with bins bins will be generated
        bins -- a user-specified binning if bins is iterable, otherwise the number of bins
        '''
        if np.iterable(bins):
            nbins = np.array(bins)
        else:
            if mode=='logit':
                nbins = np.exp(0.9*np.log(self.n)*np.linspace(-1,1,bins+1))
                nbins = nbins/(1.0+nbins)
                self.bin_center = np.sqrt(nbins[1:]*nbins[:-1])
            if mode=='linear':
                nbins=np.linspace(0,1,bins+1)
                self.bin_center = 0.5*(nbins[1:]+nbins[:-1])
            if mode=='log':
                nbins = np.exp(0.9*np.log(self.n)*np.linspace(-1,0,bins+1))
                self.bin_center = np.sqrt(nbins[1:]*nbins[:-1])

        self.bin_width = nbins[1:]-nbins[:-1]

        self.binned_sfs, tmp = np.histogram(np.linspace(0,1,self.n+1), 
                                            weights = self.sfs, bins = nbins)
        self.binned_sfs/=self.bin_width


    def saveSFS(self,fname):
        '''
        writes and existing sfs to the file fname. if fname ends on gz, the file 
        will be gzipped. if no sfs has been calculated, it will call getSFS()
        '''
        if self.sfs is None:
            self.getSFS()        
        np.savetxt(fname, self.sfs)


    def loadSFS(self,fname, alpha=None):
        '''
        load a previously saved SFS from file. checks for one-d vector, uses length of the 
        vector as SFS as sample size
        '''
        if os.path.isfile(fname):
            tmp = np.loadtxt(fname)
        else:
            print("file",fname,"does not exist")

        if len(tmp.shape)==1:
            self.n = tmp.shape[0]-1
            self.sfs = tmp
            self.alpha=alpha
        else:
            self.sfs=None
            print("expect a one dimensional vector, got object with shape",tmp.shape)




if __name__=='__main__':
    import matplotlib.pyplot as plt
    file_ending = '.dat.gz'
    calc = False
    # for alpha in [2,1.5, 1]:
    alpha = 2
    for length_of_selection in [i/10 for i in range(11)]:
        n = 500
        mySFS = SFS(n, alpha=alpha, length_of_selection=length_of_selection)
        if calc:
            print("calculating spectra for length_of_selection = ", length_of_selection)
            mySFS.getSFS(ntrees=1000)
            mySFS.saveSFS('../example_SFS/sfs_'+'_'.join(map(str, ['ls', length_of_selection, 'n', n]))
                          +file_ending)
        else:
            print("loading spectra for alpha =",alpha)
            mySFS.loadSFS('../example_SFS/sfs_'+'_'.join(map(str, ['ls', length_of_selection, 'n', n]))
                          +file_ending)

        color = str(length_of_selection)
        if length_of_selection == 0:
            color = "red"
        elif length_of_selection == 1:
            color = "green"

        plt.figure('regular')
        mySFS.binSFS(mode='logit', bins=20)
        plt.plot(mySFS.bin_center, mySFS.binned_sfs, label='t='+str(length_of_selection), lw=2, color=color)

        plt.figure('log')
        mySFS.binSFS(mode='log', bins=20)
        plt.plot(mySFS.bin_center, mySFS.binned_sfs, label='t='+str(length_of_selection), lw=2, color=color)

        plt.figure('logit')
        mySFS.binSFS(mode='logit', bins=20)
        plt.plot(logit(mySFS.bin_center), mySFS.binned_sfs, label='t='+str(length_of_selection), lw=2, color=color)

    #make linear binning figure
    plt.figure('regular')
    plt.title('linear binning')
    # plt.plot(mySFS.bin_center, 2.0/mySFS.bin_center, label=r'$1/x$', c='k', ls='--', lw=2)
    plt.yscale('log')
    plt.xlim(0,1)
    # plt.legend() #loc=9
    plt.xlabel('derived allele frequency')
    plt.xlabel('site frequency spectrum')
    plt.savefig("regular.png", dpi=200)

    plt.figure('log')
    plt.title('exponential binning')
    # plt.plot(mySFS.bin_center, 2.0/mySFS.bin_center, label=r'$1/x$', c='k', ls='--', lw=2)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(1.0/mySFS.n,1)
    # plt.legend(loc=9)
    plt.xlabel('derived allele frequency')
    plt.xlabel('site frequency spectrum')
    # plt.savefig('../example_SFS/sfs_logit_binning.pdf')
    plt.show()
    plt.savefig("log.png")

    plt.figure('logit')
    plt.title('logistic binning')
    plt.yscale('log')
    # plt.plot(logit(mySFS.bin_center), 2.0/mySFS.bin_center, label=r'$1/x$', c='k', ls='--', lw=2)
    plt.xlim(-np.log(n),np.log(n))
    tick_locs = np.array([0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999])
    plt.xticks(logit(tick_locs), map(str,tick_locs))
    plt.legend(loc=9)
    plt.xlabel('derived allele frequency')
    plt.ylabel('site frequency spectrum')
    # plt.savefig('../example_SFS/sfs_logit_binning.pdf')
    # plt.savefig('../example_SFS/sfs_logit_binning.png')
    plt.savefig("logit.png")