import numpy as np
import numpy.ma as ma
from matplotlib import pyplot as plt
import scipy
import scipy.spatial
import scipy.cluster
import pickle
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from PDBClean import pdbclean_io as pcio
from PDBClean import pdbclean_cifutils as cifutils
from PDBClean import pdbclean_viz as viz
from PDBClean import biputils as bip
from PDBClean import svdutils as svd

def reduce_feature_keep_samples(feature, sample, verbose=False, show=False):
    """
    """
    map = cifutils.keychain_to_atommap(feature, sample, verbose=verbose)
    if show:
        viz.show_atommap(map)
    feature_score   = np.sum(map, axis=1)
    feature_mask    = np.where(feature_score < map.shape[1], True, False)
    feature_reduced = ma.masked_array(feature, mask=feature_mask).compressed()
    if show:
        map = cifutils.keychain_to_atommap(feature_reduced, sample, verbose=verbose)
        viz.show_atommap(map)
    return feature_reduced

def reduce_optimized(feature, sample, verbose=False, show=False):
    """
    """
    map = cifutils.keychain_to_atommap(feature, sample, verbose=verbose)
    if show:
        viz.show_atommap(map)
    status, feature_index, sample_index, map_new = bip.bip_ilp(map)
    if verbose:
        print('Solution is {0}. Score went from {1} to {2}. We keep {3} / {4} samples and {5} / {6} features.'.format(status, 
                                                                         np.sum(map), 
                                                                         np.sum(map_new), 
                                                                         np.sum(sample_index), 
                                                                         map.shape[1],
                                                                         np.sum(feature_index), 
                                                                         map.shape[0]))
    feature_mask    = np.where(feature_index < 1, True, False)
    feature_reduced = ma.masked_array(feature, mask=feature_mask).compressed()
    sample_mask     = np.where(sample_index < 1, True, False)
    sample_array = np.array(sample)
    sample_reduced  = ma.masked_array(sample_array, mask=sample_mask).compressed()
    if show:
        viz.show_atommap(map_new)
        map_reduced = cifutils.keychain_to_atommap(feature_reduced, sample_reduced.tolist(), verbose=verbose)
        viz.show_atommap(map_reduced)
    return feature_reduced, sample_reduced.tolist()

def assign_clusters(keychain, input_list, n_clusters=-1, cutoff=-1, path=None, verbose=False, show=False):
    """
    assign_clusters
    """
    if path is None:
        print("Warning! provide path...")
    #
    file_clusters = open(path+'/clusters.pkl', 'rb')
    clusters = pickle.load(file_clusters)
    if(n_clusters<=0):
        if(cutoff<=0):
            file_cutoff = open(path+'/cutoff.pkl', 'rb')
            cutoff = pickle.load(file_cutoff)['cutoff']
        n_clusters = get_nclusters(clusters, cutoff)
    assignment = get_assignment(clusters, n_clusters)
    figname=path+'/assignment_stats.png'
    if verbose:
        print('... writing {}'.format(figname))
    unique, counts = plot_assignment_stats(assignment, figname=figname)
    for i in unique:
        pcio.check_project(projdir=path, level='cluster'+str(i), action='create', verbose=True)
    return assignment

def cluster(feature, sample, path=None, verbose=False, show=False):
    """
    cluster
    """
    if path is None:
        path='.'
    #
    map = cifutils.keychain_to_atommap(feature, sample, verbose=verbose)
    if show:
        viz.show_atommap(map)
    map_denoised = svd.denoise(map)
    #
    distance = scipy.spatial.distance.pdist(map.T)
    cutoff_dict = {"cutoff": np.amax(scipy.spatial.distance.squareform(distance))}
    file_cutoff = open(path+'/cutoff.pkl', 'wb')
    if verbose:
        print('... writing {}'.format(file_cutoff))
    pickle.dump(cutoff_dict, file_cutoff)
    file_cutoff.close()
    #
    clusters = scipy.cluster.hierarchy.linkage(distance, method='ward')
    #
    figname=path+'/dendogram.png'
    if verbose:
        print('... writing {}'.format(figname))
    plot_clusters(clusters, figname=figname)
    #
    file_clusters = open(path+'/clusters.pkl', 'wb')
    if verbose:
        print('... writing {}'.format(file_clusters))
    pickle.dump(clusters, file_clusters)
    file_clusters.close()

def get_nclusters(clusters,cutoff):
    """ get_nclusters
    """
    n_clusters=1
    keepongoing=True
    n = clusters.shape[0]
    for i in np.arange(1,n): #np.arange(0,n):
        score = clusters[i,2] - clusters[i-1,2]
        if(score > cutoff and keepongoing):
            n_clusters = n-i + 1
            print("Number of domains: ",n_clusters," (",score,")")
            keepongoing=False
    return n_clusters

def get_assignment(l,n_clusters):
    return scipy.cluster.hierarchy.fcluster(l, t=n_clusters, criterion='maxclust')

def plot_clusters(clusters,figsize=12,figname=''):
    fig = plt.figure(figsize=(figsize, figsize), dpi= 160, facecolor='w', edgecolor='k')
    nrow=2
    ncol=1
    # look at the number of natural clusters using the linkage object
    variance = clusters[:,2][::-1]
    plt.subplot(nrow,ncol,1)
    plt.scatter(np.arange(variance.shape[0]), variance)
    plt.title('Objective function change', fontsize=15)
    plt.ylabel('Variance', fontsize=13)
    plt.xlabel('Number of macrostates', fontsize=13)
    # see dendogram
    plt.subplot(nrow,ncol,2)
    scipy.cluster.hierarchy.dendrogram(clusters)
    plt.xticks(fontsize=13)
    plt.title('Ward linkage', fontsize=15)
    #
    plt.tight_layout()
    plt.show()
    if(figname):
        fig.savefig(figname)

def plot_assignment_stats(assignment, figsize=4, figname=''):
    """
    """
    unique, counts = np.unique(assignment, return_counts=True)
    fig = plt.figure(figsize=(figsize, figsize), dpi= 160, facecolor='w', edgecolor='k')
    plt.title('fraction of total population in cluster (N={0})'.format(assignment.shape[0]))
    plt.plot(unique, counts/assignment.shape[0], 'o')
    plt.show()
    if(figname):
        fig.savefig(figname)
    return unique, counts
