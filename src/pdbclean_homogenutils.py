import numpy as np
import numpy.ma as ma
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from PDBClean import pdbclean_cifutils as cifutils
from PDBClean import pdbclean_viz as viz
import cvxopt
from cvxopt import glpk
from cvxopt import matrix

def reduce_feature_keep_samples(feature, sample, verbose=False, show=False):
    """
    """
    map = cifutils.keychain_to_atommap(feature, sample, verbose=verbose)
    if show:
        viz.show_atommap(map)
    feature_score   = np.sum(map, axis=1)
    feature_mask    = np.where(feature_score < map.shape[1], True, False)
    feature_reduced = ma.masked_array(feature, mask=feature_mask).compressed()
    map = cifutils.keychain_to_atommap(feature_reduced, sample, verbose=verbose)
    if show:
        viz.show_atommap(map)
    return feature_reduced

def reduce_optimized(feature, sample, verbose=False, show=False):
    """
    """
    map = cifutils.keychain_to_atommap(feature, sample, verbose=verbose)
    if show:
        viz.show_atommap(map)
    status, feature_index, sample_index, map_new = bip_ilp(map)
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

def bip_ilp(A):
    """
    """
    print('Setting up and solving Binary linear program.')
    m,n = A.shape
    G1 = build_G(A, constraint=1)
    G2 = build_G(A, constraint=2)
    G3 = build_G(A, constraint=3)
    G4 = build_G(A, constraint=4)
    h1 = build_h(A, constraint=1)
    h2 = build_h(A, constraint=2)
    h3 = build_h(A, constraint=3)
    h4 = build_h(A, constraint=4)
    M = assemble_G(A,G1,G2,G3,G4)
    b = assemble_h(A,h1,h2,h3,h4)
    w = assemble_w(A)
    c, G, h = matrix(w), matrix(M.T), matrix(b)
    status, sol = glpk.ilp(c, G.T, h, B=set(range(len(c))))
    solution = np.array(sol).astype('int')
    solution_row = solution[0:m]
    solution_col = solution[m:m+n]
    solution_mat = solution[m+n:].reshape((m,n))
    return status, 1-solution_row, 1-solution_col, 1-solution_mat

##################### 
# BIP ILP UTILITIES #
#####################

def assemble_w(A, mode='e'):
    """
    """
    m, n = A.shape
    if(mode=='all'):
        wr = 1
        wc = 1
        we = 1
    else:
        wr = 0
        if(mode=='r'):
            wr=1
        wc = 0
        if(mode=='c'):
            wc=1
        we = 0
        if(mode=='e'):
            we=1
    w = np.empty(n+m+n*m)
    w[0:m]   = wr
    w[m:m+n] = wc
    w[n+m:n+m+n*m] = we
    print('w.shape = {0}'.format(w.shape))
    return w.astype('double')

def assemble_h(A,h1,h2,h3,h4):
    """
    """
    m, n = A.shape
    nrows = m*n
    h = np.empty(4*nrows)
    h[0:nrows]         = h1
    h[nrows:2*nrows]   = h2
    h[2*nrows:3*nrows] = h3
    h[3*nrows:4*nrows] = h4
    print('h.shape = {0}'.format(h.shape))
    return -h.astype('double')

def assemble_G(A,G1,G2,G3,G4):
    """
    """
    m, n = A.shape
    nrows = m*n
    ncols = m + n + m*n
    G = np.empty((4*nrows,ncols))
    G[0:nrows,:]         = G1
    G[nrows:2*nrows,:]   = G2
    G[2*nrows:3*nrows,:] = G3
    G[3*nrows:4*nrows,:] = G4
    print('G.shape = {0}'.format(G.shape))
    return -G.astype('double')

def build_h(A, constraint=1): #n_sample=5, n_feature=20):
    """
    """
    m, n = A.shape
    nrows = m*n
    if(constraint==1):
        h = np.invert(A.reshape(nrows).astype('bool')).astype('int')
    else:
        h = np.zeros(nrows)
    #print('h.shape = {0}'.format(h.shape))
    return h

def build_G(A, constraint=1):
    """
    """
    m, n = A.shape
    nrows = m*n
    ncols = m + n + m*n
    Gr = build_Gr(m,n,constraint=constraint)
    Gc = build_Gc(m,n,constraint=constraint)
    Ge = build_Ge(m,n,constraint=constraint)
    G  = np.zeros((nrows,ncols))
    G[:,0:m]   = Gr
    G[:,m:m+n] = Gc
    G[:,m+n:]  = Ge
    #print('G.shape = {0}'.format(G.shape))
    return G

def build_Gr(m,n,constraint=1):
    nrows = m*n
    Gr = np.zeros((nrows,m))
    if(constraint==2):
        for i in np.arange(m):
            Gr[i*n:(i+1)*n,i] = 1
    elif(constraint==3):
        for i in np.arange(m):
            Gr[i*n:(i+1)*n,i] = -1
    #print('... Gr built')
    return Gr
        
def build_Gc(m,n,constraint=1):
    nrows = m*n
    Gc = np.zeros((nrows,n))
    if(constraint==2):
        for i in np.arange(m): 
            Gc[i*n:(i+1)*n,:] = np.identity(n)
    elif(constraint==4):
        for i in np.arange(m):
            Gc[i*n:(i+1)*n,:] = -np.identity(n)
    #print('... Gc built')
    return Gc

def build_Ge(m,n,constraint=1):
    nrows = m*n
    Ge = np.identity(nrows)
    if(constraint==2):
        Ge = -Ge
    #print('... Ge built')
    return Ge
