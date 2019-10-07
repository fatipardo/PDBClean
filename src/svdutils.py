import numpy as np
from scipy import linalg

def denoise(map, n_truncated=-1):
    """
    denoise
    """
    U, L, Vh = get_svd(map)
    if(n_truncated<=0):
        n_truncated = get_truncate_order(L)
    map_denoised = reassemble(U,L,Vh,n_truncated)
    return map_denoised
    
def center(map):
    """
    center: assuming map is [features]x[samples]
    """
    return (map.T - np.mean(map, axis=1)).T
    
def get_svd(map, centered=False):
    if centered:
        map_centered = map
    else:
        map_centered = center(map)
    return linalg.svd(map_centered)

def get_truncate_order(L,threshold=0.9):
    var = L**2
    var /= np.sum(var)
    n_components = 1
    if(threshold==1.0):
        n_components = len(var)
    else:
        for i in np.arange(0,len(var)-1,1):
            var_current = np.cumsum(var)[i]
            var_next = np.cumsum(var)[i+1]
            if(var_current < threshold and var_next > threshold):
                n_components=i+1
    return n_components

def reassemble(U,L,Vh,n=-1):
    if(n==-1):
        n=L.shape[0]
    return np.dot(U[:,0:n], np.dot(np.diag(L[0:n]), Vh[0:n,:]))

