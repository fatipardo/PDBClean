#
import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
#
def filelist_to_crdarray(filelist):
    """
    filelist_to_crd_array
    """
    crds = None
    for ciffile in filelist:
        cifdict = MMCIF2Dict(ciffile)
        x_list = cifdict["_atom_site.Cartn_x"]
        y_list = cifdict["_atom_site.Cartn_y"]
        z_list = cifdict["_atom_site.Cartn_z"]
        new_crd = np.column_stack((np.array(x_list),np.array(y_list),np.array(z_list)))[np.newaxis,...]
        if crds is None:
            crds = new_crd
        else:
            print(crds.shape, new_crd.shape)
            crds = np.concatenate((crds, new_crd), axis=0)
    return crds
