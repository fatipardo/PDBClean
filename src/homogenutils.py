#
import os
import glob
import numpy as np
import numpy.ma as ma
import itertools
from matplotlib import pyplot as plt
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.mmcifio import MMCIFIO
#
def process(projdir=None, source='final_bank', target='homogmax_bank',
            pdbformat='.cif', verbose=False, show=False, dry_run=False,
            mode='void'):
    """
    process
    """
    source_dir, target_dir = define_dirs(project_dir=projdir, source=source, target=target)
    input_list = list_files_in_dir(path=source_dir, ext=pdbformat)
    keychain   = pdbs_to_keychain(input_list, verbose=verbose)
    if(mode=='keep_all_samples'):
        keychain_reduced = reduce_feature_keep_samples(keychain, input_list, verbose=verbose, show=show)
        reduce_and_save(input_list, target_dir, keychain_reduced, verbose=verbose)

def reduce_and_save(input_files, output_dir, keychain, verbose=False):
    """
    """
    i=0
    for input_file in input_files:
        cifdict      = get_dict_of_cif_file(input_file, keychain=keychain)
        output_file  = new_filepath(input_file, output_dir)
        if verbose:
            i+=1
            print('[{0}/{1}]: {2}'.format(i,len(input_files),output_file))
        write_dict_to_cif(cifdict, output_file)

###############
# HOMOGENEIZE #
###############

def write_dict_to_cif(cifdict, ciffile):
    """
    write_dict_to_cif:
    write input dictionary into ciffile
    """
    io = MMCIFIO()
    io.set_dict(cifdict)
    io.save(ciffile)

def get_dict_of_cif_file(input_file, keychain=None):
    """
    get_dict_of_cif_file:
    reads input_file and outputs its dictionary cifdict.
    If optional argument keychain is given, only the atom
    information of keys found in keychain will be returned.
    """
    cifdict = MMCIF2Dict(input_file)
    if keychain is not None:
        mask = cifdict_to_mask(cifdict, keychain)
        cifdict = reduce_cifdict(cifdict, mask)
    return cifdict

def cifdict_to_mask(cifdict, keychain):
    """
    """
    cifkeys = cifdict_to_cifkeys(cifdict, sort=False)
    mask    = np.in1d(cifkeys, keychain)
    return mask

def reduce_cifdict(cifdict, mask):
    """
    """
    cifdict_reduced = {}
    for key in cifdict.keys():
        if(key.split('.')[0]=='_atom_site'):
            cifdict_reduced[key] = list(itertools.compress(cifdict[key], mask))
        else:
            cifdict_reduced[key] = cifdict[key]
        #print(key, cifdict_reduced[key])
    return cifdict_reduced

def reduce_feature_keep_samples(feature, sample, verbose=False, show=False):
    """
    """
    map = keychain_to_atommap(feature, sample, verbose=verbose, show=show)
    feature_score   = np.sum(map, axis=1)
    feature_mask    = np.where(feature_score < map.shape[1], True, False)
    feature_reduced = ma.masked_array(feature, mask=feature_mask).compressed()
    map = keychain_to_atommap(feature_reduced, sample, verbose=verbose, show=show)
    return feature_reduced

############
# KEYCHAIN #
############

def pdbs_to_keychain(file_list, verbose=False):
    """
    pdbs_to_keychain
    each atom is identified by a unique key. We list them in a keychain
    
    INPUT  :
        - file_list: the list of input filename
    OUTPUT :
        - keychain: sorted list of unique atoms found in filelist, identified with an atom key.
    """
    if verbose:
        print('Building keychain from all files. Length of keychain: ')
    keychain = None
    for ciffile in file_list:
        cifdict = MMCIF2Dict(ciffile)
        cifkeys = cifdict_to_cifkeys(cifdict)
        if keychain is None:
            keychain = cifkeys
        else:
            keychain = np.unique(np.concatenate((keychain, cifkeys)))
        if verbose:
            print('{0}<'.format(keychain.shape[0]), end='')
    if verbose:
        print('3')
    return keychain

def cifdict_to_cifkeys(cifdict, sort=True):
    """
    cifdict_to_cifkeys
    read the dictionary of a CIF file and return its keys

    INPUT  :
        - a CIF dictionary
    OUTPUT :
        - the sorted list of all atom keys. KEY=SEQ_ID+ASYM_ID+ATOM_ID
    """
    seq_id  = cifdict["_atom_site.auth_seq_id"]
    asym_id = cifdict["_atom_site.auth_asym_id"]
    atom_id = cifdict["_atom_site.auth_atom_id"]
    arrkeys = np.column_stack((np.array(asym_id),
                               np.array(seq_id),
                               np.array(atom_id)))
    cifkeys = np.empty(arrkeys.shape[0], dtype='U25')
    for key in np.arange(cifkeys.shape[0]):
        cifkeys[key] = arrkeys[key,0]+'-'+str(arrkeys[key,1])+'-'+arrkeys[key,2]
    if sort:
        cifkeys = np.sort(cifkeys)
    return cifkeys

def cifkeys_to_cifdict(cifkeys, cifdict=None):
    """
    cifkeys_to_cifdict
    """
    seq_ids  = []
    asym_ids = []
    atom_ids = []
    for key in np.arange(cifkeys.shape[0]):
        seq_id, asym_id, atom_id = cifkeys[key].split("-")
        seq_ids.append(seq_id)
        asym_ids.append(asym_id)
        atom_ids.append(atom_id)
    if cifdict is None:
        cifdict = {}
    cifdict["_atom_site.auth_seq_id"] = seq_ids
    cifdict["_atom_site.auth_asym_id"] = asym_ids
    cifdict["_atom_site.auth_atom_id"] = atom_ids
    return cidfict

def keychain_to_atommap(keychain, input_list, verbose=False, show=False):
    """
    keychain_to_atommap
    """
    if verbose:
        print('Building atom map from keychain.')
    atommap = np.zeros((len(keychain), len(input_list)))
    for i in np.arange(len(input_list)):
        cifkeys = cifdict_to_cifkeys(MMCIF2Dict(input_list[i]))
        atommap[:,i] = np.in1d(keychain, cifkeys).astype(int)
    if show:
        show_atommap(atommap)
    return atommap

######
# IO #
######

def new_filepath(input_file, output_dir):
    """
    """
    filename = os.path.basename(input_file)
    return output_dir+'/'+filename

def list_files_in_dir(path, ext):
    """
    list_files_in_dir
    """
    filelist = ()
    if os.path.exists(path):
        filelist = glob.glob(path+'/*'+ext)
    return filelist

def define_dirs(project_dir=None, source=None, target=None):
    """
    define_dirs
    """
    if project_dir is None:
        source_dir = None
        target_dir = None
    else:
        if source is None:
            source_dir = None
        else:
            source_dir = project_dir+'/'+source
        if target is None:
            target_dir = None
        else:
            target_dir = project_dir+'/'+target
    return source_dir, target_dir

####### 
# VIZ #
#######

def show_atommap(atommap):
    """
    show_atommap
    """
    coverage_per_sample = np.sum(atommap,axis=0)
    coverage_per_feature = np.sum(atommap,axis=1)
    fig = plt.figure(figsize=(12,4))
    ax = fig.add_subplot(1,3,1)
    ax.set_title('Atom map')
    ax.set_xlabel('list of CIF files')
    ax.set_ylabel('list of atoms')
    plt.imshow(atommap, aspect='auto', cmap='Greys')
    ax = fig.add_subplot(1,3,2)
    ax.set_title('coverage per atom')
    ax.set_xlabel('list of atoms')
    plt.plot(coverage_per_feature, '.')
    ax = fig.add_subplot(1,3,3)
    ax.set_title('coverage per structure')
    ax.set_xlabel('list of CIF files')
    plt.plot(coverage_per_sample, '.')
    plt.tight_layout()


