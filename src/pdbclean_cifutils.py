import numpy as np
import itertools
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import PDBIO
#
def cif2pdb(ciffile, pdbfile):
    """
    cif2pdb
    """
    parser = MMCIFParser()
    structure = parser.get_structure('traj', ciffile)
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdbfile)

def write_dict_to_cif(cifdict, ciffile):
    """
    write_dict_to_cif:
    write input dictionary into ciffile
    """
    io = MMCIFIO()
    io.set_dict(cifdict)
    io.save(ciffile)

def read_dict_from_file(input_file, keychain=None):
    """
    read_dict_from_file:
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
    return cifdict_reduced

def cifdict_to_cifkeys(cifdict, sort=True, **kwargs):
    """
    cifdict_to_cifkeys
    read the dictionary of a CIF file and return its keys

    INPUT  :
        - a CIF dictionary
    OUTPUT :
        - the sorted list of all atom keys. KEY=SEQ_ID+ASYM_ID+ATOM_ID
    """
    arrkeys = cifdict_to_arrkeys(cifdict, **kwargs)
    cifkeys = np.empty(arrkeys.shape[0], dtype='U25')
    for key in np.arange(cifkeys.shape[0]):
        cifkeys[key] = arrkeys[key,0]+'-'+str(arrkeys[key,1])+'-'+arrkeys[key,2]
    if sort:
        cifkeys = np.sort(cifkeys)
    return cifkeys

def cifdict_to_arrkeys(cifdict, **kwargs):
    """
    """
    asym_id = np.array(cifdict["_atom_site.auth_asym_id"])
    seq_id  = np.array(cifdict["_atom_site.auth_seq_id"])
    atom_id = np.array(cifdict["_atom_site.auth_atom_id"])
    #
    chains = False
    resid  = False
    atoms  = False
    for arg in kwargs:
        if(arg=='chains'):
            chains= True
        if(arg=='resid'):
            resid = True
        if(arg=='atoms'):
            atoms = True
    select = bool(chains+resid+atoms)
    mask = np.full(atom_id.shape[0], False, dtype=bool)
    if chains:
        mask_chain = np.copy(mask)
        for chain in kwargs['chains']:
            mask_chain += np.ma.masked_where(asym_id==chain, asym_id).mask
    if atoms:
        mask_atom = np.copy(mask)
        for atom in kwargs['atoms']:
            mask_atom += np.ma.masked_where(atom_id==atom, atom_id).mask
    if select:
        mask = mask_chain*mask_atom
        mask = np.invert(mask.astype(bool))
    arrkeys = np.column_stack((np.ma.masked_array(asym_id, mask=mask).compressed(),
                               np.ma.masked_array(seq_id,  mask=mask).compressed(),
                               np.ma.masked_array(atom_id, mask=mask).compressed()))
    return arrkeys



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

def pdbs_to_keychain(file_list, verbose=False, **kwargs):
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
        cifkeys = cifdict_to_cifkeys(cifdict, **kwargs)
        if keychain is None:
            keychain = cifkeys
        else:
            keychain = np.unique(np.concatenate((keychain, cifkeys)))
        if verbose:
            print('{0}<'.format(keychain.shape[0]), end='')
    if verbose:
        print('3')
    return keychain

def keychain_to_atommap(keychain, input_list, verbose=False):
    """
    keychain_to_atommap
    """
    if verbose:
        print('Building atom map from keychain.')
    atommap = np.zeros((len(keychain), len(input_list)))
    for i in np.arange(len(input_list)):
        cifkeys = cifdict_to_cifkeys(MMCIF2Dict(input_list[i]))
        atommap[:,i] = np.in1d(keychain, cifkeys).astype(int)
    return atommap

