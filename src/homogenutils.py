#
import glob
import numpy as np
import numpy.ma as ma
from matplotlib import pyplot as plt
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
#
def process(projdir=None, mode='show_atommap', source='final_bank', target='homogmax_bank', pdbformat='.cif', verbose=True, dry_run=False):
    """
    process
    """
    if projdir is not None:
        source_dir = projdir+'/'+source
        target_dir = projdir+'/'+target
        input_list = glob.glob(source_dir+'/*'+pdbformat)
        atommap, keychain = filelist_to_atommap(input_list)
        if(mode=='show_atommap'):
            show_atommap(atommap)
        else
            if(mode=='max_sample'):
                atommap, keychain, input_list = maximize_samplesize(atommap, keychain, input_list)
            if not dry_run:
                i=0
                for input_cif in input_list:
                    if verbose:
                        i+=1
                        print('[{0}/{1}]: {2}'.format(i,len(input_list),cif_name))
                    output_cif=target_dir+'/'+cif_name
                    keychain_to_cif(input_cif, output_cif, keychain)
    return atommap, keychain, input_list
#
def keychain_to_cif(oldfile, newfile, keychain):
    """
    keychain_to_cif
    """
    #
    # Split keychain here
    # 
    with open(oldfile) as myfile:
        newciffile = open(newfile, 'w')
        for line in myfile:
            line_split = line.strip()
            line_split = line.split()
            if (line_split[0] == "ATOM") or (line_split[0] == "HETATM"):
                #
                # check if line_slit[...] is in split keychain. Do not write if not.
                # A BETTER APPROACH WOULD BE TO MASK LINES BASED ON KEYCHAIN AND ONLY WRITE THE TRUE ONES...
                #
                newciffile.write(line_split[0] + " " + line_split[1] + " " + line_split[2] + " " + line_split[3] + " " + line_split[4] + " " + line_split[5] + " " + line_split[17] + " " + line_split[7] + " " + line_split[8] + " " + line_split[9] + " " + line_split[10] + " " + line_split[11] + " " + line_split[12] + " " + line_split[13] + " " + line_split[14] + " " + line_split[15] + " " + line_split[16] + " " + line_split[17] + " " + line_split[18] + " " + line_split[19] + "\n")
            else:
                newciffile.write(line)

def maximize_samplesize(atommap, keychain, input_list):
    """
    maximize_samplesize
    remove every atom that does not belong to all samples
    """
    coverage_per_feature = np.sum(atommap,axis=1)
    mask     = np.where(coverage_per_feature < atommap.shape[1], True, False)
    mask_map = np.repeat(mask[:, np.newaxis], atommap.shape[1], axis=1)
    new_keychain = ma.masked_array(keychain, mask=mask).compressed()
    new_atommap  = ma.masked_array(atommap,  mask=mask_map).compressed().reshape(new_keychain.shape[0],atommap.shape[1])
    return new_atommap, new_keychain, input_list

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

def filelist_to_atommap(filelist, verbose=False):
    """
    filelist_to_atommap
    """
    keychain = filelist_to_keychain(filelist)
    atommap = np.zeros((len(keychain), len(filelist)))
    for i in np.arange(len(filelist)):
        if verbose:
            print(filelist[i])
        cifkeys = cifdict_to_cifkeys(MMCIF2Dict(filelist[i]))
        atommap[:,i] = np.in1d(keychain, cifkeys).astype(int)
    return atommap, keychain

def filelist_to_keychain(filelist, verbose=False):
    """
    filelist_to_keychain:
    each atom is identified by a unique key. We list them in a keychain
    """
    keychain = None
    for ciffile in filelist:
        cifdict = MMCIF2Dict(ciffile)
        cifkeys = cifdict_to_cifkeys(cifdict)
        if keychain is None:
            keychain = cifkeys
        else:
            keychain = np.unique(np.concatenate((keychain, cifkeys)))
        if verbose:
            print('Length of keychain: {0}'.format(keychain.shape[0]))
    return keychain

def cifdict_to_cifkeys(cifdict):
    """
    cifdict_to_cifkeys
    read the dictionary of a CIF file and return its keys
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
    return np.sort(cifkeys)
