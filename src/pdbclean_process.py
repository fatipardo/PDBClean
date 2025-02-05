import os, glob
import re
import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from PDBClean import pdbclean_io as pcio
from PDBClean import pdbclean_homogenutils as homogen
from PDBClean import pdbclean_cifutils as cif
#
def process(projdir=None, source='raw_bank', target='clean_bank', 
            pdbformat='.cif', verbose=True, show=False, dry_run=False,
            step='clean', **kwargs):
    """
    process
    """
    source_dir, target_dir = pcio.define_dirs(project_dir=projdir, source=source, target=target)
    input_list = pcio.list_files_in_dir(path=source_dir, ext=pdbformat)
    process_inputlist(input_list, target_dir, pdbformat=pdbformat,
                      verbose=verbose, show=show, dry_run=dry_run,
                      step=step, **kwargs)

def process_inputlist(input_list, target_dir, pdbformat='.cif',
                      verbose=True, show=False, dry_run=False,
                      step='clean', **kwargs):
    """
    process_inputlist
    """
    #print(kwargs)
    # - initialize if needed
    keychain, input_list, assignment = init_process(input_list, target_dir=target_dir, step=step, 
                                                    verbose=verbose, show=show, **kwargs )
    # - iteratively process each file in the list
    i=0
    for input_file in input_list:
        subdir='.'
        if(step=='homogenize'):
            if(kwargs['mode']=='cluster_out'):
                subdir='cluster'+str(assignment[i])
        output_file = pcio.new_filepath(input_file, target_dir, subdir=subdir)
        if verbose:
            print('[{0}/{1}]: {2}'.format(i+1,len(input_list), output_file))
        if(step=='clean'):
            clean(input_file, output_file)
        elif(step=='simplify'):
            simplify(input_file, output_file, pdbformat)
        elif(step=='fixhet'):
            fixatom(input_file, output_file, atom="HETATM")
        elif(step=='fixatm'):
            fixatom(input_file, output_file, atom="ATOM")
        elif(step=='finalize'):
            finalize(input_file, output_file)
        elif(step=='homogenize'):
            reduce_to_keychain(input_file, output_file, keychain)
        elif(step=='select'):
            reduce_to_keychain(input_file, output_file, keychain)
        i+=1

####

def init_process(input_list, target_dir=None, step='clean', verbose=True, show=False, **kwargs):
    """
    init_process
    """
    keychain=[]
    assignment=[]
    if(step=='select'):
        keychain = cif.pdbs_to_keychain(input_list, verbose=verbose, **kwargs)
    elif(step=='homogenize'):
        if(kwargs['mode']=='cluster_out'):
            keychain = np.load(target_dir+'/keychain.npy')
        else:
            keychain = cif.pdbs_to_keychain(input_list, verbose=verbose, **kwargs)
        if(kwargs['mode']=='keep_all_samples'):
            keychain = homogen.reduce_feature_keep_samples(keychain, input_list,
                                                           verbose=verbose, show=show)
        elif(kwargs['mode']=='optimize'):
            keychain, input_list = homogen.reduce_optimized(keychain, input_list,
                                                            verbose=verbose, show=show)
        elif(kwargs['mode']=='see_cluster_hierarchy'):
            homogen.cluster(keychain, input_list, path=target_dir,
                            verbose=verbose, show=show)
            np.save(target_dir+'/keychain.npy', keychain)
            input_list = ()
        elif(kwargs['mode']=='cluster_out'):
            n_clusters=-1
            if 'n_clusters' in kwargs:
                n_clusters=kwargs['n_clusters']
            cutoff=-1
            if 'cutoff' in kwargs:
                cutoff=kwargs['cutoff']
            assignment = homogen.assign_clusters(keychain, input_list, 
                                                 n_clusters=n_clusters, cutoff=cutoff,
                                                 path=target_dir, 
                                                 verbose=verbose, show=show)
    return keychain, input_list, assignment

def reduce_to_keychain(input_file, output_file, keychain):
    """
    reduce_to_keychain
    """
    cifdict = cif.read_dict_from_file(input_file, keychain=keychain)
    cif.write_dict_to_cif(cifdict, output_file)

#
def finalize(oldfile, newfile):
    """
    finalize
    """
    with open(oldfile) as myfile:
        newciffile = open(newfile, 'w')
        for line in myfile:
            line_split = line.strip()
            line_split = line.split()
            if (line_split[0] == "ATOM") or (line_split[0] == "HETATM"):
                newciffile.write(line_split[0] + " " + line_split[1] + " " + line_split[2] + " " + line_split[3] + " " + line_split[4] + " " + line_split[5] + " " + line_split[17] + " " + line_split[7] + " " + line_split[8] + " " + line_split[9] + " " + line_split[10] + " " + line_split[11] + " " + line_split[12] + " " + line_split[13] + " " + line_split[14] + " " + line_split[15] + " " + line_split[16] + " " + line_split[17] + " " + line_split[18] + " " + line_split[19] + "\n")
            else:
                newciffile.write(line)
#
def fixatom(oldfile, newfile, atom="HETATM"):
    """
    fixatom
    """
    with open(oldfile) as myfile:
        hetch_count_map = {}
        newciffile = open(newfile, 'w')
        last_res = 0
        last_res = {}
        usage = {}
        for line in myfile:
            line_split = line.strip()
            line_split = line_split.split()
            # Can Update so that it looks at all atoms?
            if (line_split[0] == atom):
                key = str(line_split[15]) + "_" + str(line_split[16]) + "_" + str(line_split[17]) + "_" + str(line_split[18])
                if line_split[17] in hetch_count_map:
                    if (line_split[15] != last_res[line_split[17]]) or (key in usage):
                        if key not in usage:
                            usage[key] = 1
                        hetch_count_map[line_split[17]] += 1
                        newciffile.write(line_split[0] + " " + line_split[1] + " " + line_split[2] + " " + line_split[3] + " " + line_split[4] + " " + line_split[5] + " " + line_split[6] + " " + line_split[7] + " " + str(hetch_count_map[line_split[17]]) + " " + line_split[9] + " " + line_split[10] + " " + line_split[11] + " " + line_split[12] + " " + line_split[13] + " " + line_split[14] + " " + str(hetch_count_map[line_split[17]]) + " " + line_split[16] + " " + line_split[17] + " " + line_split[18] + " " + line_split[19] + "\n")
                        last_res[line_split[17]] = line_split[15]
                    else:
                        usage[key] = 1
                        newciffile.write(line_split[0] + " " + line_split[1] + " " + line_split[2] + " " + line_split[3] + " " + line_split[4] + " " + line_split[5] + " " + line_split[6] + " " + line_split[7] + " " + str(hetch_count_map[line_split[17]]) + " " + line_split[9] + " " + line_split[10] + " " + line_split[11] + " " + line_split[12] + " " + line_split[13] + " " + line_split[14] + " " + str(hetch_count_map[line_split[17]]) + " " + line_split[16] + " " + line_split[17] + " " + line_split[18] + " " + line_split[19] + "\n")
                else:
                    usage[key] = 1
                    hetch_count_map[line_split[17]] = 1
                    last_res[line_split[17]] = line_split[15]
                    newciffile.write(line_split[0] + " " + line_split[1] + " " + line_split[2] + " " + line_split[3] + " " + line_split[4] + " " + line_split[5] + " " + line_split[6] + " " + line_split[7] + " " + str(hetch_count_map[line_split[17]]) + " " + line_split[9] + " " + line_split[10] + " " + line_split[11] + " " + line_split[12] + " " + line_split[13] + " " + line_split[14] + " " + str(hetch_count_map[line_split[17]]) + " " + line_split[16] + " " + line_split[17] + " " + line_split[18] + " " + line_split[19] + "\n")
            else:
                newciffile.write(line)
#
def simplify(oldfile, newfile, pdbformat):
    """
    simplify
    """
    mmcif_dict = MMCIF2Dict(oldfile)
    #
    # Create map from asym_id to assembly_id
    # assembly_id information may be either str or list type, so I'm going to force
    # to list type
    asym_assembly_map = {}
    if '_pdbx_struct_assembly_gen.assembly_id' in mmcif_dict:
        assembly_id = mmcif_dict['_pdbx_struct_assembly_gen.assembly_id']
    else:
        assembly_id = '1'
    if '_pdbx_struct_assembly_gen.asym_id_list' in mmcif_dict:
        asym_id = mmcif_dict['_pdbx_struct_assembly_gen.asym_id_list']
    else:
        asym_id = 'A'
    #
    if not isinstance(assembly_id, list):
        assembly_id_list = []
        asym_id_list = []
        assembly_id_list.append(assembly_id)
        asym_id_list.append(asym_id)
    else:
        assembly_id_list = []
        assembly_id_list = assembly_id
        asym_id_list = asym_id
    # Only keep one asymmetric unit
    id_remove = []
    if(len(asym_id_list)>1):
        for i in np.arange(len(asym_id_list)-1):
            for j in np.arange(i+1,len(asym_id_list)):
                if(asym_id_list[j] == asym_id_list[i]):
                    id_remove.append(j)
    if(len(id_remove)>0):
        for i in id_remove:
            del asym_id_list[i]
            del assembly_id_list[i]
    # Convert asym_id entry into a list of asym_ids
    for i in range(len(assembly_id_list)):
        asym_id = asym_id_list[i]
        asym_id = asym_id.strip()
        asym_id = re.sub(' ', '', asym_id)
        asym_id = asym_id.split(',')
        for ident in asym_id:
            asym_assembly_map[ident] = assembly_id_list[i]
    #
    # Create entity_id -> [asym_id] correspondence map.
    # This is needed in order to determine whether or not a chain's molID info
    # should be printed to the file
    if '_atom_site.label_entity_id' in mmcif_dict:
        entity_id_list = mmcif_dict['_atom_site.label_entity_id']
    else:
        entity_id_list = ['1']
    if '_atom_site.label_asym_id' in mmcif_dict:
        asym_id_list = mmcif_dict['_atom_site.label_asym_id']
    else:
        asym_id_list = ['1']

    entity_asym_assembly_map = {}
    for i in range(len(entity_id_list)):
        if entity_id_list[i] in entity_asym_assembly_map:
            if asym_id_list[i] not in entity_asym_assembly_map[entity_id_list[i]]:
                entity_asym_assembly_map[entity_id_list[i]].append(asym_id_list[i])
        else:
            entity_asym_assembly_map[entity_id_list[i]] = [asym_id_list[i]]
    #
    # Start writing to file for each biological assembly unit (different structure)
    for assembly in assembly_id_list:
        if (len(assembly_id_list)==1):
            newciffilename = str(re.sub(pdbformat, '', newfile))+"+00"
        else:
            newciffilename = str(re.sub(pdbformat, '', newfile))+"+0"+str(assembly)
        newciffile = open(newciffilename+pdbformat, 'w')
        newciffile.write("data_"+newciffilename+"\n")
        #Write entry.id
        newciffile.write("#\n")
        L = mmcif_dict['_entry.id']
        entryid = '_entry.id   ' + L
        newciffile.write(entryid + "\n")
        # Write Audit category
        newciffile.write("#\n")
        newciffile.write("loop_\n")
        newciffile.write("_citation_author.name\n")
        L = mmcif_dict['_citation_author.name']
        if isinstance(L, list):
            for i in L:
                newciffile.write("'" + re.sub("'", "", i) + "'" + "\n")
        else:
            newciffile.write("'" + re.sub("'", "", L) + "'" + "\n")
        # Write Citation category
        newciffile.write("#" + "\n")
        newciffile.write("loop_" + "\n")
        newciffile.write("_citation.title" + "\n")
        newciffile.write("_citation.year" + "\n")
        newciffile.write("_citation.pdbx_database_id_DOI" + "\n")
        L1 = mmcif_dict['_citation.title']
        L2 = mmcif_dict['_citation.year']
        L3 = mmcif_dict['_citation.pdbx_database_id_DOI']
        if isinstance(L1, list):
            for i in range(len(L1)):
                #newciffile.write("'" + re.sub("'", "", L1[i]) + "' " + L2[i] + " " + L3[i] + "\n")
                newciffile.write("'" + re.sub("'|\n|\r", "", L1[i]) + "' " + L2[i] + " " + L3[i] + "\n")
        else:
            #newciffile.write("'" + re.sub("'", "", L1) + "' " + L2 + " " + L3 + "\n")
            newciffile.write("'" + re.sub("'|\n|\r", "", L1) + "' " + L2 + " " + L3 + "\n")
        # Write Resolution category
        newciffile.write("#" + "\n")
        newciffile.write("loop_" + "\n")
        newciffile.write("_exptl.method" + "\n")
        newciffile.write("_exptl.resolution" + "\n")
        if '_exptl.method' in mmcif_dict:
            L1 = mmcif_dict['_exptl.method']
        elif '_refine_hist.pdbx_refine_id' in mmcif_dict:
            L1 = mmcif_dict['_refine_hist.pdbx_refine_id']
        elif '_refine.pdbx_refine_id' in mmcif_dict:
            L1 = mmcif_dict['_refine.pdbx_refine_id']
        else:
            L1 = 'X'
        if '_refine.ls_d_res_high' in mmcif_dict:
            L2 = mmcif_dict['_refine.ls_d_res_high']
        elif '_em_3d_reconstruction.resolution' in mmcif_dict:
            L2 = mmcif_dict['_em_3d_reconstruction.resolution']
        elif '_refine_hist.d_res_high' in mmcif_dict:
            L2 = mmcif_dict['_refine_hist.d_res_high']
        else:
            L2 = 'X'
        if isinstance(L1, list) and isinstance(L2, list):
            for i in range(len(L1)):
                newciffile.write("'" + L1[i] + "' " + L2[i] + " " + "\n")
        elif isinstance(L1, list) and not isinstance(L2,list):
            newciffile.write("'" + L1[0] + "' " + L2 + " " + "\n")
        elif not isinstance(L1,list) and isinstance(L2,list):
            newciffile.write("'" + L1 + "' " + L2[0] + " " + "\n")
        else:
            newciffile.write("'" + L1 + "' " + L2 + " " + "\n")
        # Write Entity category
        # Note: Special check included to make sure molID pertains to chain in assembly
        newciffile.write("#" + "\n")
        newciffile.write("loop_" + "\n")
        newciffile.write("_entity.id" + "\n")
        newciffile.write("_entity.pdbx_description" + "\n")
        L1 = mmcif_dict['_entity.id']
        L2 = mmcif_dict['_entity.pdbx_description']
        if isinstance(L1, list):  # Have to check list or not because files with only one assembly unit won't be list
            for i in range(len(L1)):
                if L1[i] in entity_asym_assembly_map:
                    for asym in entity_asym_assembly_map[L1[i]]:
                        if asym in asym_assembly_map:
                            if (assembly == asym_assembly_map[asym]): # Only print molIDs pertaining to particular assembly unit
                                L2[i] = L2[i].upper()
                                L2[i] = L2[i].replace(":", "")
                                newciffile.write(L1[i] + " '" + L2[i].replace("'", "") + "'\n")
                                break  # Because multiple asym per assembly and don't want to print multiple molID lines
        else:
            for asym in entity_asym_assembly_map[L1]:
                if asym in asym_assembly_map:
                    if (asym_assembly_map[asym] == assembly):
                        L2 = L2.upper()
                        L2 = L2.replace(":", "")
                        newciffile.write(L1 + " '" + L2.replace("'", "") + "'\n")
                        break  # Because multiple asym per assembly and don't want to print multiple molID lines
        # Now write the coordinate portion of the file
        newciffile.write("#" + "\n")
        newciffile.write("loop_" + "\n")
        newciffile.write("_atom_site.group_PDB" + "\n")
        newciffile.write("_atom_site.id" + "\n")
        newciffile.write("_atom_site.type_symbol" + "\n")
        newciffile.write("_atom_site.label_atom_id" + "\n")
        newciffile.write("_atom_site.label_alt_id" + "\n")
        newciffile.write("_atom_site.label_comp_id" + "\n")
        newciffile.write("_atom_site.label_asym_id" + "\n")
        newciffile.write("_atom_site.label_entity_id" + "\n")
        newciffile.write("_atom_site.label_seq_id" + "\n")
        newciffile.write("_atom_site.pdbx_PDB_ins_code" + "\n")
        newciffile.write("_atom_site.Cartn_x" + "\n")
        newciffile.write("_atom_site.Cartn_y" + "\n")
        newciffile.write("_atom_site.Cartn_z" + "\n")
        newciffile.write("_atom_site.occupancy" + "\n")
        newciffile.write("_atom_site.B_iso_or_equiv" + "\n")
        newciffile.write("_atom_site.auth_seq_id" + "\n")
        newciffile.write("_atom_site.auth_comp_id" + "\n")
        newciffile.write("_atom_site.auth_asym_id" + "\n")
        newciffile.write("_atom_site.auth_atom_id" + "\n")
        newciffile.write("_atom_site.pdbx_PDB_model_num" + "\n")
        L1 = mmcif_dict['_atom_site.group_PDB']
        L2 = mmcif_dict['_atom_site.id']
        L3 = mmcif_dict['_atom_site.type_symbol']
        L4 = mmcif_dict['_atom_site.label_atom_id']
        L5 = mmcif_dict['_atom_site.label_alt_id']
        L6 = mmcif_dict['_atom_site.label_comp_id']
        L7 = mmcif_dict['_atom_site.label_asym_id']
        L8 = mmcif_dict['_atom_site.label_entity_id']
        L9 = mmcif_dict['_atom_site.label_seq_id']
        L10 = mmcif_dict['_atom_site.pdbx_PDB_ins_code']
        L11 = mmcif_dict['_atom_site.Cartn_x']
        L12 = mmcif_dict['_atom_site.Cartn_y']
        L13 = mmcif_dict['_atom_site.Cartn_z']
        L14 = mmcif_dict['_atom_site.occupancy']
        L15 = mmcif_dict['_atom_site.B_iso_or_equiv']
        L16 = mmcif_dict['_atom_site.auth_seq_id']
        L17 = mmcif_dict['_atom_site.auth_comp_id']
        L18 = mmcif_dict['_atom_site.auth_asym_id']
        L19 = mmcif_dict['_atom_site.auth_atom_id']
        L20 = mmcif_dict['_atom_site.pdbx_PDB_model_num']
        # This is strictly to clean up a common problem found in ribosome structures
        # HETATMs are defined mid chain, ie in an asymmetric unit of non-HETATMs
        # Want to treat these atoms different from other HETATMs
        asymatom_list = []
        for i in range(len(L1)):
            if L1[i] == "ATOM":
                if L7[i] not in asymatom_list:
                    asymatom_list.append(L7[i])
        #
        hetatm_map = {}
        hetatm_count = 0
        for i in range(len(L1)):
            if L7[i] in asym_assembly_map:
                if (assembly == asym_assembly_map[L7[i]]): # Only print molIDs pertaining to particular assembly unit
                    if (L1[i]=="ATOM"):
                        newciffile.write(L1[i] + " " + L2[i] + " " + L3[i] + ' "' + L4[i] + '" ' + L5[i] + " " + L6[i] + " " + L7[i] + " " + L8[i] + " " + L9[i] + " " + L10[i] + " " + L11[i] + " " + L12[i] + " " + L13[i] + " " + L14[i] + " " + L15[i] + " " + L16[i] + " " + L17[i] + " " + L18[i] + ' "' + L19[i] + '" ' + L20[i] + "\n")
                    elif (L1[i]=="HETATM"):
                        if L7[i] in asymatom_list:
                            newciffile.write("ATOM" + " " + L2[i] + " " + L3[i] + ' "' + L4[i] + '" ' + L5[i] + " " + L6[i] + " " + L7[i] + " " + L8[i] + " " + L9[i] + " " + L10[i] + " " + L11[i] + " " + L12[i] + " " + L13[i] + " " + L14[i] + " " + L15[i] + " " + L16[i] + " " + L17[i] + " " + L18[i] + ' "' + L19[i] + '" ' + L20[i] + "\n")
                        else:
                            newciffile.write(L1[i] + " " + L2[i] + " " + L3[i] + ' "' + L4[i] + '" ' + L5[i] + " " + L6[i] + " " + L7[i] + " " + L8[i] + " " + L9[i] + " " + L10[i] + " " + L11[i] + " " + L12[i] + " " + L13[i] + " " + L14[i] + " " + L15[i] + " " + L16[i] + " " + L17[i] + " " + "het" + L6[i] + L18[i] + ' "' + L19[i] + '" ' + L20[i] + "\n")
        newciffile.write("#" + "\n")
#
def clean(oldfile, newfile):
    """
    clean
    """
    entry_list = ['_entry.id',
                  '_atom_site.group_PDB',
                  '_citation_author.name',
                  '_citation.title',
                  '_pdbx_struct_assembly_gen.assembly_id',
                  '_entity.pdbx_description',
                  '_exptl.method',
                  '_em_3d_reconstruction.resolution',
                  '_refine_hist.pdbx_refine_id',
                  '_refine.pdbx_refine_id']
    with open(oldfile) as old_file:
        alllines = []
        linecount = 0
        poundline = 0
        flag = 0
        for line in old_file:
            alllines.append(line)
            if linecount == 0:
                with open(newfile, 'a') as new_file:
                    new_file.write(alllines[0])
            for entry in entry_list:
                flag = check_and_write_entry(entry, line, alllines, line[0:len(entry)], flag, range(poundline, linecount), newfile)
            if '#' in line[0]:
                poundline = linecount
            linecount += 1
        with open(newfile, 'a') as new_file:
            new_file.write('#\n')
#
def check_and_write_entry(entry, line, alllines, key, flag, linerange, newfile):
    """
    check_and_write_entry
    """
    if entry in key:
        flag = 1
    elif (flag==1) and '#' in line[0]:
        with open(newfile, 'a') as new_file:
            for i in linerange:
                new_file.write(alllines[i])
        flag=0
    return flag

