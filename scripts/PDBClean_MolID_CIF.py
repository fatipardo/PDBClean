#!/usr/bin/env python
# coding: utf-8
#
# Comments from Nicholas Corsepius:
## ! ! ! master_molID_class_list is very important
# This is the list that contains every file's MolID class
# !! molIDConversion_list is important . . . it contains the objects
# MolIDConversion and is what is going to be updated by the user and evaulated
# to determine when the next step in the program is unlocked
# Create list of MolIDConversion objects using unique_molID_occur_map

from __future__ import print_function
import os, sys
import re
import numpy as np
import csv
import copy
import glob
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from PDBClean import pdbcleanmolidcifutils as molidutils


########################
# READ INPUT ARGUMENTS #
########################
n_arg = len(sys.argv)
if(n_arg<3):
    print('Usage error: {0} <source directory> <target directory>'.format(sys.argv[0]))
    sys.exit()
source_dir=sys.argv[1]
target_dir=sys.argv[2]


#########################################
# READ PDB FILES AND DEFINE MolID LISTS #
#########################################

filelist=glob.glob(source_dir+'/*.cif')
master_molID_class_list = molidutils.pdb_to_masterlist(filelist)
unique_molID_occur_map  = molidutils.CreateMasterUniqueMolIDMap(master_molID_class_list)
molIDConversion_list    = molidutils.uniquelist_to_conversionlist(unique_molID_occur_map)


#####################################
# INTERACTIVE MOLID CONVERSION MENU #
#####################################
# Goal: 
# users complete their molID conversion templates by ensuring that each member of
# molIDConversion_list has status complete = True
input_menu = ""
input_menu_complete = ""
# For use in the next section
concat_menu = ""

 # ADD Help Option and Quit options
while(input_menu != "QUIT"):
    if (input_menu_complete == "1"):
        print("""Congratulations! You have successfully constructed your
    conversion templates.You can proceed to the next section by selection option
    7 or, continue to edit your conversion template through this menu""")
    print("""PDBClean MolID Conversion Build Menu
             Select one of the following options to proceed:
             1) Show full conversion
             2) Show only unassigned conversions
             3) Enter input file
             4) Search MolID to add chain ID conversion
             5) Go entry by entry to add chain ID conversion
             6) Remove a chain ID conversion
          """)
    if (input_menu_complete == "1"):
        print("    7) Continue to next step of curation")

    input_menu = input('Option Number: ')
    if (input_menu == "1"):
        molidutils.show_full_conversion(molIDConversion_list)
    elif (input_menu == "2"):
        molidutils.show_unassigned_conversion(molIDConversion_list)
    elif (input_menu == "3"):
        input_cnv_file = input('Conversion File: ')
        user_molID_chID_map =  molidutils.read_input_file(input_cnv_file)
        molIDConversion_list = molidutils.add_user_conversion(user_molID_chID_map, molIDConversion_list)
    elif (input_menu == "4"):
        search_term = input('MolID search term: ')
        molIDConversion_list, search_molIDConversion_list = molidutils.search_conversion(molIDConversion_list, search_term)
        input_submenu_4_1 = 0
        while (input_submenu_4_1 != "DONE"):
            print("    1) Further narrow down search results")
            print("    2) Add chain ID to conversion templates")
            input_submenu_4_1 = input('Option Number: ')
            if (input_submenu_4_1 == "QUIT"):
                input_submenu_4_1 = "DONE"
            if (input_submenu_4_1 == "1"):
                search_term = input('MolID search term: ')
                molIDConversion_list, search_molIDConversion_list = molidutils.search_again_conversion(molIDConversion_list, search_molIDConversion_list, search_term)
            elif (input_submenu_4_1 == "2"):
                print("""Enter new chain IDs, comma separated, no spaces""")
                chID_list = input('Chain IDs: ')
                if (chID_list == "") or (chID_list == "QUIT"):
                    pass
                else:
                    molIDConversion_list, search_molIDConversion_list = molidutils.edit_chain_conversion(molIDConversion_list, 
                                                                                                         search_molIDConversion_list, 
                                                                                                         chID_list,
                                                                                                         action='add')
                input_submenu_4_1 = "DONE"
    elif (input_menu == "5"):
        print("Enter chain IDs for each of the following MolID.")
        print("Comma separated, no spaces")
        for molIDConversion in molIDConversion_list:
            chID_list = input(molIDConversion.molID+":")
            if (chID_list == "") or (chID_list == "QUIT"):
                pass
            else:
                chID_list = chID_list.split(',')
            molIDConversion.add_chID_list(chID_list)
            molIDConversion.check_for_completeness()
    elif (input_menu == "6"):
        search_term = input('MolID search term: ')
        molIDConversion_list, search_molIDConversion_list = molidutils.search_conversion(molIDConversion_list, search_term)
        input_submenu_6_1 = 0
        while (input_submenu_6_1 != "DONE"):
            print("    1) Further narrow down search results")
            print("    2) Remove chain ID from conversion templates")
            input_submenu_6_1 = input('Option Number: ')
            if (input_submenu_6_1 == "1"):
                search_term = input('MolID search term: ')
                molIDConversion_list, search_molIDConversion_list = molidutils.search_again_conversion(molIDConversion_list, search_molIDConversion_list, search_term)
            elif (input_submenu_6_1 == "2"):
                print("""Enter chain ID to be removed, comma separated, no spaces""")
                chID_list = input('Chain IDs: ')
                if (chID_list == "") or (chID_list == "QUIT"):
                    pass
                else:
                    molIDConversion_list, search_molIDConversion_list = molidutils.edit_chain_conversion(molIDConversion_list, 
                                                                                                         search_molIDConversion_list, 
                                                                                                         chID_list,
                                                                                                         action='remove')
                input_submenu_6_1 = "DONE"
    elif (input_menu == "7"):
        if (input_menu_complete == "1"):
            input_menu = "QUIT"
            concat_menu = "START"
    #
    input_menu_complete = molidutils.check_complete(molIDConversion_list)

########################################
# INTERACTIVE MOLID CONCATENATION MENU #
########################################
# Goal: 

if (concat_menu == "START"):
    # Prepare for concatenation step
    # We now have to take the information contained in the MolIDConversion objects
    # in molIDConversion_list to update the MolID objects in master_molID_class_list
    # We then need to mine these updated MolID objects to figure out which ones
    # contain concatenated chains. These will be presented to the user in another
    # interactive menu section where they can update the planned conversion on
    # a file by file basis
    
    master_molID_class_list = molidutils.update_masterlist(master_molID_class_list, molIDConversion_list)

    concat_menu = ""
    concat_menu_complete = ""

    while(concat_menu != "QUIT"):

        count_problems = molidutils.problem_counter(master_molID_class_list)
        if (count_problems == 0):
            concat_menu_complete = "1"

        if (concat_menu_complete == "1"):
            print("""Congratulations! You have successfully constructed your
                     conversion templates.You can proceed to the next section by selection option
                     7 or, continue to edit your conversion template through this menu""")
        print("""PDBClean Concatenations Menu
                 Note: All proposed concatenations must be accepted before the curation can
                 be completed.
                 Select one of the following options to proceed:
                 1) Show all conversions
                 2) Show only unaccepted concatenations
                 3) Search and modify destination chainIDs of proposed concatenations
                 4) Search and modify order of proposed concatenations
                 5) Search and accept proposed concatenations
              """)
        if (concat_menu_complete == "1"):
            print("    6) Finalize Curation")

        concat_menu = input('Option Number: ')

        if (concat_menu == "1"):
            molidutils.show_full_conversion(master_molID_class_list, step='concatenation')
        elif (concat_menu == "2"):
            molidutils.show_unassigned_conversion(master_molID_class_list, step='concatenation')
        elif (concat_menu == "3"):
            concat_submenu_3_1 = 0
            while(concat_submenu_3_1 != "QUIT"):
                search_term, concat_submenu_3_1 = molidutils.get_search_term(concat_submenu_3_1)
                if concat_submenu_3_1 != "QUIT":
                    found_molID_class_chID_map, molID_class_been_copied = molidutils.search_chains(master_molID_class_list, search_term)
                while (concat_submenu_3_1 != "QUIT"):
                    print("""Select one of the following options to proceed:
                             1) Perform new search
                             2) Update new chain ID
                          """)
                    concat_submenu_3_1 = input('Option Number: ')
                    if (concat_submenu_3_1 == "1"):
                        break
                    elif (concat_submenu_3_1 == "2"):
                        newchID = input('New Chain ID: ')
                        found_molID_class_chID_map = molidutils.edit_chain_order(found_molID_class_chID_map, newchID, action='try')
                        molidutils.print_conflicts(found_molID_class_chID_map)
                        print("Would you like to accept or deny these changes?")
                        while (concat_submenu_3_1 != "QUIT"):
                            concat_submenu_3_1 = input('Enter ACCEPT or DENY: ')
                            if (concat_submenu_3_1 == "ACCEPT"):
                                master_molID_class_list = molidutils.accept_newchain(master_molID_class_list, found_molID_class_chID_map)
                            concat_submenu_3_1 = "QUIT"
                            concat_menu = ""
        elif (concat_menu == "4"):
            concat_submenu_4_1 = 0
            while(concat_submenu_4_1 != "QUIT"):
                search_term, concat_submenu_4_1 = molidutils.get_search_term(concat_submenu_4_1)
                if concat_submenu_4_1 != "QUIT":
                    found_molID_class_chID_map, molID_class_been_copied = molidutils.search_chains(master_molID_class_list, search_term)
                while (concat_submenu_4_1 != "QUIT"):
                    print("""Select one of the following options to proceed:
                             1) Perform new search
                             2) Update concatenation order
                          """)
                    concat_submenu_4_1 = input('Option Number: ')
                    if (concat_submenu_4_1 == "1"):
                        break
                    elif (concat_submenu_4_1 == "2"):
                        new_order = input('New concatenation order: ')
                        found_molID_class_chID_map = molidutils.edit_chain_order(found_molID_class_chID_map, int(new_order), action='update')
                        master_molID_class_list = molidutils.accept_newchain(master_molID_class_list, found_molID_class_chID_map)
                        concat_submenu_4_1 = "QUIT"
        elif (concat_menu == "5"):
            concat_submenu_5_1 = 0
            while(concat_submenu_5_1 != "QUIT"):
                search_term, concat_submenu_5_1 = molidutils.get_search_term(concat_submenu_5_1)
                if concat_submenu_5_1 != "QUIT":
                    found_molID_class_chID_map, molID_class_been_copied = molidutils.search_chains(master_molID_class_list, search_term)
                while (concat_submenu_5_1 != "QUIT"):
                    print("""Select one of the following options to proceed:
                             1) Perform new search
                             2) Accept planned concatenation
                          """)
                    concat_submenu_5_1 = input('Option Number: ')
                     if (concat_submenu_5_1 == "1"):
                        break
                    elif (concat_submenu_5_1 == "2"):
                        found_molID_class_chID_map = molidutils.edit_chain_order(found_molID_class_chID_map, int(new_order), action='accept')
                        master_molID_class_list = molidutils.accept_newchain(master_molID_class_list, found_molID_class_chID_map)
                        concat_submenu_5_1 = "QUIT"

        # Option 6) Finalize Curation.
        elif (concat_menu == "6"):
            print("Finalizing Curation ...")

            for my_files in filelist:
                with open(my_files) as myfile:
                    with open(target_dir+"/%s" % (myfile.name), 'w') as newciffile:
                        for molID_class in master_molID_class_list:
                            if (molID_class.file_name == myfile.name):
                                for line in myfile:
                                    if (line[0:4] == "ATOM") or (line[0:6]=="HETATM"):
                                        # Chains outside map should not exist but just in case
                                        line_split = line.strip()
                                        line_split = line.split()
                                        if line_split[17] in molID_class.chID_newchID_map:
                                            # Residues have to be renumbered due to concatenations
                                            if line_split[17] in molID_class.concat_order:
                                                residue_offset = (molID_class.concat_order[line_split[17]] - 1) * 1000
                                                new_resinum = int(line_split[15]) + int(residue_offset)
                                                heresthenew_resinum = int(new_resinum)
                                                newline = line_split[0] + " " + line_split[1] + " " + line_split[2] + " " + line_split[3] + " " + line_split[4] + " " + line_split[5] + " " + molID_class.chID_newchID_map[line_split[17]] + " " + line_split[7] + " " + line_split[8] + " " + line_split[9] + " " + line_split[10] + " " + line_split[11] + " " + line_split[12] + " " + line_split[13] + " " + line_split[14] + " " + str(heresthenew_resinum) + " " + line_split[16] + " " + molID_class.chID_newchID_map[line_split[17]] + " " + line_split[18] + " " + line_split[19] + "\n"
                                                newciffile.write(newline)
                                            else:
                                                newline = line_split[0] + " " + line_split[1] + " " + line_split[2] + " " + line_split[3] + " " + line_split[4] + " " + line_split[5] + " " + molID_class.chID_newchID_map[line_split[17]] + " " + line_split[7] + " " + line_split[8] + " " + line_split[9] + " " + line_split[10] + " " + line_split[11] + " " + line_split[12] + " " + line_split[13] + " " + line_split[14] + " " + line_split[15] + " " + line_split[16] + " " + molID_class.chID_newchID_map[line_split[17]] + " " + line_split[18] + " " + line_split[19] + "\n"
                                                newciffile.write(newline)

                                        else:
                                            newciffile.write(line)
                                    # elif (line[0:6] == "COMPND"):
                                    #     if "CHAIN:" in line:
                                    #         newline = line[0:17] + molID_class.chID_newchID_map[line[17]] + line[18:]
                                    #         newciffile.write(newline)
                                    #     else:
                                    #         newciffile.write(line)
                                    else:
                                        newciffile.write(line)
            concat_menu = "QUIT"
        # Option 6) Finalize Curation.
