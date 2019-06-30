#!/usr/bin/env python
# coding: utf-8

from __future__ import print_function
import sys
import re
import numpy as np
import csv
import os
import copy
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


def pdb_to_masterlist(filelist):
    """
    pdb_to_masterlist
    """
    master_molID_class_list = []
    N=0
    for my_files in filelist:
        N += 1
        print("Reading:"+' '+my_files+"  ("+str(N)+" of "+str(len(filelist))+")")
        with open(my_files) as myfile:
            this_molID_class = make_MolID_cif(myfile)
            master_molID_class_list.append(this_molID_class)
    return master_molID_class_list

def uniquelist_to_conversionlist(unique_map):
    """
    uniquelist_to_conversionlist
    """
    molIDConversion_list = []
    for key in unique_map:
        molIDConversion_list.append(initial_MolIDConversion(key, unique_map[key]))
    return molIDConversion_list

def update_masterlist(master_molID_class_list, molIDConversion_list):
    """
    update_masterlist
    """
    molIDConvert_molID_chID_map = {}
    for molIDConversion in molIDConversion_list:
        molIDConvert_molID_chID_map[molIDConversion.molID] = molIDConversion.chID_list
    # Update master_molID_class_list
    for molID_class in master_molID_class_list:
        for molID in molID_class.molID_chID:
            molID_class.add_chID_newchID_map(molID, molIDConvert_molID_chID_map[molID])
    # Determine which MolID objects now contain conflicts (concatenations)
    for molID_class in master_molID_class_list:
        molID_class.check_for_concatenations()
    return master_molID_class_list

# The class containing all the information about each file neccessary to build
# a conversion template
class MolID(object):
        # file_name = ""
        # chID_newchID_map = {}
        # molID_chID = {}
        # concat_order_map = {}
        # complete_order = {}

    def __init__(self, file_name, chID_newchID_map, molID_chID,
                 concat_order, complete_order):
        self.file_name = file_name
        self.chID_newchID_map = chID_newchID_map
        self.molID_chID = molID_chID
        self.concat_order = concat_order
        self.complete_order = complete_order

    def add_chID_newchID_map(self, molID, newchID_list):
        N = 0
        for oldchID in self.molID_chID[molID]:
            self.chID_newchID_map[oldchID] = newchID_list[N]
            N += 1

    def check_for_concatenations(self):
        self.concat_order = {}
        # Figure out which newchIDs appear more than once which will be
        # concatenated
        usage = {}
        duplicates = {}
        for chID in self.chID_newchID_map:
            if self.chID_newchID_map[chID] not in usage:
                usage[self.chID_newchID_map[chID]] = 1
            else:
                duplicates[self.chID_newchID_map[chID]] = 1
        # Create a map of oldchID to newchID for those being concatenated
        usage = {}
        for chID in self.chID_newchID_map:
            if self.chID_newchID_map[chID] in duplicates:
                if self.chID_newchID_map[chID] in usage:
                    usage[self.chID_newchID_map[chID]] += 1
                else:
                    usage[self.chID_newchID_map[chID]] = 1
                self.concat_order[chID] = usage[self.chID_newchID_map[chID]]
        # Update complete_order
        for chID in self.chID_newchID_map:
            if chID in self.concat_order:
                self.complete_order[chID] = False
            else:
                self.complete_order[chID] = True
    # END check_for_concatenations

    # This will force the complete_order[chID] to the input "complete" which is either True or False
    def force_complete_order(self, chID, complete):
        self.complete_order[chID] = complete

    def update_concat_order(self, chID, neworder):
        otheroldchID = ""
        if chID in self.concat_order:
            newchID = self.chID_newchID_map[chID]
            for oldchID in self.chID_newchID_map:
                if (self.chID_newchID_map[oldchID] == newchID) and (oldchID != chID) and (self.concat_order[oldchID] == neworder):
                    otheroldchID = oldchID
            if (otheroldchID != ""):
                self.concat_order[otheroldchID] = self.concat_order[chID]
                self.concat_order[chID] = neworder

# END MolID class


# make_molID is the function that will grab all relevant information from each
# file input
def make_MolID(myfile):

    file_name = myfile.name
    molID_list = []
    chID_list = []
    molID_chID = {}
    chID_newchID_map = {}
    concat_order = {}
    complete_order = {}
    complete = False
    for line in myfile:
        if "COMPND" in line:
            if "MOLECULE:" in line or "MOLID:" in line:
                line_list = line.rsplit(':')
                molID = line_list[1].strip()
                # Remove leading and trailing spaces
                molID = re.sub(' \{0,\}$|^ \{0,\}', '', molID)
                # Remove any special characters
                molID = re.sub(':|;|’|“|”', '', molID)
                # Add single leading and trailing space and bookmark with ;
                molID = re.sub('^', '; ', molID)
                molID = re.sub('$', ' ;', molID)
                molID_list.append(molID.upper())
            if "CHAIN:" in line:
                line_list = line.rsplit(':')
                chID = line_list[1].strip()
                # Remove all spaces and ;
                chID = re.sub(' |;', '', chID)
                chID = chID.split(',')
                chID_list.append(chID)
    for i in range(len(chID_list)):
        if molID_list[i] not in molID_chID:
            molID_chID[molID_list[i]] = chID_list[i]
        else:
            for chid in chID_list[i]:
                    molID_chID[molID_list[i]].append(chid)
    my_molID_class = MolID(file_name, chID_newchID_map,
                           molID_chID, concat_order, complete_order)
    my_molID_class.check_for_concatenations()
    return my_molID_class
# END make_MolID


# make_molID is the function that will grab all relevant information from each
# file input
def make_MolID_cif(myfile):

    file_name = myfile.name
    chID_newchID_map = {}
    molID_chID = {}
    concat_order = {}
    complete_order = {}
    complete = False

    mmcif_dict = MMCIF2Dict(myfile)

    # CIF files contain entity_id's which are used to link molID and chID
    # Need the entity_id and auth_asym_id correspondence
    # Put into entity_chIDlist_map
    entity_list = mmcif_dict['_atom_site.label_entity_id']
    chID_list = mmcif_dict['_atom_site.auth_asym_id']
    entity_chIDlist_map = {}

    for i in range(len(entity_list)):
        if entity_list[i] not in entity_chIDlist_map:
            entity_chIDlist_map[entity_list[i]] = [chID_list[i]]
        else:
            if chID_list[i] not in entity_chIDlist_map[entity_list[i]]:
                entity_chIDlist_map[entity_list[i]].append(chID_list[i])

    entity_list = mmcif_dict['_entity.id']
    molID_list = mmcif_dict['_entity.pdbx_description']
    molID_list = [i.upper() for i in molID_list]

    for i in range(len(entity_list)):
        # Need this step because some MolIDs may be present whose chains have been removed
        if entity_list[i] in entity_chIDlist_map:
            for chID in entity_chIDlist_map[entity_list[i]]:
                if molID_list[i] not in molID_chID:
                    molID_chID[molID_list[i]] = [chID]
                else:
                    molID_chID[molID_list[i]].append(chID)

    my_molID_class = MolID(file_name, chID_newchID_map,
                           molID_chID, concat_order, complete_order)
    my_molID_class.check_for_concatenations()
    return my_molID_class
# END make_MolID_cif


# MolIDConversion Class. The idea here is to create an object that
# will be modified by the user in the first step to map a MolID to a given
# chain ID with the correct number of occurances. This is the master list of
# molID's
class MolIDConversion(object):
    molID = ""
    chID_list = []
    occur = 0
    complete = bool

    def __init__(self, molID, chID_list, occur, complete):
        self.molID = molID
        self.chID_list = chID_list
        self.occur = occur
        self.complete = complete

    def check_for_completeness(self):
        if (len(self.chID_list) >= self.occur):
            self.complete = True
        else:
            self.complete = False

    def add_chID(self, chID):
        # if chID not in chID_list:
        #    self.chID_list.append(chID)
        self.chID_list.append(chID)

    def add_chID_list(self, chID_list):
        for chID in chID_list:
            # if chID not in self.chID_list:
            #    self.chID_list.append(chID)
            self.chID_list.append(chID)

    def remove_chID_list(self, chID_list):
        for chID in chID_list:
            if chID in self.chID_list:
                self.chID_list.remove(chID)


def initial_MolIDConversion(molID, occur):
    chID_list = []
    complete = False
    my_molID_conversion = MolIDConversion(molID, chID_list, occur, complete)
    return my_molID_conversion


def add_chID_list(MolIDConversion, chID_list):
    for chID in chID_list:
        if chID not in MolIDConversion.chID_list:
            MolIDConversion.chID_list.append(chID)
    return MolIDConversion
#

# Compile information from each file to create a master map of molID to max
# number of occurances
def CreateMasterUniqueMolIDMap(molID_class_list):
    unique_molID_map = {}
    for my_molID_class in molID_class_list:
        for molID in my_molID_class.molID_chID:
            if (unique_molID_map.get(molID) is None):
                unique_molID_map[molID] = \
                    len(my_molID_class.molID_chID[molID])
            elif (unique_molID_map.get(molID) <
                    len(my_molID_class.molID_chID[molID])):
                unique_molID_map[molID] = len(my_molID_class.molID_chID[molID])
    return unique_molID_map

#
#
# Read input file function
def read_input_file(input_cnv_file):
    user_molID_chID_map = {}
    user_molID_list_v = []
    user_chID_list_v = []
    # Read user input file and create user_molID_chID_map
    if (input_cnv_file != ""):
        if (os.path.isfile(input_cnv_file) is True):
            # Reading of input file
            with open(input_cnv_file) as mycnv:
                for line in mycnv:
                    line = line.strip()
                    my_line = line.split(':')
                    molID = my_line[0]
                    chID = my_line[1].split(',')
                    user_molID_list_v.append(molID)
                    user_chID_list_v.append(chID)
            # Mining of input file for information

            # for i in range(len(user_molID_list_v)):
            #     if (user_molID_chID_map.get(user_molID_list_v[i]) is None):
            #         user_molID_chID_map[user_molID_list_v[i]] = user_chID_list_v[i]
            #     elif (user_molID_chID_map.get(user_molID_list_v[i]) is not None):
            #         for this_chID in user_chID_list_v[i]:
            #             if this_chID not in user_molID_chID_map[user_molID_list_v[i]]:
            #                 user_molID_chID_map[user_molID_list_v[i]].append(this_chID)

            for molID_i, chID_i in zip(user_molID_list_v, user_chID_list_v):
                if molID_i not in user_molID_chID_map:
                    user_molID_chID_map[molID_i] = chID_i
                else:
                    for this_chID in chID_i:
                        if this_chID not in user_molID_chID_map[molID_i]:
                            user_molID_chID_map[molID_i].append(this_chID)

        return user_molID_chID_map
# End Read user input file function


####################################
# INTERACTIVE CONVERSION FUNCTIONS #
####################################

def show_full_conversion(current_list, step='conversion'):
    """
    show_full_conversion
    """
    if(step=='conversion'):
        for molIDConversion in current_list:
            molIDCon_chID_list_forPrint = re.sub('\[|\]| |\'', '', str(molIDConversion.chID_list))
            print(molIDConversion.molID+":"+molIDCon_chID_list_forPrint)
    elif(step=='concatenation'):
        for molID_class in current_list:
            for molID in molID_class.molID_chID:
                for chID in molID_class.molID_chID[molID]:
                    if chID in molID_class.concat_order:
                        print(molID_class.file_name + ":" + molID + ":"
                              + chID + ":" +
                              molID_class.chID_newchID_map[chID] + ":"
                              + str(molID_class.concat_order[chID]))
                    else:
                        print(molID_class.file_name + ":" + molID + ":"
                              + chID + ":" +
                              molID_class.chID_newchID_map[chID] + ":" +
                              str(0))

def show_unassigned_conversion(current_list, step='conversion'):
    """
    show_unassigned_conversion
    """
    if(step=='conversion'):
        for molIDConversion in current_list:
            molIDConversion.check_for_completeness()
            if molIDConversion.complete is False:
                molIDCon_chID_list_forPrint = re.sub('\[|\]| |\'', '', str(molIDConversion.chID_list))
                print(str(molIDConversion.occur)+":"+str(molIDConversion.molID)+":"+molIDCon_chID_list_forPrint)
    elif(step=='concatenation'):
        for molID_class in current_list:
            for molID in molID_class.molID_chID:
                for chID in molID_class.molID_chID[molID]:
                    if molID_class.complete_order[chID] is False:
                        print(molID_class.file_name + ":" + molID + ":"
                              + chID + ":" +
                              molID_class.chID_newchID_map[chID] + ":"
                              + str(molID_class.concat_order[chID]))

def add_user_conversion(user_molID_chID_map, molIDConversion_list):
    """
    add_user_conversion:
    Add contents of user's input file to your MolIDConversion class and check for completeness
    """
    for molIDConversion in molIDConversion_list:
        # This is currently strict inclusion but perhaps should be except out
        # of convenience to the user
        for key in user_molID_chID_map:
            if (str(key)==str(molIDConversion.molID)):
                molIDConversion.add_chID_list(user_molID_chID_map[key])
        molIDConversion.check_for_completeness()
    # !! molIDConversion_list has been updated
    return molIDConversion_list

def search_conversion(molIDConversion_list, search_term):
    """
    search_conversion
    """
    search_molIDConversion_list = []
    for molIDConversion in molIDConversion_list:
        if search_term in molIDConversion.molID:
            search_molIDConversion_list.append(molIDConversion)
            molIDCon_chID_list_forPrint = re.sub('\[|\]| |\'', '', str(molIDConversion.chID_list))
            print(molIDConversion.molID+":"+molIDCon_chID_list_forPrint)
    for molIDConversion in search_molIDConversion_list:
        molIDConversion_list.remove(molIDConversion)
    return molIDConversion_list, search_molIDConversion_list

def search_again_conversion(molIDConversion_list, search_molIDConversion_list, search_term):
    """
    search_again_conversion
    """
    new_search_molIDConversion_list = []
    for molIDConversion in search_molIDConversion_list:
        # If search term is found on search_molIDConversion_list
        # then add to new temporary list, if not, add back to
        # molIDConversion_list
        if search_term in molIDConversion.molID:
            new_search_molIDConversion_list.append(molIDConversion)
            molIDCon_chID_list_forPrint = re.sub('\[|\]| |\'', '', str(molIDConversion.chID_list))
            print(molIDConversion.molID+":"+molIDCon_chID_list_forPrint)
        else:
            molIDConversion_list.append(molIDConversion)
    search_molIDConversion_list = new_search_molIDConversion_list
    return molIDConversion_list, search_molIDConversion_list

def edit_chain_conversion(molIDConversion_list,search_molIDConversion_list, chID_list, action='add'):
    """
    edit_chain_conversion
    """
    chID_list = chID_list.split(',')
    # Performing chain addition to searched terms
    for molIDConversion in search_molIDConversion_list:
        if(action=='add'):
            molIDConversion.add_chID_list(chID_list)
        elif(action=='remove'):
            molIDConversion.remove_chID_list(chID_list)
    # Adding modified molIDConversions back to molIDConversion_list
    for molIDConversion in search_molIDConversion_list:
        molIDConversion.check_for_completeness()
        molIDConversion_list.append(molIDConversion)
    return molIDConversion_list, search_molIDConversion_list

def check_complete(molIDConversion_list):
    """
    check_complete
    """
    num_unassigned = 0
    for molIDConversion in molIDConversion_list:
        if molIDConversion.complete is False:
            num_unassigned += 1
    if (num_unassigned > 0):
        input_menu_complete = "0"
    else:
        input_menu_complete = "1"
    return input_menu_complete

#######################################
# INTERACTIVE CONCATENATION FUNCTIONS #
#######################################

def problem_counter(master_molID_class_list):
    """
    problem_counter
    """
    count_problems = 0
    for molID_class in master_molID_class_list:
        for chID in molID_class.complete_order:
            if molID_class.complete_order[chID] is False:
                count_problems += 1
    return count_problems
