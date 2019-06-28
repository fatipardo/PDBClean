#!/usr/bin/env python
# coding: utf-8

from __future__ import print_function
from __future__ import division
import sys
import argparse
import re
import numpy as np
import csv
import os
import copy
from Bio.PDB.MMCIFParser import FastMMCIFParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC
from Bio import pairwise2
from urllib.request import urlopen


ap = argparse.ArgumentParser(description="A random script")
ap.add_argument('file', type=str, nargs='+', help='Input file to read')
filelist = ap.parse_args()

#First file contains a list of sequences to search
# seqres file is hardcoded in
# Will look for previous curation summary and read curated list and ignore list
# This is currently hardcoded and will need to be updated when the final decision
# is made on the curation summary file.



##
## urlseqres = pdb_seqres.txt location
urlseqres = 'ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt'
## urlcif = prefix to the location of the cif file. Needs key + '.cif' ending
urlcif = 'https://files.rcsb.org/download/'


Seq = []
for my_file in filelist.file:
    with open(my_file) as myfile:
        # Save sequences to search
        for line in myfile:
            Seq.append(line.strip())

MatchingSeq = []
SeqKeyMap = {}

print('Downloading latest pdb_seqres.txt file')
seqresfile = urlopen(urlseqres)
#with open("pdb_seqres.txt") as seqresfile:

print('Matching to input sequence')
N = 0
for line in seqresfile:
    line = str(line)
    line = re.sub("b'", "", line)
    line = line[0:len(line)-3]

    if (line[0]==">"):
        key = line[1:5]
    for seq in Seq:
        if seq in str(line):
            if seq in SeqKeyMap:
                SeqKeyMap[seq].append(key)
            else:
                SeqKeyMap[seq] = [key]

MasterKey = SeqKeyMap[Seq[0]]
for seq in SeqKeyMap:
    MasterKey = sorted(list(set(MasterKey) & set(SeqKeyMap[seq])))

#curcode is the RCSB code for the structure to be downloaded
curcode_list = []
if os.path.exists("CuratedFiles.txt"):
    print('Exluding previously curated files')
    with open('CuratedFiles.txt') as myfile:
        prevcurcode_list = []
        for line in myfile:
            prevcurcode_list.append(line.strip())

    for key in MasterKey:
        if key not in prevcurcode_list:
            curcode_list.append(key)

else:
    for key in MasterKey:
        curcode_list.append(key)


print(len(curcode_list))
print('Printing new cif files')
if not os.path.exists("NewCIF"):
    os.mkdir("NewCIF")
else:
    os.popen('rm -f ./NewCIF/*')

with open('NewCIFList.txt', 'w') as curcodefile:
    for curcode in curcode_list:
        curcodefile.write(curcode + "\n")


for curcode in curcode_list:
    print("Downloading PDB Entry: " + curcode)
    urlcifnow = urlcif + curcode.upper() + '.cif'
    cifstruct = urlopen(urlcifnow)

    newciffilename = 'NewCIF/' + curcode + '.cif'
    with open(newciffilename, 'w') as newciffile:
        for line in cifstruct:
            line = str(line)
            #line = re.sub("b'", "", line)
            #line = line[0:len(line)-3]
            newciffile.write(str(line[2:len(line)-2]) + "\n")

print('Complete')
