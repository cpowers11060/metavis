'''
Author: Christopher Powers
Institution: University of Rhode Island

This script will add transport and export
reactions based on a two column table in
the format of:
compound_id	compartment_from
and add a transport reaction from that compartment
to the extracellular compartment. Additionally, it
will add an exchange reaction allowing that to leave
the model as well
'''


import argparse
import re
import sys
import os

# Generate the argparse object
parser = argparse.ArgumentParser(description='Overwrite model gene assocations')
parser.add_argument('--cpds', type=str,
    help='path to a file of compounds of interest')
parser.add_argument('--out', type=str,
    help='path to the output')
parser.add_argument('--ext', type=str,
    help='extracellular compartment')
args = parser.parse_args()

# generate a dictionary of reactions to gene associations
with open(args.cpds, 'r') as infile:
    with open('{}/reactions.yaml'.format(args.out), 'a') as outfile:
        for line in infile:
            line=line.rstrip()
            listall=re.split('\t', line)
            outfile.write('- id: {}_TP_{}{}\n'.format(listall[0], listall[1], listall[2]))
            outfile.write('  equation: {}[{}] <=> {}[{}]\n'.format(listall[0], listall[1], listall[0], listall[2]))
            if listall[1] == args.ext or listall[2] == args.ext:
                outfile.write('- id: {}_EX\n'.format(listall[0]))
                outfile.write('  equation: {}[{}] <=>\n'.format(listall[0], args.ext))            

with open(args.cpds, 'r') as infile:
    with open('{}/model_def.tsv'.format(args.out), 'a') as outfile:
        for line in infile:
            line=line.rstrip()
            listall=re.split('\t', line)
            outfile.write("{}_TP_{}{}\n".format(listall[0], listall[1], listall[2]))
            if listall[1] == args.ext or listall[2] == args.ext:
                outfile.write("{}_EX\n".format(listall[0]))
