#!/bin/python
"""
Generate peptide mol and gaussian job file(gjf) file
"""
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from itertools import product
import argparse
import datetime
import time
import os

def parse_args():
    parser = argparse.ArgumentParser(description="Main script for generate 2D /3D mol file for specific length of peptide")
    parser.add_argument('--dir', action='store', dest='dir', default="./3d", help='path for save file.')
    parser.add_argument('--num', action='store', dest='num', default=2, type=int, help='length of peptide')
    parser.add_argument('--AA', action='store', dest='aa_list', type=str,default="#", help='specific amino acid list for generate peptide')
    parser.add_argument('--3d', action='store', dest='3d', default=False, type=bool, help='generate 3d structure for peptide')
    parser.add_argument('--type', action='store', dest='type', default='gjf', type=str, help='generate file type for peptide')
    args = vars(parser.parse_args())
    return args


def mol2gjf(m):
    f = open(args['dir']+"/"+m.GetProp('_Name')+".gjf",'w')
    #keywords = '#HF/6-31G* SCF=tight Test Pop=MK iop(6/33=2) iop(6/42=6)'
    keywords = '#B3LYP/6-31G* opt(maxcyc=300) scf(maxcyc=300)'
    charge,multiplicity = 0,1
    name = m.GetProp('_Name')
    line = "%%nproc=%d\n%%chk=%s\n%%mem=%dMB\n%s\n\n%s\n\n%d %d\n" % (24,name,4096,keywords,name,charge,multiplicity)
    f.write(line)
    atoms_num = m.GetNumAtoms()
    for i in range(atoms_num):
        x,y,z = m.GetConformer().GetAtomPosition(i)
        sybol = m.GetAtomWithIdx(i).GetSymbol()
        atom_line = "%s%20.5f%20.5f%20.5f\n" % (sybol,x,y,z)
        f.write(atom_line)
    f.write('\n')
    f.close() 

def gen_peptide(args):
    aa_all_list = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    if args['aa_list'] == "#": aa_list = aa_all_list 
    else: aa_list = list(set(args['aa_list'].split(','))) 
    pep_str_list = [''.join(x) for x in product(aa_list,repeat=args['num'])]
    print(len(pep_str_list))
    print(pep_str_list)
    total_cnt = len(pep_str_list)
    for i,aa in enumerate(pep_str_list):
        print("processing ",i+1,"/",total_cnt," ", aa)
        m = Chem.MolFromFASTA(aa)
        m.SetProp('_Name',aa)
        if args['3d']:
            m = Chem.AddHs(m)
            AllChem.EmbedMolecule(m, randomSeed=10)
        molblock = Chem.MolToMolBlock(m)
        with open (args['dir']+"/"+aa+".mol",'w') as f: f.write(molblock)
        if args['type'] == 'gjf': mol2gjf(m)

if __name__ == "__main__":
    starttime = datetime.datetime.now()
    args = parse_args()
    gen_peptide(args)
    endtime = datetime.datetime.now()
    print ('Congratulations! run time:{} seconds'.format((endtime - starttime).seconds))


