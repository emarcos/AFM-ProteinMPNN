#!/usr/bin/env python3
import os
import sys
import matplotlib
import re
import glob
import numpy as np
from argparse import ArgumentParser
from prody import *
import json
import subprocess
import matplotlib.pyplot as plt
from matplotlib import cm
import seaborn as sns
import pandas as pd
from scipy.stats import linregress
import math
from matplotlib.ticker import FuncFormatter

### Some functions needed ###

def tmalign_result( soln, query, norm=True ):

        # path to TMalign executable
        TMALIGN_BIN = '/data/Software/MMalign/MMalign'
        rmsd = 9999.99
        tmscore = 0.0

        # run TMalign with option '-a' to normalize by average length of two chains
        if norm:
           out, err = subprocess.Popen( (TMALIGN_BIN, soln, query, '-a T'), stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8' ).communicate()
        else:
           out, err = subprocess.Popen( (TMALIGN_BIN, soln, query), stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8' ).communicate()

        # if TMalign error, report to stdout
        if err:
                print('Error during TMalign of {0} and {1}:\n{2}'.format( soln, query, err ))
                return tmscore, rmsd
        # parse TMalign output
        for outline in out.split( '\n' ):
                if outline.startswith( 'Aligned' ):
                        rmsd_str = re.split( ', |=', outline)[3]
                        rmsd = float( rmsd_str )
                if outline.startswith( 'TM-score' ):
                        tmscore_str = re.split( ' ', outline)[1]
                        tmscore = float( tmscore_str )

        return tmscore, rmsd


def AAContentPDB(pdb):
    aa31 = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'}

    filein = open(pdb,'r')
    fasta=''
    for line in filein:
        if line.startswith('ATOM') and 'CA' in line:
            res = line[17:20]
            aa = aa31[res]
            fasta+=aa
    return fasta


def AAContentPDB_ligand(pdb):
    aa31 = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'}

    filein = open(pdb,'r')
    fasta_A=''
    fasta_B=''
    fasta_C=''
    for line in filein:
        if line.startswith('ATOM') and 'CA' in line:
            if line[21] == 'A':
                res = line[17:20]
                aa = aa31[res]
                fasta_A+=aa
            elif line[21] == 'B':
                res = line[17:20]
                aa = aa31[res]
                fasta_B+=aa
            elif line[21] == 'C':
                res = line[17:20]
                aa = aa31[res]
                fasta_C+=aa

    return fasta_A, fasta_B, fasta_C


def interface_PAE_model_competition(json_file, pdb):
    f = open(json_file)
    data = json.load(f)

    res1 = 0

    pack_list = []

    for dist_list in data['pae']:
        res1 +=1
        res2 = 1
        for dist in dist_list:
            pack = (res1, res2, dist)
            pack_list.append(pack)
            res2 +=1

    pdb_file = parsePDB(pdb)
    n_chainA = len(pdb_file.select("protein and chain A and ca"))
    n_chainB = len(pdb_file.select("protein and chain B and ca"))
    n_chainC = len(pdb_file.select("protein and chain C and ca"))

    PAE_A_list=[]
    PAE_B_list=[]

    for r1,r2,d in pack_list:
        if (r1<=n_chainA and r2>(n_chainA+n_chainB) or (r2<=n_chainA and r1>(n_chainA+n_chainB))):
            PAE_A_list.append(d)
        elif ((n_chainA<r1<=(n_chainA+n_chainB)) and r2>(n_chainA+n_chainB) or (n_chainA<r2<=(n_chainA+n_chainB)) and r1>(n_chainA+n_chainB)):
            PAE_B_list.append(d)

    interface_AC_PAE=np.mean(PAE_A_list)
    interface_BC_PAE=np.mean(PAE_B_list)

    return interface_AC_PAE, interface_BC_PAE


def GetChainResidues(pdbfile):
    chain_residues={'A':[],'B':[]}
    for line in open(pdbfile):
        if line[:4]=="ATOM":
            ch=line.split()[4]
            respos=int(line.split()[5])
            chain_residues[ch].append(respos)
    return chain_residues['A'][0],chain_residues['A'][-1],chain_residues['B'][0],chain_residues['B'][-1]


def GetChainLenght(pdb):
    pdb_file = parsePDB(pdb)
    n_ligand = len(pdb_file.select("protein and chain A and ca"))
    n_receptor = len(pdb_file.select("protein and chain B and ca"))
    return n_ligand, n_receptor


def CalcRmsd(original_pdb, model_pdb):
    #ORIGINAL OR REFERENCE SEQUENCE
    ref_TOTAL = parsePDB(original_pdb)
    ref_fasta = AAContentPDB_ligand(original_pdb)[0]
    ref_chainA = ref_TOTAL.select("protein and chain A").ca
    ref_chainB = ref_TOTAL.select("protein and chain B").ca

    #MODEL SEQUENCE
    mod_TOTAL = parsePDB(model_pdb)
    mod_fasta_A = AAContentPDB_ligand(model_pdb)[0]
    mod_fasta_B = AAContentPDB_ligand(model_pdb)[1]
    mod_fasta_C = AAContentPDB_ligand(model_pdb)[2]
    mod_chainA = mod_TOTAL.select("protein and chain A").ca
    mod_chainB = mod_TOTAL.select("protein and chain B").ca
    mod_chainC = mod_TOTAL.select("protein and chain C").ca

    #RMSD AA, RMSD AB
    ## chain A rmsd while alinning in B (interchain rmsd)
    t = calcTransformation(mod_chainC, ref_chainB) #calculate the transformation needed to align chain B of the model to chain B of reference pdb
    sup = applyTransformation(t, mod_TOTAL) # apply this transformation to the model pdb
    sup_chainA = sup.select("protein and chain A").ca #select ca atoms of chain A of the model (already aligned)
    sup_chainB = sup.select("protein and chain B").ca #select ca atoms of chain A of the model (already aligned)
    if ref_fasta == mod_fasta_A:
        rmsd_AB = calcRMSD(ref_chainA, sup_chainA)
        t_a = calcTransformation(mod_chainA, ref_chainA)
        sup_a = applyTransformation(t_a, mod_TOTAL)
        sup_a_chainA = sup_a.select("protein and chain A").ca
        rmsd_AA = calcRMSD (ref_chainA, sup_a_chainA)
        return rmsd_AB
    elif ref_fasta == mod_fasta_B:
        rmsd_BC = calcRMSD(ref_chainA, sup_chainB)
        t_b = calcTransformation(mod_chainB, ref_chainA)
        sup_b = applyTransformation(t_b, mod_TOTAL)
        sup_b_chainB = sup_b.select("protein and chain B").ca
        rmsd_CC = calcRMSD (ref_chainA, sup_b_chainB)
        return rmsd_BC


def SelectInterfaceOriginal(original_pdb):
    original = parsePDB(original_pdb)
    # select the receptor and ligand
    receptor = original.select('protein and chain B and name CB')
    ligand = original.select('protein and chain A and name CB')

    # calculate the contacts of both of them
    receptor_contacts = Contacts(receptor)
    ligand_contacts = Contacts(ligand)

    # Atoms at the interface
    int_receptor = receptor_contacts.select(8, ligand)
    int_ligand = ligand_contacts.select(8, receptor)

    # Residues at the interface (nums)
    int_receptor_res = int_receptor.getResnums()
    int_ligand_res = int_ligand.getResnums()

    return int_receptor_res, int_ligand_res


def SelectInterface(model_pdb):
    model = parsePDB(model_pdb)
    # select the receptor and ligand
    receptor = model.select('protein and chain C and name CB')
    ligand_A = model.select('protein and chain A and name CB')
    ligand_B = model.select('protein and chain B and name CB')

    # calculate the contacts of both of them
    receptor_contacts = Contacts(receptor)
    ligand_A_contacts = Contacts(ligand_A)
    ligand_B_contacts = Contacts(ligand_B)

    # Atoms at the interface
    #int_receptor = receptor_contacts.select(8, ligand)
    int_ligand_A = ligand_A_contacts.select(8, receptor)
    int_ligand_B = ligand_B_contacts.select(8, receptor)

    # Residues at the interface (nums)
    #int_receptor_res = int_receptor.getResnums()
    if int_ligand_A:
        int_ligand_A_res = int_ligand_A.getResnums()
    else:
        int_ligand_A_res = []
    if int_ligand_B:
        int_ligand_B_res = int_ligand_B.getResnums()
    else:
        int_ligand_B_res = []
    return int_ligand_A_res, int_ligand_B_res


def exponential_normalize_value(value, a):
    normalized_value = 0 + math.exp(-a * value)
    return normalized_value


colab_multimer_competition = pd.DataFrame(columns=['name_A','name_B','rank','model','IntOriginal_PLDDT_A','IntOriginal_PLDDT_B','IntModel_PLDDT_A','IntModel_PLDDT_B','rmsd_AC','rmsd_BC', 'PAE_AC_model', 'PAE_BC_model'])

#You are supposed to have n folders (one per design) with the AFM predictions and original pdb there

lines = glob.glob('6m0j*')
for line in lines:
    name_A = line.split('_vs_')[0].split('_and_')[1].strip()
    name_B = line.split('_vs_')[1].strip()
    os.chdir(f'{line.strip()}')
    original_name_A = f'{name_A}.pdb'
    original_name_B = f'{name_B}.pdb'
    n_A_ligand, n_A_receptor = GetChainLenght(original_name_A)
    n_B_ligand, n_B_receptor = GetChainLenght(original_name_B)
    #interface_selection
    oriA_ligand_interface = SelectInterfaceOriginal(original_name_A)[1]
    oriB_ligand_interface = SelectInterfaceOriginal(original_name_B)[1]
    #Calculate the num resi of the original interface at the model pdb
    oriA_ligand_interface_mod = [x for x in oriA_ligand_interface]
    oriB_ligand_interface_mod = [y for y in oriB_ligand_interface]
    #models
    models = glob.glob('*_relaxed_rank_*_model_*.pdb')
    for model in models:
        model_TOTAL = parsePDB(model)
        if not model == []:
            rank_short=f'rank_{model.split("_rank_")[1][:3]}'
            model_short=f'model_{model.split("_model_")[1][0]}'
            #calculate rmsd for both
            rmsd_A = CalcRmsd(original_name_A, model)
            rmsd_B = CalcRmsd(original_name_B, model)
            #calculate interface PAE for both
            json_file_model = f'{name_A}-{name_B}_mpnn_custom_msa_scores_{rank_short}_alphafold2_multimer_v2_{model_short}_seed_000.json'
            AC_PAE_model = interface_PAE_model_competition(json_file_model,model)[0]
            BC_PAE_model = interface_PAE_model_competition(json_file_model,model)[1]
            #Select residues at the interfaces
            ligand_A_interface, ligand_B_interface = SelectInterface(model)
            #calculate PLDDT at the ORIGINAL interfaces
            IntOriginal_plddt_A = np.mean((model_TOTAL.select(f"protein and chain A and resnum {' '.join(map(str, oriA_ligand_interface_mod))}").getBetas()))
            IntOriginal_plddt_B = np.mean((model_TOTAL.select(f"protein and chain B and resnum {' '.join(map(str, oriB_ligand_interface_mod))}").getBetas()))
            #calculate PLDDT at the MODEL interfaces
            if len(ligand_A_interface) > 0:
                IntModel_plddt_A = np.mean((model_TOTAL.select(f"protein and chain A and resnum {' '.join(map(str, ligand_A_interface))}").getBetas()))
            else:
                IntModel_plddt_A = 0.001
            if len(ligand_B_interface) > 0:
                IntModel_plddt_B = np.mean((model_TOTAL.select(f"protein and chain B and resnum {' '.join(map(str, ligand_B_interface))}").getBetas()))
            else:
                IntModel_plddt_B = 0.001
            #Append all data to the pandas data frame
            colab_multimer_competition = colab_multimer_competition.append({'name_A':name_A,'name_B':name_B,'rank':rank_short,'model':model_short,'IntOriginal_PLDDT_A':IntOriginal_plddt_A,'IntOriginal_PLDDT_B':IntOriginal_plddt_B,'IntModel_PLDDT_A':IntModel_plddt_A,'IntModel_PLDDT_B':IntModel_plddt_B,'rmsd_AC':rmsd_A,'rmsd_BC':rmsd_B, 'PAE_AC_model':AC_PAE_model, 'PAE_BC_model':BC_PAE_model}, ignore_index=True)
        else:
            print ('There is an error in '+line)

### Save the data in a csv file
colab_multimer_competition.to_csv(f'./competition_all_data.csv')

### Processing the data and getting scores ###
# Counting win or lose
for i in colab_multimer_competition.index:
    if colab_multimer_competition.loc[i,'rmsd_AC']<=2.5 and colab_multimer_competition.loc[i,'rmsd_BC']>2.5:
        colab_multimer_competition.loc[i,'counts_A']=1
        colab_multimer_competition.loc[i,'counts_B']=0
    elif colab_multimer_competition.loc[i,'rmsd_AC']>2.5 and colab_multimer_competition.loc[i,'rmsd_BC']<=2.5:
        colab_multimer_competition.loc[i,'counts_B']=1
        colab_multimer_competition.loc[i,'counts_A']=0
    elif colab_multimer_competition.loc[i,'rmsd_AC']>2.5 and colab_multimer_competition.loc[i,'rmsd_BC']>2.5:
        colab_multimer_competition.loc[i,'counts_A']=0
        colab_multimer_competition.loc[i,'counts_B']=0
#plddt values between 0 and 1
colab_multimer_competition['IntOriginal_PLDDT_A']=colab_multimer_competition['IntOriginal_PLDDT_A']/100
colab_multimer_competition['IntOriginal_PLDDT_B']=colab_multimer_competition['IntOriginal_PLDDT_B']/100
colab_multimer_competition['IntModel_PLDDT_A']=colab_multimer_competition['IntModel_PLDDT_A']/100
colab_multimer_competition['IntModel_PLDDT_B']=colab_multimer_competition['IntModel_PLDDT_B']/100
#normalization of rmsd (get a value between 0 and 1)
colab_multimer_competition['normalized_rmsd_AC']=colab_multimer_competition['rmsd_AC'].apply(lambda x: exponential_normalize_value(x, 0.05))
colab_multimer_competition['normalized_rmsd_BC']=colab_multimer_competition['rmsd_BC'].apply(lambda x: exponential_normalize_value(x, 0.05))
#normalization of PAE (get a value between 0 and 1)
colab_multimer_competition['normalized_PAE_AC_model']=colab_multimer_competition['PAE_AC_model'].apply(lambda x: exponential_normalize_value(x, 0.03))
colab_multimer_competition['normalized_PAE_BC_model']=colab_multimer_competition['PAE_BC_model'].apply(lambda x: exponential_normalize_value(x, 0.03))
#Calculate competition score(i)
colab_multimer_competition['score_A']=colab_multimer_competition['IntModel_PLDDT_A']*colab_multimer_competition['normalized_rmsd_AC']*colab_multimer_competition['normalized_PAE_AC_model']
colab_multimer_competition['score_B']=colab_multimer_competition['IntModel_PLDDT_B']*colab_multimer_competition['normalized_rmsd_BC']*colab_multimer_competition['normalized_PAE_BC_model']

### Getting Competitive Binding Scores(i-j)
# Select relevant columns for counts and scores
counts = colab_multimer_competition[['name_A', 'name_B', 'counts_A', 'counts_B']]
scores = colab_multimer_competition[['name_A', 'name_B', 'score_A', 'score_B']]
# Group counts and scores by name_A and name_B, aggregating appropriately
counts_grouped = counts.groupby(['name_A', 'name_B'], as_index=False).sum()
scores_grouped = scores.groupby(['name_A', 'name_B'], as_index=False).mean()
# Merge grouped counts and scores on name_A and name_B
data_values = pd.merge(scores_grouped, counts_grouped, on=['name_A', 'name_B'])
# Calculate comparative metrics
data_values['comparative_score'] = data_values['score_A'] - data_values['score_B']
data_values['comparative_counts'] = data_values['counts_A'] - data_values['counts_B']
# Reset the index to ensure a clean DataFrame
data_values.reset_index(drop=True, inplace=True)
data_values.to_csv('./competition_matrix.csv')

### Getting Global Competition Score
pdbs = set(data_values['name_A'].unique()) | set(data_values['name_B'].unique())
pdbs = sorted(list(pdbs))
matrix = pd.DataFrame(index=pdbs, columns=pdbs)
for i in range(len(data_values)):
    ligand_1 = data_values.iloc[i]['name_A']
    ligand_2 = data_values.iloc[i]['name_B']
    score = data_values.iloc[i]['comparative_score']
    matrix.at[ligand_1, ligand_2] = score
    matrix.at[ligand_2, ligand_1] = score*-1
# Fill diagonal with zeros
matrix.values[[range(len(pdbs))]*2] = 0
pd_matrix = pd.DataFrame(matrix)
# Creating the dictionary
row_sums_dict = {index: row.sum() for index, row in pd_matrix.iterrows()}
# Sorting the dictionary based on values
global_score = dict(sorted(row_sums_dict.items(), key=lambda x: x[1], reverse=True))
# get a sorted list of the designs based on their global score
sorted_names = list(global_score.keys())
# Reorder matrix2 based on the sorted row and column labels
sorted_matrix = matrix.loc[sorted_names, sorted_names]

# Convert the dictionary to a list of tuples
data = list(global_score.items())
# Create a DataFrame from the list of tuples
df_global_score_competition = pd.DataFrame(data, columns=['name', 'Global_score_competition'])
df_global_score_competition.to_csv('./competition_global_score.csv')
