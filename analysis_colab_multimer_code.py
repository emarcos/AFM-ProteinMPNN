import os, glob, re, sys
from prody import *
import re
import json
import subprocess
import numpy as np
import pandas as pd
from tqdm import tqdm
import re
import ast
import math
from pyrosetta import init, pose_from_pdb
from pyrosetta import *
from pyrosetta.rosetta import *
import pyrosetta.rosetta.core.select.residue_selector as selector
import warnings
warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)

# START OF FUNCTIONS # 

# metrics generation functions
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

    pdb_dir = os.path.abspath(os.path.join('..', '..', pdb))
    filein = open(pdb_dir,'r')
    fasta=''
    for line in filein:
        if line.startswith('ATOM') and 'CA' in line:
            res = line[17:20]
            aa = aa31[res]
            fasta+=aa
    return fasta

def get_chain_resindices(pdb, third_chain=None):

       ## reading pdb chains
    input_pdb=parsePDB(pdb)
    chainA_pdb=input_pdb.select("protein and chain A")
    chainB_pdb=input_pdb.select("protein and chain B")

    ## obtaining list of residue indices for each chain
    chainA_resindices = [i+1 for i in set(chainA_pdb.getResindices())] # we sum 1 to each index as prody starts indexing at 0
    chainB_resindices = [i+1 for i in set(chainB_pdb.getResindices())] # same as above
    
    if third_chain != None:
        chainCDE_pdb=input_pdb.select(f"protein and chain {third_chain}")
        chainCDE_resindices = [i+1 for i in set(chainCDE_pdb.getResindices())]
        return chainA_resindices, chainB_resindices, chainCDE_resindices
    else:
        return chainA_resindices, chainB_resindices
 
def interface_PAE(json_file, pdb, third_chain=None):

    # READING JSON
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

    ''' REMOVED FROM ORIGINAL
    res1 = data[0]["residue1"] # list of residue indexes being compared with all other residues
    res2 = data[0]["residue2"] # second list of residue indexes being compared with res1
    distance = data[0]["distance"] # list of distances for each residue pair if we zipped both previous lists

    #seq = AAContentPDB(pdb)
    #you may change the aa that "separate" your both ligand and receptor chain"
    #cha_nres = len(seq[:seq.index('AFTVTVPKDL')])
    '''
    # CALCULATING INTERFACE PAE
    PAE_list=[]
    
    # obtaining each chain residues and calculating PAE
    if third_chain != None:
        chainA_resindices, chainB_resindices, chainCDE_resindices = get_chain_resindices(pdb, third_chain)
        for r1, r2, d in pack_list:
            if r1 in chainA_resindices and (r2 in chainB_resindices or r2 in chainCDE_resindices):
                PAE_list.append(d)
    else: 
        chainA_resindices, chainB_resindices = get_chain_resindices(pdb)
        for r1, r2, d in pack_list:
            if r1 in chainA_resindices and r2 in chainB_resindices:
                PAE_list.append(d)

    interface_PAE=np.mean(PAE_list)
    return interface_PAE

def calc_secondary_structure(pdbfile, type='beta', subset_residues_A=None):
    pose = pose_from_pdb(pdbfile)

    # Define Chain Selectors
    chain_A_selector = selector.ChainSelector()
    chain_A_vector_string = pyrosetta.rosetta.utility.vector1_std_string()
    chain_A_vector_string.extend(['A'])
    chain_A_selector.set_chain_strings(chain_A_vector_string)
    chain_A_resis = list(pyrosetta.rosetta.core.select.get_residues_from_subset(chain_A_selector.apply(pose)))
    print(f"Chain A residues: {chain_A_resis}")

    # Generate secondary structure
    DSSP = protocols.moves.DsspMover()
    DSSP.apply(pose)
    ss = pose.secstruct()

    print(ss)

    beta_sheet_symbols = ['B', 'E']
    alpha_helix_symbols = ['G', 'H', 'I']

    # Calculate residues in secondary structure from chain A
    if type == 'beta':
        beta_A_resis = [resi for resi in chain_A_resis if ss[resi-1] in beta_sheet_symbols]
        print(f"Beta residues in Chain A: {beta_A_resis}")

        return beta_A_resis
    
    elif type == 'all':
        beta_A_resis = [resi for resi in chain_A_resis if ss[resi-1] in beta_sheet_symbols]
        alpha_A_resis = [resi for resi in chain_A_resis if ss[resi-1] in alpha_helix_symbols]
        print(f"Beta residues in Chain A: {beta_A_resis}")
        print(f"Alpha residues in Chain A: {alpha_A_resis}")

        return alpha_A_resis, beta_A_resis

def calc_interface_pyrosetta(pdbfile, third_chain=None):

    init('-mute all')

    pose = pose_from_pdb(pdbfile)

    ## DEFINING CHAIN SELECTORS ##
    # A
     # Define Chain Selectors
    print("Defining chain selectors for chain A and receptor chains")
    chain_A_selector = selector.ChainSelector()
    chain_A_vector_string = pyrosetta.rosetta.utility.vector1_std_string()
    chain_A_vector_string.extend(['A'])
    chain_A_selector.set_chain_strings(chain_A_vector_string)
    chain_A_resis = list(pyrosetta.rosetta.core.select.get_residues_from_subset(chain_A_selector.apply(pose)))
    print(f"Chain A residues: {chain_A_resis}")

    # B
    receptor_selector = selector.OrResidueSelector()
    chain_B_selector = selector.ChainSelector()
    chain_B_vector_string = pyrosetta.rosetta.utility.vector1_std_string()
    chain_B_vector_string.extend(['B'])
    chain_B_selector.set_chain_strings(chain_B_vector_string)
    chain_B_resis = list(pyrosetta.rosetta.core.select.get_residues_from_subset(chain_B_selector.apply(pose)))
    print(f"Chain B residues: {chain_B_resis}")

    chain_A_positions_dict = {r:i+1 for i,r in enumerate(chain_A_resis)}
    chain_B_positions_dict = {r:i+1 for i,r in enumerate(chain_B_resis)}

    receptor_selector.add_residue_selector(chain_B_selector)

    # C
    if third_chain is not None:
        chain_C_selector = selector.ChainSelector()
        chain_C_vector_string = pyrosetta.rosetta.utility.vector1_std_string()
        chain_C_vector_string.extend(['C'])
        chain_C_selector.set_chain_strings(chain_C_vector_string)
        chain_C_resis = list(pyrosetta.rosetta.core.select.get_residues_from_subset(chain_C_selector.apply(pose)))
        print(f"Chain C residues: {chain_C_resis}")

        chain_C_positions_dict = {r:i+1 for i,r in enumerate(chain_C_resis)}

        receptor_selector.add_residue_selector(chain_C_selector)

    print("Calculating residues at interface for multimer predictions")
    interface_by_vector_selector = selector.InterGroupInterfaceByVectorSelector()
    interface_by_vector_selector.group1_selector(chain_A_selector)
    interface_by_vector_selector.group2_selector(receptor_selector)
    interface_by_vector_selector.cb_dist_cut(11.0)
    interface_by_vector_selector.nearby_atom_cut(5.5)
    interface_resis = list(pyrosetta.rosetta.core.select.get_residues_from_subset(interface_by_vector_selector.apply(pose)))
    print(f"Interface residues: {interface_resis}")

    # Obtaining residues from binder and receptor
    chain_A_resis = list(pyrosetta.rosetta.core.select.get_residues_from_subset(chain_A_selector.apply(pose)))
    print(f"Chain A residues: {chain_A_resis}")
    chain_B_resis = list(pyrosetta.rosetta.core.select.get_residues_from_subset(chain_B_selector.apply(pose)))
    if third_chain is not None:
        chain_C_resis = list(pyrosetta.rosetta.core.select.get_residues_from_subset(chain_C_selector.apply(pose)))
    # receptor_resis = list(pyrosetta.rosetta.core.select.get_residues_from_subset(receptor_selector.apply(pose)))
    # print(f"Receptor residues: {receptor_resis}")

    interface_resis_A = []
    interface_resis_B = []
    interface_resis_C = []
    for residue in interface_resis:
        if residue in chain_A_resis:
            interface_resis_A.append(residue)
        if residue in chain_B_resis:
            interface_resis_B.append(residue)
        if third_chain is not None:
            if residue in chain_C_resis:
                interface_resis_C.append(residue)

    interface_positions_A = [chain_A_positions_dict[r] for r in interface_resis_A]
    interface_positions_B = [chain_B_positions_dict[r] for r in interface_resis_B]

    print(f"Interface residues in Chain A: {interface_resis_A}")
    print(f"Interface residues in Chain B: {interface_resis_B}")
    print(f"Interface residues in Chain A starting from 1: {interface_positions_A}")
    print(f"Interface residues in Chain B starting from 1: {interface_positions_B}")
    
    if third_chain is not None:
        interface_positions_C = [chain_C_positions_dict[r] for r in interface_resis_C]
        print(f"Interface residues in Chain C: {interface_resis_C}")
        print(f"Interface residues in Chain C starting from 1: {interface_positions_C}")

        return interface_positions_A, interface_positions_B, interface_positions_C
    else:
        return interface_positions_A, interface_positions_B
    
def get_consecutive_numbers(numbers):
    if not numbers:
        return []

    sublists = [[numbers[0]]]

    for i in range(1, len(numbers)):
        # Check if the current number is consecutive to the last one in the last sublist
        if numbers[i] == numbers[i - 1] + 1:
            sublists[-1].append(numbers[i])
        else:
            sublists.append([numbers[i]])

    return sublists

def fuse_overlapping_sublists(sublists):
    # Convert sublists to a list of sets for easier merging
    fused_list = [set(sublist) for sublist in sublists]

    # Keep merging until no further changes
    merged = True
    while merged:
        merged = False
        new_fused_list = []
        
        while fused_list:
            # Take the first set
            current_set = fused_list.pop(0)
            
            # Check if it overlaps with any of the remaining sets
            overlapping = False
            for i, s in enumerate(fused_list):
                if current_set & s:  # If there is an intersection
                    # Merge the sets
                    current_set |= s
                    # Remove the merged set from the list
                    fused_list.pop(i)
                    overlapping = True
                    merged = True
                    break
            
            # Add the (potentially merged) set to the new fused list
            new_fused_list.append(current_set)
        
        # Update fused list with the merged sets
        fused_list = new_fused_list

    # Convert back to a list of lists
    final_list_of_lists = [list(s) for s in fused_list]
    return final_list_of_lists

def calculate_interface_loops(ref_structure, ref_int_resis_A):

    ### SELECTING INTERFACE LOOPS FROM CHAIN A ###
    ## Calculating interface loop residues
    # calculate secondary structure residues from chain A
    _, ref_beta_resis_A = calc_secondary_structure(ref_structure, type='all')
    # retrieving all chain A residues
    ref_chainA_resis = [r.getResnum() for r in ref_chainA]
    # selecting chain A loop residues by substracting beta residues from chain A residues (alpha residues are considered as loops)
    ref_loopA_resis = [r for r in ref_chainA_resis if r not in ref_beta_resis_A]

    ## Mapping interface residues to select interface loops and b-strands
    ref_A_loops = get_consecutive_numbers(ref_loopA_resis)
    ref_beta_strands_A = get_consecutive_numbers(ref_beta_resis_A)
    print(f"Loops in reference structure: {ref_A_loops}")
    print(f"Initial beta strands in reference structure: {ref_beta_strands_A}")

    # Generating position constraints to curate loop selection
    if pattern == 'cab': # nanobodies
        loop_a_core = [i for i in range(105,110)]
        false_beta_strands_A = [strand for strand in ref_beta_strands_A if len(strand) <= 3 or any(r in strand for r in loop_a_core)]
        print(f"Initial false beta strands: {false_beta_strands_A}")
        beta_strand_cores = [[57, 58, 59, 60], [69, 70, 71]]
        beta_strand_cores_flat = [r for strand in beta_strand_cores for r in strand]
    elif pattern == 'AB': # monobodies
        false_beta_strands_A = [strand for strand in ref_beta_strands_A if len(strand) <= 3]
        beta_strand_cores = [[59, 60, 61], [19,20]]
        beta_strand_cores_flat = [r for strand in beta_strand_cores for r in strand]

    recovered_strands = [strand for strand in ref_beta_strands_A if len(strand) <= 4 and any(r in strand for r in beta_strand_cores_flat)]
    if len(recovered_strands) > 0:
        print(f"Recovered true strands: {recovered_strands}")
        for recovered_strand in recovered_strands:
            if recovered_strand in false_beta_strands_A:
                false_beta_strands_A.remove(recovered_strand)
        print(f"Final false beta strands: {false_beta_strands_A}")

    false_beta_resis_A = [r for strand in false_beta_strands_A for r in strand]
    print(f"False beta residues: {false_beta_resis_A}")

    ## Reorganizing loop sections due to incorporation of false beta strand residues
    ref_loops_resis = [r for loop in ref_A_loops for r in loop]
    ref_loops_resis.extend(false_beta_resis_A)
    ref_loops_resis.sort()
    ref_A_loops = get_consecutive_numbers(ref_loops_resis)
    print(f"Initial chain A loops after including false beta residues: {ref_A_loops}")

    ## Processing loops
    # Excluding "loops" that are part of the n or c terminal
    nterm_resis, cterm_resis = [1,2], [len(ref_chainA_resis)-1, len(ref_chainA_resis)]
    ref_loops_final1 = [l for l in ref_A_loops if not (any(nr in l for nr in nterm_resis) or any(cr in l for cr in cterm_resis))]
    print(f"Chain A loop residues after removing N/C-term loops: {ref_loops_final1}")

    # We fuse together the loops that are at less or equal than 5 residues of distance (meaning, we are connecting again broken loops)
    loops_to_connect = []
    connected_loops = []
    length_cutoff = 5
    for i in range(len(ref_loops_final1)):
        for j in range(len(ref_loops_final1)):
            if i != j and i < j:  # Ensure we're not comparing the same loop or repeating the same comparison inverted
                length = ref_loops_final1[j][0] - ref_loops_final1[i][-1]
                extension_residues = [int(ref_loops_final1[i][-1])+k+1 for k in range(length-1)]
                print(f"Extension residues between {ref_loops_final1[i]} and {ref_loops_final1[j]}: {extension_residues}")
                if length-1 <= length_cutoff and not any(br in extension_residues for br in beta_strand_cores_flat): # if there are 5 residues between both
                    loops_to_connect.append((ref_loops_final1[i], ref_loops_final1[j]))
                    connected_loop = ref_loops_final1[i]+extension_residues+ref_loops_final1[j]
                    connected_loop.sort()
                    connected_loops.append(connected_loop)

    print(f"Loops to connect: {loops_to_connect}")
    print(f"Connected loops: {connected_loops}")

    # Removing from the original list of loops the ones that are being connected
    loops_to_connect_flat = []
    for tup in loops_to_connect:
        for l in tup:
            if l not in loops_to_connect_flat:
                loops_to_connect_flat.append(l)
    print(f"Loops to be connected (flat non-redundant version): {loops_to_connect_flat}")

    ref_loops_final2 = [loop for loop in ref_loops_final1 if loop not in loops_to_connect_flat]
    print(f"Loop list after removing connected loops: {ref_loops_final2}")

    if len(connected_loops) > 0: # we had short loops
        if len(connected_loops) > 1: # there are enough connected loops for this step
            # We do an extra step to fuse together any overlapping connected loop
            connected_loops_final = fuse_overlapping_sublists(connected_loops)
            print(f"Final list of connected short loops after checking for overlappings: {connected_loops_final}")

            # Adding connected loops back into the loops list
            ref_loops_final2.extend(connected_loops_final)
        else:
            print("There is only 1 connected loop.")
            ref_loops_final2.extend(connected_loops)
    else: # we had short loops but they are not connected
        print("Short loops cannot be connected.")

    # finally keeping only loops that contain at least one interface residue number
    ref_intA_loops_final2 = [l for l in ref_loops_final2 if any(r in l for r in ref_int_resis_A)]
    print(f"Final interface loops in reference structure: {ref_intA_loops_final2}")

    ## Defining individual prody objects for each interface loop
    ref_intA_loop_objs = [ref_TOTAL.select(f"protein and chain A and resnum {' '.join(map(str, loop))}").ca for loop in ref_intA_loops_final2]

    print(f"Prody caught chain A interface loops in reference structure:")
    for loop in ref_intA_loop_objs:
        cur_loop_resis_prody = [r.getResnum() for r in loop]
        print(f"Current loop: {cur_loop_resis_prody}")

    ## ALL CHAIN A INTERFACE LOOPS    
    # Flattening final list of loops
    ref_intA_loops_resis_final = [r for loop in ref_intA_loops_final2 for r in loop]
    # defining a global object incorporating all interface loop residues
    if len(ref_intA_loop_objs) > 1:
        ref_intA_all_loops_obj = sum(ref_intA_loop_objs[1:], ref_intA_loop_objs[0])
    else:
        ref_intA_all_loops_obj = ref_intA_loop_objs[0]

    all_loop_resis_prody = [r.getResnum() for r in ref_intA_all_loops_obj]
    print(f"Prody caught chain A ALL interface loop residues in reference structure: {all_loop_resis_prody}")

    return ref_intA_loops_final2, ref_intA_loop_objs, all_loop_resis_prody, ref_intA_all_loops_obj, ref_intA_loops_resis_final

# data analysis functions
def exponential_normalize_value(value, a):
    normalized_value = 0 + math.exp(-a * value)
    return normalized_value

def extract_keys(dict_str, flag):
    # Replace `;` with commas to make it a valid dictionary format
    dict_str = dict_str.replace(';', ',')
    
    # Convert the string to a dictionary using ast.literal_eval
    dict_obj = ast.literal_eval(dict_str)
    
    # Return the keys as a single string separated by ;
    if flag == 'keys':
        return ';'.join(dict_obj.keys())
    elif 'values' in flag:
        residue_list = ';'.join(dict_obj.values())
        residue_list = [x for x in residue_list.split(';')]
        if flag == 'values length':
            length_list = [len(x.split(',')) for x in residue_list]
            return ';'.join(map(str, length_list))
        else:
            return ';'.join(map(str, residue_list))

def get_filtering_table(path, pred_type, top3=True, calc_loops=False):
    
    # Reading data frame
    df = pd.read_table(path, sep=' ')

    # Transforming pLDDT to 1 scale
    df['pLDDT']=df['pLDDT']/100
    df['pLDDT_A_int']=df['pLDDT_A_int']/100

    # Normalizing RMSD_AR to 0:1 scale using an exponential function
    df['RMSD_AR_norm'] = df['RMSD_AR'].apply(lambda x: exponential_normalize_value(x, 0.05))
    df['composite_score_int']= df['pLDDT_A_int']*df['iPTMS']*df['RMSD_AR_norm'] # CS using normalized RMSD_AR and pLDDT of the interface

    if top3:
        # Creating top3 data frame to only take into account ranks 1-3
        df_means = df[~df['rank'].isin(['rank_004', 'rank_005'])]
    else:
        df_means = df.copy()

    if calc_loops:
        # Calculate RMSD_AA/len column
        df_means.loc[:,'RMSD_AA_only_ind_loops'] = df_means.loc[:,'RMSD_AA_int_indiv_loops'].apply(lambda x: extract_keys(x, 'keys'))
        df_means.loc[:,'ind_loops'] = df_means.loc[:,'RMSD_AA_int_indiv_loops'].apply(lambda x: extract_keys(x, 'values'))
        df_means.loc[:,'length_ind_loops'] = df_means.loc[:,'RMSD_AA_int_indiv_loops'].apply(lambda x: extract_keys(x, 'values length'))

        # Extract set of pdbs
        pdbs = list(set(df_means['PDB']))
        temp_df2_means_list = []
        loop_pivot_final_df_list = []
        for i, pdb in enumerate(pdbs):
            
            # Creating a temp df for each PDB with our target columns
            temp_df = df_means[df_means['PDB']==pdb][['PDB', 'RMSD_AA_only_ind_loops', 'ind_loops', 'length_ind_loops']]
            
            # splitting ; separated values in new columns
            RMSD_AA_split = temp_df['RMSD_AA_only_ind_loops'].str.split(';', expand=True)
            loop_split = temp_df['ind_loops'].str.split(';', expand=True)
            looplen_split = temp_df['length_ind_loops'].str.split(';', expand=True)
            
            # Convert to numeric
            RMSD_AA_split = RMSD_AA_split.apply(pd.to_numeric)
            looplen_split = looplen_split.apply(pd.to_numeric)

            # Dynamically generate column names based on the number of columns formed
            num_cols = RMSD_AA_split.shape[1]  # Get the number of columns
            new_RMSD_AA_col_names = [f'RMSD_AA_loop_{i+1}' for i in range(num_cols)]
            new_looplen_col_names = [f'Length_loop_{i+1}' for i in range(num_cols)]
            new_loop_col_names = [f'Loop_{i+1}' for i in range(num_cols)]
            
            # Rename the columns
            RMSD_AA_split.columns = new_RMSD_AA_col_names
            looplen_split.columns = new_looplen_col_names
            loop_split.columns = new_loop_col_names
            # processing loop_split differently
            loop_split_pro = loop_split.iloc[0, :][new_loop_col_names].to_frame().T.reset_index().drop(columns='index')
        
            # Add the split columns in a second temp df
            temp_df2 = pd.concat([temp_df['PDB'], RMSD_AA_split], axis=1)
            temp_df2 = pd.concat([temp_df2, looplen_split], axis=1)
            
            # Calculating the mean for each numeric column
            temp_numeric_columns = temp_df2.select_dtypes(include='number').columns
            temp_df2_means = temp_df2.groupby('PDB')[temp_numeric_columns].mean().reset_index()
            temp_df2_means = pd.concat([temp_df2_means, loop_split_pro], axis=1)
            temp_df2_means_list.append(temp_df2_means)

            # Pivot the df
            temp_pivot_list = []
            for i in range(num_cols):
                loop_num = i+1
                
                RMSD_AA_col_name = f"RMSD_AA_loop_{loop_num}"
                loop_col_name = f"Loop_{loop_num}"
                looplen_col_name = f"Length_loop_{loop_num}"
                
                RMSD_AA_value = temp_df2_means.loc[0, RMSD_AA_col_name]
                loop_value = temp_df2_means.loc[0, loop_col_name]
                looplen_value = temp_df2_means.loc[0, looplen_col_name]
                
                new_row_name = pdb+'_loop_'+str(loop_num)
                
                temp_pivot_list.append([pdb, new_row_name, RMSD_AA_value, looplen_value, loop_value])
                
            temp_pivot_df = pd.DataFrame(temp_pivot_list)
            temp_pivot_df.columns = ['PDB', 'PDB_loop', 'RMSD_AA', 'Loop Length', 'Loop_resis']
            loop_pivot_final_df_list.append(temp_pivot_df)
            
        loop_pivot_final_df = pd.concat(loop_pivot_final_df_list, axis=0)
        loop_pivot_final_df['pred_type'] = pred_type
    
    # Calculating mean values accross models
    numeric_columns = df_means.select_dtypes(include='number').columns
    df_means = df_means.groupby('PDB')[numeric_columns].mean().reset_index()

    # Calculating std
    df_sd = df.groupby(['PDB'], as_index=False)[numeric_columns].std()
    df_means['std RMSD_AR'] = df_sd['RMSD_AR']

    # Adding columns to df
    df_means['pred_type'] = pred_type

    # Filtering out rows with 0 values in many columns
    empty_subset_mask = (df_means['composite_score_int'] == 0) & (df_means['TM-score'] == 0)
    if True in list(empty_subset_mask):
        # Print the warning message
        print(f"Dataset '{pred_type}' has been cut due to empty values for some elements. Please review source data.")
        
        # Retrieve the 'PDB' column values where the condition is True
        pdb_values = df_means.loc[empty_subset_mask, 'PDB'].tolist()
    
        # Print or use the pdb_values list as needed
        print("PDB values corresponding to the condition:", pdb_values)
        
        # Optionally, filter the DataFrame to exclude these rows
        df_means = df_means[~empty_subset_mask]
    
    if calc_loops:
        return df_means, loop_pivot_final_df
    else:
        return df_means

def process_loops(all_loop_df, pattern):

    if pattern == 'cab':
        loop1 = '26|27|28|29'
        loop2 = '52|53|54|55'
        loop3 = '100|101|102|103|111'

        labels = ['CDRH1', 'CDRH2', 'CDRH3']

    elif pattern == 'AB':
        loop1 = '26|27|28|29'
        loop2 = '54|55|56|57'
        loop3 = '80|81|82|83'    

        labels = ['BC loop', 'DE loop', 'FG loop']    

    # define conditions to select each cdr
    conditions = [
        all_loop_df['Loop_resis'].str.contains(loop1), # CDRH1 or BC loop
        all_loop_df['Loop_resis'].str.contains(loop2), # CDRH2 or DE loop
        all_loop_df['Loop_resis'].str.contains(loop3), # CDRH3 or FG loop
    ]

    # Add a new column with the loop labels, defaulting to 'no_loop' if no conditions are met
    all_loop_df['cdr_label'] = np.select(conditions, labels, default='no_loop')

    # Filtering to only retain CDR(-like) loops
    all_loop_df = all_loop_df[all_loop_df['cdr_label'].isin(labels)]

    return all_loop_df

# END OF FUNCTIONS # 

# START OF SCRIPT #

#Give the path where your results are
path = sys.argv[1]
pattern = sys.argv[2]

calc_loops = False
if pattern in ['cab', 'AB']:
    calc_loops = True

#Give the output path where you want the file to be saved
results_file = 'colabfold_multimer_analysis_results.txt'
outfile = open(os.path.join(path, results_file), 'w')

if calc_loops:
    outfile.write(f"PDB rank model pLDDT pLDDT_A_int TM-score RMSD RMSD_AR RMSD_AA_int_loops RMSD_AA_int_indiv_loops PTMS iPTMS PAE_interchain\n")
else:
    outfile.write(f"PDB rank model pLDDT pLDDT_A_int TM-score RMSD RMSD_AR PTMS iPTMS PAE_interchain\n")

os.chdir(path)
abs_path = os.getcwd()

# Get a list of all directories in the current directory
all_dirs = [d for d in os.listdir() if os.path.isdir(d)]
# Filter the directories to include only those with a specified patterns
# - you should have a folder for each protein complex, containing the AF2 predictions
pdb_dirs = [d for d in all_dirs if pattern in d]

print(f'A total of {len(pdb_dirs)} target directories has been found. Now iterating on them...')

loop_beta_dict = {}
for d in tqdm(pdb_dirs):

    print(f'\nCurrent PDB: {d}\n')

    print(d)
    os.chdir(os.path.join(path, d))

    ref_structures = [pdb for pdb in glob.glob('*.pdb') 
                  if 'alphafold2' not in pdb 
                  and not pdb.endswith('A.pdb') 
                  and not pdb.endswith('rec.pdb')] # UPDATE THIS TO YOUR SPECIFIC FILESYSTEM

    if len(ref_structures)!=1:
        print(f'There are {len(ref_structures)} PDBs matching the reference model pattern. Stopping script. Review path.')
        sys.exit()
    else:
        ref_structure=ref_structures[0]

    # ORIGINAL OR REFERENCE SEQUENCE
    ref_TOTAL = parsePDB(ref_structure)

    chains = sorted(list(set(chain.getChid() for chain in ref_TOTAL.iterChains())))
    if len(chains) > 2:
        third_chain = chains[2]
    else:
        third_chain = None

    ref_chainA = ref_TOTAL.select("protein and chain A").ca
    ref_chainB = ref_TOTAL.select("protein and chain B").ca
    
    if third_chain is not None:
        ref_chainCDE = ref_TOTAL.select(f"protein and chain {third_chain}").ca
        ref_receptor = ref_chainB + ref_chainCDE
        print(len(ref_chainCDE.getResnums()))
    else:
        ref_receptor = ref_chainB   

    # Calculating interface residues from chain A
    ref_int_resis_A = calc_interface_pyrosetta(ref_structure, third_chain=third_chain)[0]
    print(f"Interface residues from chain A (crystal): {ref_int_resis_A}")
    ref_chainA_interface = ref_TOTAL.select(f"protein and chain A and resnum {' '.join(map(str, ref_int_resis_A))}").ca

    # Calculating interface loops from chain A
    if calc_loops:
        ref_intA_loops_final2, ref_intA_loop_objs, all_loop_resis_prody, ref_intA_all_loops_obj, ref_intA_loops_resis_final = calculate_interface_loops(ref_structure, ref_int_resis_A)

    models = glob.glob('*_relaxed_rank_*.pdb')
    if len(models) != 5:
        print(f'\nThere are {len(models)} models for protein {d}, when there should be 5. Skipping metric calculation.\n')
        continue
    json_files = glob.glob('*_scores_*.json')
    if len(json_files) != 5:
        print(f'\nThere are {len(json_files)} json files for protein {d}, when there should be 5. Skipping metric calculation.\n')
        continue

    # Creating regexp patterns
    model_pattern = re.compile(r'model_\d+\_seed_000')
    rank_pattern = re.compile(r'_relaxed_rank_00\d+')

    # READING LOG FILE AND CALCULATING METRICS
    with open ('log.txt','r') as log_file:
        for line in log_file:
            if ' rank_' in line:
                
                model_log_string = str(model_pattern.search(line).group())
                model_log = '_'.join(model_log_string.split('_')[:2])
                
                for model in models:

                    if not model == [] and (f'model_{model[-14]}' == model_log):
                        print('####')
                        print(model)
                        print('####')

                        rank_string = str(rank_pattern.search(model).group())
                        rank_short = '_'.join(rank_string.split('_')[-2:])
                        model_string = str(model_pattern.search(model).group())
                        model_short = '_'.join(model_string.split('_')[:2])

                        # MODEL SEQUENCE
                        mod_TOTAL = parsePDB(model)
                        mod_chainA = mod_TOTAL.select("protein and chain A").ca
                        mod_chainB = mod_TOTAL.select("protein and chain B").ca

                        # if third_chain is not None:
                        #     mod_chainC = mod_TOTAL.select(f"protein and chain {third_chain}").ca
                        #     mod_receptor = mod_chainB + mod_chainC
                        # else:
                        #     mod_receptor = mod_chainB

                        # Calculating interface residues from chain A
                        # mod_chainA_interface = mod_TOTAL.select(f"protein and chain A and resnum {' '.join(map(str, ref_int_resis_A))}").ca

                        if third_chain is not None:
                            mod_chainC = mod_TOTAL.select(f"protein and chain {third_chain}").ca

                            # PDBS WHERE WE ADDED A 'G' AT THE END OF CHAIN C
                            if any(name in d for name in ['7n0i', '5a43', '5kbn']):

                                resnums = mod_chainC.getResnums()

                                # Print the three-letter code of the last residue
                                last_residue = mod_TOTAL.select(f"protein and chain {third_chain} and resnum {resnums[-1]}").getResnames()
                                print(f"Last residue in chain {third_chain} of PDB {d}: {last_residue[0]}")
                                # Exclude the last residue
                                last_resnum = resnums[-1]
                                mod_chainC = mod_TOTAL.select(f"protein and chain {third_chain} and resnum < {last_resnum}").ca
                                print(len(mod_chainC.getResnums()))

                            # PDBS WHERE WE ADDED A 'G' AT THE END OF CHAIN B
                            elif any(name in d for name in ['5nkq', '6bqo']):
                                # Get the residue numbers for the selection
                                resnums = mod_chainB.getResnums()

                                # Print the three-letter code of the last residue
                                last_residue = mod_TOTAL.select(f"protein and chain B and resnum {resnums[-1]}").getResnames()
                                print(f"Last residue in chain B of PDB {d}: {last_residue[0]}")
                                # Exclude the last residue
                                last_resnum = resnums[-1]
                                mod_chainB = mod_TOTAL.select(f"protein and chain B and resnum < {last_resnum}").ca
                                print(len(mod_chainB.getResnums()))

                            mod_receptor = mod_chainB + mod_chainC
                            mod_TOTAL = mod_chainA + mod_receptor
                        else:
                            mod_receptor = mod_chainB

                        # Calculating interface residues from chain A
                        no_interface_flag = False # we assume we will find an interface between A and the receptor
                        interface_residues = calc_interface_pyrosetta(pdbfile=model, third_chain=third_chain)
                        interface_residues_A, interface_residues_B = interface_residues[0], interface_residues[1]
                        if len(interface_residues_A) == 0:
                            print("NO INTERFACE DETECTED WITH CHAIN A!!!")
                            no_interface_flag = True
                        else:
                            mod_chainA_interface = mod_TOTAL.select(f"protein and chain A and resnum {' '.join(map(str, interface_residues_A))}").ca
                            int_selection_string_A = f"(chain A and resnum {' '.join(map(str, interface_residues_A))})"
                            int_selection_string_B = f"(chain B and resnum {' '.join(map(str, interface_residues_B))})"
                            if third_chain is not None:

                                interface_residues_C = interface_residues[2]
                                int_selection_string_C = f"(chain C and resnum {' '.join(map(str, interface_residues_C))})"

                                if len(interface_residues_B) == 0:
                                    print("NO INTERFACE DETECTED WITH CHAIN B!!!")
                                    mod_total_interface = mod_TOTAL.select(f"protein and ({int_selection_string_A} or {int_selection_string_C})").ca
                                elif len(interface_residues_C) == 0:
                                    print("NO INTERFACE DETECTED WITH CHAIN C!!!")
                                    mod_total_interface = mod_TOTAL.select(f"protein and ({int_selection_string_A} or {int_selection_string_B})").ca
                                else:
                                    mod_total_interface = mod_TOTAL.select(f"protein and ({int_selection_string_A} or {int_selection_string_B} or {int_selection_string_C})").ca

                            else:
                                mod_total_interface = mod_TOTAL.select(f"protein and ({int_selection_string_A} or {int_selection_string_B})").ca                            

                        # RMSD total, RMSD AA, RMSD AB
                        ## total
                        t_total = calcTransformation(mod_TOTAL.ca,ref_TOTAL.ca)
                        sup_total = applyTransformation(t_total, mod_TOTAL)
                        sup_ca_total = sup_total.select("protein").ca
                        rmsd_prody = calcRMSD(ref_TOTAL.ca, sup_ca_total)

                        ## INTERCHAIN RMSD - chain A (monobody) rmsd while aligning in target (B/BC/BD/BE)
                        if third_chain is not None:
                            t_ar = calcTransformation(mod_receptor, ref_receptor)
                        else:
                            t_ar = calcTransformation(mod_chainB, ref_chainB) # calculate the transformation needed to align chain B of the model to chain B of reference pdb
                        sup_ar = applyTransformation(t_ar, mod_TOTAL) # apply this transformation to the model pdb
                        # Whole chain A RMSD_AR
                        sup_chainA = sup_ar.select("protein and chain A").ca # select ca atoms of chain A of the model (already aligned)
                        rmsd_AR = calcRMSD(ref_chainA, sup_chainA) # calculate rmsd between the reference chain A and aligned model chain A

                        ### INTRACHAIN RMSD ### - chain A (monobody) rmsd while aligning in A (monobody)
                        if calc_loops:
                            ## Interface A RMSD
                            # retrieving model prody object from crystal interface residues
                            mod_intA_all_loops_obj = mod_TOTAL.select(f"protein and chain A and resnum {' '.join(map(str, ref_intA_loops_resis_final))}").ca
                            all_loop_resis_prody_mod = [r.getResnum() for r in mod_intA_all_loops_obj]
                            print(f"Prody caught chain A ALL interface loop residues in model structure: {all_loop_resis_prody_mod}")
                            # calculating RMSD
                            t_a_int = calcTransformation(mod_intA_all_loops_obj, ref_intA_all_loops_obj)
                            sup_a_int = applyTransformation(t_a_int, mod_TOTAL)
                            sup_a_chainA_int = sup_a_int.select(f"protein and chain A and resnum {' '.join(map(str, ref_intA_loops_resis_final))}").ca
                            rmsd_AA_int_loops = calcRMSD(ref_intA_all_loops_obj, sup_a_chainA_int)

                            ## Loop by loop interface A RMSD ##
                            print(f"Prody caught chain A interface loops in current model structure:")
                            rmsd_AA_loop_dict = {}
                            for ref_intA_loop_obj, ref_intA_loop_resis in zip(ref_intA_loop_objs, ref_intA_loops_final2):
                                # creating objects for loops in model
                                mod_intA_loop_obj = mod_TOTAL.select(f"protein and chain A and resnum {' '.join(map(str, ref_intA_loop_resis))}").ca
                                cur_loop_resis_mod = [r.getResnum() for r in mod_intA_loop_obj]
                                print(f"Current loop: {cur_loop_resis_mod}")
                                # calculating RMSD
                                t_a_loop = calcTransformation(mod_intA_loop_obj, ref_intA_loop_obj)
                                sup_a_loop = applyTransformation(t_a_loop, mod_TOTAL)
                                sup_a_chainA_loop = sup_a_loop.select(f"protein and chain A and resnum {' '.join(map(str, ref_intA_loop_resis))}").ca
                                rmsd_AA_loop = calcRMSD(ref_intA_loop_obj, sup_a_chainA_loop)
                                # adding to dict                  
                                rmsd_AA_loop_dict[rmsd_AA_loop]=','.join(map(str, ref_intA_loop_resis))
                            rmsd_AA_loop_dictstr = json.dumps(rmsd_AA_loop_dict, separators=(';', ':'))
                            
                        #PLDDT chain A
                        plddt_A = np.mean(mod_chainA.getBetas())
                        plddt_A_int = np.mean(mod_chainA_interface.getBetas())

                        # TMSCORE and TOTAL RMSD
                        tm_score, rmsd_mmalign = tmalign_result(ref_structure,model,norm=False)

                        # PTMSCORE
                        ptms = float(line.split()[-2].split('=')[1])

                        #IPTMSCORE
                        iptms = float(line.split()[-1].split('=')[1])

                        #PAE
                        pdb_filepath = os.path.join(os.getcwd(), ref_structure)
                        model_filepath = os.path.join(os.getcwd(), model)
                        json_file = [j for j in json_files if model_short in j][0]
                        PAE_interchain=interface_PAE(json_file, model, third_chain)
                        
                        if calc_loops:
                            outfile.write(f"{d} {rank_short} {model_short} {plddt_A} {plddt_A_int} {tm_score} {rmsd_prody} {rmsd_AR} {rmsd_AA_int_loops} {rmsd_AA_loop_dictstr} {ptms} {iptms} {PAE_interchain}\n")
                        else:
                            outfile.write(f"{d} {rank_short} {model_short} {plddt_A} {plddt_A_int} {tm_score} {rmsd_prody} {rmsd_AR} {ptms} {iptms} {PAE_interchain}")
outfile.close()

# PROCESSING DATA #
os.chdir(abs_path)
afm_data = get_filtering_table(os.path.join(abs_path, results_file), 'AF2', calc_loops=calc_loops)
afm_proteinmpnn_data = get_filtering_table(os.path.join(abs_path, results_file), 'MPNN100', calc_loops=calc_loops)

if calc_loops:
    # AF2 metrics
    afm_metrics_top3 = afm_data[0]
    afm_proteinmpnn_metrics_top3_df = afm_proteinmpnn_data[0]

    # CDR metrics
    afm_int_loop_df = afm_data[1]
    afm_proteinmpnn_int_loop_df = afm_proteinmpnn_data[1]
    # process CDR metrics
    af2_int_loop_df = process_loops(afm_int_loop_df, pattern=pattern)
    mpnn100_int_loop_df = process_loops(afm_proteinmpnn_int_loop_df, pattern=pattern)
    # saving data in csv files
    int_loop_df = pd.concat([af2_int_loop_df, mpnn100_int_loop_df])
    int_loop_df.to_csv('int_loop_df.csv', index=False)
else:
    afm_metrics_top3 = afm_data
    afm_proteinmpnn_metrics_top3_df = afm_proteinmpnn_data

# saving data in csv files
metrics_df = pd.concat([afm_metrics_top3, afm_proteinmpnn_metrics_top3_df])
metrics_df_final = metrics_df[['PDB', 'pLDDT', 'pLDDT_A_int', 'TM-score', 'RMSD', 'RMSD_AR', 'RMSD_AR_norm', 'RMSD_AA_int_loops', 'PTMS', 'iPTMS', 'PAE_interchain', 'composite_score_int',
                                        'std RMSD_AR', 'pred_type']]
metrics_df_final.columns = ['PDB', 'binder pLDDT', 'binder interface pLDDT', 'TM-score', 'RMSD', 'RMSD_AB(C)', 'RMSD_AB(C)_norm', 'RMSD_int_loops', 'PTM', 'ipTM', 'interchain_PAE', 'Composite Score',
                                  'std RMSD_AB(C)', 'pred_type']
metrics_df_final.to_csv('metrics_df.csv', index=False)

# END OF SCRIPT #