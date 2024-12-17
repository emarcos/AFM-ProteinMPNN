#!/usr/bin/env python3
from prody import *
from argparse import *
from pyrosetta import init, pose_from_pdb
#from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.core.select.residue_selector import *
from argparse import ArgumentParser

def process_mpnn_output(input_file, output_file, chains, verbose=False):

    """
    Process the ProteinMPNN output file to 
    update sequence titles.
    """

    if verbose:
        print(f"Processing ProteinMPNN paired MSA file: {os.path.basename(input_file)} -> {os.path.basename(output_file)}")

    nseq = 0
    with open(input_file, 'r') as f, open(output_file, 'w') as o:
        for line_number, line in enumerate(f):
            if line_number < 1:
                o.write(line)
                continue
            if not line.startswith('>'):
                line = line.replace('/', '')
            else:
                if chains == 3:
                    line = f'>seqA{nseq:03d}\tseqB{nseq:03d}\tseqC{nseq:03d}\n'
                elif chains == 2:
                    line = f'>seqA{nseq:03d} \t seqB{nseq:03d}\n'
                nseq += 1
            o.write(line)

    if verbose:
        print(f"ProteinMPNN paired MSA processing completed.")

def write_a3m_file(fastaA, fastaB, processed_mpnnfile, multimer_msa_a3m, output_a3m_file, fastaC=None, verbose=False):

    """
    Write the a3m file.
    """

    if verbose:
        print(f"Generating a3m file: {output_a3m_file}")

    with open(output_a3m_file, 'w') as fileout:

        # 1. Write header and input sequence into a3m file
        if fastaC:
            header = f'#{len(fastaA)},{len(fastaB)},{len(fastaC)}\t1,1,1\n'
            input_seq_title = '>101\t102\t103\n'
            input_seq = f'{fastaA}{fastaB}{fastaC}\n'
        else:
            header = f'#{len(fastaA)},{len(fastaB)}\t1,1\n'
            input_seq_title = '>101\t102\n'
            input_seq = f'{fastaA}{fastaB}\n'

        fileout.write(header)
        fileout.write(input_seq_title)
        fileout.write(input_seq)

        # 2. Appending paired afm-proteinmpnn MSA to a3m
        with open(processed_mpnnfile, 'r') as mf:
            for line_number, line in enumerate(mf):
                if line_number >= 2:
                    fileout.write(line)

        # 3. Appending unpaired MSA a3m file content (monobody and receptor)
        flag = False
        with open(multimer_msa_a3m, 'r') as f:
            next(f); next(f)  # skipping first two lines
            for line in f:
                line = line.strip('\n')
                if fastaC and line in ['>101', '>102', '>103']:
                    fileout.write(f'{line}\n')
                    flag = True
                elif not fastaC and line in ['>101', '>102']:
                    fileout.write(f'{line}\n')
                    flag = True
                elif flag:
                    fileout.write(f'{line}\n')

    if verbose:
        print(f"a3m file generation completed, saved in: {output_a3m_file}")

if __name__ == "__main__":

    parser = ArgumentParser(description=f"""Generate a3m file for AFM-ProteinMPNN prediction.
Example of usage: ./a3m_generator.py -pdb <input_pdb> -fasta <input_fasta> -afm_a3m <input_afm_a3m>""")
    parser.add_argument("-pdb", "--pdbfile", required=True, help="Path to the input PDB file.")
    parser.add_argument("-fasta", "--fastafile", required=True, help="Path to the input FASTA file.")
    parser.add_argument("-afm_a3m", "--multimer_msa_a3m", required=True, help="Path to the input default AlphaFold-Multimer MSA a3m file.")
    parser.add_argument("-o", "--output", help="Path to the desired output file. Defaults to current directory.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose mode.")

    args = parser.parse_args()

    if args.verbose:
        print(f"\nRunning a3m_generator.py ...\n")

    # 1. Analyzing PDB
    pdbfilename = '.'.join(os.path.basename(args.pdbfile).split('.')[:-1])
    structure = parsePDB(args.pdbfile)
    input_chains = sorted(list(set(chain.getChid() for chain in structure.iterChains())))

    if len(input_chains) not in  [2,3]:
        sys.exit(f"Input PDB {args.pdbfile} contains {len(input_chains)} and therefore cannot be processed. Please review and try again.")

    # 2. Analyzing FASTA
    with open(args.fastafile, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                fasta_all = line

    fastaA = fasta_all.split(':')[0]
    fastaB = fasta_all.split(':')[1]
    fastaC = None

    num_fasta_seqs = len(fasta_all.split(':'))
    if num_fasta_seqs == 3:
        fastaC = fasta_all.split(':')[2]

    # 3. Quality check before proceeding
    if num_fasta_seqs not in [2,3] or len(input_chains) not in [2,3] or num_fasta_seqs != len(input_chains):
        sys.exit(f"""Input FASTA and PDB files must have the same number of chains (2 or 3). 
        - Your FASTA ({args.fastafile}) contains {num_fasta_seqs} chains.
        - Your PDB ({args.pdbfile}) contains {len(input_chains)} chains.
        Please review and try again.""")       

    # 4. Processing mpnn output
    mpnnfile = os.path.join(os.getcwd(), f'mpnn/seqs/{pdbfilename}.fa')
    processed_mpnnfile = os.path.join(os.getcwd(), f'mpnn/seqs/processed_{pdbfilename}.fa')
    process_mpnn_output(mpnnfile, processed_mpnnfile, chains=num_fasta_seqs, verbose=args.verbose)

    # 5. Generating a3m file
    output_a3m_file = os.path.join(args.output, f"{pdbfilename}_mpnn_custom_msa.a3m") if args.output else os.path.join(os.getcwd(), f"{pdbfilename}_mpnn_custom_msa.a3m")
    write_a3m_file(fastaA, fastaB, processed_mpnnfile, args.multimer_msa_a3m, output_a3m_file, fastaC=fastaC, verbose=args.verbose)

# end of script