#!/usr/bin/env python3
from prody import *
from argparse import *
from pyrosetta import init, pose_from_pdb
#from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.core.select.residue_selector import *
from argparse import ArgumentParser

def process_mpnn_output(input_file, output_file, chain, fastaA, fastaB, verbose=False):

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
                if chain == 'A':
                    line = line.split('/')[0]+fastaB+line.split('/')[1]
                if chain == 'B':
                    line = fastaA+line.split('/')[0]+line.split('/')[1]
            else:
                if chain == 'A':
                    line = f'>seqA_{nseq:03d} \t seqB_{nseq:03d} \t seqC_{nseq:03d}\n'
                if chain == 'B':
                    line = f'>seqA_{(nseq+100):03d} \t seqB_{(nseq+100):03d} \t seqC_{(nseq+100):03d}\n'
                nseq += 1
            o.write(line)

    if verbose:
        print(f"ProteinMPNN paired MSA processing completed.")


def write_a3m_file(fastaA, fastaB, fastaC, processed_mpnnfile_AC, processed_mpnnfile_BC, receptor_msa, output_a3m_file, verbose=False):

    """
    Write the a3m file.
    """

    if verbose:
        print(f"Generating a3m file: {output_a3m_file}")

    with open(output_a3m_file, 'w') as fileout:

        # 1. Write header and input sequence into a3m file
        fileout.write(f'#{len(fastaA)},{len(fastaB)},{len(fastaC)}\t1,1,1\n')
        fileout.write('>101\t102\t103\n')
        fileout.write(fastaA+fastaB+fastaC)

        # 2. Appending paired afm-proteinmpnn MSA to a3m
        # complex_AC
        with open(processed_mpnnfile_AC, 'r') as mf_AC:
            for line_number, line in enumerate(mf_AC):
                if line_number >= 2:
                    fileout.write(line)
        # complex BC
        with open(processed_mpnnfile_BC, 'r') as mf_BC:
            for l_n, l in enumerate(mf_BC):
                if l_n >= 2:
                    fileout.write(l)

        # 3. Appending individual filename seqs for A and B (binders)
        #Binder A
        fileout.write('>101\n')
        fileout.write(fastaA+'-'*(len(fastaB)+len(fastaC)))
        fileout.write(os.linesep)
        #Binder B
        fileout.write('>102\n')
        fileout.write('-'*len(fastaA)+fastaB+'-'*len(fastaC))
        fileout.write(os.linesep)

        # 4. Appending unpaired MSA a3m file content (receptor)
        with open(receptor_msa, 'r') as f:
            for line in f:
                line = line.strip('\n')
                if line.startswith('>'):
                    fileout.write(f'{line}\n')
                else:
                    blank_spaces = '-'*(len(fastaA)+len(fastaB))
                    fileout.write(f'{blank_spaces+line}\n')

    if verbose:
        print(f"a3m file generation completed, saved in: {output_a3m_file}")

if __name__ == "__main__":

    parser = ArgumentParser(description=f"""Generate a3m file for AFM-ProteinMPNN prediction.
Example of usage: ./a3m_generator.py -pdbA <input_pdbA> -pdbB <input_pdbB> -fasta <input_fasta> -receptor_a3m <input_receptor_a3m>""")
    parser.add_argument("-pdbA", "--pdbfileA", required=True, help="Path to the input PDB file.")
    parser.add_argument("-pdbB", "--pdbfileB", required=True)
    parser.add_argument("-fasta", "--fastafile", required=True, help="Path to the input FASTA file.")
    parser.add_argument("-receptor_a3m", "--receptor_msa", required=True, help="Path to the input default AlphaFold-Multimer MSA a3m file.")
    parser.add_argument("-o", "--output", help="Path to the desired output file. Defaults to current directory.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose mode.")

    args = parser.parse_args()

    if args.verbose:
        print(f"\nRunning a3m_generator.py ...\n")

    # 1. Analyzing PDB
    pdbfilename_A = args.fastafile.split('_vs_')[0].split('_and_')[1]
    pdbfilename_B = args.fastafile.split('_vs_')[1][:-6]
    structureA = parsePDB(args.pdbfileA)
    input_chainsA = sorted(list(set(chain.getChid() for chain in structureA.iterChains())))

    if len(input_chainsA) not in  [2,3]:
        sys.exit(f"Input PDB {args.pdbfileA} contains {len(input_chainsA)} and therefore cannot be processed. Please review and try again.")

    structureB = parsePDB(args.pdbfileB)
    input_chainsB = sorted(list(set(chain.getChid() for chain in structureB.iterChains())))

    if len(input_chainsB) not in  [2,3]:
        sys.exit(f"Input PDB {args.pdbfileB} contains {len(input_chainsB)} and therefore cannot be processed. Please review and try again.")

    # 2. Analyzing FASTA
    with open(args.fastafile, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                fasta_all = line

    fastaA = fasta_all.split(':')[0]
    fastaB = fasta_all.split(':')[1]
    fastaC = fasta_all.split(':')[2]

    # 3. Processing mpnn output
    mpnnfile_AC = os.path.join(os.getcwd(), f'{pdbfilename_A}_mpnn/seqs/{pdbfilename_A}.fa')
    mpnnfile_BC = os.path.join(os.getcwd(), f'{pdbfilename_B}_mpnn/seqs/{pdbfilename_B}.fa')
    processed_mpnnfile_AC = os.path.join(os.getcwd(), f'{pdbfilename_A}_mpnn/seqs/processed_{pdbfilename_A}.fa')
    processed_mpnnfile_BC = os.path.join(os.getcwd(), f'{pdbfilename_B}_mpnn/seqs/processed_{pdbfilename_B}.fa')
    process_mpnn_output(mpnnfile_AC, processed_mpnnfile_AC, 'A', fastaA, fastaB, verbose=args.verbose)
    process_mpnn_output(mpnnfile_BC, processed_mpnnfile_BC, 'B', fastaA, fastaB, verbose=args.verbose)

    # 4. Generating a3m file
    output_a3m_file = os.path.join(args.output, f"{pdbfilename_A}-{pdbfilename_B}_mpnn_custom_msa.a3m") if args.output else os.path.join(os.getcwd(), f"{pdbfilename_A}-{pdbfilename_B}_mpnn_custom_msa.a3m")
    write_a3m_file(fastaA, fastaB, fastaC, processed_mpnnfile_AC, processed_mpnnfile_BC, args.receptor_msa, output_a3m_file, verbose=args.verbose)

# end of script
