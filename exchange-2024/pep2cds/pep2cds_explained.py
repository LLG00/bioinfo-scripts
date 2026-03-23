#!/usr/bin/env python3
# Author: Luiza Lima Galli
# Date: July 2025
# Description: Converts a protein alignment (PEP) into a codon-based CDS alignment
#              using mapped codons from CDS and PEP FASTA files.
# Note: Some parts of this script were developed with the help of AI tools (e.g., ChatGPT),
#       and were reviewed and tested by the author.
# License: MIT


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def map_codons(cds, pep):
    """
    Maps each amino acid to its corresponding codon(s) based on paired CDS and PEP sequences.

    Performs sanity checks:
      - Warns if a PEP sequence has no matching CDS entry
      - Warns if CDS length is not exactly 3x the PEP length
      - Warns if any CDS sequence does not start with ATG

    Args:
        cds (dict): {sequence_name: nucleotide_sequence}
        pep (dict): {sequence_name: amino_acid_sequence}

    Returns:
        dict: {amino_acid: [list of corresponding codons]}
    """
    # Verifies if every PEP sequence is present in the CDS file
    for seq_name in pep:
        if seq_name not in cds:
            print(f"Warning: Sequence name {seq_name} in PEP is not found in CDS.")

    # Verifies if the length of the CDS is 3 times the length of the corresponding PEP sequence
    for seq_name in cds:
        if seq_name in pep:
            if len(cds[seq_name]) != len(pep[seq_name]) * 3:
                print(f"Warning: Length of CDS sequence {seq_name} is not three times longer than the corresponding PEP sequence.")

    # Verifies if the CDS sequences start with the start codon "ATG"
    for seq_name, seq in cds.items():
        if not seq.startswith("ATG"):
            print(f"Warning: Sequence {seq_name} does not start with a start codon (ATG).")

    codons = {}     # Dictionary to store codons corresponding to each amino acid
    aminoacids = set()  # Set to store unique amino acids in PEP sequences

    for seq_name, seq in pep.items():
        for i in range(len(seq)):
            aminoacids.add(seq[i])  # Adds unique amino acids

    # Initialises empty lists for each amino acid in codons dictionary
    for aminoacid in aminoacids:
        codons[aminoacid] = []

    # For each CDS sequence, extracts codons and associates them to the corresponding amino acid in the PEP sequence
    for seq_name, seq in cds.items():
        for i in range(0, len(seq), 3):
            codon = seq[i:i+3]
            aminoacid = pep[seq_name][i//3]
            codons[aminoacid].append(codon)

    # Removes duplicated codons for each amino acid
    for aminoacid in codons:
        codons[aminoacid] = list(set(codons[aminoacid]))

    return codons


def map_codons_sequence(cds, codons_dict):
    """
    Produces an ordered list of (amino_acid, codon_index) pairs for a single CDS sequence.
    Used to maintain codon identity when rewriting a PEP alignment as CDS.

    Args:
        cds (str): nucleotide sequence (CDS)
        codons_dict (dict): {amino_acid: [codons]} from map_codons()

    Returns:
        list of (str, int) tuples: [(amino_acid, index_in_codon_list), ...]
    """
    codons_sequence = []  # List to store the ordered sequence of codons as pairs (amino_acid, index)

    for i in range(0, len(cds), 3):
        codon = cds[i:i+3]
        for aminoacid, codons in codons_dict.items():  # Searches the amino acid that corresponds to this codon in the codons dictionary
            if codon in codons:
                aminoacid_idx_in_dict = codons_dict[aminoacid].index(codon)
                codons_sequence.append((aminoacid, aminoacid_idx_in_dict))

    return codons_sequence


def rewrite_pep_alignment(codons_order, codons_dict, pep_alignment_dict):
    """
    Reconstructs a CDS alignment by replacing each amino acid in the PEP alignment
    with its corresponding codon. Gap characters ('-') are expanded to '---'.

    Args:
        codons_order (list): output of map_codons_sequence()
        codons_dict (dict): {amino_acid: [codons]}
        pep_alignment_dict (dict): {sequence_name: aligned_pep_sequence} (single entry)

    Returns:
        str: the resulting codon-based alignment string
    """
    # Gets the aligned PEP sequence (assuming one key only is in the dictionary)
    pep_alignment = pep_alignment_dict[list(pep_alignment_dict.keys())[0]]
    cds_alignment = ''  # String to store the final alignment in codons

    for i in range(len(pep_alignment)):
        if pep_alignment[i] == '-':  # If there is a gap in the PEP sequence, adds triple gaps in the CDS sequence
            cds_alignment += '---'
        elif codons_order[0][0] != pep_alignment[i]:  # If the expected amino acid and the actual one do not match, removes entries until they match (to synchronise)
            while codons_order[0][0] != pep_alignment[i]:
                codons_order.pop(0)
            cds_alignment += codons_dict[codons_order[0][0]][codons_order[0][1]]  # Adds the corresponding codon and removes the pair
            codons_order.pop(0)
        else:  # If the amino acid matches, adds the corresponding codon and removes the pair
            cds_alignment += codons_dict[codons_order[0][0]][codons_order[0][1]]
            codons_order.pop(0)

    # Warns if any amino acid in the alignment does not have a corresponding codon
    for aminoacid in pep_alignment:
        if aminoacid != '-' and aminoacid not in codons_dict:
            print(f"Warning: Amino acid {aminoacid} in PEP alignment does not correspond to any codon.")

    return cds_alignment


# ---------------------------------------------------------------------------
# Open input files
# ---------------------------------------------------------------------------

# Open the CDS FASTA file
dict_cds = {}
try:
    with open("example_data/example_CDS.fasta") as file:
        for line in file:
            if line.startswith('>'):  # Extracts the sequence name, removing '>' and possible indentation
                seq_name = line.strip().split()[0][1:]
                dict_cds[seq_name] = ''
            else:  # Concatenates lines of the sequence
                dict_cds[seq_name] += line.strip()
except FileNotFoundError:
    print("Warning: CDS file not found or format is incorrect.")

# Open the PEP FASTA file
dict_pep = {}
try:
    with open("example_data/example_PEP.fasta") as file:
        for line in file:
            if line.startswith('>'):
                seq_name = line.strip().split()[0][1:]
                dict_pep[seq_name] = ''
            else:
                dict_pep[seq_name] += line.strip()
except FileNotFoundError:
    print("Warning: PEP file not found or format is incorrect.")

# Open the PEP alignment file
dict_pep_alignment = {}
try:
    with open("example_data/example_PEP_align.fasta") as file:
        for line in file:
            if line.startswith('>'):
                seq_name = line.strip().split()[0][1:]
                dict_pep_alignment[seq_name] = ''
            else:
                dict_pep_alignment[seq_name] += line.strip()
except FileNotFoundError:
    print("Warning: PEP Alignment file not found or format is incorrect.")


# ---------------------------------------------------------------------------
# Main pipeline — create CDS alignment based on PEP alignment
# ---------------------------------------------------------------------------

with open("output_CDS_align.fasta", "w") as out:
    for seq_name, seq in dict_cds.items():

        # Creates auxiliary dictionaries for the current CDS sequence,
        # its corresponding PEP sequence, and the PEP alignment
        cds_aux           = {seq_name: seq}
        pep_aux           = {seq_name: dict_pep[seq_name]}
        pep_alignment_aux = {seq_name: dict_pep_alignment[seq_name]}

        codons_dict   = map_codons(cds_aux, pep_aux)                              # Gets the mapped codons dictionary for the sequence
        codons_order  = map_codons_sequence(seq, codons_dict)                     # Maps the CDS sequence in the correct codon order
        cds_alignment = rewrite_pep_alignment(codons_order, codons_dict, pep_alignment_aux)  # Creates the CDS alignment based on the PEP alignment

        out.write(f">{seq_name}\n{cds_alignment}\n")  # Writes output file

        # Prints the result
        print(f">{seq_name}")
        print(cds_alignment)

print("\nDone! CDS alignment saved to output_CDS_align.fasta")
