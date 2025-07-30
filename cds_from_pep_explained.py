# Author: Luiza Lima Galli
# Date: July 2025
# Description: Converts a protein alignment (PEP) into a codon-based CDS alignment using mapped codons from CDS and PEP FASTA files.
# Note: Some parts of this script were developed with the help of AI tools (e.g., ChatGPT), and were reviewed and tested by the author.
# License: MIT

# Map codons function

def map_codons(cds, pep):

    # Verifies if every PEP sequnce is present in the CDS file
    
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
    
    codons = {} # Dictionary to store codons corresponding to each aminoacid
    
    aminoacids = set() # Set to store unique aminoacids in PEP sequences
    for seq_name, seq in pep.items():
        for i in range(len(seq)):
            aminoacids.add(seq[i]) # Adds unique aminoacids
            
    # Initializes empty lists to each aminoacid in codons dictionary
    
    for aminoacid in aminoacids: 
        codons[aminoacid] = [] 

    # For each CDS sequence, extracts codons and associates it to the aminoacid corresponding to the ´PEP sequence
    
    for seq_name, seq in cds.items(): 
        for i in range(0, len(seq), 3):
            codon = seq[i:i+3]
            aminoacid = pep[seq_name][i//3] 
            codons[aminoacid].append(codon)
            
    # Removes duplicated codons to each aminoacid 
    
    for aminoacid in codons:
        codons[aminoacid] = list(set(codons[aminoacid]))        
    
    return codons

# Map codons sequence

def map_codons_sequence(cds, codons_dict):
    codons_sequence = [] # List to store the ordered sequence of codons as pairs (aminoacid, index)
    for i in range(0, len(cds), 3):
        codon = cds[i:i+3]
        for aminoacid, codons in codons_dict.items(): # Searches the aminoacid that corresponds to these codons at the codons dictionary
            if codon in codons:
                aminoacid_idx_in_dict = codons_dict[aminoacid].index(codon)
                codons_sequence.append((aminoacid, aminoacid_idx_in_dict))
    
    return codons_sequence

# Rewrite PEP-alignment with codon-sequence

# Gets the aligned PEP sequence (assuming one key only is in the dictionary)

def rewrite_pep_alignment(codons_order, codons_dict, pep_alignment_dict):
    pep_alignment = pep_alignment_dict[list(pep_alignment_dict.keys())[0]]
    cds_alignment = '' # String to store the final alignment in codons
    
    for i in range(len(pep_alignment)):
        if pep_alignment[i] == '-': # If there is a gap in the PEP sequence, it adds triple gaps in the CDS sequence
            cds_alignment += '---' 
        elif codons_order[0][0] != pep_alignment[i]: # If the expected aminoacid and the actual one does not match, removes it until they match (to synchronize)
            while codons_order[0][0] != pep_alignment[i]: # Adds the corresponding codon and removes the pair
                codons_order.pop(0) 
            cds_alignment += codons_dict[codons_order[0][0]][codons_order[0][1]]
            codons_order.pop(0) 
        else: # If the aminoacid matches, adds the corresponding codon and removes the pair
            cds_alignment += codons_dict[codons_order[0][0]][codons_order[0][1]]
            codons_order.pop(0)

    # Warns if any aminoacid in the alignment does not have a corresponding codon
    
    for aminoacid in pep_alignment:
        if aminoacid != '-' and aminoacid not in codons_dict:
            print(f"Warning: Amino acid {aminoacid} in PEP alignment does not correspond to any codon.")
    
    return cds_alignment

# Open the CDS fasta file

dict_cds = {}
try:
    with open("example_CDS.fasta") as file:
        for line in file:
            if line.startswith('>'): # Extracts the name of the sequence, removing '>' and possible indentations
                seq_name = line.strip().split()[0][1:]
                dict_cds[seq_name] = ''
            else: # Concatenates lines of the sequence
                dict_cds[seq_name] += line.strip()
except FileNotFoundError:
    print("Warning: CDS file not found or format is incorrect.")

# Open the PEP fasta file

dict_pep = {}
try:
    with open("example_PEP.fasta") as file:
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
    with open("example_PEP_align.txt") as file:
        for line in file:
            if line.startswith('>'):
                seq_name = line.strip().split()[0][1:]
                dict_pep_alignment[seq_name] = ''
            else:
                dict_pep_alignment[seq_name] += line.strip()
except FileNotFoundError:
    print("Warning: PEP Alignment file not found or format is incorrect.")

# Loop to create cds alignment based on PEP alignment

with open("example_CDS_align.txt", "w") as out:
    for seq_name, seq in dict_cds.items(): # Creates auxiliary dictionaries for the current CDS sequence, its corresponding PEP sequence and the PEP alignment
        cds_aux = {}
        cds_aux[seq_name] = seq
    
        pep_aux = {}
        pep_aux[seq_name] = dict_pep[seq_name]
    
        pep_alignment_aux = {}
        pep_alignment_aux[seq_name] = dict_pep_alignment[seq_name]
    
        codons_dict = map_codons(cds_aux, pep_aux) # Gets the mapped codons dictionary to the sequence
        codons_order = map_codons_sequence(seq, codons_dict) # Maps the CDS sequence in the codon correct order according to the codon dictionary
        cds_alignment = rewrite_pep_alignment(codons_order, codons_dict, pep_alignment_aux) # Creates the CDS alignment based on the PEP alignment and codon order

        out.write(f">{seq_name}\n{cds_alignment}\n") # Writes an output file

        # Prints the result 

        print(f">{seq_name}")
        print(cds_alignment)
