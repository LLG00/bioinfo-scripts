# Author: Luiza Lima Galli
# Date: July 2025
# Description: Script to calculate query and subject coverage by merging HSP intervals from BLAST XML output.
# Note: Some code components were generated with the help of AI tools (e.g., ChatGPT) and adapted/tested by me.
# License: MIT

from Bio.Blast import NCBIXML
from Bio import SeqIO

query_lengths = {record.id: len(record.seq) for record in SeqIO.parse("query_fasta_example.fasta", "fasta")}

hit_lengths = {record.id: len(record.seq) for record in SeqIO.parse("subject_fasta_example.fasta", "fasta")}

def coverage_calculate(intervals):
    if not intervals:
        return 0
    intervals.sort(key=lambda x: x[0])  
    coverage = 0
    start, end = intervals[0]
    for actual_start, actual_end in intervals[1:]:
        if actual_start <= end + 1:
            end = max(end, actual_end)
        else:
            coverage += end - start + 1
            start, end = actual_start, actual_end
    coverage += end - start + 1  
    return coverage

with open("example_output.txt", "w") as output:
    output.write(
        "query_id\tquery_length\thit_id\thit_length\tidentity\tbit_score\te_value\tquery_start\tquery_end\thit_start\thit_end\tquery_coverage(%)\tsubject_coverage(%)\n"
    )

    with open("example_blast.xml") as handle:
        for record in NCBIXML.parse(handle):
            query_id = record.query.split()[0] if record.query else "unknown_query"
            query_len = query_lengths.get(query_id, 0)

            for alignment in record.alignments:
                hit_id = alignment.hit_id.split('|')[-2] if '|' in alignment.hit_id else alignment.hit_id
                hit_len = hit_lengths.get(hit_id, 0)

                query_intervals = []
                subject_intervals = []

                for hsp in alignment.hsps:
                    align_len = hsp.align_length
                    identity = hsp.identities / align_len if align_len else 0
                    bit_score = hsp.bits
                    e_value = hsp.expect
                    query_start = hsp.query_start
                    query_end = hsp.query_end
                    hit_start = hsp.sbjct_start
                    hit_end = hsp.sbjct_end

                    query_intervals.append((min(query_start, query_end), max(query_start, query_end)))
                    subject_intervals.append((min(hit_start, hit_end), max(hit_start, hit_end)))

                query_coverage = coverage_calculate(query_intervals) / query_len * 100 if query_len else 0
                subject_coverage = coverage_calculate(subject_intervals) / hit_len * 100 if hit_len else 0

                output.write(
                    f"{query_id}\t{query_len}\t{hit_id}\t{hit_len}\t{identity:.6f}\t{bit_score:.2f}\t{e_value:.2e}\t"
