from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import os

def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    return arg

parser = argparse.ArgumentParser(description= 'Script for translating DNA sequence FASTA file into six-frame translated amino acid sequence FASTA file')

parser.add_argument('--input', '-i', metavar='input.fasta',
        help='Input FASTA file',
        type=lambda x: is_valid_file(parser, x))

parser.add_argument('--output', '-o', metavar='output.pep.fasta',
        help='Output translated FASTA file')

args = parser.parse_args()

with open(args.output, "w") as output_handle:
    for record in SeqIO.parse(args.input, "fasta"):
        dna = record.seq

        for i in range(3):
            sliced_dna = dna[i:]
            sliced_dna = sliced_dna[:(len(sliced_dna)//3)*3]
            pep = sliced_dna.translate()
            pep_record = SeqRecord(
                pep,
                id=record.id + "_rframe" + str(i+1),
                name=record.name,
                description=record.description
            )
            
            SeqIO.write(pep_record, output_handle, "fasta")

        revcom_dna = record.seq.reverse_complement()

        for i in range(3):
            sliced_dna = revcom_dna[i:]
            sliced_dna = sliced_dna[:(len(sliced_dna)//3)*3]
            pep = sliced_dna.translate()
            pep_record = SeqRecord(
                pep,
                id=record.id + "_rframe-" + str(i+1),
                name=record.name,
                description=record.description
            )
            
            SeqIO.write(pep_record, output_handle, "fasta")