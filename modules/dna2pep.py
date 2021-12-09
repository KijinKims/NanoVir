from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

input_fasta = "example.fasta"
output_fasta = "example_output.fasta"

with open(output_fasta, "w") as output_handle:
    for record in SeqIO.parse(input_fasta, "fasta"):
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
            sliced_dna = dna[i:]
            sliced_dna = sliced_dna[:(len(sliced_dna)//3)*3]
            pep = sliced_dna.translate()
            pep_record = SeqRecord(
                pep,
                id=record.id + "_rframe-" + str(i+1),
                name=record.name,
                description=record.description
            )
            
            SeqIO.write(pep_record, output_handle, "fasta")