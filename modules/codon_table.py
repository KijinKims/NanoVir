from types import MappingProxyType

# Standard genetic code:
codon = {}
# nTn
codon["TTT"] = codon["TTC"] = "F"
codon["TTA"] = codon["TTG"] = codon["CTT"] = codon["CTC"] = codon["CTA"] = codon["CTG"] = "L"
codon["ATT"] = codon["ATC"] = codon["ATA"] = "I"
codon["ATG"] = "M"
codon["GTT"] = codon["GTC"] = codon["GTA"] = codon["GTG"] = "V"

# nCn
codon["TCT"] = codon["TCC"] = codon["TCA"] = codon["TCG"] = "S"
codon["CCT"] = codon["CCC"] = codon["CCA"] = codon["CCG"] = "P"
codon["ACT"] = codon["ACC"] = codon["ACA"] = codon["ACG"] = "T"
codon["GCT"] = codon["GCC"] = codon["GCA"] = codon["GCG"] = "A"

# nAn
codon["TAT"] = codon["TAC"] = "Y"
codon["TAA"] = codon["TAG"] = "*" 	#Stop
codon["CAT"] = codon["CAC"] = "H"
codon["CAA"] = codon["CAG"] = "Q"
codon["AAT"] = codon["AAC"] = "N"
codon["AAA"] = codon["AAG"] = "K"
codon["GAT"] = codon["GAC"] = "D"
codon["GAA"] = codon["GAG"] = "E"

# nGn
codon["TGT"] = codon["TGC"] = "C"
codon["TGA"] = "*"  			#Stop
codon["TGG"] = "W"
codon["CGT"] = codon["CGC"] = codon["CGA"] = codon["CGG"] = "R"
codon["AGT"] = codon["AGC"] = "S"
codon["AGA"] = codon["AGG"] = "R"
codon["GGT"] = codon["GGC"] = codon["GGA"] = codon["GGG"] = "G"

codon_read_only = MappingProxyType(codon)