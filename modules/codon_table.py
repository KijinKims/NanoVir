from types import MappingProxyType

# Standard genetic code:
standard_codon = {}
# nTn
standard_codon["TTT"] = standard_codon["TTC"] = "F"
standard_codon["TTA"] = standard_codon["TTG"] = standard_codon["CTT"] = standard_codon["CTC"] = standard_codon["CTA"] = standard_codon["CTG"] = "L"
standard_codon["ATT"] = standard_codon["ATC"] = standard_codon["ATA"] = "I"
standard_codon["ATG"] = "M"
standard_codon["GTT"] = standard_codon["GTC"] = standard_codon["GTA"] = standard_codon["GTG"] = "V"

# nCn
standard_codon["TCT"] = standard_codon["TCC"] = standard_codon["TCA"] = standard_codon["TCG"] = "S"
standard_codon["CCT"] = standard_codon["CCC"] = standard_codon["CCA"] = standard_codon["CCG"] = "P"
standard_codon["ACT"] = standard_codon["ACC"] = standard_codon["ACA"] = standard_codon["ACG"] = "T"
standard_codon["GCT"] = standard_codon["GCC"] = standard_codon["GCA"] = standard_codon["GCG"] = "A"

# nAn
standard_codon["TAT"] = standard_codon["TAC"] = "Y"
standard_codon["TAA"] = standard_codon["TAG"] = "*" 	#Stop
standard_codon["CAT"] = standard_codon["CAC"] = "H"
standard_codon["CAA"] = standard_codon["CAG"] = "Q"
standard_codon["AAT"] = standard_codon["AAC"] = "N"
standard_codon["AAA"] = standard_codon["AAG"] = "K"
standard_codon["GAT"] = standard_codon["GAC"] = "D"
standard_codon["GAA"] = standard_codon["GAG"] = "E"

# nGn
standard_codon["TGT"] = standard_codon["TGC"] = "C"
standard_codon["TGA"] = "*"  			#Stop
standard_codon["TGG"] = "W"
standard_codon["CGT"] = standard_codon["CGC"] = standard_codon["CGA"] = standard_codon["CGG"] = "R"
standard_codon["AGT"] = standard_codon["AGC"] = "S"
standard_codon["AGA"] = standard_codon["AGG"] = "R"
standard_codon["GGT"] = standard_codon["GGC"] = standard_codon["GGA"] = standard_codon["GGG"] = "G"

standard_codon = MappingProxyType(standard_codon)