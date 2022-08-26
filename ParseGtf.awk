#!/usr/bin/awk

# inputs .gtf
# outputs two fasta files: XXX.bit_1.fasta and XXX.bit_2.fasta


BEGIN{
	FS="\t"
	OFS="\t"
}



$3 = "gene"{

	print a[1]
}


