import os
import sys
import re


genes = dict()
geneStrands = dict()
with sys.stdin as lines:
	for line in lines:

		if line[0] == "#":
			continue 
			
		cols = line.rstrip().split("\t")
		feature = cols[2]

		if feature != "exon":
			continue

		chrom, c1, c2, strand, data = cols[0], cols[3], cols[4], cols[6], cols[-1]

		gid = re.search('gene_name "([^"]+)";', data).group(1)
		tid = re.search('transcript_id "([^"]+)";', data).group(1)
		
		if chrom not in genes:
			genes[chrom] = dict()

		if gid not in genes[chrom]:
			genes[chrom][gid] = dict()

		if tid not in genes[chrom][gid]:
			genes[chrom][gid][tid] = list()

		geneStrands[gid] = strand
		genes[chrom][gid][tid].append((int(c1),int(c2)))

uniqueDict = dict()
for chrom, allgenes in genes.items():
	for gene, txns in allgenes.items():
		for tid, coords in txns.items():

			if geneStrands[gene] == "-":
				coords = coords[::-1]

			for num, coord in enumerate(coords,0):
				if num == 0 :
					continue

				j2, null = coord
				null, j1 = coords[num-1]
				j1 , j2 = j1-1, j2-1
				jcnSet = (j1, j2, gene)

				if jcnSet not in uniqueDict:
					print(chrom, j1, j2, gene, ".", geneStrands[gene], sep="\t")
					uniqueDict[jcnSet] = 1

