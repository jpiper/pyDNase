#!/usr/bin/env python

# Copyright (C) 2015 Jason Piper - j.piper@warwick.ac.uk
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os, argparse, tempfile, operator, pickle
from collections import defaultdict
import pybedtools, pysam

def rev_comp(str):
	"""
	Given a DNA string, returns the reverse complement
	
	>>> rev_comp("AATTGGCC")
	'GGCCAATT'
	"""
	rev_dic = {'A':'T','G':'C','C':'G','T':'A'}
	return ''.join([rev_dic[i] for i in str[::-1]])

def get_kmers(file_path):
	"""
	Counts the number of occurences of each sequence in a FASTA file 
	and returns a dictionary with the format {sequence:count}, e.g. {'AATTCC':5}
	
	Note: Ignores anything with Ns in, and treats all DNA as uppercase (i.e. not repeatmasked)
	"""
	kmers = defaultdict(int)
	with open(file_path,'r') as fasta:
		for i in fasta:
			if i[0] != ">":
				if 'N' not in i:
					kmers[i.strip().upper()] +=1
	return kmers

def generate_6mer_bed(bam_file, gdict):
	"""
	Loads data from a BAM file, and prints the 6-mers from the 5' ends of the reads
	as a BED file to a temp file, returns the file path of this temp file
	"""
	outfile = tempfile.NamedTemporaryFile(delete=False)
	samfile = pysam.AlignmentFile(bam_file, "rb")
	for i in samfile:
		# Ignore unmapped reads
		if not i.is_unmapped:
			chrom = samfile.getrname(i.reference_id)
			if chrom in gdict.keys():
				# Determine which end of the read is the 5' end
				if i.is_reverse:
					strand = "-"
					startbp, endbp = i.reference_end - 3, i.reference_end + 3
				else:
					strand  = "+"
					startbp, endbp = i.reference_start - 3, i.reference_start + 3
				if startbp > 0 and endbp < gdict[chrom]:
					print >> outfile, "\t".join((str(i) for i in (chrom, startbp, endbp, 0, 0, strand)))
	outfile.close()
	return outfile.name

def genome_dic(g_file):
	"""
	Make a dictionary of chromosome sizes from a .chrom.sizes file
	e.g. {chrom_name:chrom_size}
	"""
	gdict = {}
	with open(g_file) as ifile:
		for i in ifile:
			i = i.split()
			gdict[i[0]] = int(i[1])
	return gdict

if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description='Calculates the 6-mer 5\' insertion bias for a NGS dataset')
	parser.add_argument("regions", help="BED file of the regions you want to exclude from calculating the bias. This is usually the DHSs.")
	parser.add_argument("reads", help="The sorted, indexed BAM file containing the DNase-seq data")
	parser.add_argument("genome_sequence", help="The sorted, indexed FASTA file containing the genome sequence")
	parser.add_argument("genomesize", help="The .chrom.sizes file containing chromosome sizes generated using something like  \"mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \"select chrom, size from hg19.chromInfo\"  > hg19.chrom.sizes\"")
	parser.add_argument("output", help="output file prefix to write the observed/expected ratios to (will append .txt and .pickle)")
	args  = parser.parse_args()
	
	test_bam        = args.reads
	test_bed        = args.regions
	genome_sequence = args.genome_sequence
	genome          = args.genome_size
	outfile         = args.output
	
	# First, pull all the 6mers surrounding 5' ends
	print("Determining transposition sites (roughly 60s per 1E6 reads)...")
	bed_file_for_6mers = generate_6mer_bed(test_bam, genome_dic(genome))
	all_6mers = pybedtools.BedTool(bed_file_for_6mers)
	
	# Can't guarantee peaks are sorted, so sort them
	peaks = pybedtools.BedTool(test_bed)	
	peaks = peaks.sort()
	
	# Get the transposition sites outside of the peaks only
	print("Filtering for those outside peaks")
	bg_6mers      = all_6mers.intersect(peaks,v=True)
	
	print("Generating shuffled background")
	shuf_bg_6mers = bg_6mers.shuffle(g=genome, noOverlapping=True, excl=peaks.fn)
		
	print("Generating FASTA file for 6mers...")
	genome_fasta  = pybedtools.BedTool(genome_sequence)
	observed_cuts = bg_6mers.sequence(fi=genome_fasta)
	
	print("Generating FASTA file for shuffled 6mers...")
	shuffled_cuts = shuf_bg_6mers.sequence(fi=genome_fasta)

	print("Getting 6mers for observed...")
	observed = get_kmers(observed_cuts.seqfn)

	print("Getting 6mers for shuffled...")
	expected = get_kmers(shuffled_cuts.seqfn)
	
	print("Calculating...")
	enriched = {i:observed[i]/float(expected[i]) for i in observed.keys()}

	print("Dumping bias txt file...")
	with open(outfile + ".txt", 'w') as ofile:
		for i in sorted(enriched.items(), key=operator.itemgetter(1)):
			print >> ofile, "\t".join(map(str,i))

	print("Writing bias pickle file...")
	totalsum = float(sum(enriched.values()))
	whatdic = {key:{'forward':val/totalsum,'reverse':enriched[rev_comp(key)]/totalsum} for key,val in enriched.iteritems()}
	with open(outfile + ".pickle", "w") as bias_file:
	    pickle.dump(whatdic,bias_file)
		
	os.remove(bed_file_for_6mers)
