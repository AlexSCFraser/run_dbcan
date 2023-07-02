#!/usr/bin/env python3
##########################################################
# hmmscan parser for dbCAN meta server
#
# Based off the hmmscan parser used in the dbCAN server,
# written by Dr. Yin
#
# Written by Tanner Yohe under the supervision
# of Dr. Yin in the YinLab at NIU.
#
# Updated by Le Huang from tips the contributor WATSON Mick <mick.watson@roslin.ed.ac.uk>,
# Thank you!
#
# Modified by Alex Fraser to have a run() method that can be called and returns data for better integration with other
# scripts. This script also retains the ability to be called from shell and output to pipe redirection.
# This file had to be renamed from "hmmscan-parser.py" to "hmmscan_parser.py" because of python module import conventions.
# Modified on 07/06/22
#
# Modified by Alex Fraser to be a pure python implementation, completely removing the UNIX system requirement.
# This makes it portable across UNIX and windows, more readable, and allows for interactive debugging.
# 02/09/22
#
# INPUT
# python hmmscan-parser-dbCANmeta.py [inputFile] [eval] [coverage]
# eval and coverage are optional, inputFile is required
# -updating info:
# -adds pid for every subprocess to make codes robust.
# Last updated: 1/10/19
###########################################################
import re
from collections import defaultdict
from subprocess import call
import subprocess
import sys
import os


def parse_hmmer_file(filename):
	# OLD SHELL PIPE CODE, Replaced with python file parsing into memory and sorting below
	# cat "+input_file+"  | grep -v '^#' | awk '{print $1,$3,$4,$6,$13,$16,$17,$18,$19}' | sed 's/ /\t/g' | sort -k 3,3 -k 8n -k 9n
	space_re = re.compile(" +")
	data = []
	with open(filename, 'r') as f:
		for line in f:
			if line[0] != '#':
				line_array = space_re.sub(" ", line).split(' ')
				if len(line_array) < 22:
					raise Exception(f"hmmscan output file {filename} seems to have truncated output. One possible "
									f"cause of this is that your system ran out of memory for HMMER hmmscan process.")
				if line_array[17] == line_array[18]:  # Same as the "next if $a[-1] == $a[-2];" line in old perl code
					pass  # skip modules with same start and end index, they are erroneous hits
				else:
					data.append(
						[line_array[0], int(line_array[2]), line_array[3], int(line_array[5]), float(line_array[12]),
						 int(line_array[15]), int(line_array[16]), int(line_array[17]), int(line_array[18])]
					)

	sorted_data = sorted(data, key=lambda x: (x[2], x[7], x[8]))

	# OLD PERL CODE, this has been replaced by the python dict creation below.
	# 	while ( <> ) {
	#     chomp;
	#     @a = split;
	#     next if $a[-1] == $a[-2];
	#     push(@ {
	#         $b {
	#             $a[2]
	#         }
	#     }, $_);
	# }

	# group the entries by accession using a dict with key: accession and value: list of entries of that accession
	data_dict = defaultdict(list)
	for entry in sorted_data:
		data_dict[entry[2]].append(entry)

	# OLD PERL CODE, this has been replaced by the python loops below.
	# foreach(sort keys %b) {
	#         @a = @ {
	#             $b {
	#                 $_
	#             }
	#         };
	#         for ($i = 0; $i < $#a; $i++) {
	#             @b = split(/\t/, $a[$i]);
	#             @c = split(/\t/, $a[$i + 1]);
	#             $len1 = $b[-1] - $b[-2];
	#             $len2 = $c[-1] - $c[-2];
	#             $len3 = $b[-1] - $c[-2];
	#             if ($len3 > 0 and($len3 / $len1 > 0.5 or $len3 / $len2 > 0.5)) {
	#                 if ($b[4] < $c[4]) {
	#                     splice(@a, $i + 1, 1);
	#                 } else {
	#                     splice(@a, $i, 1);
	#                 }
	#                 $i = $i - 1;
	#             }
	#         }
	#         foreach(@a) {
	#                 print $_.\"\n\";}}

# 	Outer loop runs through each accession grouped list
	for entry_list in data_dict.values():
		# Inner loop compares consecutive entries to remove modules from the same accession which overlap, keeping the
		# one with higher e-value
		i = 0
		while i < len(entry_list)-1:
			b = entry_list[i]
			c = entry_list[i+1]
			len1 = b[-1] - b[-2]
			len2 = c[-1] - c[-2]
			len3 = b[-1] - c[-2]
			if len3 > 0 and (len3 / len1 > 0.5 or len3 / len2 > 0.5):  # checks for effective overlap of modules
				if b[4] < c[4]:
					# c has higher probability of being erroneous, therefore is removed
					entry_list.pop(i+1)  # splice(@a, $i + 1, 1)
				else:
					# b has higher probability of being erroneous, therefore is removed
					entry_list.pop(i)  # splice(@a, $i, 1)

				i -= 1
			i += 1

	valid_entries = []
	for accession_group in data_dict.values():
		for accession_entry in accession_group:
			valid_entries.append(accession_entry)

	return valid_entries


def run(input_file, eval_num=1e-15, coverage=0.35, verbose=False):

	# bash pipeline dumping data into tempfile replaced with parse_hmmer_file() for clarity and better cross-platform
	# compatibility
	# tmpfile = "temp." + str(os.getpid())
	# 	call("cat "+input_file+"  | grep -v '^#' | awk '{print $1,$3,$4,$6,$13,$16,$17,$18,$19}' | sed 's/ /\t/g' | sort -k 3,3 -k 8n -k 9n | perl -e 'while(<>){chomp;@a=split;next if $a[-1]==$a[-2];push(@{$b{$a[2]}},$_);}foreach(sort keys %b){@a=@{$b{$_}};for($i=0;$i<$#a;$i++){@b=split(/\t/,$a[$i]);@c=split(/\t/,$a[$i+1]);$len1=$b[-1]-$b[-2];$len2=$c[-1]-$c[-2];$len3=$b[-1]-$c[-2];if($len3>0 and ($len3/$len1>0.5 or $len3/$len2>0.5)){if($b[4]<$c[4]){splice(@a,$i+1,1);}else{splice(@a,$i,1);}$i=$i-1;}}foreach(@a){print $_.\"\n\";}}' > " + tmpfile, shell=True)

	hmmer_data = parse_hmmer_file(input_file)

	intermediate_data = ""

	for row in hmmer_data:
		row.append(float(int(row[6])-int(row[5]))/int(row[1]))
		if float(row[4]) <= float(eval_num) and float(row[-1]) >= float(coverage):
			if verbose:
				print('\t'.join([str(x) for x in row]))
			intermediate_data += '\t'.join([str(x) for x in row]) + '\n'

	# with open(tmpfile) as f:
	# 	for line in f:
	# 		row = line.rstrip().split('\t')
	# 		row.append(float(int(row[6])-int(row[5]))/int(row[1]))
	# 		if float(row[4]) <= eval_num and float(row[-1]) >= coverage:
	# 			if verbose:
	# 				print('\t'.join([str(x) for x in row]))
	# 			intermediate_data += '\t'.join([str(x) for x in row]) + '\n'
	# # call(['rm', tmpfile])
	# os.remove(tmpfile)

	return intermediate_data


if __name__ == "__main__":
	if len(sys.argv) > 3:
		file = sys.argv[1]
		eval_arg = float(sys.argv[2])
		coverage_arg = float(sys.argv[3])
		run(file, eval_arg, coverage_arg, verbose=True)
	if len(sys.argv) > 1:
		file = sys.argv[1]
		run(file, verbose=True)
	else:
		print("Please give a hmmscan output file as the first command")
		exit()
