### convert_NCBI_DMfasta_for_gff.py

import sys
import os
import getopt

try:
	opts, args = getopt.getopt(sys.argv[1:], 'o:f:g:h')
																					 	
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)


NCBI_fasta_file_name = None
old_annot_file_name  = None
output_file_name = None

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h', '--help'):
		print("\n**** convert_NCBI_DMfasta_for_gff.py | Written by DJP, 21/10/18 in Python 3.5 in Basel, Swiss ****\n")
		print("converts the NCBI fasta headers to match the gff (dryad)") 
		print("\n**** USAGE ****\n")
		print("python3 convert_NCBI_DMfasta_for_gff.py -f [NCBI fasta file] -g [gff file (dryad)] -o [output file name]") 		
		
		sys.exit(2)
		
	elif opt in ('-f'):
		NCBI_fasta_file_name = arg
	elif opt in ('-g'):
		old_annot_file_name = arg
	elif opt in ('-o'):		
		output_file_name = arg
	else:
		print("i dont know")
		sys.exit(2)


### NCBI seqs
##### FIRST unwrap fasta - precautionary will be necessary for some files 
### note making a temp unwrapped fasta file  - removed at end

seqF1 = NCBI_fasta_file_name
output_fasta_name = seqF1 + ".TEMP_extract_fasta_file" 

output_file = open(output_fasta_name, "w")
print("Unwrapping fasta file\n")
count = 0
in_file = open(seqF1)
for line in in_file:
	count = count + 1
	line = line.rstrip("\n")
	if line.startswith(">") and count == 1:
		output_file.write(line + "\n")
	elif line.startswith(">") and count > 1:
		output_file.write("\n" + line + "\n")
	else: 
		output_file.write(line)	

output_file.close()

										
### add seqs to dictionary
name_list = []
seq_list  = []
seq_dict_NCBI_fasta  = {}

all_conv_name = set()
all_conv_name_list = []

done = 0
seq_file_1 = open(output_fasta_name)
for line in seq_file_1:
	lineA = line.rstrip("\n")
	if lineA.startswith(">"):
		lineB = lineA.replace(">", "")
		name_list.append(lineB)
	else:
		seq_list.append(lineA)
		done = done + 1

for element in range(0,len(name_list)):
	name1 = name_list[element].split(",")[0].split(" ")[-1]
	#print(name1)
	all_conv_name.add(name1)
	all_conv_name_list.append(name1)
	seq1 = seq_list[element]
	seq_dict_NCBI_fasta[name1] = seq1

seq_N = str(len(name_list))

print ("\nLoaded " + seq_N + " seqeunces from " + seqF1 + ".\n")

## tidyup
seq_file_1.close()
os.remove(output_fasta_name)


###### check seq names match to the gff


# read gff

gff_file = open(old_annot_file_name)
for line in gff_file:
	line = line.rstrip("\n").split("\t")[0]
	if line not in all_conv_name:
		print("Name conversion failed, exiting\n\n")
		sys.exit(2)
gff_file.close()

print("All seq names in " + NCBI_fasta_file_name + " successfully converted to match those in the gff (" + old_annot_file_name + ")\n")


#### output converted fasta

output_file = open(output_file_name, "w")

for s in all_conv_name_list:
	seq_out = seq_dict_NCBI_fasta.get(s)
	output_file.write(">" + s + "\n" + seq_out + "\n")

print("converted fasta file outputted to " + output_file_name + "\n\n\nFinished, William Black\n\n\n" )















