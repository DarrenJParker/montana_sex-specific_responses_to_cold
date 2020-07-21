##### CT_SB_edgeR_tidier.py


import sys
import os
import getopt
import decimal

try:
	opts, args = getopt.getopt(sys.argv[1:], 'a:b:f:m:g:i:F:o:v:h')
																						
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)

SB_6_filename        = None
SB_19_filename       = None

CT_F_filename        = None
CT_M_filename        = None

interaction_filename = None
gff_filename         = None
out_base = "Ishouldhavenamedthis"
vir_SBGE_filename = None

FDR_cutoff = None

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h', '--help'):
		print("\n**** CT_SB_edgeR_tidier.py ****\n")
		print("\ntakes edgeR output files from CT_edgeR.R and sticks them together with Dvir orth info from gff")
		print("see readme for details on use")

		sys.exit(2)

	elif opt in ('-F'):
		FDR_cutoff = arg
	elif opt in ('-a'):
		SB_6_filename = arg
	elif opt in ('-b'):
		SB_19_filename = arg
	elif opt in ('-f'):
		CT_F_filename = arg
	elif opt in ('-m'):
		CT_M_filename = arg
	elif opt in ('-i'):
		interaction_filename = arg
	elif opt in ('-o'):
		out_base = arg
	elif opt in ('-g'):
		gff_filename = arg
	elif opt in ('-v'):
		vir_SBGE_filename = arg
	else:
		print("i dont know")
		sys.exit(2)


if FDR_cutoff == None:
	print("\nFDR cutoff for calling sex-bias set to 0.05. Use -F to alter this\n")
	FDR_cutoff = 0.05
else:
	print("\nFDR cutoff set to " + str(FDR_cutoff) +"\n")


### all file info to dicts

all_gene_names_from_TTT = set()

### get SB_6

SB_dict_6 = {}

SB_6_file = open(SB_6_filename)

line_N = 0
for line in SB_6_file:
	line_N = line_N + 1
	if line_N > 1:
		line = line.rstrip("\n").split(",")
		# print(line)
		gene_name    = line[0]
		all_gene_names_from_TTT.add(gene_name)
		log2_FC      = decimal.Decimal(line[1])
		FDR_val      = decimal.Decimal(line[5])
		bias = ""
		
		if FDR_val < FDR_cutoff and log2_FC > 0: #### female biased
			bias = "female_biased"
		elif FDR_val < FDR_cutoff and log2_FC < 0: #### male biased
			bias = "male_biased"
		elif FDR_val < FDR_cutoff and log2_FC == 0: #### just in case
			bias = "ERROR?"
		else:
			bias = "Unbiased"
		
		
		### add to dict
		
		new_rec = str(log2_FC) + "," + str(FDR_val) + "," + bias
		SB_dict_6[gene_name] = new_rec 
		# print(gene_name)
		# print(new_rec)
SB_6_file.close()


### get SB_19

SB_dict_19 = {}

SB_19_file = open(SB_19_filename)

line_N = 0
for line in SB_19_file:
	line_N = line_N + 1
	if line_N > 1:
		line = line.rstrip("\n").split(",")
		# print(line)
		gene_name    = line[0]
		all_gene_names_from_TTT.add(gene_name)
		log2_FC      = decimal.Decimal(line[1])
		FDR_val      = decimal.Decimal(line[5])
		bias = ""
		
		if FDR_val < FDR_cutoff and log2_FC > 0: #### female biased
			bias = "female_biased"
		elif FDR_val < FDR_cutoff and log2_FC < 0: #### male biased
			bias = "male_biased"
		elif FDR_val < FDR_cutoff and log2_FC == 0: #### just in case
			bias = "ERROR?"
		else:
			bias = "Unbiased"
		
		
		### add to dict
		
		new_rec = str(log2_FC) + "," + str(FDR_val) + "," + bias
		SB_dict_19[gene_name] = new_rec 
		# print(gene_name)
		# print(new_rec)
SB_19_file.close()


### get CT_F

CT_dict_F = {}

CT_F_file = open(CT_F_filename)

line_N = 0
for line in CT_F_file:
	line_N = line_N + 1
	if line_N > 1:
		line = line.rstrip("\n").split(",")
		# print(line)
		gene_name    = line[0]
		all_gene_names_from_TTT.add(gene_name)
		log2_FC      = decimal.Decimal(line[1])
		FDR_val      = decimal.Decimal(line[5])
		
		### add to dict
		
		new_rec = str(log2_FC) + "," + str(FDR_val) 
		CT_dict_F[gene_name] = new_rec 
		# print(gene_name)
		# print(new_rec)
CT_F_file.close()


### get CT_M

CT_dict_M = {}

CT_M_file = open(CT_M_filename)

line_N = 0
for line in CT_M_file:
	line_N = line_N + 1
	if line_N > 1:
		line = line.rstrip("\n").split(",")
		# print(line)
		gene_name    = line[0]
		all_gene_names_from_TTT.add(gene_name)
		log2_FC      = decimal.Decimal(line[1])
		FDR_val      = decimal.Decimal(line[5])
		
		### add to dict
		
		new_rec = str(log2_FC) + "," + str(FDR_val) 
		CT_dict_M[gene_name] = new_rec 
		# print(gene_name)
		# print(new_rec)
CT_M_file.close()


##### interaction

dict_interaction = {}

interaction_file = open(interaction_filename)

line_N = 0
for line in interaction_file:
	line_N = line_N + 1
	if line_N > 1:
		line = line.rstrip("\n").split(",")
		# print(line)
		gene_name    = line[0]
		all_gene_names_from_TTT.add(gene_name)
		log2_FC      = decimal.Decimal(line[1])
		FDR_val      = decimal.Decimal(line[5])
		
		### add to dict
		
		new_rec = str(log2_FC) + "," + str(FDR_val) 
		dict_interaction[gene_name] = new_rec 
		# print(gene_name)
		# print(new_rec)
interaction_file.close()



################################################################################################
### get Dvir RBBH from gff


RBBH_dict = {}

gff_file = open(gff_filename)
for line in gff_file:
	line = line.rstrip("\n").split("\t")
	if line[2] == "gene":
		desc = line[8]
		gene_name = desc.split("ID=")[1].split(";")[0]
		RBBH_vir_FB = ""
		try:
			RBBH_vir_FB = desc.split("_RBBH_Vir_FB_N=")[1].split(";")[0]
		except:
			RBBH_vir_FB = "NA"

		RBBH_vir_FBsym = ""
		try:
			RBBH_vir_FBsym = desc.split("_RBBH_Vir_FB_sym=")[1].split(";")[0]
		except:
			RBBH_vir_FBsym = "NA"
		
		
		RBBH_dict[gene_name] = RBBH_vir_FB + "," + RBBH_vir_FBsym

		

######################################################################################################################################
### vir SB info

vir_SB_dict = {}

vir_SBGE_file = open(vir_SBGE_filename)

vir_SB_head = ""

line_N = 0
for line in vir_SBGE_file:
	line_N = line_N + 1
	line = line.strip().split(",")
	if line_N == 1:
		vir_SB_head = "vir_" + line[2] + ",vir_" + line[6]
	else:
		
		bias = ""
		
		if line[6] != "NA":

			### note values are opp dir to mine (female biased = -ve )
			if decimal.Decimal(line[6]) < FDR_cutoff and decimal.Decimal(line[2]) < 0: #### female biased
				bias = "female_biased"
			elif decimal.Decimal(line[6]) < FDR_cutoff and decimal.Decimal(line[2]) > 0: #### male biased
				bias = "male_biased"
			elif decimal.Decimal(line[6]) < FDR_cutoff and decimal.Decimal(line[2]) == 0: #### just in case
				bias = "ERROR?"
			else:
				bias = "Unbiased"
		else:
			bias = "NA"
			
		
		
		vir_SB_dict[line[0]] = line[2] + "," + line[6] + "," + bias



######################################################################################################################################
###### output

out_gene_list = sorted(list(all_gene_names_from_TTT))
out_file = open(out_base + "_CTSB_tidied.csv", "w")

out_file.write("gene_name" + "," +
			   "CT_F_logFC,CT_F_FDR," +
			   "CT_M_logFC,CT_M_FDR," +
			   "SB_6_logFC,SB_6_FDR,SB_6_sexbias," +
			   "SB_19_logFC,SB_19_FDR,SB_19_sexbias," +
			   "CT_SB_inter_logFC,CT_SB_inter_FDR,"
			   "Dvir_FBN,Dvir_FBsym," + vir_SB_head + ",vir_sexbias\n")

for gene in out_gene_list:
	

	SB_6_info        = SB_dict_6.get(gene)
	SB_19_info       = SB_dict_19.get(gene)
	
	CT_F_info         = CT_dict_F.get(gene)
	CT_M_info         = CT_dict_M.get(gene)
	
	interaction_info  = dict_interaction.get(gene)
	
	rec_list = [SB_6_info, SB_19_info, CT_F_info, CT_M_info, interaction_info]
	for i in rec_list:
		if i == None:
			print("\n\nMissing genes - sort this out. Exiting\n\n")
			sys.exit(2)
	
	
	RBBH_info = RBBH_dict.get(gene).strip()
	vir_SB_info = vir_SB_dict.get(RBBH_info.split(",")[0])
	if vir_SB_info == None:
		vir_SB_info = "NA,NA,NA"
	
	print(RBBH_info )
	#print(vir_SB_info)
	# 
	out_file.write(gene + "," + CT_F_info + "," + CT_M_info + "," + SB_6_info + "," + SB_19_info + ","  + interaction_info + "," + RBBH_info + "," + vir_SB_info + "\n")
		
	# print(gene)
	# print(rec_list)



print("\n\n\nFinished\n\n")
	



