import os
import sys
import math
import numpy

# make a file that will be used as input in R for the y-axis (diff-seq data)		
def prepareDS(args):
	files = args["<data>"]
	files1 = files.split(',')
	ref = args["--reference"]
	dirIN = 'DataFolder/Coverage'+ref+'/'
	os.system('mkdir DataFolder/NXDS_comparisonGraphs') # consider the diff-seq coverage
	dirOUT = 'DataFolder/NXDS_comparisonGraphs/'
	mut = args["--SNPfile"]
	
	base=['A','C','G','T']	
	
	NR={} # store the SNP positions as keys, and the non reference nucleotide as value
	Fs = open('Polymorphisms/' + mut, 'rU')
	for line in Fs:
		Ls = line.split()
		if Ls[2]=='SNP':
			NR[Ls[0]]=Ls[3][1]
	Fs.close()	
	
	for i in files1:
		D={} # store per position the non-reference allele count
		Total=0.0 # store total coverage to derive frequencies
		
		f=open(dirIN+i+'_Coverage.txt','rU')
		for line in f:
			if line.startswith("Position"):
				continue
			else:
				L1=line.split()
				if L1[0] in NR:
					k = base.index(NR[L1[0]]) # find the index of the non-ref allele in base
					D[L1[0]] = float(L1[k+12]) + float(L1[k+16])
				else:
					Alleles = [float(L1[12])+float(L1[16]), float(L1[13])+float(L1[17]), float(L1[14])+float(L1[18]), float(L1[15])+float(L1[19])]
					Alleles[base.index(L1[23])]=0
					maxAll=max(Alleles)		
					D[L1[0]] = maxAll		
				Total+=float(L1[22])
		f.close()
		F1=open(dirOUT + i + '_nonref.txt', 'w')
		F1.write('Position\tNonRef\n')
		for i in range (1, len(D)+1):
			key=str(i)
			F1.write('{}\t{}\n'.format(key,D[key]/Total))
		F1.close()
		
# make a file that will be used as input in R for the x-axis (nextera data)		
def prepareNX(args):
	files = args["<data>"]
	files1 = files.split(',')
	dirIN = 'DataFolderNX/Coverage/'
	os.system('mkdir DataFolder/NXDS_comparisonGraphs')
	dirOUT = 'DataFolder/NXDS_comparisonGraphs/'
	mut = args["--SNPfile"]
	
	base=['A','C','G','T']	
	
	NR={} # store the SNP positions as keys, and the non reference nucleotide as value
	Fs = open('Polymorphisms/' + mut, 'rU')
	for line in Fs:
		Ls = line.split()
		if Ls[2]=='SNP':
			NR[Ls[0]]=Ls[3][1]
	Fs.close()	
	
	for i in files1:
		D={} # store per position the non-reference allele count
		Total=0.0 # store total coverage to derive frequencies
		
		f=open(dirIN+i+'_Coverage.txt','rU')
		for line in f:
			if line.startswith("Position"):
				continue
			else:
				L1=line.split()
				if L1[0] in NR:
					k = base.index(NR[L1[0]]) # find the index of the non-ref allele in base
					D[L1[0]] = float(L1[k+1])
				else:
					Alleles = [float(L1[1]), float(L1[2]), float(L1[3]), float(L1[4])]
					Alleles[base.index(L1[5])]=0
					maxAll=max(Alleles)	
					D[L1[0]] = maxAll				
				for a in range (1,5):
					Total+=float(L1[a])
		f.close()
		F1=open(dirOUT + i + '_nonref.txt', 'w')
		F1.write('Position\tNonRef\n')
		for i in range (1, len(D)+1):
			key=str(i)
			F1.write('{}\t{}\n'.format(key,D[key]/Total))
		F1.close()
					
# Figure 4D
def dilutioncomparison(args):
	DS_files=['V1_S1','V5_S2','V25_S3','V125_S4','L50_S27_noPhiX','L10_S28_noPhiX','L5_S29_noPhiX','L1_S30_noPhiX','L0-5_S31_noPhiX']
	NX_files=['v1Nex1_S7','v1Nex2_S12','v4Nex1_S9','v4Nex2_S14','v9Nex1_S8','v9Nex2_S13','v19Nex1_DA1_S2','v19Nex2_DA2_S3','v24Nex1_DA7_S8','v24Nex2_DA8_S9','v99Nex1_DA3_S4','v99Nex2_DA4_S5','v124Nex1_DA9_S10','v124Nex2_DA10_S11','v199Nex1_DA5_S6','v199Nex2_DA6_S7']
	DS_dilutions = [2,5,25,125,2,10,20,100,200]
	NX_dilutions = [2,2,5,5,10,10,20,20,25,25,100,100,125,125,200,200]
	mut = args["--SNPfile"]
	outfile = args["--out"]
	dirIN ='DataFolder/NXDS_comparisonGraphs/'
	dirOUT = 'DataFolder/NXDS_comparisonGraphs/'
	
	mmPos={}
	Fs = open('Polymorphisms/' + mut, 'rU')
	for line in Fs:
		Ls = line.split()
		if Ls[2]=='SNP':
			mmPos[Ls[0]]='SNP' # dictionary to store the SNP positions as keys
	Fs.close()
	
	NonRef_Freq={} # store diffseq frequencies, nextera and averaged for each per dilution
	for i in range (0, len(DS_files)):
		f=open(dirIN+DS_files[i]+'_nonref.txt','rU')
		NonRef_Freq[DS_files[i]] = []
		for line in f:
			L1=line.split()
			if L1[0] in mmPos:
				NonRef_Freq[DS_files[i]].append(float(L1[1])) # for each file make an entry and store the SNP NR frequencies
		f.close()
		
	for i in range (0, len(NX_files)):
		f=open(dirIN+NX_files[i]+'_nonref.txt','rU')
		NonRef_Freq[NX_files[i]] = []
		for line in f:
			L1=line.split()
			if L1[0] in mmPos:
				NonRef_Freq[NX_files[i]].append(float(L1[1])) # for each file make an entry and store the SNP NR frequencies
		f.close()
	
	for key in NonRef_Freq:
		q=numpy.mean(NonRef_Freq[key])
		NonRef_Freq[key].append(q) # last entry for each key is average
	
	DilutionD={} # correspond each diffseq dataset to the dilution factor and the nx data with the same dilution factor
	for i in range (0, len(DS_files)):
		DilutionD[DS_files[i]] = [DS_dilutions[i],NonRef_Freq[DS_files[i]][len(NonRef_Freq[DS_files[i]])-1]] # dilution factor, average frequencies for dataset
		nxlist=[]
		for u in range (0,len(NX_files)):
			if NX_dilutions[u]==DS_dilutions[i]:
				nxlist.append(NonRef_Freq[NX_files[u]][len(NonRef_Freq[NX_files[u]])-1])

		DilutionD[DS_files[i]].append(numpy.mean(nxlist))

	F1=open(dirOUT + outfile+'.txt', 'w')
	F1.write('DilutionFactor\tNX\tDS\n')
	for key in DilutionD:
		F1.write('{}\t{}\t{}\n'.format(DilutionD[key][0], DilutionD[key][2], DilutionD[key][1]))
	F1.close()
	
	os.system('Rscript AllDilutionsNXDS.R ' + dirOUT + outfile+'.txt '+ dirOUT + outfile+'.pdf')

