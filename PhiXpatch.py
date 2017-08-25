import os
import sys

def preprocess(args):
    files = args["<data>"]
    files1 = files.split(',')
    os.system('mkdir DataFolder/MaskPhiX')
    dirIN = 'DataFolder/RawData/'
    dirOUT = 'DataFolder/MaskPhiX/'
    qualityCutoff = str(30)
    DAo82rc='AGCATGAGCGCTCGTCTCTGAAG' # auxiliary sequence
    R2rc='CTGTCTCTTATACACATCTCCGAGCCCACGAGAC' # Nextera read 2 RC
    ref = "PhiX"
	
    for f in files1:
		# merge reads    
    	os.system('software/FLASH-1.2.11/flash -M 100 -O -o '+f+' -d '+dirOUT+' '+dirIN+f+'_L001_R1_001.fastq '+dirIN+f+'_L001_R2_001.fastq')
    	# trim adaptors
    	os.system("cutadapt -q " + qualityCutoff  + " -g "+ DAo82rc +' -o '+dirOUT+f+'_trim82rc.fastq '+dirOUT+f+'.extendedFrags.fastq')
    	os.system("cutadapt -q " + qualityCutoff  + ' -a '+R2rc+' -o '+dirOUT+f+'_trimR2rc.fastq '+dirOUT+f+'_trim82rc.fastq')
    	# align
    	os.system('software/bowtie2-2.2.6/bowtie2 -x Reference/'+ref+'/bowtie2/'+ref+' '+dirOUT+f+'_trimR2rc.fastq -S '+dirOUT+f+'.sam')
    	F1=open(dirOUT+f+'.sam','rU')
    	mapped={} # dictionary where ID, read and q are stored
    	for line in F1:
    		L1=line.split()
    		if L1[1]!='4': # find those that mapped
    			mapped[L1[0]]=[] # hold the identifier read and q to find in fastq file
    	F1.close()
    	# eliminate reads that mapped to PhiX from the raw data files	
    	F3=open(dirIN+f+'_noPhiX_L001_R1_001.fastq','w')
    	F2=open(dirIN+f+'_L001_R1_001.fastq','rU')
    	for line in F2:
			if line.startswith ('@M00'):
				L1=line.split()
				L2=L1[0].replace('@','')
				if L2 not in mapped:
					F3.write(line)
					p=1 
				else:
					p=0
			else:
				if p==1: 
					F3.write(line) 
				else:
					continue
    	F3.close()
    	F2.close()
    	F3=open(dirIN+f+'_noPhiX_L001_R2_001.fastq','w')
    	F2=open(dirIN+f+'_L001_R2_001.fastq','rU') 
    	for line in F2:
			if line.startswith ('@M00'):
				L1=line.split()
				L2=L1[0].replace('@','')
				if L2 not in mapped:
					F3.write(line)
					p=1 
				else:
					p=0
			else:
				if p==1: 
					F3.write(line) 
				else:
					continue
    	F3.close()
    	F2.close()