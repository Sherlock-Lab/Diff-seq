import os
import numpy as np
import pysam

def trim(args):
	read1Adaptor = 'CTGTCTCTTATACACATCTCCGAGCCCACGAGAC'
	read2Adaptor = 'CTGTCTCTTATACACATCTGACGCTGCCGACGA'
	
	files = args["<data>"]
	files1 = files.split(',')
	os.system('mkdir DataFolderNX/Trimmed')
	dirIN = 'DataFolderNX/rawData/'
	dirOUT = 'DataFolderNX/Trimmed/'
	for f in files1:
		os.system("cutadapt -q 30 --minimum-length 12 -u 10 -a " + read1Adaptor + " -o " + dirOUT + f + "_tmp.1.fastq -p " + dirOUT + f + "_tmp.2.fastq " + dirIN + f + '_L001_R1_001.fastq' + " " + dirIN + f + '_L001_R2_001.fastq' + " >> " + dirOUT + f +"_output_data.txt")
		os.system("cutadapt -q 30 --minimum-length 12 -u 10 -a " + read2Adaptor + " -o " + dirOUT + f + "_trimmedR2.fastq -p " + dirOUT + f + "_trimmedR1.fastq " + dirOUT + f + "_tmp.2.fastq " + dirOUT + f + "_tmp.1.fastq >> " + dirOUT + f + "_output_data.txt")

def align(args):
	files = args["<data>"]
	files1 = files.split(',')
	ref = args["--reference"]
	os.system('mkdir DataFolderNX/Aligned')
	dirIN = 'DataFolderNX/Trimmed/'
	dirOUT = 'DataFolderNX/Aligned/'
	for f in files1:
		readGroupInfo = '@RG\\tID:group1\\tSM:' + f + '\\tPL:illumina\\tLB:12_08_2014\\tPU:' + f
		BWArefIndex = 'Reference/' + ref + '/bwa/' + ref + '.fsa'
		os.system('software/bwa-0.7.15/bwa mem -M -R "' + readGroupInfo + '" ' + BWArefIndex + ' ' + dirIN + f + '_trimmedR1.fastq ' + dirIN + f + '_trimmedR2.fastq > ' + dirOUT + f + '.sam')
		#Sort SAM file (saves as bam file)
		os.system('java -Xmx2g -jar software/picard.jar SortSam INPUT=' + dirOUT + f + '.sam OUTPUT=' + dirOUT + f + '.bam SORT_ORDER=coordinate >> ' + dirOUT + f + "_output_data.txt")
		#Build BAM index
		os.system('java -Xmx2g -jar software/picard.jar BuildBamIndex INPUT= ' + dirOUT + f + '.bam >> ' + dirOUT + f + "_output_data.txt")
	
def coverage(args):
	files = args["<data>"]
	files1 = files.split(',')
	ref = args["--reference"]
	os.system('mkdir DataFolderNX/Coverage')
	dirIN = 'DataFolderNX/Aligned/'
	dirOUT = 'DataFolderNX/Coverage/'
	
	txt=open('Reference/'+ref+'/'+ref+'.txt','rU')
	S=txt.read()
	txt.close()
 	base=['A','C','G','T','N']
 	
	for f in files1:
		samfile = pysam.AlignmentFile(dirIN + f + '.bam', "rb")

		REF_D={} # store coverage
		for i in range (0,len(S)):
			REF_D[i]=[]
			for k in range (0,5): # store A, C, G, T, N
				REF_D[i].append(0.0)
		
		for pileupcolumn in samfile.pileup(ref):
			for pileupread in pileupcolumn.pileups:
				if not pileupread.is_del and not pileupread.is_refskip:
					x=base.index(pileupread.alignment.query_sequence[pileupread.query_position])
					REF_D[pileupcolumn.pos][x]+=1

		samfile.close()
							
		F3=open(dirOUT + f + '_Coverage.txt','w')
		F3.write('Position\tA\tC\tG\tT\tREF\n')
		
		for key in REF_D:
			F3.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(key+1,REF_D[key][0],REF_D[key][1],REF_D[key][2],REF_D[key][3],S[key]))	

		F3.close()