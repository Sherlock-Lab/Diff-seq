import os
import sys
import pysam

# merges paired-end reads
def merge(args):
    files = args["<data>"]
    files1 = files.split(',')
    os.system('mkdir DataFolder/Merged')
    dirIN = 'DataFolder/RawData/'
    dirOUT = 'DataFolder/Merged/'
    for f in files1:
    	print f
    	os.system('flash -M 100 -O -o '+f+' -d '+dirOUT+' '+dirIN+f+'_L001_R1_001.fastq '+dirIN+f+'_L001_R2_001.fastq')

# based on the whole read keeps one per multiplicate
def dedup(args):
	files = args["<data>"]
	files1 = files.split(',')
	dirIN='DataFolder/Merged/'
	os.system('mkdir DataFolder/Dedup')
	dirOUT='DataFolder/Dedup/'
	for f in files1:
		print f
		F1=open(dirIN+f+'.extendedFrags.fastq','rU')
		F2=open(dirOUT+f+'_Dedup.fastq','w')
		F3=open(dirOUT+f+'duplicateCount.txt','w')
		F4=open(dirOUT+f+'duplicate_hist.txt','w')
	# 	Initiate a dictionary to count how many times each read appears
		RD={}
		unique=0
		total=0
	# 	Read the fastq file
		S1=F1.read()
		L1=S1.split('\n')
		for i in range (1,len(L1),4):
			if L1[i] not in RD:
				RD[L1[i]]=1
				F2.write('{}\n{}\n{}\n{}\n'.format(L1[i-1],L1[i],L1[i+1],L1[i+2]))
				unique+=1
				total+=1
			else:
				total+=1
				RD[L1[i]]+=1
		fraction=unique*1.0/total
		F3.write('Total='+str(total)+'\nUnique='+ str(unique)+'\nUnique fraction='+ str(fraction))

		F4.write('Read\t# it appears\n')
		for key in RD:
			F4.write('{}\t{}\n'.format(key,RD[key]))
		F1.close()
		F2.close()
		F3.close()
		F4.close()

# trims adaptors, UMIs and adjacent sequences
def trimadaptors(args):
	files = args["<data>"]
	files1 = files.split(',')
	qualityCutoff = str(30)
	DAo82rc='AGCATGAGCGCTCGTCTCTGAAG' # auxiliary sequence
	R2rc='CTGTCTCTTATACACATCTCCGAGCCCACGAGAC' # Nextera read 2 RC
	dirIN='DataFolder/Dedup/'
	os.system('mkdir DataFolder/AdaptorsTrimmed')
	dirOUT='DataFolder/AdaptorsTrimmed/'
	for f in files1:
		print f
		os.system("cutadapt -q " + qualityCutoff  + " -g "+ DAo82rc +' -o '+dirOUT+f+'_trim82rc.fastq '+dirIN+f+'_Dedup.fastq')
		os.system("cutadapt -q " + qualityCutoff  + ' -a '+R2rc+' -o '+dirOUT+f+'_trimR2rc.fastq '+dirOUT+f+'_trim82rc.fastq')

# Filters out reads that do not start with a G and are reads less than 8 nucleotides long
def filter(args):
	files = args["<data>"]
	files1 = files.split(',')
	dirIN='DataFolder/AdaptorsTrimmed/'
	os.system('mkdir DataFolder/Filtered')
	dirOUT='DataFolder/Filtered/'
	for f in files1:
		print f
		F1=open(dirIN+f+'_trimR2rc.fastq','rU')
		F2=open(dirOUT+f+'_Filtered.fastq','w')
		Sf=F1.read()
		Lf=Sf.split('\n')
		for i in range (1,len(Lf),4): 
			l=len(Lf[i])
			if l>7:
				if Lf[i][0]=='G': # keep and trim the G
					F2.write(Lf[i-1]+'\n'+Lf[i][1:]+'\n'+Lf[i+1]+'\n'+Lf[i+2][1:]+'\n')
				else: 
					continue
			else:
				continue
		F1.close()
		F2.close()

# Aligns to reference	
def align(args):
	files = args["<data>"]
	files1 = files.split(',')
	ref = args["--reference"]
	dirIN='DataFolder/Filtered/'
	os.system('mkdir DataFolder/Aligned'+ref+'')
	dirOUT='DataFolder/Aligned'+ref+'/'
	for f in files1:
		print f
		os.system('software/bowtie2-2.2.6/bowtie2 -x Reference/'+ref+'/bowtie2/'+ref+' '+dirIN+f+'_Filtered.fastq -S '+dirOUT+f+'_Filtered.sam')
		
# split sam file in 2 for forward and reverse reads
		F1=open(dirOUT+f+'_Filtered.sam','rU')
		F2=open(dirOUT+f+'_Forward.sam','w')
		F3=open(dirOUT+f+'_Reverse.sam','w')
		#CHANGE THE M00653 to if it does not start with @
		for line in F1:
			if line.startswith('@'):
				F2.write(line)
				F3.write(line)
			else:
				Ls=line.split()
				if Ls[1]=='0':
					F2.write(line)
				elif Ls[1]=='16':
					F3.write(line)
				else:
					continue

		F1.close()
		F2.close()
		F3.close()
				
		#Sort SAM file (saves as bam file)
		os.system('java -Xmx2g -jar software/picard.jar SortSam INPUT=' + dirOUT + f + '_Forward.sam OUTPUT=' + dirOUT + f + '_Forward.bam SORT_ORDER=coordinate >> ' + dirOUT + f + "_Forward_output_data.txt")
		os.system('java -Xmx2g -jar software/picard.jar SortSam INPUT=' + dirOUT + f + '_Reverse.sam OUTPUT=' + dirOUT + f + '_Reverse.bam SORT_ORDER=coordinate >> ' + dirOUT + f + "_Reverse_output_data.txt")
		os.system('java -Xmx2g -jar software/picard.jar SortSam INPUT=' + dirOUT + f + '_Filtered.sam OUTPUT=' + dirOUT + f + '_Filtered.bam SORT_ORDER=coordinate >> ' + dirOUT + f + "_Filtered_output_data.txt")

		#Build BAM index
		os.system('java -Xmx2g -jar software/picard.jar BuildBamIndex INPUT= ' + dirOUT + f + '_Forward.bam >> ' + dirOUT + f + "_Forward_output_data.txt")
		os.system('java -Xmx2g -jar software/picard.jar BuildBamIndex INPUT= ' + dirOUT + f + '_Reverse.bam >> ' + dirOUT + f + "_Reverse_output_data.txt")
		os.system('java -Xmx2g -jar software/picard.jar BuildBamIndex INPUT= ' + dirOUT + f + '_Filtered.bam >> ' + dirOUT + f + "_Filtered_output_data.txt")

def coverage(args):
	files = args["<data>"]
	files1 = files.split(',')
	ref = args["--reference"]
	os.system('mkdir DataFolder/Coverage'+ref+'')
	dirIN = 'DataFolder/Aligned'+ref+'/'
	dirOUT = 'DataFolder/Coverage'+ref+'/'
	
	txt=open('Reference/'+ref+'/'+ref+'.txt','rU')
	S=txt.read()
	txt.close()
 	base=['A','C','G','T','N']
 	base2=['A','C','G','T']
 	
	for f in files1:
		print f
		REF_D={} # store coverage
		dsREF_D={} # store diffseq coverage
		for i in range (0,len(S)):
			REF_D[i]=[]
			dsREF_D[i]=[]
			for k in range (0,10): # store A, C, G, T, N fwd and rvs
				REF_D[i].append(0.0)
				dsREF_D[i].append(0.0)		
		
		samfilef = pysam.AlignmentFile(dirIN + f + '_Forward.bam', "rb")
		samfiler = pysam.AlignmentFile(dirIN + f + '_Reverse.bam', "rb")
		
# 		Make the REF_D for the total coverage
		for pileupcolumn in samfilef.pileup(ref):
			for pileupread in pileupcolumn.pileups:
				if not pileupread.is_del and not pileupread.is_refskip:
					x=base.index(pileupread.alignment.query_sequence[pileupread.query_position])
					REF_D[pileupcolumn.pos][x]+=1			
		for pileupcolumn in samfiler.pileup(ref):
			for pileupread in pileupcolumn.pileups:
				if not pileupread.is_del and not pileupread.is_refskip:
					x=base.index(pileupread.alignment.query_sequence[pileupread.query_position])
					REF_D[pileupcolumn.pos][x+5]+=1
		
		samfilef = pysam.AlignmentFile(dirIN + f + '_Forward.bam', "rb")
		samfiler = pysam.AlignmentFile(dirIN + f + '_Reverse.bam', "rb")
# 		Make dsREF_D for Diffseq coverage			
		for read in samfilef:
			if read.query_alignment_start==0:
				u=[0,1,2]
			elif read.query_alignment_start==1:
				u=[-1,0,1]
			elif read.query_alignment_start==2:
				u=[-2,-1,0]
			elif read.query_alignment_start==3:
				u=[-3,-2,-1]
			else:
				print read
				
			for i in [0,1,2]:
				if read.seq[i] in base2:
					B=base.index(read.seq[i])
					dsREF_D[read.reference_start+u[i]][B]+=1
					if read.seq[i]!='G': 
						break
						
		for read in samfiler:
			if read.query_alignment_end==(read.query_length-read.query_alignment_start):
				u=[-1,-2,-3]
			elif read.query_alignment_end==(read.query_length-read.query_alignment_start-1):
				u=[0,-1,-2]
			elif read.query_alignment_end==(read.query_length-read.query_alignment_start-2):
				u=[1,0,-1]
			elif read.query_alignment_end==(read.query_length-read.query_alignment_start-3):
				u=[2,1,0]
			else:
				print read
			
			iter=[read.query_length-1,read.query_length-2,read.query_length-3]	
			for i in [0,1,2]:
				if read.reference_end+u[i] in dsREF_D:
					if read.seq[iter[i]] in base2:
						B=base.index(read.seq[iter[i]])
						dsREF_D[read.reference_end+u[i]][B+5]+=1
						if read.seq[iter[i]]!='C': 
							break
						
		samfilef.close()
		samfiler.close()
		
		F4=open(dirOUT+f+'_Coverage.txt','w')
		F4.write('Position\tTotAFWD\tTotCFWD\tTotGFWD\tTotTFWD\tTotARVS\tTotCRVS\tTotGRVS\tTotTRVS\tTotFWD\tTotRVS\tTot\tdsAFWD\tdsCFWD\tdsGFWD\tdsTFWD\tdsARVS\tdsCRVS\tdsGRVS\tdsTRVS\tdsFWD\tdsRVS\tds\tREF\n')

		for key in REF_D:
			F4.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(key+1,REF_D[key][0],REF_D[key][1],REF_D[key][2],REF_D[key][3],REF_D[key][5],REF_D[key][6],REF_D[key][7],REF_D[key][8], REF_D[key][0]+REF_D[key][1]+REF_D[key][2]+REF_D[key][3],REF_D[key][5]+REF_D[key][6]+REF_D[key][7]+REF_D[key][8], REF_D[key][0]+REF_D[key][1]+REF_D[key][2]+REF_D[key][3]+REF_D[key][5]+REF_D[key][6]+REF_D[key][7]+REF_D[key][8],dsREF_D[key][0],dsREF_D[key][1],dsREF_D[key][2],dsREF_D[key][3],dsREF_D[key][5],dsREF_D[key][6],dsREF_D[key][7],dsREF_D[key][8], dsREF_D[key][0]+dsREF_D[key][1]+dsREF_D[key][2]+dsREF_D[key][3],dsREF_D[key][5]+dsREF_D[key][6]+dsREF_D[key][7]+dsREF_D[key][8], dsREF_D[key][0]+dsREF_D[key][1]+dsREF_D[key][2]+dsREF_D[key][3]+dsREF_D[key][5]+dsREF_D[key][6]+dsREF_D[key][7]+dsREF_D[key][8], S[key]))

		F4.close()

def coverageplots(args):
	files = args["<data>"]
	files1 = files.split(',')
	ref = args["--reference"]
	dirIN = 'DataFolder/Coverage'+ref+'/'
	
	for f in files1:
		os.system('Rscript CoveragePlots.R ' + dirIN + f + '_Coverage.txt ' + dirIN + f +'_TotalCoverageFrequencies.pdf ' + dirIN + f + '_DiffSeqCoverageFrequencies.pdf')
		
