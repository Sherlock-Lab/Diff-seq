import os
import sys
import math
import numpy

# The following makes a file with diff-seq coverage for all libraries to be used for stacked graphs (figure 3A)
def stacked(args):
	files = args["<data>"]
	files1 = files.split(',')
	outFile = args["--out"]
	ref = args["--reference"]
	SNPfile = args["--SNPfile"]
	dirIN='DataFolder/Coverage'+ref+'/'
	os.system('mkdir DataFolder/StackedGraphs')
	dirOUT = 'DataFolder/StackedGraphs/'
	DS={}
	txt=open('Reference/'+ref+'/'+ref+'.txt','rU')
	S=txt.read()
	txt.close()
	RefLen=len(S)
	F5=open(dirOUT+outFile+'_DiffSeqStacked.txt','w')
	F5.write('Position\tFrequency\tLibrary\n')
	
	Tot={} # Here store per library the total coverage, so as to derive frequencies for the common file
	for f in files1:
		Tot[f]=0.0
		F0=open(dirIN+f+'_Coverage.txt','rU')
		for line in F0:
			if line.startswith('Position'):
				continue
			else:
				L1=line.split()
				if L1[0] not in DS:
					DS[L1[0]]=[]
				DS[L1[0]].append(float(L1[22]))
				Tot[f]+=float(L1[22])
								

	for i in range (1,RefLen+1):
		key=str(i)
		for u in range (0, len(files1)):
			F5.write('{}\t{}\t{}\n'.format(key,DS[key][u]/Tot[files1[u]],str(u+1)+'_'+files1[u]))
	F5.close()
	os.system('Rscript StackedGraphs.R Polymorphisms/'+SNPfile+' '+dirOUT+outFile+'_DiffSeqStacked.txt '+dirOUT+outFile+'_DSstackedGraphs.pdf ' + str(len(files1)-1) + ' 2668')

# figure 3B
def heatmap(args):
	files = args["<data>"]
	files1 = files.split(',')
	ref = args["--reference"]
	dirIN='DataFolder/Coverage'+ref+'/'
	os.system('mkdir DataFolder/Heatmaps')
	dirOUT = 'DataFolder/Heatmaps/'
	mut = args["--SNPfile"]
	mmPos={}

	Fs = open('Polymorphisms/' + mut, 'rU')
	for line in Fs:
		Ls = line.split()
		if Ls[2]=='SNP':
			mmPos[int(Ls[0])]=Ls[3]
	Fs.close()

	for f in files1:
		RfAr={} # store per position the ratio (reference fwd coverage + alternative reverse coverage)/total coverage at position for mismatch positions for all libraries
		Af={} # store per position the ratio alternative fwd coverage/total fwd coverage at position for mismatch positions for all libraries
		Ar={} # store per position the ratio alternative rvs coverage/total rvs coverage at position for mismatch positions for all libraries
		Afr={} # store per position the AVERAGE of the ratio alternative rvs coverage/total rvs coverage and alternative fwd coverage/total fwd coverage at position for mismatch positions for all libraries
		TAfr={} # store per position the weighted AVERAGE (instead of Af/2+Ar/2 do (Altfwd+Altrvs)/(Ref+Alt)
		maxAfr={} # store the maximum value
		Coverage={} # store the total coverage at position for mismatch positions for all libraries
		ID={} # store the Ref fwd, alternative forward and trinucleotide per position
	
		txt=open('Reference/'+ref+'/'+ref+'.txt','rU')
		S=txt.read()
		txt.close()
	
		for i in range (2, len(S)):
			if i in mmPos:
				ID[i]=[]
				ID[i].append(mmPos[i][0]) # reference allele
				ID[i].append(mmPos[i][1]) # alternative allele
				ID[i].append(S[i-2:i+1]) # trinucleotide
	
		REF_D={}
		for i in mmPos:
			REF_D[i]=[]
			for k in range (0,8): # store A, C, G, T fwd and rvs
				REF_D[i].append(0.0)
		
		Total=0.0	
		base=['A','C','G','T']
	
		F1=open(dirIN+f+'_Coverage.txt','rU')
		
		for line in F1:
			if line.startswith('Position'):
				continue
			else:
				L1=line.split()
				for i in range (12,20):
					Total+=float(L1[i])
				i=int(L1[0])
				if i in REF_D:
					for i2 in range (12,20):
						REF_D[i][i2-12]=float(L1[i2])
	
		for key in REF_D:
			indexRfwd=base.index(ID[key][0]) # base index of the reference allele
			indexAfwd=base.index(ID[key][1]) # base index of the alt allele
			if REF_D[key][indexRfwd]+REF_D[key][indexAfwd]>0 and REF_D[key][indexRfwd+4]+REF_D[key][indexAfwd+4]>0:
				RfAr[key]=(REF_D[key][indexRfwd]+REF_D[key][indexAfwd+4])/(REF_D[key][indexRfwd]+REF_D[key][indexRfwd+4]+REF_D[key][indexAfwd]+REF_D[key][indexAfwd+4])
				Af[key]=REF_D[key][indexAfwd]/(REF_D[key][indexRfwd]+REF_D[key][indexAfwd])
				Ar[key]=REF_D[key][indexAfwd+4]/(REF_D[key][indexRfwd+4]+REF_D[key][indexAfwd+4])
				Afr[key]=(Af[key]+Ar[key])/2
				TAfr[key]=(REF_D[key][indexAfwd]+REF_D[key][indexAfwd+4])/(REF_D[key][indexRfwd]+REF_D[key][indexAfwd]+REF_D[key][indexRfwd+4]+REF_D[key][indexAfwd+4])
				maxAfr[key]=max(Af[key],Ar[key])
			elif REF_D[key][indexRfwd]+REF_D[key][indexAfwd]>0 and REF_D[key][indexRfwd+4]+REF_D[key][indexAfwd+4]==0:
				RfAr[key]=(REF_D[key][indexRfwd]+REF_D[key][indexAfwd+4])/(REF_D[key][indexRfwd]+REF_D[key][indexRfwd+4]+REF_D[key][indexAfwd]+REF_D[key][indexAfwd+4])
				Af[key]=REF_D[key][indexAfwd]/(REF_D[key][indexRfwd]+REF_D[key][indexAfwd])
				Ar[key]='NA'
				Afr[key]=Af[key]/2
				TAfr[key]=Af[key]
				maxAfr[key]=Af[key]		
			elif REF_D[key][indexRfwd]+REF_D[key][indexAfwd]==0 and REF_D[key][indexRfwd+4]+REF_D[key][indexAfwd+4]>0:
				RfAr[key]=(REF_D[key][indexRfwd]+REF_D[key][indexAfwd+4])/(REF_D[key][indexRfwd]+REF_D[key][indexRfwd+4]+REF_D[key][indexAfwd]+REF_D[key][indexAfwd+4])
				Af[key]='NA'
				Ar[key]=REF_D[key][indexAfwd+4]/(REF_D[key][indexRfwd+4]+REF_D[key][indexAfwd+4])
				Afr[key]=Ar[key]/2
				TAfr[key]=Ar[key]
				maxAfr[key]=Ar[key]
			else:
				RfAr[key]='NA'
				Af[key]='NA'
				Ar[key]='NA'
				Afr[key]='NA'
				TAfr[key]='NA'
				maxAfr[key]='NA'
			Coverage[key]=(REF_D[key][indexRfwd]+REF_D[key][indexRfwd+4]+REF_D[key][indexAfwd]+REF_D[key][indexAfwd+4])/Total

		Fm=open(dirOUT+f+'Heatmap.txt','w')
		Fm.write('Position\tRefAlt\ttriN\tMismatchBias\tAltFWD\tAltRVS\tAlterAvgBias\tAlterTotBias\tAlterMaxBias\tLog2CoverageFreq1\tLog2CoverageFreq2\tLog2CoverageFreq3\tLog2CoverageFreq4\tLog2CoverageFreq5\tLog2CoverageFreq6\n')
		for key in ID:
			Fm.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(key,ID[key][0]+ID[key][1],ID[key][2],RfAr[key],Af[key],Ar[key],Afr[key],TAfr[key],maxAfr[key],math.log(Coverage[key],2),math.log(Coverage[key],2),math.log(Coverage[key],2),math.log(Coverage[key],2),math.log(Coverage[key],2),math.log(Coverage[key],2)))
		Fm.close()
		
		os.system('Rscript Heatmaps.R '+dirOUT+f+'Heatmap.txt '+dirOUT+f+'_Heatmaps.pdf')
		
# Figure 3C, D	
def diffsignal(args):
	files = args["<data>"]
	files1 = files.split(',')
	f0 = args["--control"]
	ref = args["--reference"]
	SNPs = args["--SNPfile"]
	percentage = args["--percentage"]
	per1 = percentage.split(',')
	dirIN='DataFolder/Coverage'+ref+'/'
	os.system('mkdir DataFolder/DiffSignal')
	dirOUT = 'DataFolder/DiffSignal/'
	minDS={} # store per position minor allele ds coverage for library and control library
	DS={} # store per position total ds coverage for library and control library
	Total = [0.0,0.0] # store total ds coverage for library and control library to derive frequencies
	
	txt=open('Reference/'+ref+'/'+ref+'.txt','rU')
	S=txt.read()
	txt.close()
	for f in files1:
		F1=open(dirIN+f+'_Coverage.txt','rU')
		for line in F1:
			if line.startswith('Position'):
				continue
			else:
				L1=line.split()
				DS[L1[0]]=[float(L1[22])]
				Total[0]+=float(L1[22])
				Alleles = [float(L1[12])+float(L1[16]), float(L1[13])+float(L1[17]), float(L1[14])+float(L1[18]), float(L1[15])+float(L1[19])]
				maxAll=max(Alleles)
				Alleles[Alleles.index(maxAll)] = 0
				minDS[L1[0]] = [max(Alleles)]
		F1.close()
		
		Fc=open(dirIN+f0+'_Coverage.txt','rU')
		for line in Fc:
			if line.startswith('Position'):
				continue
			else:
				L1=line.split()
				DS[L1[0]].append(float(L1[22]))
				Total[1]+=float(L1[22])
				Alleles = [float(L1[12])+float(L1[16]), float(L1[13])+float(L1[17]), float(L1[14])+float(L1[18]), float(L1[15])+float(L1[19])]
				maxAll=max(Alleles)
				Alleles[Alleles.index(maxAll)] = 0
				minDS[L1[0]].append(max(Alleles))
		Fc.close()
			
		F6=open(dirOUT+f+'diffDiffSeqFreq.txt','w')
		F6.write('Position\tminDS\tDS\tminDS0\tDS0\n')
		for i in range (1, len(S)+1):
			key=str(i)
			F6.write('{}\t{}\t{}\t{}\t{}\n'.format(key,minDS[key][0]/Total[0],DS[key][0]/Total[0],minDS[key][1]/Total[1],DS[key][1]/Total[1]))
		F6.close()
		
		os.system('Rscript DifferentialSignal.R ' + dirOUT + f + 'diffDiffSeqFreq.txt Polymorphisms/' + SNPs + ' ' + dirOUT + f + '_DiffSignal.pdf "log2 ("' + per1[files1.index(f)] + '"% rare variant/no rare variant)"')
		
# Figure 3E-G, S5
def model(args):
	files = args["<data>"]
	files1 = files.split(',')
	ref = args["--reference"]
	mut = args["--SNPfile"]
	dirIN='DataFolder/Coverage'+ref+'/'
	os.system('mkdir DataFolder/Model')
	dirOUT = 'DataFolder/Model/'
	
	mmPos={}
	Fs = open('Polymorphisms/' + mut, 'rU')
	for line in Fs:
		Ls = line.split()
		if Ls[2]=='SNP':
			mmPos[int(Ls[0])]='SNP'
	Fs.close()
	
	txt=open('Reference/'+ref+'/'+ref+'.txt','rU')
	S=txt.read()
	txt.close()
		
	for f in files1:
		F1=open(dirIN+f+'_Coverage.txt','rU')
		REF_D={} # store diff-seq coverage
		for line in F1:
			if line.startswith('Position'):
				continue
			else:
				L1=line.split()
				REF_D[int(L1[0])]=[]
				for i in range (12,20):
					REF_D[int(L1[0])].append(float(L1[i]))
		F1.close()		
	
		base=['A','C','G','T']
	
	# Find per position the major and minor allele
		Alleles={}
		for key in REF_D:
			TotalPos=REF_D[key][0]+REF_D[key][4]+REF_D[key][1]+REF_D[key][5]+REF_D[key][2]+REF_D[key][6]+REF_D[key][3]+REF_D[key][7]
			if TotalPos>0:
				# FIND the major and minor alleles
				AllelesDS=[REF_D[key][0]+REF_D[key][4],REF_D[key][1]+REF_D[key][5],REF_D[key][2]+REF_D[key][6],REF_D[key][3]+REF_D[key][7]]
				major=max(AllelesDS) # Major Coverage of the sum FWD and RVS
				MAJi=AllelesDS.index(major) 
				majorAllele=base[MAJi] # Major Allele	
				AllelesDS[MAJi]=0
				minor=max(AllelesDS) # Minor Coverage FWD and RVS
				if minor!=0:
					MINi=AllelesDS.index(minor)
					minorAllele=base[MINi] # Minor Allele	
				else:
					minorAllele='-'	
			else:
				majorAllele='-' 
				minorAllele='-'
			Alleles[key]=[]
			Alleles[key].append(majorAllele) 
			Alleles[key].append(minorAllele) 

	#  z-scores derived when the whole diff-seq coverage is considered
		TotalD={}
		TotalD['fwd']=[]
		TotalD['rvs']=[]
		for key in REF_D:
			TotalPosF=REF_D[key][0]+REF_D[key][1]+REF_D[key][2]+REF_D[key][3]
			TotalPosR=REF_D[key][4]+REF_D[key][5]+REF_D[key][6]+REF_D[key][7]
			TotalD['fwd'].append(TotalPosF)
			TotalD['rvs'].append(TotalPosR)
		MeanT={}
		StDevT={}
		for key in TotalD:
			MeanT[key]=numpy.mean(TotalD[key])
			StDevT[key]=numpy.std(TotalD[key])	
		
		TotalStat={}
		for key in REF_D:
			TotalStat[key]=[]
			TotalPosF=REF_D[key][0]+REF_D[key][1]+REF_D[key][2]+REF_D[key][3]
			TotalPosR=REF_D[key][4]+REF_D[key][5]+REF_D[key][6]+REF_D[key][7]
			bf=(TotalPosF-MeanT['fwd'])/StDevT['fwd']
			TotalStat[key].append(bf) # z-score
			br=(TotalPosR-MeanT['rvs'])/StDevT['rvs']
			TotalStat[key].append(br)

	#  z-scores derived when the minor allele diff-seq coverage is considered	
		MinTD={}
		MinTD['fwd']=[]
		MinTD['rvs']=[]
		for key in REF_D:
			if Alleles[key][1] in base:
				TotalPosF=REF_D[key][base.index(Alleles[key][1])]
				TotalPosR=REF_D[key][base.index(Alleles[key][1])+4]
				MinTD['fwd'].append(TotalPosF)
				MinTD['rvs'].append(TotalPosR)
		MeanT={}
		StDevT={}
		for key in MinTD:
			MeanT[key]=numpy.mean(MinTD[key])
			StDevT[key]=numpy.std(MinTD[key])
		MinorTotalStat={}
		for key in REF_D:
			MinorTotalStat[key]=[]
			if Alleles[key][1] in base:
				a=base.index(Alleles[key][1])
				MinorF=REF_D[key][a]
				MinorR=REF_D[key][a+4]
			else:
				MinorF=0.0
				MinorR=0.0
			bf=(MinorF-MeanT['fwd'])/StDevT['fwd']
			MinorTotalStat[key].append(bf) # z-score
			br=(MinorR-MeanT['rvs'])/StDevT['rvs']
			MinorTotalStat[key].append(br)

	#  z-scores derived when the major allele diff-seq coverage is considered per major allele ID
		Base1={}
		for key in REF_D:
			if Alleles[key][0]!='-':
				if Alleles[key][0]+'f' not in Base1:
					Base1[Alleles[key][0]+'f']=[]
					Base1[Alleles[key][0]+'r']=[]
				Base1[Alleles[key][0]+'f'].append(REF_D[key][base.index(Alleles[key][0])])
				Base1[Alleles[key][0]+'r'].append(REF_D[key][base.index(Alleles[key][0])+4])
		MeanT={}
		StDevT={}
		for key in Base1:
			MeanT[key]=numpy.mean(Base1[key])
			StDevT[key]=numpy.std(Base1[key])	
		MajorStat={}
		for key in REF_D:
			MajorStat[key]=[]
			for i in range (0,4):
				b=(REF_D[key][i]-MeanT[base[i]+'f'])/StDevT[base[i]+'f']
				MajorStat[key].append(b)
			for i in range (4,8):
				b=(REF_D[key][i]-MeanT[base[i-4]+'r'])/StDevT[base[i-4]+'r']
				MajorStat[key].append(b)
	
	# z-scores derived when the minor allele diff-seq coverage is considered per major and minor allele ID
		Base2={}
		for key in REF_D:
			if Alleles[key][0]!='-':
				for n in base:
					if Alleles[key][0]!=n:
						if Alleles[key][0]+n+'f' not in Base2:
							Base2[Alleles[key][0]+n+'f']=[]
							Base2[Alleles[key][0]+n+'r']=[]
						Base2[Alleles[key][0]+n+'f'].append(REF_D[key][base.index(Alleles[key][0])])
						Base2[Alleles[key][0]+n+'r'].append(REF_D[key][base.index(Alleles[key][0])+4])
		MeanT={}
		StDevT={}
		for key in Base2:
			MeanT[key]=numpy.mean(Base2[key])
			StDevT[key]=numpy.std(Base2[key])
		MinorStat={}
		for key in REF_D:
			MinorStat[key]=[]
			for i in range (0,4):
				if Alleles[key][0]!=base[i] and Alleles[key][0]!='-':
					b=(REF_D[key][i]-MeanT[Alleles[key][0]+base[i]+'f'])/StDevT[Alleles[key][0]+base[i]+'f']
					MinorStat[key].append(b)
				else: 
					MinorStat[key].append(-10000)
			for i in range (4,8):
				if Alleles[key][0]!=base[i-4] and Alleles[key][0]!='-':
					b=(REF_D[key][i]-MeanT[Alleles[key][0]+base[i-4]+'r'])/StDevT[Alleles[key][0]+base[i-4]+'r']
					MinorStat[key].append(b)
				else: 
					MinorStat[key].append(-10000)
				
	# z-scores derived when the trinucleotide ID is considered
		TriN={}
		for key in range (2,len(REF_D)):
			if Alleles[key-1][0]!='-':
				U=Alleles[key-1][0]
			if Alleles[key-1][0]=='-':
				U=S[key-2]
			if Alleles[key+1][0]!='-':
				D=Alleles[key+1][0]
			if Alleles[key+1][0]=='-':
				D=S[key]
			if Alleles[key][0]!='-':
				P=Alleles[key][0]
			if Alleles[key][0]=='-':
				P=S[key-1]
			if U+P+D+'f' not in TriN:
				TriN[U+P+D+'f']=[]
				TriN[U+P+D+'r']=[]
			TriN[U+P+D+'f'].append(REF_D[key][base.index(P)])
			TriN[U+P+D+'r'].append(REF_D[key][base.index(P)+4])
		MeanT={}
		StDevT={}
		for key in TriN:
			if len(TriN[key])>1:
				MeanT[key]=numpy.mean(TriN[key])
				StDevT[key]=numpy.std(TriN[key])		
		TriNStat={}
		for key in range (2,len(REF_D)):
			TriNStat[key]=[]
			if Alleles[key-1][0]!='-':
				U=Alleles[key-1][0]
			if Alleles[key-1][0]=='-':
				U=S[key-2]
			if Alleles[key+1][0]!='-':
				D=Alleles[key+1][0]
			if Alleles[key+1][0]=='-':
				D=S[key]
			for i in range (0,4):
				if U+base[i]+D+'f' in MeanT:
					b=(REF_D[key][i]-MeanT[U+base[i]+D+'f'])/StDevT[U+base[i]+D+'f']
					TriNStat[key].append(b)
				else:
					TriNStat[key].append(-10000)		
			for i in range (4,8):
				if U+base[i-4]+D+'r' in MeanT:
					b=(REF_D[key][i]-MeanT[U+base[i-4]+D+'r'])/StDevT[U+base[i-4]+D+'r']
					TriNStat[key].append(b)
				else:
					TriNStat[key].append(-10000)
				
		F3=open(dirOUT+f+'CandidateParameters.txt','w')
		F3.write('Position\tTotfZ\tTotrZ\tTotMinfZ\tTotMinrZ\tMajfZ\tMajrZ\tMinfZ\tMinrZ\tTriN1f\tTriN1r\tTriN2f\tTriN2r\tStatus\n')
		for key in range (2,len(REF_D)):
			F3.write('{}\t{}\t{}\t{}\t{}\t'.format(key,TotalStat[key][0],TotalStat[key][1],MinorTotalStat[key][0],MinorTotalStat[key][1]))			
			Plist=[MajorStat[key][0],MajorStat[key][1],MajorStat[key][2],MajorStat[key][3],MajorStat[key][4],MajorStat[key][5],MajorStat[key][6],MajorStat[key][7]]
			K=max(Plist)
			if Plist.index(K) in range (0,4):
				F3.write('{}\t{}\t'.format(K,MajorStat[key][Plist.index(K)+4]))
			if Plist.index(K) in range (4,8):
				F3.write('{}\t{}\t'.format(MajorStat[key][Plist.index(K)-4],K))			
			Plist=[MinorStat[key][0],MinorStat[key][1],MinorStat[key][2],MinorStat[key][3],MinorStat[key][4],MinorStat[key][5],MinorStat[key][6],MinorStat[key][7]]
			K=max(Plist)
			if K!=-10000:
				if Plist.index(K) in range (0,4):
					F3.write('{}\t{}\t'.format(K,MinorStat[key][Plist.index(K)+4]))
				if Plist.index(K) in range (4,8):
					F3.write('{}\t{}\t'.format(MinorStat[key][Plist.index(K)-4],K))
			else:
				F3.write('NA\tNA\t')
			Plist=[TriNStat[key][0],TriNStat[key][1],TriNStat[key][2],TriNStat[key][3],TriNStat[key][4],TriNStat[key][5],TriNStat[key][6],TriNStat[key][7]]
			TriNMax=max(Plist)
			if TriNMax!=-10000:
				if Plist.index(TriNMax) in range (0,4):
					a=Plist.index(TriNMax)
					b=Plist.index(TriNMax)+4
					F3.write('{}\t{}\t'.format(TriNMax,TriNStat[key][b]))
				if Plist.index(TriNMax) in range (4,8):
					a=Plist.index(TriNMax)-4
					b=Plist.index(TriNMax)
					F3.write('{}\t{}\t'.format(TriNStat[key][a],TriNMax))
				for u in [a,b]:
					Plist[u]=-10000
				TriNMin=max(Plist)
				if TriNMin!=-10000:
					if Plist.index(TriNMin) in range (0,4):
						F3.write('{}\t{}\t'.format(TriNMin,TriNStat[key][Plist.index(TriNMin)+4]))
					if Plist.index(TriNMin) in range (4,8):
						F3.write('{}\t{}\t'.format(TriNStat[key][Plist.index(TriNMin)-4],TriNMin))
				else:
					F3.write('NA\tNA\t')
			else:
				F3.write('NA\tNA\tNA\tNA\t')
			if key in mmPos:
				F3.write('SNP\n')
			else:
				F3.write('match\n')
		F3.close()		

# Table 1
def posfreq(args):
	files = args["<data>"]
	files1 = files.split(',')
	mut = args["--SNPfile"]
	ref = args["--reference"]
	dirIN='DataFolder/Coverage'+ref+'/'
	os.system('mkdir DataFolder/Table1')
	dirOUT = 'DataFolder/Table1/'
	
	mmPos={}
	Fs = open('Polymorphisms/' + mut, 'rU')
	for line in Fs:
		Ls = line.split()
		if Ls[2]=='SNP':
			mmPos[Ls[0]]=Ls[3] # correspond the SNP positions to their ref-alt
	Fs.close()
	
	for f in files1:
		F1=open(dirIN+f+'_Coverage.txt','rU')
		REF_D={} # store diff-seq coverage
		for line in F1:
			ds_coverage=0.0
			L1=line.split()
			if L1[0] in mmPos:
				REF_D[L1[0]]=[]
				for i in range (12,20):
					REF_D[L1[0]].append(float(L1[i]))
					ds_coverage+=float(L1[i])
				REF_D[L1[0]].append(ds_coverage)
		F1.close()
		F2=open(dirOUT+f+'_Table1freq.txt','w')
		F2.write('Position\tAf\tCf\tGf\tTf\tAr\tCr\tGr\tTr\tREFALT\n')
		for key in REF_D:
			F2.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(key,REF_D[key][0]/REF_D[key][8],REF_D[key][1]/REF_D[key][8],REF_D[key][2]/REF_D[key][8],REF_D[key][3]/REF_D[key][8],REF_D[key][4]/REF_D[key][8],REF_D[key][5]/REF_D[key][8],REF_D[key][6]/REF_D[key][8],REF_D[key][7]/REF_D[key][8],mmPos[key][0]+'/'+mmPos[key][1]))
		F2.close()

		os.system('Rscript Table1.R ' + dirOUT + f + '_Table1freq.txt '+ dirOUT + f + '_Table1freq.pdf')
	
		F3=open(dirOUT+f+'_Table1counts.txt','w')
		F3.write('Position\tAf\tCf\tGf\tTf\tAr\tCr\tGr\tTr\tREFALT\n')
		for key in REF_D:
			F3.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(key,REF_D[key][0],REF_D[key][1],REF_D[key][2],REF_D[key][3],REF_D[key][4],REF_D[key][5],REF_D[key][6],REF_D[key][7],mmPos[key][0]+'|'+mmPos[key][1]))
		F3.close()

		os.system('Rscript Table1.R ' + dirOUT + f + '_Table1counts.txt '+ dirOUT + f + '_Table1counts.pdf')
