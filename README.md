Diffseq version 1.0 7.5.2017

software:
bowtie2-2.2.6
bwa-0.7.15
cutadapt-1.10
flash-1.2.11
picard-2.7.1
R

python dependencies:
math
numpy
docopt
pysam

from the directory that the diffseq.py is in run:

To eliminate sequences that derive from PhiX intending for internal control
$ python diffseq.py removePhiX L0_S26,L50_S27,L10_S28,L5_S29,L1_S30,L0-5_S31,L0-1_S32,L0-05_S33,L0-01_S34

To merge paired-end reads
$ python diffseq.py merge 1kb-1_S21,1kb-2_S22,1kb-3_S35,1kb-4_S23,1kb-5_S24,1kb-6_S36,1kb-7_S25,1kb-8_S26,L0_S26_noPhiX,L50_S27_noPhiX,L10_S28_noPhiX,L5_S29_noPhiX,L1_S30_noPhiX,L0-5_S31_noPhiX,L0-1_S32_noPhiX,L0-05_S33_noPhiX,L0-01_S34_noPhiX,V1_S1,V5_S2,V25_S3,V125_S4,M4150_S4,M4150_S11,V4240_S3,V4240_S10,v5248_S1,v5248V4240_S2

To consider one read per multiplicate
$ python diffseq.py dedup 1kb-1_S21,1kb-2_S22,1kb-3_S35,1kb-4_S23,1kb-5_S24,1kb-6_S36,1kb-7_S25,1kb-8_S26,L0_S26_noPhiX,L50_S27_noPhiX,L10_S28_noPhiX,L5_S29_noPhiX,L1_S30_noPhiX,L0-5_S31_noPhiX,L0-1_S32_noPhiX,L0-05_S33_noPhiX,L0-01_S34_noPhiX,V1_S1,V5_S2,V25_S3,V125_S4,M4150_S4,M4150_S11,V4240_S3,V4240_S10,v5248_S1,v5248V4240_S2

To trim adaptors and other sequences
$ python diffseq.py trimadaptors 1kb-1_S21,1kb-2_S22,1kb-3_S35,1kb-4_S23,1kb-5_S24,1kb-6_S36,1kb-7_S25,1kb-8_S26,L0_S26_noPhiX,L50_S27_noPhiX,L10_S28_noPhiX,L5_S29_noPhiX,L1_S30_noPhiX,L0-5_S31_noPhiX,L0-1_S32_noPhiX,L0-05_S33_noPhiX,L0-01_S34_noPhiX,V1_S1,V5_S2,V25_S3,V125_S4,M4150_S4,M4150_S11,V4240_S3,V4240_S10,v5248_S1,v5248V4240_S2

To filter out reads that do not have the expected structure according to the molecular biology
$ python diffseq.py filter 1kb-1_S21,1kb-2_S22,1kb-3_S35,1kb-4_S23,1kb-5_S24,1kb-6_S36,1kb-7_S25,1kb-8_S26,L0_S26_noPhiX,L50_S27_noPhiX,L10_S28_noPhiX,L5_S29_noPhiX,L1_S30_noPhiX,L0-5_S31_noPhiX,L0-1_S32_noPhiX,L0-05_S33_noPhiX,L0-01_S34_noPhiX,V1_S1,V5_S2,V25_S3,V125_S4,M4150_S4,M4150_S11,V4240_S3,V4240_S10,v5248_S1,v5248V4240_S2

To map to the reference genome
$ python diffseq.py align 1kb-1_S21,1kb-2_S22,1kb-3_S35,1kb-5_S24,1kb-6_S36,1kb-7_S25,1kb-8_S26 --reference PET17
$ python diffseq.py align 1kb-4_S23 --reference PET17_1
$ python diffseq.py align L0_S26_noPhiX,L50_S27_noPhiX,L10_S28_noPhiX,L5_S29_noPhiX,L1_S30_noPhiX,L0-5_S31_noPhiX,L0-1_S32_noPhiX,L0-05_S33_noPhiX,L0-01_S34_noPhiX,V1_S1,V5_S2,V25_S3,V125_S4 --reference v36-13
$ python diffseq.py align M4150_S4,M4150_S11 --reference M4150
$ python diffseq.py align V4240_S3,V4240_S10 --reference V4240
$ python diffseq.py align v5248_S1,v5248V4240_S2 --reference v5248

To make coverage files to be used in plots and further analysis
$ python diffseq.py coverage 1kb-1_S21,1kb-2_S22,1kb-3_S35,1kb-5_S24,1kb-6_S36,1kb-7_S25,1kb-8_S26 --reference PET17
$ python diffseq.py coverage 1kb-4_S23 --reference PET17_1
$ python diffseq.py coverage L0_S26_noPhiX,L50_S27_noPhiX,L10_S28_noPhiX,L5_S29_noPhiX,L1_S30_noPhiX,L0-5_S31_noPhiX,L0-1_S32_noPhiX,L0-05_S33_noPhiX,L0-01_S34_noPhiX,V1_S1,V5_S2,V25_S3,V125_S4 --reference v36-13
$ python diffseq.py coverage M4150_S4,M4150_S11 --reference M4150
$ python diffseq.py coverage V4240_S3,V4240_S10 --reference V4240
$ python diffseq.py coverage v5248_S1,v5248V4240_S2 --reference v5248

To make coverage and diff-seq plots for 1 kb (figures 2, S2 and S3)
$ python diffseq.py coverageplots 1kb-1_S21,1kb-2_S22,1kb-3_S35,1kb-5_S24,1kb-6_S36,1kb-7_S25,1kb-8_S26 --reference PET17
$ python diffseq.py coverageplots 1kb-4_S23 --reference PET17_1

To make stacked diff-seq coverage plots (figure 3A)
$ python diffseq.py stacked L50_S27_noPhiX,L10_S28_noPhiX,L5_S29_noPhiX,L1_S30_noPhiX,L0-5_S31_noPhiX,L0-1_S32_noPhiX,L0-05_S33_noPhiX,L0-01_S34_noPhiX,L0_S26_noPhiX --out Clonal_vDNA_seriesL --reference v36-13 --SNPfile viralClones.txt

To make biases heatmap (figures 3B, S4)
$ python diffseq.py heatmap L5_S29_noPhiX --reference v36-13 --SNPfile viralClones.txt

To make differential diffseq coverage plots with control as denominator (Figure 3C and D)
$ python diffseq.py diffsignal L5_S29_noPhiX --control L0_S26_noPhiX --reference v36-13 --SNPfile viralClones.txt --percentage 5

To make ROCs of different parameters and their ranks (figures 3E, S5)
$ python diffseq.py model L5_S29_noPhiX --reference v36-13 --SNPfile viralClones.txt
$ Rscript ModelParameters.R DataFolder/Model/L5_S29_noPhiXCandidateParameters.txt DataFolder/Model/L5_S29_noPhiXCandidateParameters.pdf
$ Rscript rankZscores.R DataFolder/Model/L5_S29_noPhiXCandidateParameters.txt DataFolder/Model/L5_S29_noPhiX_rankedZscores.pdf Polymorphisms/viralClones.txt

To make ROCs of the same parameter for different data (figure 3F,G)
$ python diffseq.py model M4150_S4,M4150_S11 --reference M4150 --SNPfile M4150.txt
$ python diffseq.py model V4240_S3,V4240_S10 --reference V4240 --SNPfile V4240.txt
$ python diffseq.py model v5248_S1,v5248V4240_S2 --reference v5248 --SNPfile V4240.txt
$ Rscript Model.R DataFolder/Model/V4240population.pdf DataFolder/Model/V4240_S3CandidateParameters.txt 8 9 12 DataFolder/Model/V4240_S10CandidateParameters.txt DataFolder/Model/v5248V4240_S2CandidateParameters.txt DataFolder/Model/v5248_S1CandidateParameters.txt replicate1 replicate2 diluted control
$ Rscript Model.R DataFolder/Model/M4150population.pdf DataFolder/Model/M4150_S4CandidateParameters.txt 6 7 8 DataFolder/Model/M4150_S11CandidateParameters.txt replicate1 replicate2

from the directory that the nextera.py is in run:
To trim sequencing adaptors
$ python nextera.py trim v1Nex1_S7,v1Nex2_S12,v4Nex1_S9,v4Nex2_S14,v9Nex1_S8,v9Nex2_S13,v19Nex1_DA1_S2,v19Nex2_DA2_S3,v24Nex1_DA7_S8,v24Nex2_DA8_S9,v99Nex1_DA3_S4,v99Nex2_DA4_S5,v124Nex1_DA9_S10,v124Nex2_DA10_S11,v199Nex1_DA5_S6,v199Nex2_DA6_S7
$ python nextera.py trim M4150Nex1_S6,M4150Nex2_S11,V4240Nex1_S5,V4240Nex2_S10

To map to the reference
$ python nextera.py align v1Nex1_S7,v1Nex2_S12,v4Nex1_S9,v4Nex2_S14,v9Nex1_S8,v9Nex2_S13,v19Nex1_DA1_S2,v19Nex2_DA2_S3,v24Nex1_DA7_S8,v24Nex2_DA8_S9,v99Nex1_DA3_S4,v99Nex2_DA4_S5,v124Nex1_DA9_S10,v124Nex2_DA10_S11,v199Nex1_DA5_S6,v199Nex2_DA6_S7 --reference v36-13
$ python nextera.py align M4150Nex1_S6,M4150Nex2_S11 --reference M4150
$ python nextera.py align V4240Nex1_S5,V4240Nex2_S10 --reference V4240

To make coverage files
$ python nextera.py coverage v1Nex1_S7,v1Nex2_S12,v4Nex1_S9,v4Nex2_S14,v9Nex1_S8,v9Nex2_S13,v19Nex1_DA1_S2,v19Nex2_DA2_S3,v24Nex1_DA7_S8,v24Nex2_DA8_S9,v99Nex1_DA3_S4,v99Nex2_DA4_S5,v124Nex1_DA9_S10,v124Nex2_DA10_S11,v199Nex1_DA5_S6,v199Nex2_DA6_S7 --reference v36-13

To compare Diffseq to Nextera data (figures 4A-C, S6)
$ python diffseq.py prepareDS L50_S27_noPhiX,L10_S28_noPhiX,L5_S29_noPhiX,L1_S30_noPhiX,L0-5_S31_noPhiX,V1_S1,V5_S2,V25_S3,V125_S4 --reference v36-13 --SNPfile viralClones.txt
$ python diffseq.py prepareNX v1Nex1_S7,v1Nex2_S12,v4Nex1_S9,v4Nex2_S14,v9Nex1_S8,v9Nex2_S13,v19Nex1_DA1_S2,v19Nex2_DA2_S3,v24Nex1_DA7_S8,v24Nex2_DA8_S9,v99Nex1_DA3_S4,v99Nex2_DA4_S5,v124Nex1_DA9_S10,v124Nex2_DA10_S11,v199Nex1_DA5_S6,v199Nex2_DA6_S7 --SNPfile viralClones.txt
$ Rscript DiffSeqNextera.R DataFolder/NXDS_comparisonGraphs/V1_S1_nonref.txt DataFolder/NXDS_comparisonGraphs/L50_S27_noPhiX_nonref.txt DataFolder/NXDS_comparisonGraphs/v1Nex1_S7_nonref.txt  DataFolder/NXDS_comparisonGraphs/v1Nex2_S12_nonref.txt Polymorphisms/viralClones.txt DataFolder/NXDS_comparisonGraphs/1to1.pdf 1:1
$ Rscript DiffSeqNextera.R DataFolder/NXDS_comparisonGraphs/V5_S2_nonref.txt DataFolder/NXDS_comparisonGraphs/V5_S2_nonref.txt DataFolder/NXDS_comparisonGraphs/v4Nex1_S9_nonref.txt  DataFolder/NXDS_comparisonGraphs/v4Nex2_S14_nonref.txt Polymorphisms/viralClones.txt DataFolder/NXDS_comparisonGraphs/1to4.pdf 1:4
$ Rscript DiffSeqNextera.R DataFolder/NXDS_comparisonGraphs/L10_S28_noPhiX_nonref.txt DataFolder/NXDS_comparisonGraphs/L10_S28_noPhiX_nonref.txt DataFolder/NXDS_comparisonGraphs/v9Nex1_S8_nonref.txt  DataFolder/NXDS_comparisonGraphs/v9Nex2_S13_nonref.txt Polymorphisms/viralClones.txt DataFolder/NXDS_comparisonGraphs/1to9.pdf 1:9
$ Rscript DiffSeqNextera.R DataFolder/NXDS_comparisonGraphs/L5_S29_noPhiX_nonref.txt DataFolder/NXDS_comparisonGraphs/L5_S29_noPhiX_nonref.txt DataFolder/NXDS_comparisonGraphs/v19Nex1_DA1_S2_nonref.txt  DataFolder/NXDS_comparisonGraphs/v19Nex2_DA2_S3_nonref.txt Polymorphisms/viralClones.txt DataFolder/NXDS_comparisonGraphs/1to19.pdf 1:19
$ Rscript DiffSeqNextera.R DataFolder/NXDS_comparisonGraphs/V25_S3_nonref.txt DataFolder/NXDS_comparisonGraphs/V25_S3_nonref.txt DataFolder/NXDS_comparisonGraphs/v24Nex1_DA7_S8_nonref.txt  DataFolder/NXDS_comparisonGraphs/v24Nex2_DA8_S9_nonref.txt Polymorphisms/viralClones.txt DataFolder/NXDS_comparisonGraphs/1to24.pdf 1:24
$ Rscript DiffSeqNextera.R DataFolder/NXDS_comparisonGraphs/L1_S30_noPhiX_nonref.txt DataFolder/NXDS_comparisonGraphs/L1_S30_noPhiX_nonref.txt DataFolder/NXDS_comparisonGraphs/v99Nex1_DA3_S4_nonref.txt  DataFolder/NXDS_comparisonGraphs/v99Nex2_DA4_S5_nonref.txt Polymorphisms/viralClones.txt DataFolder/NXDS_comparisonGraphs/1to99.pdf 1:99
$ Rscript DiffSeqNextera.R DataFolder/NXDS_comparisonGraphs/V125_S4_nonref.txt DataFolder/NXDS_comparisonGraphs/V125_S4_nonref.txt DataFolder/NXDS_comparisonGraphs/v124Nex1_DA9_S10_nonref.txt  DataFolder/NXDS_comparisonGraphs/v124Nex2_DA10_S11_nonref.txt Polymorphisms/viralClones.txt DataFolder/NXDS_comparisonGraphs/1to124.pdf 1:124
$ Rscript DiffSeqNextera.R DataFolder/NXDS_comparisonGraphs/L0-5_S31_noPhiX_nonref.txt DataFolder/NXDS_comparisonGraphs/L0-5_S31_noPhiX_nonref.txt DataFolder/NXDS_comparisonGraphs/v199Nex1_DA5_S6_nonref.txt  DataFolder/NXDS_comparisonGraphs/v199Nex2_DA6_S7_nonref.txt Polymorphisms/viralClones.txt DataFolder/NXDS_comparisonGraphs/1to199.pdf 1:199

To make the last panel of figure 4 (4D)
$ python diffseq.py dilutioncomparison --out NXDS_comparison --SNPfile viralClones.txt

To organize data for table 1
$ python diffseq.py posfreq --data L50_S27_noPhiX --SNPfile viralClones.txt --reference v36-13
