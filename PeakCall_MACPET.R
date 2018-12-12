#############################################################################################
###############################      run all analysis for MACPET   ##########################
#############################################################################################
library(MACPET)
# Important: One has to specify a save directory and a data directory, the current module only gives information about how the models were ran.
#########################################################################
##################      TF data         #################################
#########################################################################
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------Data A1
# ------------------------------------------------------------------------------------
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000BZZ"
Repetition="REP1"
SA_AnalysisDir=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACPET_results")
if(!dir.exists(SA_AnalysisDir)) dir.create(SA_AnalysisDir)

snow <- BiocParallel::SnowParam(workers = 4, type = "SOCK", progressbar=FALSE)
BiocParallel::register(snow, default=TRUE)

MACPETUlt(SA_AnalysisDir=SA_AnalysisDir,
          SA_stages=c(0:3),
          SA_prefix="MACPET",
          S0_fastq1="YOURPATH/ChIA_PET_Data/Experiment_ENCSR000BZZ/REP1/fastq.files/ENCFF000KZA1.fastq.gz",
          S0_fastq2="YOURPATH/ChIA_PET_Data/Experiment_ENCSR000BZZ/REP1/fastq.files/ENCFF000KZI2.fastq.gz",
          S0_LinkerA="GTTGGATCAT",
          S0_LinkerB="GTTGGATCCG",
          S0_MinReadLength=18,
          S0_MaxReadLength=50,
          S0_LinkerOccurence=3,#The specific dataset has already many of its linkers removed, so this parameter has to be set at 3 (just like MANGO does/suggests)
          S0_image=TRUE,
          S0_fastqStream=2000000,
          S1_fastq1_usable_dir="",
          S1_fastq2_usable_dir="",
          S1_image=TRUE,
          S1_BAMStream=2000000,
          S1_makeSam=TRUE,
          S1_genome="hg19",
          S1_RbowtieIndexBuild=FALSE,
          S1_RbowtieIndexDir="YOURPATH/hg19.ebwt",
          S1_RbowtieIndexPrefix="hg19",
          S1_RbowtieRefDir="",
          S2_PairedEndBAMpath="",
          S2_image=TRUE,
          S2_BlackList=TRUE,
          S3_fileSelfDir="YOURPATH/ChIA_PET_Data/Experiment_ENCSR000BZZ/REP1/MACPET_results/S2_results/MACPET_pselfData",
          S3_image=TRUE,
          S3_method="BH")

rm(list=ls())
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------Data A2
#------------------------------------------------------------------------------------
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000CAD"
Repetition="REP1"
SA_AnalysisDir=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACPET_results")
if(!dir.exists(SA_AnalysisDir)) dir.create(SA_AnalysisDir)

snow <- BiocParallel::SnowParam(workers = 4, type = "SOCK", progressbar=FALSE)
BiocParallel::register(snow, default=TRUE)

MACPETUlt(SA_AnalysisDir=SA_AnalysisDir,
          SA_stages=c(0:3),
          SA_prefix="MACPET",
          S0_fastq1="YOURPATH/ChIA_PET_Data/Experiment_ENCSR000CAD/REP1/fastq.files/ENCFF000KYZ1.fastq.gz",
          S0_fastq2="YOURPATH/ChIA_PET_Data/Experiment_ENCSR000CAD/REP1/fastq.files/ENCFF000KYY2.fastq.gz",
          S0_LinkerA="GTTGGATAAG",
          S0_LinkerB="GTTGGAATGT",
          S0_MinReadLength=18,
          S0_MaxReadLength=50,
          S0_LinkerOccurence=0,
          S0_image=TRUE,
          S0_fastqStream=2000000,
          S1_fastq1_usable_dir="",
          S1_fastq2_usable_dir="",
          S1_image=TRUE,
          S1_BAMStream=2000000,
          S1_makeSam=TRUE,
          S1_genome="hg19",
          S1_RbowtieIndexBuild=FALSE,
          S1_RbowtieIndexDir="YOURPATH/hg19.ebwt",
          S1_RbowtieIndexPrefix="hg19",
          S1_RbowtieRefDir="",
          S2_PairedEndBAMpath="",
          S2_image=TRUE,
          S2_BlackList=TRUE,
          S3_fileSelfDir="YOURPATH/ChIA_PET_Data/Experiment_ENCSR000CAD/REP1/MACPET_results/S2_results/MACPET_pselfData",
          S3_image=TRUE,
          S3_method="BH")

rm(list=ls())
# #------------------------------------------------------------------------------------
# #------------------------------------------------------------------------------------Data A3
# #------------------------------------------------------------------------------------
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000CAC"
Repetition="REP1"
SA_AnalysisDir=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACPET_results")
if(!dir.exists(SA_AnalysisDir)) dir.create(SA_AnalysisDir)

snow <- BiocParallel::SnowParam(workers = 4, type = "SOCK", progressbar=FALSE)
BiocParallel::register(snow, default=TRUE)

MACPETUlt(SA_AnalysisDir=SA_AnalysisDir,
          SA_stages=c(0:3),
          SA_prefix="MACPET",
          S0_fastq1="YOURPATH/ChIA_PET_Data/Experiment_ENCSR000CAC/REP1/fastq.files/ENCFF000KYB1.fastq.gz",
          S0_fastq2="YOURPATH/ChIA_PET_Data/Experiment_ENCSR000CAC/REP1/fastq.files/ENCFF000KYQ2.fastq.gz",
          S0_LinkerA="GTTGGATAAG",
          S0_LinkerB="GTTGGAATGT",
          S0_MinReadLength=18,
          S0_MaxReadLength=50,
          S0_LinkerOccurence=0,
          S0_image=TRUE,
          S0_fastqStream=2000000,
          S1_fastq1_usable_dir="",
          S1_fastq2_usable_dir="",
          S1_image=TRUE,
          S1_BAMStream=2000000,
          S1_makeSam=TRUE,
          S1_genome="hg19",
          S1_RbowtieIndexBuild=FALSE,
          S1_RbowtieIndexDir="YOURPATH/hg19.ebwt",
          S1_RbowtieIndexPrefix="hg19",
          S1_RbowtieRefDir="",
          S2_PairedEndBAMpath="",
          S2_image=TRUE,
          S2_BlackList=TRUE,
          S3_fileSelfDir="YOURPATH/ChIA_PET_Data/Experiment_ENCSR000CAC/REP1/MACPET_results/S2_results/MACPET_pselfData",
          S3_image=TRUE,
          S3_method="BH")

rm(list=ls())
#########################################################################
##################      Histones data   #################################
#########################################################################

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------ Data 4
#------------------------------------------------------------------------------------
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000FDD"
Repetition="REP1"
SA_AnalysisDir=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACPET_results")
if(!dir.exists(SA_AnalysisDir)) dir.create(SA_AnalysisDir)

snow <- BiocParallel::SnowParam(workers = 4, type = "SOCK", progressbar=FALSE)
BiocParallel::register(snow, default=TRUE)

MACPETUlt(SA_AnalysisDir=SA_AnalysisDir,
          SA_stages=c(0:3),
          SA_prefix="MACPET",
          S0_fastq1="YOURPATH/ChIA_PET_Data/Experiment_ENCSR000FDD/REP1/fastq.files/ENCFF002ACB1_correct.fastq.gz",
          S0_fastq2="YOURPATH/ChIA_PET_Data/Experiment_ENCSR000FDD/REP1/fastq.files/ENCFF002ACC2_correct.fastq.gz",
          S0_LinkerA="GTTGGATAAG",
          S0_LinkerB="GTTGGAATGT",
          S0_MinReadLength=18,
          S0_MaxReadLength=50,
          S0_LinkerOccurence=0,
          S0_image=TRUE,
          S0_fastqStream=2000000,
          S1_fastq1_usable_dir="",
          S1_fastq2_usable_dir="",
          S1_image=TRUE,
          S1_BAMStream=2000000,
          S1_makeSam=TRUE,
          S1_genome="hg19",
          S1_RbowtieIndexBuild=FALSE,
          S1_RbowtieIndexDir="YOURPATH/hg19.ebwt",
          S1_RbowtieIndexPrefix="hg19",
          S1_RbowtieRefDir="",
          S2_PairedEndBAMpath="",
          S2_image=TRUE,
          S2_BlackList=TRUE,
          S3_fileSelfDir="YOURPATH/ChIA_PET_Data/Experiment_ENCSR000FDD/REP1/MACPET_results/S2_results/MACPET_pselfData",
          S3_image=TRUE,
          S3_method="BH")

rm(list=ls())
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------ Data 5
#------------------------------------------------------------------------------------
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000FDG"
Repetition="REP1"
SA_AnalysisDir=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACPET_results")
if(!dir.exists(SA_AnalysisDir)) dir.create(SA_AnalysisDir)

snow <- BiocParallel::SnowParam(workers = 4, type = "SOCK", progressbar=FALSE)
BiocParallel::register(snow, default=TRUE)

MACPETUlt(SA_AnalysisDir=SA_AnalysisDir,
          SA_stages=c(0:3),
          SA_prefix="MACPET",
          S0_fastq1="YOURPATH/ChIA_PET_Data/Experiment_ENCSR000FDG/REP1/fastq.files/ENCFF002ACQ1_correct.fastq.gz",
          S0_fastq2="YOURPATH/ChIA_PET_Data/Experiment_ENCSR000FDG/REP1/fastq.files/ENCFF002ACR2_correct.fastq.gz",
          S0_LinkerA="GTTGGATAAG",
          S0_LinkerB="GTTGGAATGT",
          S0_MinReadLength=18,
          S0_MaxReadLength=50,
          S0_LinkerOccurence=0,
          S0_image=TRUE,
          S0_fastqStream=2000000,
          S1_fastq1_usable_dir="",
          S1_fastq2_usable_dir="",
          S1_image=TRUE,
          S1_BAMStream=2000000,
          S1_makeSam=TRUE,
          S1_genome="hg19",
          S1_RbowtieIndexBuild=FALSE,
          S1_RbowtieIndexDir="YOURPATH/hg19.ebwt",
          S1_RbowtieIndexPrefix="hg19",
          S1_RbowtieRefDir="",
          S2_PairedEndBAMpath="",
          S2_image=TRUE,
          S2_BlackList=TRUE,
          S3_fileSelfDir="YOURPATH/ChIA_PET_Data/Experiment_ENCSR000FDG/REP1/MACPET_results/S2_results/MACPET_pselfData",
          S3_image=TRUE,
          S3_method="BH")

rm(list=ls())
#########################################################################
##################      POL 2 data      #################################
#########################################################################
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------ Data 6
#------------------------------------------------------------------------------------
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000BZY"
Repetition="REP1"
SA_AnalysisDir=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACPET_results")
if(!dir.exists(SA_AnalysisDir)) dir.create(SA_AnalysisDir)

snow <- BiocParallel::SnowParam(workers = 4, type = "SOCK", progressbar=FALSE)
BiocParallel::register(snow, default=TRUE)

MACPETUlt(SA_AnalysisDir=SA_AnalysisDir,
          SA_stages=c(0:3),
          SA_prefix="MACPET",
          S0_fastq1="YOURPATH/ChIA_PET_Data/Experiment_ENCSR000BZY/REP1/fastq.files/ENCFF000KYG1.fastq.gz",
          S0_fastq2="YOURPATH/ChIA_PET_Data/Experiment_ENCSR000BZY/REP1/fastq.files/ENCFF000KYK2.fastq.gz",
          S0_LinkerA="GTTGGATAAG",
          S0_LinkerB="GTTGGAATGT",
          S0_MinReadLength=18,
          S0_MaxReadLength=50,
          S0_LinkerOccurence=0,
          S0_image=TRUE,
          S0_fastqStream=2000000,
          S1_fastq1_usable_dir="",
          S1_fastq2_usable_dir="",
          S1_image=TRUE,
          S1_BAMStream=2000000,
          S1_makeSam=TRUE,
          S1_genome="hg19",
          S1_RbowtieIndexBuild=FALSE,
          S1_RbowtieIndexDir="YOURPATH/hg19.ebwt",
          S1_RbowtieIndexPrefix="hg19",
          S1_RbowtieRefDir="",
          S2_PairedEndBAMpath="",
          S2_image=TRUE,
          S2_BlackList=TRUE,
          S3_fileSelfDir="YOURPATH/ChIA_PET_Data/Experiment_ENCSR000BZY/REP1/MACPET_results/S2_results/MACPET_pselfData",
          S3_image=TRUE,
          S3_method="BH")

rm(list=ls())
#########################################################################
################## long range chia pet        ###########################
# IMPORTANT: this data is huge, so a machine with a lot of RAM has to be used or else R will crush, note: this is NOT a problem of MACPET rather than a problem of memory handling in R.
#########################################################################
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------ Data 7
#------------------------------------------------------------------------------------
DATA.type="LR_ChIA-PET_Data"
Experiment="Experiment_GSM1872886"
Repetition="REP1"
SA_AnalysisDir=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACPET_results")
if(!dir.exists(SA_AnalysisDir)) dir.create(SA_AnalysisDir)

snow <- BiocParallel::SnowParam(workers = 4, type = "SOCK", progressbar=FALSE)
BiocParallel::register(snow, default=TRUE)

MACPETUlt(SA_AnalysisDir=SA_AnalysisDir,
          SA_stages=c(0:3),
          SA_prefix="MACPET",
          S0_fastq1="YOURPATH/LR_ChIA-PET_Data/Experiment_GSM1872886/REP1/fastq.files/SRR2312566_pass_1_correct.fastq.gz",
          S0_fastq2="YOURPATH/LR_ChIA-PET_Data/Experiment_GSM1872886/REP1/fastq.files/SRR2312566_pass_2_correct.fastq.gz",
          S0_LinkerA="ACGCGATATCTTATC",
          S0_LinkerB="AGTCAGATAAGATAT",
          S0_MinReadLength=18,
          S0_MaxReadLength=1000,
          S0_LinkerOccurence=3,# since the data is with tegmentation, this should be set at 3 (like MANGO suggests)
          S0_image=TRUE,
          S0_fastqStream=2000000,
          S1_fastq1_usable_dir="",
          S1_fastq2_usable_dir="",
          S1_image=TRUE,
          S1_BAMStream=2000000,
          S1_makeSam=TRUE,
          S1_genome="hg19",
          S1_RbowtieIndexBuild=FALSE,
          S1_RbowtieIndexDir="YOURPATH/hg19.ebwt",
          S1_RbowtieIndexPrefix="hg19",
          S1_RbowtieRefDir="",
          S2_PairedEndBAMpath="",
          S2_image=TRUE,
          S2_BlackList=TRUE,
          S3_fileSelfDir="",
          S3_image=TRUE,
          S3_method="BH")

rm(list=ls())