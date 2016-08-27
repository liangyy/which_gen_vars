##
## The following makefile controls the pipeline described in 
## Moyerbrailean et al (http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005875)
## Wayne State University, 2016
##

## Required software
## * samtools
## * bedtools
## * UCSC utils (http://hgdownload.soe.ucsc.edu/admin/exe/)

## Custom utils
## * scanPwm
## * scanPwmVar
## * makeXbFile
## * bedXbFileXtract

# Folders
dnaseFolder  = dnase/
seedFolder   = pwms/seed/
t5kBedFolder = $(seedFolder)pwmScan.sites/
t5kSeqFolder = $(seedFolder)pwmScan.seq/
pwmFolder    = pwms/recal/pwmFiles.new/
selFolder    = pwms/selPwms/
initCenti    = initial/
tmp          = /dev/shm/ # Any temporary storage directory

# Files
1kGenomes     = genome/variants.txt.gz # Tab separated, chr / pos / ref allele / alt allele
genome        = genome/hg19.fa    # Genome in fasta format
genome2bit    = genome/hg19.2bit  # Genome in 2bit format
genometmp2bit = $(tmp)/hg19.2bit
chromSizes    = genome/chromSizes.hg19.txt # From UCSC

# Scripts
lambdaScript    = src/parseCentipedeTop5K.Lambda.v2.gmb.R
pwmScript       = src/parseCentipedeTop5K.Pwm.v2.gmb.R
betaScript      = src/parseCentipedeTop5K.AllBeta.v5.gmb.R
initCentiScript = src/runCentipede.R
fullCentiScript = src/runCentipedeOnAll.R

# Misc
awkDnaseStr = '{if(and($$2,0x10)==0){print $$3,$$4-1,"+"};if(and($$2,0x10)==16){print $$3,$$4-1+length($$10)-1,"-"}}'

# Parameters
ncpus       = 10
zThresh     = 5
numJobs     = 60
numJobsFull = 42

# Helper commands
.SECONDARY:

%.Folder:
	mkdir -p $*;

$(genometmp2bit): $(genome2bit)
	cp $(genome2bit) $(genometmp2bit)

%.Qsub64:
	touch $@
	echo "cd $${PWD}; make $*" | qsub -q masxq -l nodes=1:ppn=64 -o $@ -e $@.e

%.Qsub12:
	touch $@
	echo "cd $${PWD}; make $*" | qsub -q mmtxq -l nodes=1:ppn=12 -o $@ -e $@.e


######################
## Initial PWM scan ##
######################

## Using the seed motifs, generate the training sites using 'scanPwm'
## 1) scan pwm human, choose top-scoring matches
## 2) scan pwm chimp, liftover top-scoring matches (not in [1])
## 3) scan pwm macaque, liftover top-scoring matches (not in [1,2])
## 4) combine sequences from [1,2,3], use for next step

# Extract the sequences at the Top5K sites
$(seedFolder)pwmScan.seq/%.seq.gz: $(t5kBedFolder)/%.bed.gz $(genometmp2bit)
	less $< | seqBedFor2bit -noMask -window=20 $(genometmp2bit) stdin stdout | gzip > $@
all: ${patsubst $(pwmFolder)/%.bed.gz,%.seq.gz,${wildcard $(pwmFolder)/*.bed.gz}}


############################
## Initial CENTIPEDE scan ##
############################

# Convert DNase-seq data to x8b files
$(dnaseFolder)%.x8b.bz2: $(dnaseFolder)%.bam
	samtools view -q10 $^ | awk -F"\t" $(awkDnaseStr) | makeXbFile $(chromSizes) -readSize=0 stdin $(tmp)/$*.x8b -verbose=2
	pbzip2 -b1m2000vp$${NCPUS:-1} $(tmp)/$*.x8b
	mv $(tmp)/$*.x8b.bz2 $@
allDnase: $(patsubst $(dnaseFolder)%.bam, $(dnaseFolder)%.x8b.bz2.Qsub12, $(wildcard $(dnaseFolder)*.bam))

# Get the coverage for each sample
$(dnaseFolder)numReads.txt:
	grep '... completed with ' *.Qsub12.e | sed 's/:.*with /\t/g' | sed 's/ reads//g' > $@

# Run CENTIPEDE for a Tissue-TF pair (see below)
$(initCenti)$(expName)/%.out.gz:
	-R --vanilla --args $* $(expName) ./$(initCenti) < $(initCentiScript) | gzip > $@ 2> $@.err
$(expName).mk: $(tmp)/$(expName).x8b $(patsubst $(t5kBedFolder)/%.bed.gz, $(initCenti)$(expName)/%.out.gz, $(wildcard $(t5kBedFolder)/*.bed.gz))

# Submit a tissue to run on its own node (preferred method)
$(initCenti)%.exp: $(dnaseFolder)/%.x8b.bz2
	mkdir -p $(initCenti)$*
	cp $^ $(tmp)
	pbzip2 -dm1000vp$(numJobs) $(tmp)/$*.x8b.bz2
	make -j $(numJobs) $*.mk expName=$*
	rm -f $(tmp)$*.x8b
initQsub: $(patsubst $(dnaseFolder)/%.x8b.bz2, $(initCenti)%.exp.Qsub64, $(wildcard $(dnaseFolder)/*.x8b.bz2))

#############################
## PWM model Recalibration ##
#############################

# Parse out the CENTIPEDE data necessary for PWM recalibration
$(initCenti)combo/%.lambdaParse.gz:
	-R --vanilla --args $* < $(lambdaScript) | gzip > $@ 2> $@.err
comboParse: $(initCenti)combo.Folder $(patsubst $(t5kBedFolder)/%.bed.gz, $(initCenti)combo/%.lambdaParse.gz, $(wildcard $(t5kBedFolder)/*.bed.gz))

# Recalibration of PWM motif models
$(pwmFolder)%.txt.gz: $(pwmFolder).Folder
	-R --vanilla --args $* $(seedFolder) $(pwmFolder) < $(pwmScript) | gzip > $@ 2> $@.err &

# Identification of active tissue-TF pairs
selNames.txt: selPwms.Folder
	-R --vanilla --args $(zThresh) < $(betaScript) > selPwms/output.txt 2> selPwms/output.txt.e;
	cat selPwms/*.pwms.txt | sort -k 1 | uniq > $@


##########################
## Genome-wide PWM scan ##
##########################

# Run the scans for motifs with active footprints in >=1 tissues
all: selNames.txt $(genome)
	cat $< | while read f; do \
		${MAKE} pwms/recal/pwmRescan.sites/$${f}.bed.gz; \
		${MAKE} pwms/recal/pwmSNPs.sites/$${f}.bed.gz; \
	done
	rm $(genome)

# Scan the genome for all matches above a threshold
pwms/recal/pwmRescan.sites/%.bed.gz: $(pwmFolder)/%.pwm
	T=`less $^ | grep -w '#p.Tcut2' | cut -f5`; \
	scanPwm $^ $(genome) -t=$$T -base=2 -omp=$(ncpus) | awk -v OFS='\t' '{print $$1,$$2,$$3,"$*",$$5,$$4}' | bedSort stdin stdout | gzip > $@

# Scan the genome again, this time considering 1KG variants
pwms/recal/pwmSNPs.sites/%.bed.gz: $(pwmFolder)/%.pwm
	T=`less $^ | grep -w '#p.Tcut2' | cut -f5`; \
	scanPwmVar $^ $(genome) $(1kGenomes) -base=2 -t=$$T -omp=$(ncpus) | bedSort stdin stdout | gzip > $@ 2>$@.e



#########################
## Full CENTIPEDE scan ##
#########################

# Basic job for a given tissue-TF pair
full/res/$(expName)/%.out.gz: full/res/$(expName).Folder
	-nice -n10 ionice -c2 -n5 R --vanilla --args $* $(expName) ./ < $(fullCentiScript) 2> $@.err | gzip > $@

# Submit a tissue to run on its own node
full/res/%.exp: $(dnaseFolder)/%.x8b.bz2
	hostname
	mkdir -p full/res/$*
	cp $^ $(tmp)/
	pbzip2 -dm1000fvp$(numJobsFull) $(tmp)/$*.x8b.bz2
	ionice -c2 -n5 make -j $(numJobsFull) full/res/$*.mk expName=$*
	rm -f $(tmp)/$*.x8b
fullQsub: $(patsubst $(dnaseFolder)/%.x8b.bz2, full/res/%.exp.Qsub64, $(wildcard $(dnaseFolder)/*.x8b.bz2))

