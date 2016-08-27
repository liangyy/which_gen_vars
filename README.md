# Which Genetics Variants in DNase-Seq Footprints Are More Likely to Alter Binding?
The code in this repository constitutes the pipeline described in [Moyerbrailean *et al.* (2016)](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005875)

For more detailed instructions on how to run the code, see the [makefile](Makefile).

### Required software
* [samtools](http://samtools.sourceforge.net/)
* [bedtools](http://bedtools.readthedocs.io/en/latest/)
* [UCSC command-line utils](http://hgdownload.soe.ucsc.edu/admin/exe/)

### Custom utils
* scanPwm
* scanPwmVar
* makeXbFile
* bedXbFileXtract

### Directory structure and supporing files

See the [example directory structure file](exampleDirectoryStructure.txt) to see how to set up your analysis directory to properly utilize the Makefile.

Several supporting files are required for the analysis, specific to the genome the DNase-seq samples were aligned to and the variants you wish to analyze:

#### chromSizes.[genomeversion].txt
A two-column file listing each chromosome and the size in bases. E.g.,

    chr1	249250621
    chr2	243199373
    ...

#### variants.txt.gz
A tab separated file with one row per variant: chromosome / position / ref allele / alt allele. E.g., 

    chr1   	10583  	G      	A
    chr1   	10611  	C      	G
    chr1   	13302  	C      	T
    ...


#### genome.fa and genome.2bit
Fasta and 2bit versions of the reference genome (e.g., hg19.fa), used for the PWM scanning steps.


## Pipeline overview

### 1) Initial PWM scan
Using the seed motifs, generate the training sites using 'scanPwm'
1. Scan Human genome for motifs, choose 5k top-scoring matches
2. Scan Chimp genome for motifs, liftover up to 5k top-scoring matches (not in [1])
3. Scan Macaque genome for motifs, liftover up to 5k top-scoring matches (not in [1,2])
4. Extract the Human DNA sequences at the positions identified in 1-3


### 2) Initial CENTIPEDE scan

1. Convert DNase-seq data to x8b files
2. Get the coverage (number of reads) for each DNase-seq sample
3. For each DNase-seq sample, run CENTIPEDE for each PWM model

Depending on your system and number of samples/PWMs, Step 3 can be parallelized such that each Sample is run on a separate node, and within a node, each Sample-PWM pair is processed on a separate core.


### 3) PWM model Recalibration

1. Parse out the CENTIPEDE data necessary for PWM recalibration
2. Recalibration of PWM motif models
3. Identification of active tissue-TF pairs

### 4) Genome-wide PWM scan

Focusing only on Sample-PWM combinations that showed footprints in Part 2,

1. Scan the genome for all matches above a match threshold
2. Scan the genome again, this time calculating the binding scores for alternate allels when overlapping 1KG variants


### 5) Full CENTIPEDE scan

Run CENTIPEDE for each Sample-PWM combination as in Part 2 step 3, but this time on the full set of PWM matches from Part 4
