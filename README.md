#VCFped
VCFped is a command line utility which quickly and confidently identifies close relationships (e.g. parent-child, trios and quartets) in multisample VCF files.
The algorithm does not rely on allele frequencies or other annotations, and is highly robust to inbreeding and poor quality variants.

#Installation and basic usage:
The simplest installation of VCFped is using pip:

    >pip install vcfped

To run VCFped on a variant file *myvariants.vcf*, the basic command is 
    
    >vcfped myvariants.vcf

In most circumstances VCFped can be run unsupervised with default parameters. For overview of the options, use the -h option:
    
    >vcfped -h
    usage: vcfped.py [-h] [--quiet] [--version] [--prefix PREFIX] [--all]
                     [--no-gender] [--no-pairwise] [--no-trio]
                     [-var VARIABLES [VARIABLES ...]]
                     [-p PERCENTILES [PERCENTILES ...]] [-e EXACTMAX]
                     [-s SAMPLESIZE] [-sAABB SAMPLESIZEAABB] [-t1 THRESHTRIO1]
                     [-t2 THRESHTRIO2] [-male THRESHMALE] [-female THRESHFEMALE]
                     file

    positional arguments:
      file                  path to VCF file

    optional arguments:
      -h, --help            show this help message and exit
      --quiet               do not print information on the screen
      --version             show program's version number and exit
      --prefix PREFIX       prefix for output files
      --all                 show results for all pairs/triples (not only the
                            inferred)
      --no-gender           skip gender analysis
      --no-pairwise         skip pairwise analysis
      --no-trio             skip trio analysis
      -var VARIABLES [VARIABLES ...]
                            quality variables to be used for filtering
      -p PERCENTILES [PERCENTILES ...]
                            filtering percentiles
      -e EXACTMAX           if approx. line count exceeds this, apply random
                            sampling
      -s SAMPLESIZE         sample at least this many variant lines (if sampling)
      -sAABB SAMPLESIZEAABB
                            sample at least this many lines where both 0/0 and 1/1
                            occur as genotypes (if sampling)
      -t1 THRESHTRIO1       threshold (%) for trio test 1 (AA + BB = AB)
      -t2 THRESHTRIO2       threshold (%) for trio test 2 (BB + BB = BB)
      -male THRESHMALE      upper limit (%) for male heterozygosity on X
      -female THRESHFEMALE  lower limit (%) for female heterozygosity on X

---

#Introduction
Correct relationships are crucial for successful analysis of family-based sequencing, e.g. detection of de novo variants in trios. Traditional methods for relatedness inference are not easily adaptable to sequencing data, since several basic assumptions of these methods are not satisfied (e.g. independent marker loci, correct allele frequencies, non-censored data). Our program sidesteps these assumptions by exploiting forced allele sharing between closely related individuals. In particular, the algorithm does not depend on allele frequencies or other variant annotations, making the program suitable as a quality control step prior to annotation and downstream analysis. 

VCFped is written in Python and has been extensively tested on variant files from whole-exome sequencing (WES), whole-genome sequencing (WGS) and the TruSightOne gene panel from Illumina.

#Filtering strategy
The relatedness testing algorithms of VCFped depend on correct genotypes for the variants. To minimize the inclusion of wrong genotypes a semi-automated filtering strategy has been implemented, using the quality variables present in the variant file (e.g. QUAL, DP, GQ). The complete distributions of the available variables are computed when loading the file, and a stepwise filtering process is created based on percentiles of these distributions (see Figure 1). The default percentile ranks are 10, 30 and 50, but these can be modified using the *-p* option.

#Relatedness inference using forced allele sharing
For a few close relationships the genotype of one individual can be determined (under Mendelian inheritance) by the genotypes of others. VCFped exploits the uniqueness of such forced allele sharing to identify trios, parent-child pairs and monozygotic twins.

###Trio relationships
To identify trios a combination of two tests are used. The first recognizes all basic trio types (regular, inverted, generational), while the second discriminates regular trios from the other two types. 

####*Trio test 1: AA + BB = AB*
What's computed: For each ordering of each triple of samples, the conditional frequency of genotype AB in one sample given AA and BB in the two others.
Default threshold: 90 %

####*Trio test 2: BB + BB = BB*
What's computed: For each ordering of each triple of samples, the conditional frequency of genotype BB in one sample given BB in both the others.
Default threshold: 95 %

###Pairwise relationships
####*Parent-child test: Identity by state = 0*
What's computed: For each pair of samples, the frequency of AA + BB given the presence of B in at least one of the individuals.
Default threshold: 5 %

####*MZ twins test: Identity by state = 2*
What's computed: For each pair of samples, the frequency of sites with identical genotypes.
Default threshold: 95 %


