# VCFped
VCFped is a command line utility which quickly and confidently identifies close relationships (e.g. parent-child, trios and quartets) in multisample VCF files.
The algorithm does not rely on allele frequencies or other annotations, and is highly robust to inbreeding and poor quality variants.

## Installation and basic usage:
The simplest installation of VCFped is using pip:

    >pip install vcfped

To run VCFped on a variant file *myvariants.vcf*, the basic command is 
    
    >vcfped myvariants.vcf

In most circumstances VCFped can be run unsupervised with default parameters. For overview of the options, use the -h option:
    
    >vcfped -h
    usage: vcfped.py [-h] [--quiet] [--version] [--no-gender] [--no-pairwise]
                     [--no-trio] [--all] [-o PREFIX]
                     [-v VARIABLES [VARIABLES ...]]
                     [-p PERCENTILES [PERCENTILES ...]] [-e EXACTMAX]
                     [-s SAMPLESIZE] [-d SAMPLESIZEAABB] [-t1 T1_THRESH]
                     [-t2 T2_THRESH] [-mz MZ_THRESH] [-po PO_THRESH]
                     [-female FEMALE_THRESH] [-male MALE_THRESH]
                     file
    

    positional arguments:
      file                  path to VCF file

    optional arguments:
      -h, --help            show this help message and exit
      --quiet               do not print information on the screen
      --version             show program's version number and exit
      --no-gender           skip gender analysis
      --no-pairwise         skip pairwise analysis
      --no-trio             skip trio analysis
      --all                 show results for all pairs/triples (not only the
                            inferred)
      -o PREFIX             prefix for output files
      -v VARIABLES [VARIABLES ...]
                            quality variables to be used for filtering
      -p PERCENTILES [PERCENTILES ...]
                            filtering percentile ranks
      -e EXACTMAX           if approx. line count exceeds this, apply random
                            sampling
      -s SAMPLESIZE         sample at least this many variant lines (if sampling)
      -d SAMPLESIZEAABB     sample at least this many lines where both 0/0 and 1/1
                            occur as genotypes (if sampling)
      -t1 T1_THRESH         threshold (%) for T1 score (AA + BB = AB)
      -t2 T2_THRESH         threshold (%) for T2 score (BB + BB = BB)
      -mz MZ_THRESH         threshold (%) for MZ score (IBS=2 | neither is AA)
      -po PO_THRESH         threshold (%) for PO score (IBS>0 | either is BB)
      -female FEMALE_THRESH
                            lower limit (%) for female heterozygosity on X
      -male MALE_THRESH     upper limit (%) for male heterozygosity on X
      
---

## Introduction
Correct relationships are crucial for successful analysis of family-based sequencing, e.g. detection of de novo variants in trios. Traditional methods for relatedness inference are not easily adaptable to sequencing data, since several basic assumptions of these methods are not satisfied (e.g. independent marker loci, correct allele frequencies, non-censored data). VCFped sidesteps these assumptions by exploiting forced allele sharing between closely related individuals. In particular, the algorithm does not depend on allele frequencies or other variant annotations, making the program well suited as a quality control step prior to annotation and downstream analysis. 

VCFped is written in Python and has been extensively tested on variant files from whole-exome sequencing (WES), whole-genome sequencing (WGS) and the TruSightOne gene panel from Illumina.

## Filtering strategy
The relatedness testing algorithms of VCFped depend on correct genotypes for the variants. To minimize the inclusion of wrong genotypes a semi-automated filtering strategy has been implemented, using the quality variables present in the variant file (e.g. QUAL, DP, GQ). The complete distributions of the available variables are computed when loading the file, and a stepwise filtering process is created based on percentiles of these distributions (see Figure 1). The default percentile ranks are 10, 30 and 50, but these can be modified using the *-p* option.

## Relatedness inference using forced allele sharing
For a few close relationships the genotype of one individual can be determined (under Mendelian inheritance) by the genotypes of others. VCFped exploits such forced allele sharing to identify trios, parent-child pairs and monozygotic twins. In the description below, *A* and *B* denote the labels of a diallelic locus, while *g<sub>i</sub>* is the genotype of individual *i*. For trio and pariwise relationships only autosomal loci are considered, while X-linked loci are used for gender inference.

### Trio relationships
To identify trios a combination of two tests are used. The first (*AA + BB = AB*) recognizes the three basic trio types: regular (parent-offspring), inverted (two siblings and one parent) and  generational (child-parent-grandparent). The second (*BB + BB = BB*) discriminates regular trios from the other two types. Both test scores are computed for every cyclic ordering of each triple of samples in the input variant file.

##### *T1 score = freq(g<sub>3</sub> = AB | g<sub>1</sub> = AA, g<sub>2</sub> = BB or vice versa)*
What's computed: The conditional frequency of genotype AB in one sample given AA and BB in the two others. The AB genotype is forced for the child of a regular trio, the parent in an inverted trio, and the middle individual in a generational trio.  
Default threshold: 90 %

##### *T2 score = freq(g<sub>3</sub> = BB | g<sub>1</sub> = g<sub>2</sub> = BB)*  
What's computed: The conditional frequency of genotype BB in one sample given BB in both the others. This is forced for the child of a regular trio, but also if individual 3 is a monozygotic twin (or a sample duplicate!) of one of the others.  
Default threshold: 95 %   (modify with option -t2)

### Pairwise relationships
For pairwise relatedness VCFped tests for monozygotic twins (MZ) and parent-offspring (PO), which both obey patterns of forced allele sharing. MZ twins always have both alleles identical by state (IBS), while parent-offspring are always IBS > 0. The conditions in the following expressions ensure that the scores are unaffected by the data censoring of VCF files.

##### *MZ score = freq(IBS=2 | neither is AA)*
What's computed: The frequency of equal genotypes among all variants where both have a least one B.  
Default threshold: 95 %

##### *PO score = freq(IBS>0 | either is BB)*
What's computed: The frequency of at least one shared allele, given that at least one is homozygous BB.    
Default threshold: 99 %

### Gender prediction
VCFped predicts the gender of each sample by using variants on X (except pseudoautosomal regions):

##### *XHET=freq(AB | AB or BB)*
What's computed: The heterozygosity on X among all non-AA variants.  
Default interpretation: *Male* if XHET < 5 %, *female* if XHET > 25 %
