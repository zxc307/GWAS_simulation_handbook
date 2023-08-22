## Table of Contents

* [About This Handbook](#About)
* [Getting Started](#getting-started)
* [Simulation Examples](#simulation-examples)
* [Parameter Effects](#parameter-effects)
* [Post Simulation Quality Control](#post-simulation-quality-control)
* [Files In This Github Repository](#files-in-this-github-repository)
* [Data and Software Resources](#data-and-software-resources)
* [Web Resources](#web-resources)
* [Getting Help](#getting-help)
* [Acknowledgements](#acknowledgements)

## About This Handbook

We aim to present a straightforward guide to the quality control, data simulation and simulation quality assessment of individual-level pseudo GWAS data. Our handbook targets on beginners of GWAS simulation. Advanced users requiring more complex modeling strategies are recommended to consult software developer’s websites to seek options that are more sophisticated. A brief list of these is provided in the [Data and Software Resources](#data-and-software-resources) section.

## Getting Started
### Install Software

### Reference Data
### A quick Simulation
## Simulation Examples
## Parameter Effects
## Post Simulation Quality Control
## Files In This Github Repository
## Data and Software Resources

### Data resources
* 1000G phase 3 reference data (2502 individuals)
You can download the open-access data from the [FTP server](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/).
A brief description of the data can be found in the [documentation](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/README_phase3_callset_20150220).
In Linux, you can download all VCF format data with FASTA format ancestral sequence using the following code:
```ruby
wget -r -np -nd ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
wget -r -np -nd ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gzip -d hs37d5.fa.gz
for chr in {1..22..1};do;samtools faidx hs37d5.fa "$chr" > hs37d5_chr$chr.fa;done
```

### Program sources (*All examples in this repository are run in Linux system*)
#### Plink and Plink2
[Plink 1.9](https://www.cog-genomics.org/plink/1.9/) and [Plink 2](https://www.cog-genomics.org/plink/2.0/) combined are used in this script.
Plink2 is an advanced version with lots of new features and most importantly runs faster.
However, it is still under development and [some modules such as fully-powered data merging are not available yet](https://www.cog-genomics.org/plink/2.0/#:~:text=its%20own%20score.-,Coming%20next,-Fully%2Dpowered%20merge).
Thus, we recommend to use Plink 1.9 for some basic data management process and only use Plink2 for new and fully-powered features.

#### SLiM
The official site of SLiM can be found [here](https://messerlab.org/).
You can also visit their [GitHub site](https://github.com/MesserLab/SLiM).
All versions of SLiM release and detailed documentations are also available [here](https://github.com/MesserLab/SLiM/releases).

#### VCFtools
Their official GitHub site is [here](https://github.com/vcftools/vcftools) with a detailed [manual](https://vcftools.github.io/man_latest.html).

#### SAMtools
Authors of SAMtools provide some standard workflows on their [homepage](http://www.htslib.org/).
You can also check their latest [manual](http://www.htslib.org/doc/samtools.html).


## Web Resources
## Getting Help
### Regarding this handbook
If you have any questions or need help with this handbook, you can either create an issue under this repository or email me at zxc307@case.edu.

### Regarding related software
Following are mailing lists for related software:
* [plink 1.9/2.0](plink2-users@googlegroups.com)
* [SLiM](slim-discuss@googlegroups.com)
* [vcftools](vcftools-help@lists.sourceforge.net)
* [SAMtools](samtools-help@lists.sourceforge.net)： You may need to subscrib to the [samtools-help list](https://sourceforge.net/projects/samtools/lists/samtools-help)


You can also get support by creating issues under GitHub repositories of the following software:
* [SLiM/issues](https://github.com/MesserLab/SLiM/issues)
* [SAMtools/issues](https://github.com/samtools/samtools/issues)
* [vcftools/issues](https://github.com/vcftools/vcftools/issues)

## Acknowledgements