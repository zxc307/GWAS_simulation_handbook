## Table of Contents

* [About This Handbook](#About)
* [Getting Started](#getting-started)
  * [Install Software](#install-software)
  * [Reference Data](#data)
  * [A quick Simulation](#a-quick-simulation)
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
#### 1000G phase 3 reference data (2502 individuals)
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

## Web Resources
## Getting Help
## Acknowledgements