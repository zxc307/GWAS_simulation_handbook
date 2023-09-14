## Table of Contents

* [About This Handbook](#About-this-handbook)
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
Visit the SLiM official GitHub repository to download the latest version of the source code: https://github.com/MesserLab/SLiM/releases.  
Here I presented an example of installation using [CMake](https://cmake.org/) in Linux, as GWAS simulations are heavy loading jobs and usually run on servers/high performance computing systems where Linux is widely used.
```ruby
wget https://github.com/MesserLab/SLiM/releases/download/v4.0.1/SLiM.zip
unzip ./SLiM.zip
cd SLiM
mkdir build
cd build
cmake ../
make slim
#The executive will be created in the "build" folder.
```
If you are using other operating systems, please refer to chapter 2 in the [manual](https://github.com/MesserLab/SLiM/releases/download/v4.0.1/SLiM_Manual.pdf). I also recommend MAC users/beginners to try [SLiM GUI](https://github.com/MesserLab/SLiM/releases/download/v4.0.1/SLiM_OSX_Installer.pkg) as it comes with a straightforward debugging system.
### Reference Data
To simulate individual-level GWAS data, we need real samples as founders. Here, we recommend using the well-established 1000-genome project (1KGP), an open-access repository of whole-genome sequencing data of diverse populations (URL: https://www.internationalgenome.org/category/ftp/). You can use the code from the [data resources](#data-resources) section to download it.  
For simplicity, we used distinct ancestral subpopulations such as British from England and Scotland (GBR, N=91) as references. We created a subset of the dataset and performed pre-simulation quality control in the following steps:
* Select a subpopulation group and exclude IN/DELs.
```ruby
for chr in {1..22..1}; do vcftools --gzvcf ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep GBR.list.txt --remove-indels --mac 1 --recode --out ./GBR.chr$chr; done
#"GBR.list.txt" contains a list of family and individual IDs of GBR subpopulation.
#"--remove-indels" option removes all IN/DELs.
#"--mac 1" option removes monomorphic sites.
```
* Exclude multi-allelic sites and transfer data to Plink format
```ruby
for chr in {1..22..1};do plink2 --vcf GBR.chr$chr.vcf --export vcf --out GBR.chr$chr.clean --rm-dup force-first --max-alleles 2;done
#"--rm-dup force-first" option removes duplicate-ID variants and only keep the first instance of each.
#"--max-alleles 2" excludes multi-allelic variants with more than one minor allele.
```
* Clean ancestral sequence/FASTA data  
The extracted reference nucleotide sequence contained many sites, particularly at the beginning and end of the chromosome, designated as N in the FASTA data, meaning that the nucleotide at that position is unknown. Besides, there are some M (bases with aMino groups) and R (bases with puRine) in chromosome 3. SLiM requires that the ancestral sequence be composed only of A/C/G/T nucleotides; N/M/R is not a legal option. For this reason, here we change all of those positions to A for SLiM to read; those sites were not involved in SNPs anyway, so this was harmless for our purposes here.
```ruby
for chr in {1..22..1};do sed -i 's/N/A/g' hs37d5_chr$chr.fa; done;sed -i 's/M/A/g' hs37d5_chr3.fa;sed -i 's/R/A/g' hs37d5_chr3.fa
#sed -i 's/N/A/g' replaced all "N"s with "A"s
#sed -i 's/M/A/g' replaced all "M"s with "A"s
#sed -i 's/R/A/g' replaced all "R"s with "A"s
```

### A quick Simulation
In this section, we performed a simple simulation using the smallest chromosome (chr22) and a small subset of GBR samples (N=30). In the mating model, we utilized Wright Fisher (WF) model and mated people randomly by ten generations. Controllers/parameters of the model are explained in the script.
```ruby
plink2 --vcf GBR.chr22.clean --keep GBR.30list.txt --mac 1 --out ./GBR.chr22.30.recode --export vcf
#"GBR.30list.txt" contains 30 family and individual IDs of GBR reference samples.
#Here, we created a subset vcf file of 30 reference GBR samples.
slim quick.txt
#Run the simulation code; "quick.txt" is provided in this repository
plink2 --vcf GBR.input30.out50.gen10.chr22.vcf
#Show quick stats of the simulated individuals
```
As we set a fixed seed for randomization of the simulation example, you shall have 50 GBR samples with 104,915 variants simulated.

## Simulation Examples
Here we present five simulation examples using different models and settings.
### Model 1 (WF model, uniform distribution of recombination, single chromosome)
In the first example, we applied WF model and uniform distribution of recombination. Neutral mutation with rate of 1e-7 and random recombination with rate of 1e-8 applied to all simulated variants. We used chromosome 22 as an example in the model. Details of the model were described in the SLiM script ([model1](./model1.txt)). Following are codes to create plink format data.
```ruby
slim model1.txt #run slim simulation
awk -F'\t' '!/^#/{$1="22";$3="chr22_"$2} 1' OFS='\t' model1.vcf > temp && mv temp model1.vcf #correct chromosome number and name SNPs by chromosome and positions
plink --vcf model1.vcf --make-bed --out model1 #convert to bi-allelic plink format files, will only keep the most frequenct ALT
```

### Model 2 (WF model, recombination based on hotspots, single chromosome)
In the second example, we applied WF model neutral mutation with rate of 1e-7. However, we simulated recombination sites based on summary of hotspots from 1KGP real data ([the genetic map](./genetic_map_chr22_combined_b37.txt)). Code details were described in the SLiM script ([model2](./model2.txt)). Following are codes to create plink format data.
```ruby
slim model2.txt #run slim simulation
awk -F'\t' '!/^#/{$1="22";$3="chr22_"$2} 1' OFS='\t' model2.vcf > temp && mv temp model2.vcf #correct chromosome number and name SNPs by chromosome and positions
plink --vcf model2.vcf --make-bed --out model2 #convert to bi-allelic plink format files, will only keep the most frequenct ALT
```
### Model 3 (WF model, recombination based on hotspots, whole genome)
SLiM uses two major tick cycles, one is "WF" and the other is "nonWF". You can specify it with option "initializeSLiMModelType()". However, only "nonWF" tick cycle supports customized reproduction model and fixed pedigree input which is required to simulate whole-genome data. Alternatively, we can build a WF model by customizing the survival model in a two-step "nonWF" tick cycle to achieve WF modeling on whole-genome simulation. Typically, we dropped all parents in the survival model in the tick cycle. All settings of recombination and mutations are the same as model 2. Pool size of each generation is 240 (30 * 2^3).  
Firstly, we simulated a fixed mating pattern. Secondly, we applied the mating pattern to chromosome 22 as an example. Code details were described in the SLiM script ([model3_1](./model3_1.txt))([model3_2](./model3_2.txt)). Following are codes to create plink format data. Additionally, you will have all survived IDs of the latest generation recorded in the "model3.mating.txt" file.
```ruby
slim model3_1.txt;slim model3_2.txt #run two-step model
awk -F'\t' '!/^#/{$1="22";$3="chr22_"$2} 1' OFS='\t' model3.vcf > temp && mv temp model3.vcf #correct chromosome number and name SNPs by chromosome and positions
plink --vcf model3.vcf --make-bed --out model3 #convert to bi-allelic plink format files, will only keep the most frequenct ALT
```
To simulate other chromosomes, simply customize ([model3_2](./model3_2.txt)) by change the chromosome number. In this model, you can apply a fixed mating pattern to all chromosomes to achieve whole-genome simulation.  

### Model 4 (non-WF model, recombination based on hotspots, single chromosome)
In model 4, we applied non-WF model to single chromosome simulation. Simulated samples end up being a mixed generation from random mating. Code details were described in the SLiM script ([model4](./model4.txt)). Following are codes to create plink format data. Samples will be simulated from a mixed generation with randomization.
```ruby
slim model4.txt 
awk -F'\t' '!/^#/{$1="22";$3="chr22_"$2} 1' OFS='\t' model4.vcf > temp && mv temp model4.vcf #correct chromosome number and name SNPs by chromosome and positions
plink --vcf model4.vcf --make-bed --out model4 #convert to bi-allelic plink format files, will only keep the most frequenct ALT
```
### Model 5 (non-WF model, recombination based on hotspots, random survival, whole genome)
In model 5, we applied non-WF model to whole-genome simulation. We used a two-step simulation. Firstly, we simulated a fixed mating pattern. Secondly, we applied the mating pattern to chromosome 22 as an example. Code details were described in the SLiM script ([model5_1](./model5_1.txt))([model5_2](./model5_2.txt)). Following are codes to create plink format data. Additionally, you will have all survived IDs of a mixed generation recorded in the "model5.mating.txt" file.
```ruby
slim model5_1.txt;slim model5_2.txt #run two-step model
awk -F'\t' '!/^#/{$1="22";$3="chr22_"$2} 1' OFS='\t' model5.vcf > temp && mv temp model5.vcf #correct chromosome number and name SNPs by chromosome and positions
plink --vcf model5.vcf --make-bed --out model5 #convert to bi-allelic plink format files, will only keep the most frequenct ALT
```
To simulate other chromosomes, simply customize ([model5_2](./model5_2.txt)) by change the chromosome number. In this model, you can apply a fixed mating pattern to all chromosomes to achieve whole-genome simulation.  

## Parameter Effects
## Post Simulation Quality Control
## Files In This Github Repository
[GBR.30list.txt](./GBR.30list.txt) contains 30 GBR family and individual IDs sampled from 1000-Genome data.  
[GBR.list.txt](./GBR.list.txt) contains 91 GBR family and individual IDs sampled from 1000-Genome data.  
[quick.txt](./quick.txt) an example SLiM script to simulate 50 GBR samples with 104,915 variants.  
[model1.txt](./model1.txt) a SLiM script to simulate samples using model 1 described in [Simulation Examples](#simulation-examples)  
[model2.txt](./model2.txt) a SLiM script to simulate samples using model 2 described in [Simulation Examples](#simulation-examples)  
[model3_1.txt](./model3_1.txt)      [model3_2.txt](./model3_2.txt) two-step scripts to simulate samples using model 3 described in [Simulation Examples](#simulation-examples)  
[model4.txt](./model4.txt) a SLiM script to simulate samples using model 4 described in [Simulation Examples](#simulation-examples)  
[model5_1.txt](./model5_1.txt)      [model5_2.txt](./model5_2.txt) two-step scripts to simulate samples using model 5 described in [Simulation Examples](#simulation-examples)  
## Data and Software Resources

### Data resources
#### 1000G phase 3 reference data (2502 individuals)
You can download the open-access data from the [FTP server](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/).
A brief description of the data can be found in the [documentation](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/README_phase3_callset_20150220).
In Linux, you can download all VCF format data with FASTA format ancestral sequence using the following code:
```ruby
wget -r -np -nd ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
#Download all VCF files and documentations
wget -r -np -nd ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
#Download the FASTA file
gzip -d hs37d5.fa.gz
#Unzip FASTA file
for chr in {1..22..1};do;samtools faidx hs37d5.fa "$chr" > hs37d5_chr$chr.fa;done
#Subset FASTA file by chromosomes
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


## Other Web Resources
* [Genetic Simulation Resources (GSR)](https://surveillance.cancer.gov/genetic-simulation-resources/)
Genetic Simulation Resources is provided by the Division of Cancer Control and Population Sciences at the NCI.
This website provides a catalogue of existing software packages that simulate genetic data of the human genome. This website is designed to help you identify simulators that meet your particular research needs, compare simulators based on attributes and provide feedback to authors of existing simulators.
* [JoniColeman's GWAS codebook](https://github.com/JoniColeman/gwas_scripts/tree/master)
The scripts in this GitHub repository provide a straight-forward guide to the quality control, imputation and analysis of genome-wide genotype data.
* [stdpopsim](https://popsim-consortium.github.io/stdpopsim-docs/stable/introduction.html)
Stdpopsim is a community-maintained standard library of population genetic models.
## Getting Help
### Regarding this handbook
If you have any questions or need help with this handbook, you can either create an issue under this repository or email me at zxc307@case.edu.

### Regarding related software
Following are mailing lists for related software:  
plink 1.9/2.0: plink2-users@googlegroups.com  
SLiM: slim-discuss@googlegroups.com  
vcftools: vcftools-help@lists.sourceforge.net  
SAMtools: samtools-help@lists.sourceforge.net  (You may need to subscrib to the [samtools-help list](https://sourceforge.net/projects/samtools/lists/samtools-help))  


You can also get support by creating issues under GitHub repositories of the following software:
* [SLiM/issues](https://github.com/MesserLab/SLiM/issues)
* [SAMtools/issues](https://github.com/samtools/samtools/issues)
* [vcftools/issues](https://github.com/vcftools/vcftools/issues)

## Acknowledgements
Great thanks to the following people for their help and input:
* [Dr. Benjamin C. Haller](http://benhaller.com/) provided invaluable assistance in coding and SLiM debugging.
* [Dr. Christopher Chang](https://www.linkedin.com/in/christopher-chang-6910a51/) provided great help with Plink-related issues and code debugging.
* [Dr. Fredrick Schumacher](https://case.edu/medicine/pqhs/about/people/primary-faculty/fredrick-r-schumacher) supervised and funded this project.
* [Dr. Wei-min Chen](https://www.chen.kingrelatedness.com/) shared his expertise and gave us great advice regarding genetic relatedness issues.


