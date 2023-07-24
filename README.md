# About This Handbook
We aim to present a straightforward guide to the quality control, data simulation and simulation quality assessment of individual-level pseudo GWAS data. Our tutorial targets on beginners of GWAS simulation. Advanced users requiring more complex modeling strategies are recommended to consult software developer’s websites to seek options that are more sophisticated. A brief list of these is provided in the **Data and Software** section.

# Methods and Codes

## Data resources
### 1000G phase 3 reference data (2502 individuals)
#### FTP server
```ruby
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
```
#### Documentation
```ruby
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/README_phase3_callset_20150220
```
#### Download vcf
```ruby
wget -r -np -nd ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
```
#### Download FASTA with the ancestral sequence
```ruby
wget -r -np -nd ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gzip -d hs37d5.fa.gz
for chr in {1..22..1}
do
samtools faidx hs37d5.fa "$chr" > hs37d5_chr$chr.fa
done
```

## Program sources (*All examples in this script are run in Linux system*)
### Plink and Plink2
Plink(version 1.9) and Plink2(alpha) combined are used in this script. Plink2 is an advanced version of the popular Plink and it has lots of new features and options in data management. However, it is a in development and only beta versions are available. Thus, we recommend to use Plink 1.9 for some basic data management process and only use Plink2 for new features that have been tested.
#### Plink and Plink2 website
```ruby
https://www.cog-genomics.org/plink/1.9/
```
#### Downloads
Both Plink and Plink2 have excutable binary versions so you can directly run them after downloads and do not need to install them into your system.
```ruby
ttps://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20210606.zip
https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20210908.zip
```

### SLiM
#### Website
```ruby
https://messerlab.org/
```
#### Documentation
```ruby
http://benhaller.com/slim/SLiM_Manual.pdf
```
#### Installation with cmake in Linux
```ruby
wget http://benhaller.com/slim/SLiM.zip
unzip ./SLiM.zip
cd SLiM
mkdir build
cd build
cmake ../
make slim
```
### VCFtools
#### GitSite and Documentation
```ruby
https://github.com/vcftools/vcftools
```
#### Installation in Linux
If you have Git:
```ruby
git clone git://github.com/vcftools/vcftools.git
cd vcftools
./autogen.sh
./configure
make
make install
```
If you do not have Git:
```ruby
wget https://cfhcable.dl.sourceforge.net/project/vcftools/vcftools_0.1.13.tar.gz
tar -xvzf vcftools_0.1.13.tar.gz
cd vcftools_0.1.13
make
```
### SAMtools
#### Website
```ruby
http://www.htslib.org/
```
#### Documentation
```ruby
http://www.htslib.org/doc/#manual-pages
```
#### Installation in Linux
```ruby
wget https://gigenet.dl.sourceforge.net/project/samtools/samtools/1.13/samtools-1.13.tar.bz2
tar -xf samtools-1.13.tar.bz2
cd samtools-1.13
make
```
### 
#### Website
#### Documentation
#### Installation in Linux

## Pre simulation quality control

### Select reference panel by ancestry and drop INDELs (use three populations as an example)
```ruby
for pop in {YRI,CHB,GBR}
do
    for chr in {1..22..1}
    do
    vcftools --gzvcf ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep $pop.90list.txt --remove-indels --mac 1 --recode --out ./$pop/$pop.chr$chr.90
    done
done
```
### Exclude tri-allelic locations
```ruby
for pop in {YRI,CHB,GBR}
do
    for chr in {1..22..1}
    do
    cat ./$pop/$pop.chr$chr.90.recode.vcf | awk '!a[$2]++{print}' > ./$pop/$pop.chr$chr.90.recode.vcf.cp
    perl -lane 'if(/^#/ or length("$F[3]$F[4]")==2){print}' ./$pop/$pop.chr$chr.90.recode.vcf.cp > ./$pop/$pop.chr$chr.90.recode.vcf
    done
done
```
### Exclude non-variant and transfer the reference to plink format
```ruby
for pop in {YRI,CHB,GBR}
do
    for chr in {1..22..1}
    do
    plink --vcf ./$pop/$pop.chr$chr.90.recode.vcf --recode --out ./$pop/$pop.chr$chr.90
    plink --file ./$pop/$pop.chr$chr.90 --make-bed --out ./$pop/$pop.chr$chr.bin.90
    done
done
```
### Clean FASTA files for SLiM
The extracted reference nucleotide sequence contained a large number of sites, particularly at the beginning and end of the chromosome, designated as N in the FASTA data, meaning that the nucleotide at that position is unknown. Besides, there are some M (bases with aMino groups) and R (bases with puRine) in chromosome 3. SLiM requires that the ancestral sequence be composed only of A/C/G/T nucleotides; N/M/R is not a legal option. For this reason, please change all of those positions to A for SLiM to read; those sites were not involved in SNPs anyway, so this was harmless for our purposes here.

```ruby
for chr in {1..22..1}
do
sed -i 's/N/A/g' hs37d5_chr$chr.fa
done
sed -i 's/M/A/g' hs37d5_chr3.fa
sed -i 's/R/A/g' hs37d5_chr3.fa
```
## WGS simulation with SLiM
### Extract reference samples (VCF) for SLiM 
```ruby
for pop in {YRI,CHB,GBR}
do
    for chr in {1..22..1}
    do
    vcftools --vcf ./$pop/$pop.chr$chr.90.recode.vcf --keep $pop.60list.txt --mac 1 --recode --out ./$pop/$pop.chr$chr.60
    vcftools --vcf ./$pop/$pop.chr$chr.60.recode.vcf --keep $pop.30list.txt --mac 1 --recode --out ./$pop/$pop.chr$chr.30
    done
done
```
### WGS simulation
#### Write simulation scripts
```ruby
for simN in {1..10..1}
do
  for popC in {YRI,CHB,GBR}
  do 
    for generationN in {30,100,300}
    do
      for inputN in {30,60,90}
      do
        for chrN in {1..22..1}
        do
        cat slim.code.template.txt | sed "s/popC/$popC/g" | sed "s/simN/$simN/g" | sed "s/generationN/$generationN/g" | sed "s/inputN/$inputN/g" | sed "s/chrN/$chrN/g" > ./$popC/$popC.POP$inputN/out2000/slim.code.$popC.POP$inputN.simN$simN.gen$generationN.chr$chrN.txt
        done
      done
    done
  done
done
```
#### Simulate WGS data
```ruby
**************************************************************************************************
**************************************************************************************************
**************************************************************************************************
**************************************************************************************************
cd /mnt/rstor/SOM_EPBI_FRS2/zxc307/1000g2021
module load slim
module load plink
module load vcftools
module load plink2

for simN in {1..10..1}
do
  for popC in {YRI,CHB,GBR}
  do 
    for generationN in {30,100,300}
    do
      for inputN in {30,60,90}
      do
        for chrN in {9..9..1}
        do
        slim ./$popC/$popC.POP$inputN/out2000/slim.code.$popC.POP$inputN.simN$simN.gen$generationN.chr$chrN.txt
        perl -lane 'if(/^#/ or length("$F[3]$F[4]")==2){print}' ./$popC/$popC.POP$inputN/out2000/sim$simN.gen$generationN.chr$chrN.vcf > ./$popC/$popC.POP$inputN/out2000/sim$simN.gen$generationN.chr$chrN.notri.vcf
        done
      done
    done
  done
done
```
### WGS data management
#### Set position as id, correct chromosome number and make files in plink2 format
```ruby
for simN in {10..10..1}
do
  for popC in {YRI,CHB,GBR}
  do 
    for generationN in {30,100,300}
    do
      for inputN in {30,60,90}
      do
        for chrN in {1..22..1}
        do
        plink --vcf ./$popC/$popC.POP$inputN/sim$simN.gen$generationN.chr$chrN.notri.vcf --make-bed --out ./$popC/$popC.POP$inputN/sim$simN.gen$generationN.chr$chrN.notri.bin
        cp ./$popC/$popC.POP$inputN/sim$simN.gen$generationN.chr$chrN.notri.bin.bim ./$popC/$popC.POP$inputN/sim$simN.gen$generationN.chr$chrN.notri.bin.bim.cp
        awk '{{print "'$chrN'","chr""'$chrN'""_"$4,$3,$4,$5,$6}}' ./$popC/$popC.POP$inputN/sim$simN.gen$generationN.chr$chrN.notri.bin.bim.cp > ./$popC/$popC.POP$inputN/sim$simN.gen$generationN.chr$chrN.notri.bin.bim
        plink2 --bfile ./$popC/$popC.POP$inputN/sim$simN.gen$generationN.chr$chrN.notri.bin --make-pgen --out ./$popC/$popC.POP$inputN/sim$simN.gen$generationN.chr$chrN.notri.bin
        done
      done
    done
  done
done
```
#### WGS data merge
##### Create mergelist
```ruby
for simN in {1..10..1}
do
  for popC in {YRI,CHB,GBR}
  do 
    for generationN in {30,100,300}
    do
      for inputN in {30,60,90}
      do
      cat mergelist.template.txt | sed "s/popC/$popC/g" | sed "s/simN/$simN/g" | sed "s/inputN/$inputN/g" | sed "s/generationN/$generationN/g" > ./$popC/$popC.POP$inputN/mergelist.$popC.POP$inputN.simN$simN.gen$generationN.txt
      done
    done
  done
done
```
##### Merge
```ruby
for simN in {1..10..1}
do
  for popC in {YRI,CHB,GBR}
  do 
    for generationN in {30,100,300}
    do
      for inputN in {30,60,90}
      do
      plink2 --pmerge-list ./$popC/$popC.POP$inputN/mergelist.$popC.POP$inputN.simN$simN.gen$generationN.txt --make-pgen --out ./simulatedWGS/simN$simN.$popC.POPin$inputN.gen$generationN.notri.bin
      done
    done
  done
done
```
## WGS to GWAS

## GWAS quality control

### Principal component analysis
#### Reference panel data management
```ruby
for chr in {1..22..1}
do  
vcftools --gzvcf ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --remove-indels --mac 1 --recode --out ./refplink/1000Gwhole.chr$chr
cat ./refplink/1000Gwhole.chr$chr.recode.vcf | awk '!a[$2]++{print}' > ./refplink/1000Gwhole.chr$chr.recode.vcf.cp
perl -lane 'if(/^#/ or length("$F[3]$F[4]")==2){print}' ./refplink/1000Gwhole.chr$chr.recode.vcf.cp > ./refplink/1000Gwhole.chr$chr.recode.vcf
plink2 --vcf ./refplink/1000Gwhole.chr$chr.recode.vcf --make-bed --out ./refplink/ref1000G.chr$chr
cp ./refplink/ref1000G.chr$chr.bim ./refplink/ref1000G.chr$chr.bim.cp
awk '{{print "'$chr'","chr""'$chr'""_"$4,$3,$4,$5,$6}}' ./refplink/ref1000G.chr$chr.bim.cp > ./refplink/ref1000G.chr$chr.bim
plink2 --bfile ./refplink/ref1000G.chr$chr --make-pgen --out ./refplink/ref1000G.chr$chr
done
plink2 --pmerge-list refmergelist.txt --merge-max-allele-ct 2 --make-pgen --out ./refplink/ref1000G.whole
plink2 --pfile ./refplink/1000G.whole --indep 50 5 2 --out ./refplink/1000G.whole
```
#### PCA data merge and prune
```ruby
for simN in {1..10..1}
do
  for popC in {YRI,CHB,GBR}
  do 
    for generationN in {30,100,300}
    do
      for inputN in {30,60,90}
      do
      sed '1 /\.\/simulatedWGS\/simN'$simN'\.'$popC'\.POPin'$inputN'\.gen'$generationN'\.notri\.bin/g;p' -i simmerge.list.txt
      done
    done
  done
done

for simN in {1..10..1}
do
  for popC in {YRI,CHB,GBR}
  do
    for generationN in {30,100,300}
    do
      for inputN in {30,60,90}
      do
      plink2 --pfile ./simulatedWGS/simN$simN.$popC.POPin$inputN.gen$generationN.notri.bin --make-bed --out ./simulatedWGS/simN$simN.$popC.POPin$inputN.gen$generationN.notri.bin
      plink --bfile ./simulatedWGS/simN$simN.$popC.POPin$inputN.gen$generationN.notri.bin --bmerge ./refplink/ref1000G.whole.prune --merge-mode 1 --geno 0.01 --make-bed --out ./simulatedWGS/simN$simN.$popC.POPin$inputN.gen$generationN.notri.bin.pca
      done
    done
  done
done

simN=1;popC=YRI;generationN=30;inputN=30;
plink2 --pfile ./simulatedWGS/simN$simN.$popC.POPin$inputN.gen$generationN.notri.bin --pmerge ./refplink/ref1000G.whole.prune --geno 0.01 --make-pgen --out ./simulatedWGS/simN$simN.$popC.POPin$inputN.gen$generationN.notri.bin.pca
      

plink2 --pfile ./refplink/ref1000G.whole.prune --make-bed --out ./refplink/ref1000G.whole.prune

plink2 --pfile ./refplink/1000G.whole --make-bed --out ./refplink/1000G.whole
      plink2 --pfile ./simulatedWGS/simN$simN.$popC.POPin$inputN.gen$generationN.notri.bin --make-bed --out ./simulatedWGS/simN$simN.$popC.POPin$inputN.gen$generationN.notri.bin
      plink --bfile ./simulatedWGS/simN$simN.$popC.POPin$inputN.gen$generationN.notri.bin --bmerge ./refplink/1000G.whole --make-bed --out ./PCA/simN$simN.$popC.POPin$inputN.gen$generationN.forPCA
      plink2 --bfile ./PCA/simN$simN.$popC.POPin$inputN.gen$generationN.forPCA --pca-aprox --out ./PCA/PCA.simN$simN.$popC.POPin$inputN.gen$generationN

      

find ./simulatedWGS/ -name "simN$simN.*.POPin$inputN.gen$generationN.notri.bin.log" | sort | sed "s/\.log//" > ./simmerge.simN$simN.POPin$inputN.gen$generationN.list.txt
      plink2 --pmerge-list ./simmerge.simN$simN.POPin$inputN.gen$generationN.list.txt --make-pgen --out ./simulatedWGS/mergeforPCA.simN$simN.POPin$inputN.gen$generationN
      


```
      cp ./simulatedWGS/simN$simN.$popC.POPin$inputN.gen$generationN.notri.bin.psam ./simulatedWGS/simN$simN.$popC.POPin$inputN.gen$generationN.notri.bin.psam.cp
      awk -v FS="\t" '{{print "'$popC'""_simN""'$simN'",$2,$3}}' ./simulatedWGS/simN$simN.$popC.POPin$inputN.gen$generationN.notri.bin.psam.cp > ./simulatedWGS/simN$simN.$popC.POPin$inputN.gen$generationN.notri.bin.psam
      sed -i "1c #FID \t IID \t SEX" ./simulatedWGS/simN$simN.$popC.POPin$inputN.gen$generationN.notri.bin.psam


*******************************************************************************************************************
*******************************************************************************************************************
*******************************************************************************************************************
*******************************************************************************************************************


### Structure analysis

### Duplicated SNPs removal

### Kinship analysis

### HWE check (common variants)

