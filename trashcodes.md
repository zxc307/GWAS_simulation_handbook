
for simN in {1..1..1}
do
  for generationN in {30,100,300}
  do
    for inputN in {30,60,90}
    do
      for chrN in {1..1..1}
      do
      let chrm=$chrN+1
      plink --bfile ./$popC/$popC.POP$inputN/sim$simN.gen$generationN.chr$chrN.notri.bin --bmerge ./$popC/$popC.POP$inputN/sim$simN.gen$generationN.chr$chrm.notri.bin --make-bed --out ./simulatedWGS/simN$simN.$popC.POPin$inputN.gen$generationN.notri.merge$chrN
      plink --bfile ./$popC/$popC.POP$inputN/sim$simN.gen$generationN.chr$chrN.notri.bin --exclude ./simulatedWGS/simN$simN.$popC.POPin$inputN.gen$generationN.notri.merge$chrN-merge.missnp --make-bed --out ./$popC/$popC.POP$inputN/sim$simN.gen$generationN.chr$chrN.notri.bin
      plink --bfile ./$popC/$popC.POP$inputN/sim$simN.gen$generationN.chr$chrm.notri.bin --exclude ./simulatedWGS/simN$simN.$popC.POPin$inputN.gen$generationN.notri.merge$chrN-merge.missnp --make-bed --out ./$popC/$popC.POP$inputN/sim$simN.gen$generationN.chr$chrm.notri.bin
      plink --bfile ./$popC/$popC.POP$inputN/sim$simN.gen$generationN.chr$chrN.notri.bin --bmerge ./$popC/$popC.POP$inputN/sim$simN.gen$generationN.chr$chrm.notri.bin --make-bed --out ./simulatedWGS/simN$simN.$popC.POPin$inputN.gen$generationN.notri.merge$chrm
      done
      for chrN in {2..21..1}
      do
      let chrm=$chrN+1
      plink --bfile ./simulatedWGS/simN$simN.$popC.POPin$inputN.gen$generationN.notri.merge$chrN --bmerge ./$popC/$popC.POP$inputN/sim$simN.gen$generationN.chr$chrm.notri.bin --make-bed --out ./simulatedWGS/simN$simN.$popC.POPin$inputN.gen$generationN.notri.merge$chrm
      plink --bfile ./simulatedWGS/simN$simN.$popC.POPin$inputN.gen$generationN.notri.merge$chrN --exclude ./simulatedWGS/simN$simN.$popC.POPin$inputN.gen$generationN.notri.merge$chrm-merge.missnp --make-bed --out ./simulatedWGS/simN$simN.$popC.POPin$inputN.gen$generationN.notri.merge$chrm
      plink --bfile ./$popC/$popC.POP$inputN/sim$simN.gen$generationN.chr$chrm.notri.bin --exclude ./simulatedWGS/simN$simN.$popC.POPin$inputN.gen$generationN.notri.merge$chrm-merge.missnp --make-bed --out ./$popC/$popC.POP$inputN/sim$simN.gen$generationN.chr$chrm.notri.bin
      done
      plink --bfile ./simulatedWGS/simN$simN.$popC.POPin$inputN.gen$generationN.notri.merge$chrm --bmerge ./$popC/$popC.POP$inputN/sim$simN.gen$generationN.chr$chrm.notri.bin --make-bed --out ./simulatedWGS/simN$simN.$popC.POPin$inputN.gen$generationN.notri.whole
      for chrN in {2..21..1}
      do 
      rm ./simulatedWGS/simN$simN.$popC.POPin$inputN.gen$generationN.notri.merge$chrN*
      done
    done
  done
done


      

      plink --vcf ./$popC/$popC.POP$inputN/sim$simN.gen$generationN.chr$chrN.notri.vcf --make-bed --out ./$popC/$popC.POP$inputN/sim$simN.gen$generationN.chr$chrN.notri.bin
      cp ./$popC/$popC.POP$inputN/sim$simN.gen$generationN.chr$chrN.notri.bin.bim ./$popC/$popC.POP$inputN/sim$simN.gen$generationN.chr$chrN.notri.bin.bim.cp
      awk '{{print $1,"chr"$1"_"$4,$3,$4,$5,$6}}' ./$popC/$popC.POP$inputN/sim$simN.gen$generationN.chr$chrN.notri.bin.bim.cp > ./$popC/$popC.POP$inputN/sim$simN.gen$generationN.chr$chrN.notri.bin.bim


for popC in {YRI,CHB,GBR}
do 
  for generationN in {30,100,300}
  do
    for inputN in {30,60,90}
    do
      for chrN in {9..10..1}
      do
      slim ./$popC/$popC.POP$inputN/slim.code.$popC.POP$inputN.simN$simN.gen$generationN.chr$chrN.txt
      perl -lane 'if(/^#/ or length("$F[3]$F[4]")==2){print}' ./$popC/$popC.POP$inputN/sim$simN.gen$generationN.chr$chrN.vcf > ./$popC/$popC.POP$inputN/sim$simN.gen$generationN.chr$chrN.notri.vcf
      plink --vcf ./$popC/$popC.POP$inputN/sim$simN.gen$generationN.chr$chrN.notri.vcf --make-bed --out ./$popC/$popC.POP$inputN/sim$simN.gen$generationN.chr$chrN.notri.bin
      done
    done
  done
done
