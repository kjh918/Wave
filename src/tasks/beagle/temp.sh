
# bcftools isec -n=2 -w1 -Oz -o overlap.vcf.gz cbNIPT_25_03_01.downsample.p0.0278_0.1x.vcf.gz /storage/home/jhkim/Projects/cbNIPT/GCX-cbNIPT_In-silico_Inputation-2025-11-11/Resources/Reference/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

for chr in {2..2}; do
  java -Xmx16g -jar /storage/home/jhkim/Apps/beagle5.5/beagle.27Feb25.75f.jar \
      gt=cbNIPT_25_03_01_0.1x.chr${chr}.vcf.gz \
      ref=/storage/home/jhkim/Projects/cbNIPT/GCX-cbNIPT_In-silico_Inputation-2025-11-11/Resources/Reference/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.chrom.vcf.gz \
      map=/storage/home/jhkim/Projects/cbNIPT/GCX-cbNIPT_In-silico_Inputation-2025-11-11/Resources/Reference/chr_in_chrom_field/plink.chrchr${chr}.GRCh38.map \
      out=imputed.chr${chr} \
      nthreads=8
done


