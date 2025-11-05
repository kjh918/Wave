# /storage/home/jhkim/Apps/Python-3.11.13/python /storage/home/jhkim/Projects/cbNIPT/GCX-cbNIPT_QC_Setup-2025-10-23/Results/ADO/classify_ado_bias_imbalance.py \
#   /storage/home/jhkim/Projects/cbNIPT/GCX-cbNIPT_QC_Setup-2025-10-23/Results/ADO/MergedVCF/merged_gvcf.bcftools.chr22.gvcf.gz \
#   --control cbNIPT_25_03_01 \
#   --min-dp-control 10 \
#   --min-dp-other 2 \
#   --min-gq 20 \
#   --ab-het-low 0.3 --ab-het-high 0.7 \
#   --only-biallelic-snp \
#   --site-pass-only \
#   --out-summary summary_per_sample.tsv \
#   --out-discordant discordant_sites.tsv
# /storage/home/jhkim/Apps/Python-3.11.13/python classify_ado_bias_imbalance.py  merged.chr22.allSites.vcf.gz \
#   --control cbNIPT_25_03_01 \
#   --min-dp-control 10 \
#   --min-dp-other 5 \
#   --only-snv \
#   --out-summary summary.tsv \
#   --out-detail details.tsv
/storage/home/jhkim/Apps/Python-3.11.13/python compute_dropout_metrics.py \
  /storage/home/jhkim/Projects/cbNIPT/GCX-cbNIPT_QC_Setup-2025-10-23/Results/ADO/CombineGVCFs.GenotypeCall/CombineGVCFs.GenotypeCall.chr11.vcf.gz \
  --control cbNIPT_25_03_01 \
  --min-dp-control 10 \
  --min-dp-sample 5 \
  --min-gq 20 \
  --only-biallelic-snp \
  --out-metrics /storage/home/jhkim/Projects/cbNIPT/GCX-cbNIPT_QC_Setup-2025-10-23/Results/ADO/CombineGVCFs.GenotypeCall/summary/metrics_per_sample.chr11.tsv \
  --out-details /storage/home/jhkim/Projects/cbNIPT/GCX-cbNIPT_QC_Setup-2025-10-23/Results/ADO/CombineGVCFs.GenotypeCall/summary/details.chr11.tsv \
  --out-genostats /storage/home/jhkim/Projects/cbNIPT/GCX-cbNIPT_QC_Setup-2025-10-23/Results/ADO/CombineGVCFs.GenotypeCall/summary/genotype_stats_per_sample.chr11.tsv \
  --prefix-plots /storage/home/jhkim/Projects/cbNIPT/GCX-cbNIPT_QC_Setup-2025-10-23/Results/ADO/CombineGVCFs.GenotypeCall/summary/dropout_plots.chr11.
