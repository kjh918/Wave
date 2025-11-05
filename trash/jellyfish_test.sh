# /storage/home/jhkim/Apps/jellyfish-2.3.1/bin/jellyfish count -m 21 -s 100M -t 8 -C -o read1.jf cbNIPT_25_03_02_R1.fasta 
/storage/home/jhkim/Apps/jellyfish-2.3.1/bin/jellyfish histo -t 16 /storage/home/jhkim/Projects/cbNIPT/GCX-cbNIPT_QC_Setup-2025-10-23/Results/cbNIPT_25_03_02/01_Trim_Fastp/mer_counts.jf > /storage/home/jhkim/Projects/cbNIPT/GCX-cbNIPT_QC_Setup-2025-10-23/Results/cbNIPT_25_03_02/01_Trim_Fastp/mer_count.histo 

D=$(awk '$1>1{if($2>m){m=$2;d=$1}}END{print d}' /storage/home/jhkim/Projects/cbNIPT/GCX-cbNIPT_QC_Setup-2025-10-23/Results/cbNIPT_25_03_02/01_Trim_Fastp/mer_count.histo 
) 
THR=$(awk -v D=$D 'BEGIN{print int(D/2)}') 
S=$(awk -v T=$THR '$1>=T{sum+=$1*$2}END{print sum}' /storage/home/jhkim/Projects/cbNIPT/GCX-cbNIPT_QC_Setup-2025-10-23/Results/cbNIPT_25_03_02/01_Trim_Fastp/mer_count.histo 
) 
awk -v S=$S -v D=$D 'BEGIN{printf "Genome size ≈ %.0f bp\n", S/D}'