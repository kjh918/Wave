
### 
# https://github.com/gymrek-lab/TRTools
# https://trtools.readthedocs.io/en/stable/
### 


singularity exec -B /storage,/data /storage/images/trtools-6.1.0.sif PrancSTR --vcf $vcf --out $prefix --region $region

