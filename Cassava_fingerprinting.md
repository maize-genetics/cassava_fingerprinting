# Cassava fingerprinting
bfe4@cornell.edu
May 18, 2026
DArTseq-LD data from D.Gimode@cgiar.org bevis.16@osu.edu
DArT encoding for mapping file: 0=ref homo, 1=alt homo, 2=hetero

## Step 1: Overview

Convert DArT SNP mapping files through the following steps:
1. **DArT SNP Mapping** → **VCF** (properly formatted with sample name cleaning)
2. **VCF** → **PLINK dosage** (minor allele dosage for kinship analysis)  
3. **PLINK dosage** → **Formatted dosage matrix** (samples × markers)

### Make vcf file from mapping file

```bash
Rscript code/mapping2vcf.R -m data/OrderAppendix_2_DCas22-7517/Report_DCas22-7517_SNP_mapping_2.csv -o output/Report_DCas22-7517_SNP_mapping_2.vcf -p 2
```
### Verify 

```bash
Rscript code/validateVCF_mapping.R -v output/Report_DCas22-7517_SNP_mapping_2.vcf -m data/OrderAppendix_2_DCas22-7517/Report_DCas22-7517_SNP_mapping_2.csv -s 30
```

### Sort

```bash
bgzip output/Report_DCas22-7517_SNP_mapping_2.vcf
tabix -p vcf output/Report_DCas22-7517_SNP_mapping_2.vcf.gz
bcftools sort output/Report_DCas22-7517_SNP_mapping_2.vcf.gz \
  -o output/Report_DCas22-7517_SNP_mapping_2_sorted.vcf.gz
```

### Convert VCF to PLINK dosage format
```bash
plink --vcf output/Report_DCas22-7517_SNP_mapping_2.vcf.gz \
      --allow-extra-chr \
      --export A-transpose \
      --out output/Report_DCas22-7517_SNP_mapping_2plink_dosage
```
note:  PLINK is converting to Minor Allele Dosage (not alternative allele)

### Convert to dosage matrix

```bash
Rscript code/makeDosageMartix.R -i output/Report_DCas22-7517_SNP_mapping_2plink_dosage.traw -o output/Report_DCas22-7517_SNP_mapping_2plink_dosageMatrix.txt 
```

### Validate marker
```bash
Rscript code/validateMarker.R
```

### Check the full marker info for both
```bash
echo "=== MARKER 1 (inconsistent) ==="
grep "7131698|F|0-65:A>C-65:A>C" data/OrderAppendix_2_DCas22-7517/Report_DCas22-7517_SNP_mapping_2.csv | cut -d, -f1,13

echo "=== MARKER 2 (correct) ==="  
grep "15484497|F|0-51:A>C-51:A>C" data/OrderAppendix_2_DCas22-7517/Report_DCas22-7517_SNP_mapping_2.csv | cut -d, -f1,13
```
### Check if allele frequencies are different (count is by minor allele?)
```bash
echo "=== VCF LINES ==="
grep "7131698|F|0-65:A>C-65:A>C" output/Report_DCas22-7517_SNP_mapping_2.vcf | cut -f1-9
grep "15484497|F|0-51:A>C-51:A>C" output/Report_DCas22-7517_SNP_mapping_2.vcf | cut -f1-9

echo "=== MARKER 1 GENOTYPE COUNTS ==="
grep "7131698|F|0-65:A>C-65:A>C" output/Report_DCas22-7517_SNP_mapping_2.vcf | \
  tr '\t' '\n' | tail -n +10 | sort | uniq -c

echo "=== MARKER 2 GENOTYPE COUNTS ==="  
grep "15484497|F|0-51:A>C-51:A>C" output/Report_DCas22-7517_SNP_mapping_2.vcf | \
  tr '\t' '\n' | tail -n +10 | sort | uniq -c

#=== MARKER 1 GENOTYPE COUNTS ===
#145 ./.
#612 0/0
#416 0/1
#831 1/1
#=== MARKER 2 GENOTYPE COUNTS ===
#1 ./.
#1401 0/0
#56 0/1
#546 1/1
```