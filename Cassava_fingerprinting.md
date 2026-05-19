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


## Run KING. Don't do any QC filtering yet.
```bash
conda install -c bioconda king
```

### Convert VCF to PLINK binary

```bash
grep -v "^Unknown" Report_DCas22-7517_SNP_mapping_2_sorted.vcf | \
awk '
BEGIN {OFS="\t"}
/^#/ {print; next}
/^Chromosome01/ {$1="1"; print; next}
/^Chromosome02/ {$1="2"; print; next}
/^Chromosome03/ {$1="3"; print; next}
/^Chromosome04/ {$1="4"; print; next}
/^Chromosome05/ {$1="5"; print; next}
/^Chromosome06/ {$1="6"; print; next}
/^Chromosome07/ {$1="7"; print; next}
/^Chromosome08/ {$1="8"; print; next}
/^Chromosome09/ {$1="9"; print; next}
/^Chromosome10/ {$1="10"; print; next}
/^Chromosome11/ {$1="11"; print; next}
/^Chromosome12/ {$1="12"; print; next}
/^Chromosome13/ {$1="13"; print; next}
/^Chromosome14/ {$1="14"; print; next}
/^Chromosome15/ {$1="15"; print; next}
/^Chromosome16/ {$1="16"; print; next}
/^Chromosome17/ {$1="17"; print; next}
/^Chromosome18/ {$1="18"; print; next}
/^Scaffold/ {gsub(/^Scaffold/, "100", $1); print; next}
' > Report_DCas22-7517_SNP_mapping_2_sorted_Names.vcf
```

Get plink formatted files for KING:
```bash
plink --vcf output/Report_DCas22-7517_SNP_mapping_2_sorted_Names.vcf \
  --make-bed \
  --out output/Report_DCas22-7517_SNP_mapping_2_sorted_Names \
  --allow-extra-chr \
  --const-fid 0
```

# Run KING robust kinship
```bash
king -b output/Report_DCas22-7517_SNP_mapping_2_sorted_Names.bed --kinship --prefix output/Report_DCas22-7517_SNP_mapping_2_sorted_Names_KING
```
```bash
king -b output/Report_DCas22-7517_SNP_mapping_2_sorted_Names.bed --ibdseg --prefix output/Report_DCas22-7517_SNP_mapping_2_sorted_Names_KING_ibdseg
```
Total length of 4 chromosomal segments usable for IBD segment analysis is 88.6 Mb.
  Information of these chromosomal segments can be found in file Report_DCas22-7517_SNP_mapping_2_sorted_Names_KING_ibdsegallsegs.txt

Segments too short.

# Examine *kin file
```bash
awk '$9 > 0.35 {print}' Report_DCas22-7517_SNP_mapping_2_sorted_Names_KING.kin | wc -l

#60013

awk '$9 > 0.17 && $9 < 0.35 {print}' Report_DCas22-7517_SNP_mapping_2_sorted_Names_KING.kin | wc -l

#9213

awk '$9 > 0.08 && $9 < 0.17 {print}' Report_DCas22-7517_SNP_mapping_2_sorted_Names_KING.kin | wc -l

#12687

awk '$10 == 1 {print}' Report_DCas22-7517_SNP_mapping_2_sorted_Names_KING.kin | head -5

#0       BALICol1998     MERCURY 5991    1.000   0.0000  0.1482  0.0030  0.4320  1
#0       BALICol1998     Nase5   5832    1.000   0.0000  0.0746  0.0271  0.0556  1
#0       BALICol1998     NASE8   5858    1.000   0.0000  0.0775  0.0280  0.0583  1
#0       BALICol1998     UG110052        5905    1.000   0.0000  0.1429  0.0056  0.4071  1
#0       BALICol1998     600021  5858    1.000   0.0000  0.0848  0.0347  0.0450  1
```
