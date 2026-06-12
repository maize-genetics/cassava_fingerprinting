# Cassava fingerprinting analyses
##### June 2026
##### bfe4@cornell.edu

## Raw files
Starts with DArTseq-LD data in SNP_mapping_2.csv format
- Report_DCas22-7517_SNP_mapping_2.csv

## Code
- code/01.mapping2vcf.R: csv needs to be converted to a vcf file
- code/02.ibs0vskinship.R: analyze PLINK and King results
- code/03.knowledge_graph.R: make knowledge graph 
- code/04.makeDosageMartix.R: convert PLINK dosage file to matrix (for other software program testing)

    Example useage: 

### Convert to VCF format
```bash
Rscript code/mapping2vcf.R -m data/OrderAppendix_2_DCas22-7517/Report_DCas22-7517_SNP_mapping_2.csv -o output/Report_DCas22-7517_SNP_mapping_2.vcf -p 2
```

### Sort
```bash
bgzip output/Report_DCas22-7517_SNP_mapping_2.vcf
tabix -p vcf output/Report_DCas22-7517_SNP_mapping_2.vcf.gz
bcftools sort output/Report_DCas22-7517_SNP_mapping_2.vcf.gz \
-o output/Report_DCas22-7517_SNP_mapping_2_sorted.vcf.gz
```

### Convert VCF to PLINK binary

```bash (deal with contigs for PLINK)
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
/^Scaffold/ {
    if (!seen[$1]) {  # If scaffold hasnt been seen before
        scaffold_map[$1] = "contig" contig_counter
        seen[$1] = 1
        contig_counter++
    }
    $1 = scaffold_map[$1]
    print
    next
}
' output/Report_DCas22-7517_SNP_mapping_2_sorted.vcf > output/Report_DCas22-7517_SNP_mapping_2_sorted_Names.vcf
```

### Get filtered vcf:
```bash
plink --vcf output/Report_DCas22-7517_SNP_mapping_2_sorted_Names.vcf \
      --allow-extra-chr \
      --mind 0.2 \
      --geno 0.2 \
      --maf 0.01 \
      --recode vcf \
      --out output/Report_DCas22-7517_SNP_mapping_2_sorted_Names_geno0.2_mind0.2_maf0.01
      --genome 
```

### Make King table in PLINK
```bash
./plink2 --allow-extra-chr --bfile output/Report_DCas22-7517_SNP_mapping_2_sorted_Names_geno0.2_mind0.2_maf0.01 --make-king-table --out output/plinkAndKing_geno0.2_mind0.2_maf0.01
```

### Use R to analyze PLINK and KING results
- code/02.ibs0vskinship.R
- code/03.knowledge_graph.R

### Extras:

### To get PLINK dosage file
```bash
plink --vcf output/Report_DCas22-7517_SNP_mapping_2_sorted_Names_geno0.2_mind0.2_maf0.01.vcf \
      --allow-extra-chr \
      --export A-transpose \
      --out output/Report_DCas22-7517_SNP_mapping_2_sorted_Names_geno0.2_mind0.2_maf0.01_dosage
```

### Convert to dosage matrix 
```bash
Rscript code/makeDosageMartix.R -i output/Report_DCas22-7517_SNP_mapping_2_sorted_Names_geno0.2_mind0.2_maf0.01_dosage.traw -o output/Report_DCas22-7517_SNP_mapping_2_sorted_Names_geno0.2_mind0.2_maf0.01_dosageMatrix.txt 
```

## Figures 
- figures/maximum_reference_relationships_plot.png (maximum pairwise relationship between a farm x reference sample)
- figures/maximum_unconnected_relationships_plot.png (relationships among samples that did not have a strong relationship with a reference sample in the first analysis)
- figures/reference_strong_relationships.png (relationships amoung reference samples)

## Network visualization
[View knowledge graph of strongest farm x reference and reference x reference pairwise relationships](https://maize-genetics.github.io/cassava_fingerprinting/cassava_knowledge_graph.html)
- Shows farm samples and reference varieties connected by genetic relationships
- Green edges: vegetative clones
- Red edges: first degree relationships  
- Orange edges: second degree relationships
- Blue edges: reference-reference relationships

## Files generated from R scripts (not in git repo)
- `complete_reference_stats.csv` - Reference variety statistics
- `max_ref_relationships.csv` - Maximum farm-reference pairwise relationships
- `ref_ref_relationships.csv` - Reference x reference relationships  
- `max_unconnected_relationships.csv` - Relationships among unconnected samples (farm samples without a clonal, 1st or 2nd degree relationship with a reference sample).

## Additional notes

Attempted to find additional reference samples to match unmatched farm samples; however, VCF files did not have overlapping positions. Can return to analysis later with imputted vcf files. 

### See if any same samples as reference from fingerprinting

### Get data into match v7 coordinates
```bash
export PYTHONPATH=/programs/CrossMap-0.7.3/lib64/python3.9/site-packages:/programs/CrossMap-0.7.3/lib/python3.9/site-packages
export PATH=/programs/CrossMap-0.7.3/bin:$PATH
sed 's/>Chromosome0*\([0-9]*\)$/>\1/' Mesculenta_520_v7.fa > Mesculenta_520_v7_Names.fa
CrossMap vcf Mesculenta_305_v6.to_v7.final.numeric.chain.gz DCas19_4459.vcf.gz Mesculenta_520_v7_Names.fa.gz DCas19_4459_v7.vcf
vcf-sort -c DCas19_4459_v7.vcf > DCas19_4459_v7_sorted.vcf
CrossMap vcf Mesculenta_305_v6.to_v7.final.numeric.chain.gz DCas19_4459_82719.vcf.gz Mesculenta_520_v7_Names.fa.gz DCas19_4459_82719_v7.vcf
vcf-sort -c DCas19_4459_82719_v7.vcf > DCas19_4459_82719_v7_sorted.vcf
```

### Check if any of the same samples between files
```bash
bcftools query -l DCas19_4459_82719_v7_sorted.vcf | grep -f ~/Documents/cassava_fingerprinting/reference_only.txt
#Mkumba
#NASE14
#UG120024
#UG120156
#UG120183
#UG120193
bcftools query -l DCas19_4459_v7_sorted.vcf | grep -f ~/Documents/cassava_fingerprinting/reference_only.txt
#Mkumba 
#NASE14 
#UG120024
#UG120156
#UG120183
#UG120193
```

### Check for overlaps
```bash
bcftools isec -p tmp_dir cleaned_DCas19_4459_v7_sorted_labeled_markerIDs_fixref.vcf.gz cleaned_Report_DCas22-7517_SNP_mapping_2_sorted_Names_labeled_markerIDs_fixref.vcf.gz
bcftools isec -p tmp_dir cleaned_Report_DCas22-7517_SNP_mapping_2_sorted_Names_labeled_markerIDs_fixref.vcf.gz DCas19_4459_v7_sorted.vcf.gz
```