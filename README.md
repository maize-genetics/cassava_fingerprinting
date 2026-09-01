# Cassava fingerprinting analyses
##### June 2026
##### bfe4@cornell.edu

[Sneak peak, described below](https://maize-genetics.github.io/cassava_fingerprinting/cassava_knowledge_graph.html)

## Raw files
Starts with DArTseq-LD data in SNP_mapping_2.csv format
- Report_DCas22-7517_SNP_mapping_2.csv

## Code
- code/01.mapping2vcf.R: csv needs to be converted to a vcf file
- code/02.ibs0vskinship.R: analyze PLINK and KING results, classify relationships, identify clones and unmatched samples
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

**Note on sample name formatting:** Several reference variety names contain underscores (e.g., `OFUMBA_CHAI`, `NJULE_RED`, `MASAKA_LOCAL-1`, `MASAKA_LOCAL-2`). PLINK/KING interpret the text before the underscore as the Family ID (FID) and the text after as the Individual ID (IID), splitting these names apart in the raw output. `code/02.ibs0vskinship.R` includes a correction step (`fix_split_names()`) to reconstruct the full sample names before downstream analysis.

## Relationship classification approach

We use **KING's robust kinship coefficient** (Φ) as the primary metric for classifying pairwise relationships, since it is robust to population structure/stratification (important for a diverse cassava germplasm panel mixing reference varieties and farmer-maintained samples):

| Kinship (Φ) | Relationship |
|---|---|
| ≥ 0.36 | Highly related / clones* |
| ≥ 0.19 | First degree (parent-offspring or full siblings) |
| ≥ 0.088 | Second degree |
| ≥ 0.044 | Third degree |
| < 0.044 | Fourth degree / unrelated |

\* We use "highly related / clones" rather than a definitive "clones" label because kinship + IBS0 alone cannot distinguish a true vegetative clone from a highly inbred relative (e.g., progeny of self-pollination). Distinguishing these would require pedigree records or additional analysis.

**IBS0** (the proportion of SNPs where two samples are opposite homozygotes, e.g., AA vs. aa) is plotted alongside kinship as a visual/diagnostic check, since it is one of the raw inputs to the kinship calculation itself and helps flag borderline cases or potential genotyping issues.

### Why not IBD segment analysis?
We attempted KING's IBD segment detection (`--ibdseg`) as a potentially more precise alternative, but our DArTseq-LD marker density was insufficient:
```
king -b output/Report_DCas22-7517_SNP_mapping_2_sorted_Names.bed --ibdseg --prefix output/cassava_ibdseg
#Total length of 4 chromosomal segments usable for IBD segment analysis is 88.4 Mb.
#Segments too short.
```
Only 88.4 Mb of the genome had sufficient marker density for segment detection (a small fraction of the ~760 Mb cassava genome), so we proceeded with the kinship coefficient + IBS0 approach instead.

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
- figures/maximum_reference_relationships_plot.png: strongest pairwise relationship for each farm sample to any reference variety (farm-reference), plus reference-reference pairs
- figures/reference_strong_relationships.png: strong (clone/1st/2nd/3rd degree) relationships among reference varieties only
- figures/combined_farm_connectivity_plot.png: for farm samples lacking a strong reference match, shows their true strongest connection (searched across the full dataset) - split into farms with a farm-farm clone/close match vs. truly isolated farms
- figures/isolated_reference_varieties_table.png: reference varieties lacking any clone/1st/2nd degree match, with their true strongest connection (table)
- figures/maximum_reference_relationships_inbreeding_plot.png: maximum farm-reference relationships, plotting KING kinship coefficient against inbreeding coefficient (F) of the farm sample
- figures/kinship_vs_deltaF_all_relationships.png: maximum farm-reference relationships, plotting kinship against ΔF (F_farm − F_reference), across all relationship categories
- figures/F_farm_by_relationship_category.png: boxplot of farm sample inbreeding coefficient (F) using the maximum farm-reference relationship, grouped by relationship category to closest reference variety
- figures/ibs0_clone_vs_firstdegree_comparison.png: comparison of IBS0 (proportion of opposite-homozygote genotypes) between "Highly related/clones" and "First degree" pairs (all pairwise relationships in the dataset, not just each sample's maximum match). Used to evaluate whether first-degree matches are consistent with direct parent-offspring relationships (genetically constrained to IBS0≈0) versus siblings sharing an unsampled parent (not similarly constrained). First-degree pairs show IBS0 roughly 4x higher than the clone-level floor (Welch's t-test, p < 2.2×10⁻¹⁶), consistent with siblings rather than direct parent-offspring.


## Network visualization
[View knowledge graph of strongest farm x reference and reference x reference pairwise relationships](https://maize-genetics.github.io/cassava_fingerprinting/figures/cassava_knowledge_graph_v2.html)
- Shows farm samples and reference varieties connected by genetic relationships
- Triangle = Reference variety, Circle = Farm sample
- Red edges: vegetative clones
- Blue edges: first degree relationships  
- Pink edges: second degree relationships

Additional network views:
- [Farm-farm clonal relationships (samples lacking a reference match)](https://maize-genetics.github.io/cassava_fingerprinting/figures/red_only_farms_knowledge_graph.html)
- [Reference x Reference relationships (clone through 3rd degree)](https://maize-genetics.github.io/cassava_fingerprinting/figures/ref_ref_clone_1st_2nd_3rd_degree_graph.html)

## Files generated from R scripts (not in git repo)
- `max_ref_relationships.csv` - Maximum farm-reference pairwise relationships (one row per farm sample)
- `ref_ref_relationships.csv` - All reference x reference pairwise relationships
- `complete_reference_stats.csv` - Summary statistics per reference variety (clone counts, relationship type breakdown, etc.)
- `true_max_unmatched_farms.csv` - True strongest connection (searched across full dataset) for farm samples lacking a strong reference match
- `farms_with_farm_farm_clones.csv` - Subset of the above: farms with a strong farm-farm match despite no reference match (likely unlabeled duplicate/sister clones)
- `truly_isolated_farms.csv` - Subset of the above: farms with no strong match to anything (reference or farm)
- `true_max_isolated_references.csv` / `isolated_reference_varieties_true_max.csv` - Reference varieties lacking a clone/1st/2nd degree match, with their true strongest connection

## Files generated from R scripts (not in git repo)
- `max_ref_relationships.csv` - Maximum farm-reference pairwise relationships (one row per farm sample)
- `ref_ref_relationships.csv` - All reference x reference pairwise relationships
- `complete_reference_stats.csv` - Summary statistics per reference variety (clone counts, relationship type breakdown, etc.)
- `max_unmatched_farms.csv` - True strongest connection (searched across full dataset) for farm samples lacking a strong reference match
- `farms_with_farm_farm_clones.csv` - Subset of the above: farms with a strong farm-farm match despite no reference match (likely unlabeled duplicate/sister clones)
- `truly_isolated_farms.csv` - Subset of the above: farms with no strong match to anything (reference or farm)
- `true_max_isolated_references.csv` / `isolated_reference_varieties_true_max.csv` - Reference varieties lacking a clone/1st/2nd degree match, with their true strongest connection
- `cytoscape_edges.csv` / `cytoscape_nodes.csv` - Node/edge lists formatted for standalone Cytoscape desktop app (not currently used for the visNetwork-based figures/graphs above). Kept for potential future use, since Cytoscape can use edge weight (e.g., kinship) to directly influence layout/branch length, unlike the force-based physics layout currently used in visNetwork.

**Important methodological note:** When identifying "unmatched"/"isolated" samples, we search each unmatched sample's strongest connection across the **entire dataset**, not just within a restricted subset of other unmatched samples. An earlier version of this analysis restricted the search to pairs where *both* samples lacked a reference match, which incorrectly excluded valid strong matches (e.g., a farm sample's true clone partner might itself have a strong reference match and would otherwise be excluded from consideration).

## Additional notes

Attempted to find additional reference samples to match unmatched farm samples; however, VCF files did not have overlapping positions. Can return to analysis later with imputed vcf files. 

### Get data into match v7 coordinates
```bash
export PYTHONPATH=/programs/CrossMap-0.7.3/lib64/python3.9/site-packages:/programs/CrossMap-0.7.3/lib/python3.9/site-packages
export PATH=/programs/CrossMap-0.7.3/bin:$PATH
sed 's/>Chromosome0*\([0-9]*\)$/>\1/' Mesculenta_520_v7.fa > Mesculenta_520_v7_Names.fa
CrossMap vcf Mesculenta_305_v6.to_v7.final.numeric.chain.gz DCas19_4459.vcf.gz M
