# comparative

# 1. Select the longest representative isoforms per gene in the protein file 

```
#!/bin/bash
# List of sample IDs
samples=(07_vicina_hic 08_vicina_longstitch)
# Loop through each sample
for sample in "${samples[@]}"
do
    echo "Processing $sample..."

    python3 <<EOF
from Bio import SeqIO
from collections import defaultdict

input_fasta = "/nesi/nobackup/uow03920/06_blowfly_assembly_jan/11_braker/${sample}/braker/braker.aa"
output_fasta = "${sample}_longest_isoforms.fa"

# Store longest isoform per gene
longest = defaultdict(lambda: ("", 0))  # gene_id → (record, length)

for record in SeqIO.parse(input_fasta, "fasta"):
    gene_id = record.id.split(".")[0]
    if len(record.seq) > longest[gene_id][1]:
        longest[gene_id] = (record, len(record.seq))

with open(output_fasta, "w") as out:
    for rec, _ in longest.values():
        SeqIO.write(rec, out, "fasta")

print(f"Finished ${sample}")
EOF
done
```

# 2. Run OrthoFinder
I provided orthofinder with my own Newick tree because OrthoFinder wasn't producing an accurate species tree. It probably wasn't accurate because blowflies are closely related and we didn't have a huge amount of species included.

Orthofinder is a fast, accurate, and comprehensive platform used to infer orthogroups (genes descended from a single gene in the last common ancestor), orthologs, rooted gene trees, and rooted species trees. It is essential for comparative genomics to identify gene duplication events, analyse gene family evolution, and compare protein sets. 

My newick tree was called tree.nwk and was in the same folder as OrthoFinder

```
(06_cuprina_longest_isoforms_modified,(((01_hilli_longest_isoforms_modified,03_stygia_longest_isoforms_modified),02_quadrimaculata_longest_isoforms_modified),04_vicina_longest_isoforms_modified));
```

Run OrthoFinder
```
#!/bin/bash -e
#SBATCH --account=uow03920
#SBATCH --job-name=orthofinder
#SBATCH --time=172:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paige.matheson14@gmail.com
#SBATCH --output orthofinder_%j.out
#SBATCH --error orthofinder_%j.err

module purge
ml OrthoFinder/2.5.2 MAFFT/7.505 IQ-TREE/2.2.2.2

# Run:
orthofinder -f proteins -M msa -T iqtree -t 8 -a 4
```

proteins is a folder containing the longest isoform protein files. 


# 3. Gene family expansion and contraction 
Convert SpeciesTree (from OrthoFinder) to Ultrametric tree (chronos method). This is done in R. 

```
library(ape)
# Load the rooted species tree
tree <- read.tree("SpeciesTree_rooted.txt")

# Check that it's rooted and binary
stopifnot(is.rooted(tree))
stopifnot(is.binary(tree))

# Convert to ultrametric with a relaxed clock
ultra_tree <- chronos(tree, lambda = 1, model = "correlated")

# Scale tree so root-to-tip distance = 17 Mya
tree_age <- max(node.depth.edgelength(ultra_tree))
ultra_tree$edge.length <- ultra_tree$edge.length * (17 / tree_age)

# Plot and save
plot(ultra_tree, main = "Ultrametric Tree (Root Scaled to 17 Mya)")
write.tree(ultra_tree, file = "SpeciesTree_ultrametric_17Mya.txt")  
```

# Filter orthogroup gene count file (from Orthofinder) and modify it for CAFE 

```
awk 'OFS="\t" {$NF=""; print}' Orthogroups.GeneCount.tsv > tmp \
&& awk '{print "(null)""\t"$0}' tmp > cafe.input.tsv \
&& sed -i '1s/(null)/Desc/g' cafe.input.tsv \
&& rm tmp

python2 cafetutorial_clade_and_size_filter.py \
    -i cafe.input.tsv -o filtered.cafe.input.tsv -s 2
```

# Run CAFE using that input
```
#!/bin/bash -e
#SBATCH --account=uow03920
#SBATCH --job-name=CAFE
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paige.matheson14@gmail.com
#SBATCH --output CAFE_%j.out    # save the output into a file
#SBATCH --error CAFE_%j.err     # save the error output into a file

#purge
module purge

ml GCC/12.3.0 

#Set variables
CAFE5_EXEC="/nesi/nobackup/uow03920/06_blowfly_assembly_jan/18_cafe/CAFE5/build/cafe5"
CAFE_DIR="/nesi/nobackup/uow03920/06_blowfly_assembly_jan/25_cafe/"
ORTHO_DIR="/nesi/nobackup/uow03920/06_blowfly_assembly_jan/15_orthofinder/proteins_feb4/OrthoFinder/Results_Feb04"
TREE="/nesi/nobackup/uow03920/06_blowfly_assembly_jan/15_orthofinder/proteins_feb4/OrthoFinder/Results_Feb04/Species_Tree/SpeciesTree_ultrametric_17Mya.txt"
COUNTS="/nesi/nobackup/uow03920/06_blowfly_assembly_jan/25_cafe/filtered.cafe.input.tsv"
K_VALS=(1 2 3 4 5 6)
REPLICATES=(1 2)
PVALUE_THRESHOLD=0.05

#Step 1: Run CAFE5 with multiple k values and replicates
for k in "${K_VALS[@]}"; do
    for rep in "${REPLICATES[@]}"; do
        OUTDIR="$CAFE_DIR/k${k}_rep${rep}"
        mkdir -p "$OUTDIR"
        echo "Running CAFE5 for k=$k, replicate=$rep..."
        $CAFE5_EXEC -i "$COUNTS" -t "$TREE" -o "$OUTDIR" -k $k > "$OUTDIR/stdout_log.txt"
    done
done

#Step 2: Extract log-likelihoods, check convergence, find best k
best_k=""
best_lnl=999999

echo -e "\n--- Log-likelihoods and convergence check ---"
for k in "${K_VALS[@]}"; do
    lnls=()
    for rep in "${REPLICATES[@]}"; do
        OUTDIR="$CAFE_DIR/k${k}_rep${rep}"
        LNL=$(grep "Final -lnL:" "$OUTDIR/stdout_log.txt" | sed -E 's/.*Final -lnL: //')
        echo "k=$k, rep=$rep: lnL=$LNL"
        lnls+=("$LNL")
done
```

# 4. Plot lineage specific gene gains and losses 

# 5 GO enrichment analysis 
Generate background gene sets using filtered.cafe.input.csv and Orthogroups.tsv

To assign GO terms, (1) prepare merged annotation files from EGGNOG and Interproscan annotations (see folder in this repository for R scripts) (2) prepare background gene sets with GO terms

Perform TOPGO enrichment analysis - use TOPGO.R script





























