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


# 6. Gene classification

I created protein files for HSP40, HSP70, and HSP90 from uniprot and NCBI. I used both Calliphoridae and Dipteran genes (specific and more broad). 

An example of my search was "HSP40" OR "DNAJ" AND "Calliphoridae" NOT partial NOT low quality". I also did "HSP40" OR "DNAJ" AND "Diptera" NOT partial NOT low quality" and combined everything (i.e., both Diptera + Calliphoridae for each database) into one fasta file per gene family. But when I ran the below script, I ran a couple of the 'classified hsps' in BLAST and they came up with some weird bacterial/fungal stuff. I think it is better to create the .hmm based on domains. For HSP40, these are J-domains. 

I did this instead:

## HMM curation ( this is to ensure everything that we predict has a J domain; I only did this for hsp40 -- hsp70 and 90 I used the full proteins)
Download the PFAM hmm library:
```
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
```

locate the domains in each of the sequences in the fasta file:
```
hmmscan \
  --cpu 6 \
  --domtblout diptera_jdomain.domtblout \
  Pfam-A.hmm \
  known_hsp40.fasta
```
this produces the file diptera_jdomain.domtblout

Then extract the protein sequence ids and remove duplicates:
```
grep -v '^#' diptera_jdomain.domtblout | awk '$2 ~ /PF00226/ {print $4}' | sort | uniq > jdomain_ids.txt
```
Then use seqkit to extract only the sequences that have PF00226 domain in the protein sequence id
```
seqkit grep -f jdomain_ids.txt known_hsp40.fasta > hsp40_jdomain.fasta
```
Remove duplicate records
```
seqkit rmdup -n hsp40_jdomain.fasta > hsp40_jdomain_dedup.fasta
```
Extract J-domains with length filter
```
grep -v '^#' diptera_jdomain.domtblout | \
awk '$2 ~ /PF00226/ {
  start=$18;
  end=$19;
  len=end-start+1;
  if(len>=60 && len<=80)
    print $4 "\t" start "\t" end
}' > jdomain_coords.tsv
```
Extract the sequences
```
seqkit subseq \
  --bed jdomain_coords.tsv \
  hsp40_jdomain_dedup.fasta \
  > jdomains_60_80.fasta
```
HPD filter
```
seqkit grep -s -p HPD jdomains_60_80.fasta > jdomains_filtered.fasta
```
Remove duplicates
```
seqkit rmdup -s jdomains_filtered.fasta > jdomains_nr.fasta
```
Multiple sequence alignment
```
mafft --auto jdomains_nr.fasta > jdomains_aligned.fasta
```
build the hmm 
```
hmmbuild hsp40.hmm jdomains_aligned.fasta
```


```
#!/bin/bash -e
#SBATCH --account=uow03920
#SBATCH --job-name=HSP40
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paige.matheson14@gmail.com
#SBATCH --output HSP40_%j.out
#SBATCH --error HSP40_%j.err

# Purge any loaded modules to avoid conflicts
module purge

# Load required modules
ml BLAST/2.16.0-GCC-12.3.0 SeqKit/2.4.0 HMMER/3.4-GCC-12.3.0 CD-HIT/4.8.1-GCC-11.3.0 BLASTDB GCC/12.3.0 

for i in 01_hilli 02_quadrimaculata 03_stygia 04_vicina 06_cuprina; do
    cd ${i}

    echo "Processing ${i} ..."

    # Step 1: HMM search directly on proteome
    hmmsearch --cpu 4 \
              --domtblout ${i}_hmmer_hsp40.domtblout \
              ../hsp40.hmm \
              ${i}_longest_isoforms_modified.faa

    # Step 2: Filter significant hits
    grep -v '^#' ${i}_hmmer_hsp40.domtblout | \
    awk '$13 <= 1e-3 {print $1}' | sort | uniq > ${i}_pfam_confirmed_ids.txt

    # Step 3: Extract sequences
    seqkit grep -f ${i}_pfam_confirmed_ids.txt ${i}_longest_isoforms_modified.faa \
        > ${i}_confirmed_HSP40.fasta

    # Step 4: Remove redundancy
    cd-hit -i ${i}_confirmed_HSP40.fasta \
           -o ${i}_HSP40_final_cdhit.faa \
           -c 0.98

    cd ../
done
```







Filtered the fasta a bit just to remove duplicates (since we used two databases, probably overlap) and by HPD motif

```
seqkit grep -r -p HPD jdomains_60_80.fasta > jdomains_filtered.fasta
```
```
seqkit rmdup -s jdomains_filtered.fasta > jdomains_nr.fasta
```




















