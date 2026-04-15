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

# Remove duplicate proteins
I had quite a number of identical sequences in my protein files. Duplicated protein sequences can sometimes arise e.g., if alternative isoforms are treated as distinct proteins (less likely because we selected for the longest isoforms above). They can also arise from fragmented genome assemblies (e.g., a single gene might be split across multiple contigs and braker may predict the same partial gene on both contigs resulting in the same sequence with different headers) or overlapping gene predictions. I removed them to reduce redundancy

```seqkit rmdup -s 06_cuprina_longest_isoforms_modified.faa > 06_cuprina_no_dup_longest_isoforms.fa```

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
#SBATCH --job-name=HSP90
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paige.matheson14@gmail.com
#SBATCH --output HSP90_%j.out
#SBATCH --error HSP90_%j.err

# =========================
# LOAD MODULES
# =========================
module purge
ml BLAST/2.16.0-GCC-12.3.0 SeqKit/2.4.0 HMMER/3.4-GCC-12.3.0 CD-HIT/4.8.1-GCC-11.3.0 BLASTDB GCC/12.3.0

# =========================
# GET REFERENCE IDS
# =========================
grep ">" hsp90.fasta | sed 's/>//' > ref_ids.txt

# =========================
# LOOP OVER SPECIES
# =========================
for i in 01_hilli 02_quadrimaculata 03_stygia 04_vicina 06_cuprina; do
    cd ${i}
    echo "Processing ${i} ..."

    # ------------------------
    # Step 1: Make BLAST DB
    # ------------------------
    makeblastdb -in ${i}_longest_isoforms_modified.faa -dbtype prot

    # ------------------------
    # Step 2: Forward BLAST
    # ------------------------
    blastp -query ../hsp90.fasta \
           -db ${i}_longest_isoforms_modified.faa \
           -out ${i}_HSP90_blastp_results.out \
           -evalue 1e-10 \
           -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"

    # ------------------------
    # Step 3: Extract BLAST hits
    # ------------------------
    bash ../extract_hsp90.sh \
         ${i}_HSP90_blastp_results.out \
         ${i}_longest_isoforms_modified.faa \
         ${i}_blast_hits.fasta

    # Safety check
    if [ ! -s ${i}_blast_hits.fasta ]; then
        echo "WARNING: No BLAST hits for ${i}, skipping..."
        cd ../
        continue
    fi

    # ------------------------
    # Step 4: Reciprocal BLAST
    # ------------------------
    blastp -query ${i}_blast_hits.fasta \
           -db ../hsp90.fasta \
           -out ${i}_reciprocal.out \
           -evalue 1e-10 \
           -outfmt "6 qseqid sseqid evalue bitscore"

    # Keep best hit per query
    sort -k1,1 -k4,4nr ${i}_reciprocal.out | awk '!seen[$1]++' > ${i}_reciprocal_best.tsv

    # FIXED RBH FILTER (this was your bug)
    grep -F -f ../ref_ids.txt ${i}_reciprocal_best.tsv | cut -f1 > ${i}_rbh_ids.txt

    # Safety check
    if [ ! -s ${i}_rbh_ids.txt ]; then
        echo "WARNING: No RBH hits found for ${i}, skipping..."
        cd ../
        continue
    fi

    # Extract RBH sequences
    seqkit grep -f ${i}_rbh_ids.txt ${i}_blast_hits.fasta > ${i}_rbh_hits.fasta

    # ------------------------
    # Step 5: HMMER validation
    # ------------------------
    hmmsearch --cpu 4 \
              --domtblout ${i}_hmmer_hsp90.domtblout \
              ../hsp90.hmm \
              ${i}_rbh_hits.fasta

    # ------------------------
    # Step 6: Filter significant hits
    # ------------------------
    grep -v '^#' ${i}_hmmer_hsp90.domtblout \
        | awk '$13 <= 1e-5 {print $1}' \
        | sort | uniq > ${i}_pfam_confirmed_ids.txt

    # Safety check
    if [ ! -s ${i}_pfam_confirmed_ids.txt ]; then
        echo "WARNING: No HMM-confirmed HSP90s for ${i}, skipping..."
        cd ../
        continue
    fi

    # ------------------------
    # Step 7: Extract confirmed sequences
    # ------------------------
    seqkit grep -f ${i}_pfam_confirmed_ids.txt ${i}_rbh_hits.fasta > ${i}_confirmed_HSP90.fasta

    # ------------------------
    # Step 8: Remove redundancy
    # ------------------------
    cd-hit -i ${i}_confirmed_HSP90.fasta \
           -o ${i}_HSP90_final_cdhit.faa \
           -c 0.99

    cd ../
done
```

## concatenate all of the species confirmed hsp90s into one fasta 
```cat *_confirmed_HSP90.fasta > all_species_hsp90.fasta```

## Add drosophila genes in 
Because HSP90 has three different homologs (paralogs), I added in drosophila genes from flybase for hsp83, trap1, and gp93 so when I build the phylogeny i can see what genes are associated with what homologs. The sequences that fall outside of these clusters probably represent isoforms (to remove) 



















