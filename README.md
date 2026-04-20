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







# figuring out the expanding/contracting gene families

# merge eggnog + interproscan annotations for each species 
```
llibrary(tidyverse)
library(dplyr)


samples <- c("01_hilli", "02_quadrimaculata", "03_stygia", "04_vicina", "06_cuprina")

for (sample in samples) {
  # File paths
  eggnog_gff <- file.path(sample, paste0(sample, "_eggnog.emapper.decorated_modified.gff"))
  interpro_tsv <- file.path(sample, paste0(sample, "_longest_isoforms_modified.faa.tsv"))
  
  cat("Processing:", sample, "\n")
  
  #---- 1. Parse EggNOG-mapper GFF3 ----#
  eggnog <- read_tsv(
    eggnog_gff,
    comment = "##",
    col_names = FALSE
  )
  
  # Filter mRNA rows where em_* annotations are present
  eggnog_mrna <- eggnog %>% filter(X3 == "mRNA")
  
  #---- 2. Extract key fields from attribute column using regex ----#
  eggnog_clean <- eggnog_mrna %>%
    transmute(Protein_ID = str_extract(X9, "(?<=ID=)[^;]+"),
              GO_terms   = str_extract(X9, "(?<=em_GOs=)[^;]+"),
              KEGG_paths = str_extract(X9, "(?<=em_KEGG_Pathway=)[^;]+"),
              KEGG_KO    = str_extract(X9, "(?<=em_KEGG_ko=)[^;]+"),
              BRITE_class= str_extract(X9, "(?<=em_BRITE=)[^;]+"),
              EggNOG_PFAMs = str_extract(X9, "(?<=em_PFAMs=)[^;]+"),
              EggNOG_PFAMs_desc = str_extract(X9, "(?<=em_desc=)[^;]+"))
  
  # Remove rows with all NA or Protein_ID missing
  eggnog_clean <- eggnog_clean %>%
    filter(!is.na(Protein_ID)) %>%
    group_by(Protein_ID) %>%
    summarise(across(everything(), ~ toString(na.omit(.))), .groups = "drop") %>%
    mutate(
      ID_number = as.numeric(str_extract(Protein_ID, "(?<=g)\\d+"))
    ) %>%
    arrange(ID_number) %>%
    select(-ID_number)
  
  #---- 2. Parse InterProScan TSV ----#
  interpro <- read_tsv(interpro_tsv, col_names = FALSE, comment = "#")
  
  # 1. Pfam domains and descriptions (only from source = "Pfam")
  pfam_data <- interpro %>%
    filter(X4 == "Pfam") %>%
    select(Protein_ID = X1, Pfam_ID = X5, Pfam_desc = X6, Pfam_Eval = X9) %>%
    mutate(Pfam_Eval = as.numeric(Pfam_Eval)) %>%
    group_by(Protein_ID, Pfam_ID, Pfam_desc) %>%
    summarise(Min_Eval = min(Pfam_Eval, na.rm = TRUE), .groups = "drop") %>%
    group_by(Protein_ID) %>%
    summarise(
      Pfam_domains = toString(Pfam_ID),
      Pfam_descriptions = toString(Pfam_desc),
      Pfam_Evals = toString(Min_Eval),
      .groups = "drop"
    )
  
  #### Panther fetch
  panther_data <- interpro %>%
    filter(X4 == "PANTHER") %>%
    select(Protein_ID = X1, Panther_ID = X5, Panther_desc = X6, Panther_Eval = X9) %>%
    mutate(Panther_Eval = as.numeric(Panther_Eval)) %>%
    group_by(Protein_ID, Panther_ID, Panther_desc) %>%
    summarise(Min_Eval = min(Panther_Eval, na.rm = TRUE), .groups = "drop") %>%
    group_by(Protein_ID) %>%
    summarise(
      Panther_ID = toString(Panther_ID),
      Panther_descriptions = toString(Panther_desc),
      Panther_Evals = toString(Min_Eval),
      .groups = "drop"
    )
  
  # 2. InterPro IDs and descriptions (from all sources)
  ipr_data <- interpro %>%
    select(Protein_ID = X1, InterPro_ID = X12, InterPro_desc = X13) %>%
    filter(!is.na(InterPro_ID)) %>%
    group_by(Protein_ID) %>%
    summarise(
      InterPro_IDs = toString(unique(na.omit(InterPro_ID))),
      InterPro_descriptions = toString(unique(na.omit(InterPro_desc))),
      .groups = "drop"
    )
  
  # 3. GO terms (column X14), clean up GO IDs
  go_data <- interpro %>%
    select(Protein_ID = X1, GO = X14) %>%
    filter(!is.na(GO)) %>%
    separate_rows(GO, sep = ",\\s*") %>%                              # Split multiple GO term blocks
    mutate(GO = str_extract_all(GO, "GO:\\d{7}")) %>%                 # Extract only GO:#######
  unnest(GO) %>%
    group_by(Protein_ID) %>%
    summarise(
      GO_terms_interpro = toString(unique(GO)),
      .groups = "drop"
    )
  
  # 4. Merge InterProScan pieces
  interpro_clean <- reduce(
    list(pfam_data, panther_data, ipr_data, go_data),
    full_join,
    by = "Protein_ID"
  ) %>%
    mutate(
      ID_number = as.numeric(str_extract(Protein_ID, "(?<=g)\\d+"))
    ) %>%
    arrange(ID_number) %>%
    select(-ID_number)
  
  # 5. Merge with EggNOG annotations (GO terms remain separate)
  merged_annotations <- full_join(eggnog_clean, interpro_clean, by = "Protein_ID") %>%
    select(
      Protein_ID,
      EggNOG_PFAMs,
      EggNOG_PFAMs_desc,
      KEGG_paths,
      KEGG_KO,
      BRITE_class,
      Pfam_domains,
      Pfam_descriptions,
      Pfam_Evals,
      Panther_ID,
      Panther_descriptions,
      Panther_Evals,
      InterPro_IDs,
      InterPro_descriptions,
      GO_terms,               # from EggNOG
      GO_terms_interpro      # cleaned InterProScan GO terms
    )
  
  # 6. Clean NA and export
  merged_annotations <- merged_annotations %>%
    mutate(across(everything(), ~str_replace_all(., "NA|^, | ,|, NA", "")))
  
  # Reorder columns in the desired order
  merged_annotations <- merged_annotations %>%
    select(
      Protein_ID,
      EggNOG_PFAMs,
      EggNOG_PFAMs_desc,
      GO_terms,             
      KEGG_KO,
      KEGG_paths,
      BRITE_class,
      Pfam_domains,
      Pfam_descriptions,
      Pfam_Evals,
      InterPro_IDs,
      InterPro_descriptions,
      GO_terms_interpro,
      Panther_ID,
      Panther_descriptions,
      Panther_Evals,
    )
  colnames(merged_annotations) <- c("Protein_ID", "EggNOG_PFAMs_domains","EggNOG_PFAMs_descriptions", "EggNOG_GO_terms", "EggNOG_KEGG_KO", "EggNOG_KEGG_paths", "EggNOG_BRITE_class", "Interpro_Pfam_Ids",
                                    "Interpro_Pfam_descriptions", "Interpro_Pfam_Evals", "InterPro_IDs", "InterPro_descriptions", "Interpro_GO_terms", "Panther_ID", "Panther_descriptions", "Panther_Evals")
  
  # ---- 4. Write to file ----
  output_file <- file.path(sample, "merged_protein_annotations.tsv")
  write_tsv(merged_annotations, output_file)
}
```

# annotate the merged annotation file with the orthogroups from orthofinder
```
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# === File paths ===
orthogroups_file <- "Orthogroups.tsv"
og_list_dir <- "Orthogroups.tsv"
annotation_dir <- "01_annotated_genes"

# === Load OrthoFinder orthogroups file ===
orthogroups <- read_tsv(orthogroups_file)

# Create a long-format OG → gene mapping
og_gene_map <- orthogroups %>%
  pivot_longer(-Orthogroup, names_to = "Species", values_to = "Genes") %>%
  mutate(Genes = strsplit(Genes, ", ")) %>%
  unnest(Genes) %>%
  distinct(Orthogroup, Species, Genes)


# === Process each species ===
ann_files <- list.files(annotation_dir, pattern = "_annotated.tsv", full.names = TRUE)

for (ann_file in ann_files) {
  species <- str_replace(basename(ann_file), "_annotated.tsv", "")
  message("Processing species: ", species)
  
  # Load OG list
  og_ids <- unique(og_gene_map$Orthogroup)
  og_subset <- og_gene_map %>%
    filter(`Orthogroup` %in% og_ids & str_detect(Genes, species))  # Filter OG + species
  
  
  # Load merged annotation
  annotation_file <- file.path(annotation_dir, paste0(species, "_annotated.tsv"))
  ann <- read_tsv(annotation_file, show_col_types = FALSE)
  
  # Clean GeneID column name
  colnames(ann)[1] <- "GeneID"
  
  
  # Merge with OG data
  annotated <- og_subset %>%
    rename(GeneID = Genes) %>%
    left_join(ann, by = "GeneID") %>%
    select(
      Orthogroup, GeneID,
      EggNOG_PFAMs_domains,
      EggNOG_PFAMs_descriptions,
      EggNOG_GO_terms,
      EggNOG_KEGG_KO,
      EggNOG_KEGG_paths,
      EggNOG_BRITE_class,
      Interpro_Pfam_Ids,
      Interpro_Pfam_descriptions,
      Interpro_Pfam_Evals,
      InterPro_IDs,
      InterPro_descriptions,
      Interpro_GO_terms,
      Panther_ID, Panther_descriptions, Panther_Evals
    )
  
  # Output file
  out_file <- file.path("02_annotated_OGs", paste0(species, "_annotated_OGs.tsv"))
  write_tsv(annotated, out_file)
  message("[✓] Output written: ", out_file)
}
```

# figure out what orthogroups are significant based on the cafe files
```
#!/usr/bin/env python3

from collections import defaultdict

input_file = "Gamma_branch_probabilities.tab"
threshold = 0.05

# Dictionary to hold column name -> list of matching FamilyIDs
results = defaultdict(list)

with open(input_file) as f:
    header = f.readline().strip().split("\t")

    for line in f:
        parts = line.strip().split("\t")
        family_id = parts[0]
        for i in range(1, len(parts)):
            try:
                val = float(parts[i])
                if val <= threshold:
                    results[header[i]].append(family_id)
            except ValueError:
                continue  # skip N/A or non-numeric values

# Output CSV-style results
for col in header[1:]:
    if results[col]:  # only print columns with at least one match
        print(f"{col}," + ",".join(results[col]))

# Output one file per column (species)
for col in header[1:]:
    if results[col]:  # only if there are matches
        filename = f"{col}.txt"
        with open(filename, "w") as out:
            out.write("\n".join(results[col]) + "\n")
```

this writes txt files that are used below:

# split into expanding contracting
```
# === Libraries ===
library(readr)
library(dplyr)
library(tidyr)
library(stringr)

# === Input files ===
gamma_file <- "Gamma_change.tab"
og_dir <- "01_species_sig_ogs"  # folder containing *_sig_ogs.txt files

# === Read CAFE Gamma_change.tab ===
gamma <- read_tsv(gamma_file)

# === Clean column names ===
colnames(gamma)[1] <- "Orthogroup"

# Clean species names
colnames(gamma) <- str_remove(colnames(gamma), "<.*?>")
colnames(gamma) <- str_replace(colnames(gamma),
                               "06_cuprina_longest_isoforms_modified", "06_cuprina")
colnames(gamma) <- str_replace(colnames(gamma),
                               "01_hilli_longest_isoforms_modified", "01_hilli")
colnames(gamma) <- str_replace(colnames(gamma),
                               "03_stygia_longest_isoforms_modified", "03_stygia")
colnames(gamma) <- str_replace(colnames(gamma),
                               "02_quadrimaculata_longest_isoforms_modified", "02_quadrimaculata")
colnames(gamma) <- str_replace(colnames(gamma),
                               "04_vicina_longest_isoforms_modified", "04_vicina")

# Remove the model <?> column names
gamma <- gamma[, !colnames(gamma) %in% c("", ".1", ".2", ".3")]


# === Convert to long format ===
gamma_long <- gamma %>%
  pivot_longer(-Orthogroup, names_to = "Species", values_to = "Change") %>%
  mutate(Change = as.numeric(Change)) %>%
  filter(Change != 0) %>%
  mutate(Direction = ifelse(Change > 0, "expanded", "contracted"))

# === List species-specific OG files ===
og_files <- list.files(og_dir, pattern = "_sig_ogs.txt", full.names = TRUE)

# === Process each species ===
for (og_file in og_files) {
  
  # Extract species name from file
  species <- str_replace(basename(og_file), "_sig_ogs.txt", "")
  message("Processing species: ", species)
  
  # Read significant OGs
  sig_ogs <- read_lines(og_file)
  
  # Match with gamma data
  changes <- gamma_long %>%
    filter(Species == species, Orthogroup %in% sig_ogs)
  
  # Split into expanded and contracted
  expanded <- changes %>% filter(Direction == "expanded") %>% pull(Orthogroup)
  contracted <- changes %>% filter(Direction == "contracted") %>% pull(Orthogroup)
  
  # Write to files
  write_lines(expanded, file.path(og_dir, paste0(species, "_expanded_OGs.txt")))
  write_lines(contracted, file.path(og_dir, paste0(species, "_contracted_OGs.txt")))
  
  message("[✓] ", species, ": ", length(expanded), " expanded, ", length(contracted), " contracted OGs")
}
```








