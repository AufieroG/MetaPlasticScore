# PlasticScore

**PlasticScore** is an R-based computational pipeline designed to assess and quantify the plastic-degrading potential of microbial taxa. By integrating HMMER `hmmsearch` outputs with normalized taxon abundance data and metadata, PlasticScore calculates degradation scores, analyzes enzyme loads, and visualizes contributions across taxonomic levels.

⚠️ **Development status:** PlasticScore is currently **under active development**. APIs, default parameters, and outputs may change in future releases.

---

## 📦 Prerequisites & Dependencies

PlasticScore runs in **R**. Before using the pipeline, ensure you have the required packages installed.

### Install CRAN and Bioconductor packages

```r
# Install CRAN packages
install.packages(c("stats", "utils", "ggplot2", "vegan", "tidyr",
                   "circlize", "grid", "ggrepel", "viridis",
                   "RColorBrewer", "magrittr", "dplyr", "devtools"))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
```

---

## 📦 Package installation

```r
# Replace USERNAME with the repository owner
devtools::install_github("AufieroG/MetaPlasticScore")
```

---

## ⚙️ Input Preparation

The pipeline requires a strict directory structure and file naming convention to link HMMER outputs with taxonomic and abundance data. All identifiers (taxon IDs, sample names) must match exactly across files (case-sensitive).

### 1. HMMER search outputs

PlasticScore requires output files from `hmmsearch` (HMMER suite) generated using enzyme-specific HMM profiles.

- **Directory:** Place all HMMER output files in a single directory (e.g. `hmmsearch/`).  
- **File formats:** `.tbl` (table) or `.domtbl` (domain table) files are required.  
- **Naming convention:** Files must be named as:

```
<taxon_id>__<enzyme_id>.tbl
or
<taxon_id>__<enzyme_id>.domtbl
```

**Important:**  
- `<taxon_id>` must exactly match the taxon identifiers used in the abundance and taxonomy CSV files.  
- `<enzyme_id>` must correspond to the HMM model name.

### Pre-generated HMM models

The package includes **pre-generated HMM models** for plastic-degrading enzymes, located in:

```
PlasticScore/inst/extdata/HMMmodel/
```

These HMM profiles were built from protein sequences identified using:

> *PlasticDB*
> [https://plasticdb.org/](https://plasticdb.org/)
> *Victor Gambarini, Olga Pantos, Joanne M Kingsbury, Louise Weaver, Kim M Handley, Gavin Lear, PlasticDB: a database of microorganisms and proteins linked to plastic biodegradation. Database*
> [https://doi.org/10.1093/database/baac008](https://doi.org/10.1093/database/baac008)

and are based on the experimental framework and annotations described in:

> *De Filippis, F., Bonelli, M., Bruno, D. et al. Plastics shape the black soldier fly larvae gut microbiome and select for biodegrading functions. Microbiome 11, 205 (2023).*
> [https://doi.org/10.1186/s40168-023-01649-0](https://doi.org/10.1186/s40168-023-01649-0)

Users may use these HMM models directly or replace them with custom HMM profiles if desired.

The HMM models are used as input to `hmmsearch` to perform profile-to-sequence searches against the protein complement of each taxon.

### Example `hmmsearch` command

```bash
for hmm in path/to/HMMmodel/*.hmm; do
  enzyme_id=$(basename "$hmm" .hmm)
  for taxon in $(cat taxa_list.txt); do
    hmmsearch \
      --cpu 1 \
      --domtblout output_dir/${taxon}__${enzyme_id}.domtbl \
      --tblout    output_dir/${taxon}__${enzyme_id}.tbl \
      -E 1e-5 \
      --incE 1e-3 \
      "$hmm" \
      path/to/protein_sequences/${taxon}.faa
  done
done
```

---

### 2. Abundance and metadata files

The following CSV files are also required and must be mutually consistent.

#### Abundance table

- File example: `Abundance.csv`
- - Required columns:
  - `Taxa` (taxon identifier)  
  - `sample IDs` (`STD_rep1`, `STD_rep2`, `PE_rep1`, …)  
- Values: normalized or absolute abundances (numeric)

Example structure:

```
Taxa,STD_rep1,STD_rep2,STD_rep3,PE_rep1,PE_rep2,PE_rep3
Bifidobacterium_pseudolongum,0.766,0.762,0.782,0,0,0
Cellulosimicrobium_sp,0.183,0.155,0.162,0.131,0.118,0.119
...
```

#### Sample metadata

- File example: `Metadata_groups.csv`
- Required columns:
  - `Samples` → sample ID (must match abundance column names)
  - `Groups` → experimental condition

Example:

```
Samples,Groups
STD_rep1,STD
STD_rep2,STD
STD_rep3,STD
PE_rep1,PE
PE_rep2,PE
PE_rep3,PE
```

#### Taxonomic metadata

- File example: `Metadata_taxa.csv`
- Required columns:
  - `Phylum` → phylum assignment
  - `taxon_id` → taxon identifier (must match `Taxa`)


Example:

```
Phylum,taxon_id
Actinobacteria,Bifidobacterium_pseudolongum
Actinobacteria,Cellulosimicrobium_sp
...
```

---

### 📂 Example datasets included in the package

The package includes **example input datasets** located in:

```
PlasticScore/inst/extdata/
```

These example files are provided **for demonstration and testing purposes only**.

In the example datasets:
- Samples labeled **PE** correspond to **larvae fed on polyethylene (PE)**.
- Samples labeled **STD** correspond to **larvae fed on a standard diet**.

---

## ▶️ Running the PlasticScore pipeline

```r
df <- run_PlasticScore(
  hmmsearch_files = "path/hmmsearch/",
  abundance_csv = "path/Abbundance_simulated_low_variance.csv",
  metadata_groups_csv = "path/Metadata_groups.csv",
  metadata_taxa_csv = "path/Metadata_taxa.csv",
  min_hmm_coverage = 0.7,
  max_dom_bias = 0.1,
  enzymes_for_plastic = c("All"),
  method = "abundance",
  normalize_abundance = FALSE
)
```

---

## 🧾 Main Output Objects

```r
# Filtered HMMER hits retained after coverage and domain-bias filtering
hits_filt <- df$hits_filt

# PlasticScore computed for each sample
plastic_score <- df$plastic_score

# Enzyme load per taxon (number or abundance-weighted count of plastic-degrading enzymes)
enzyme_load <- df$enzyme_load

# Taxon × enzyme matrix (presence / counts / abundance-weighted values)
taxon_enzyme_mat <- df$taxon_enzyme_mat

# Per-taxon summary statistics and contribution metrics
taxon_stats <- df$taxon_stats

# Abundance table actually used in the analysis (raw or normalized)
abundance <- df$abundance

# Sample metadata (e.g. experimental groups such as STD and PE)
metadata_groups <- df$metadata_groups

# Taxonomic metadata used for aggregation (e.g. Phylum)
metadata_taxa <- df$metadata_taxa

# Enzyme abundance aggregated at the sample level
enzyme_by_sample <- df$enzyme_by_sample

# Contribution of each phylum to the PlasticScore or enzyme load
contrib_phylum_mat <- df$contrib_phylum_mat

# Phylum × enzyme matrix (enzyme profiles aggregated by phylum)
phylum_enzyme_mat <- df$phylum_enzyme_mat

# Long-format table of enzyme abundances grouped by phylum (tidy format)
df_phylum_enzyme <- df$df_phylum_enzyme
```

---

## 📈 Plot generation

```r
plot_PlasticScore(
  df_phylum_enzyme = df_phylum_enzyme,
  taxon_enzyme_mat = taxon_enzyme_mat,
  phylum_enzyme_mat = phylum_enzyme_mat,
  enzyme_by_sample = enzyme_by_sample,
  plastic_score = plastic_score,
  contrib_phylum_mat = contrib_phylum_mat,
  meta = metadata_groups,
  meta_taxa = metadata_taxa,
  sample_col = "Samples",
  condition_col = "Groups"
)
```

---

## 📝 Notes

- Taxon identifiers must be consistent across all inputs.  
- Default parameters are tuned for normalized abundance data.  
- HMM models are provided for convenience and reproducibility.  
- The package is under active development; feedback is welcome.

---





