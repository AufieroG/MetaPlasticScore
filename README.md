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
devtools::install_github("USERNAME/PlasticScore")
```

---

## ⚙️ Input Preparation

The pipeline requires a strict directory structure and file naming convention to link HMMER outputs with taxonomic and abundance data. All identifiers (taxon IDs, sample names) must match exactly across files (case-sensitive).

### 1. HMMER search outputs

PlasticScore requires output files from `hmmsearch` (HMMER suite) generated using enzyme-specific HMM profiles.

- **Directory:** Place all HMMER output files in a single directory (e.g. `hmmsearch/`).  
- **File formats:** Both `.tbl` (table) and `.domtbl` (domain table) files are required.  
- **Naming convention:** Files must be named as:

```
<taxon_id>__<enzyme_id>.tbl
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

These HMM profiles were built from protein sequences identified using the **PlasticDB** database and are based on the study:

> **Zrimec et al., 2023**  
> *Plastic-degrading potential across the global microbiome*  
> Microbiome  
> https://doi.org/10.1186/s40168-023-01649-0

Users may use these models directly or replace them with custom HMM profiles if desired.

### Example `hmmsearch` command

```bash
hmmsearch   --cpu 1   --domtblout output_dir/<taxon_id>__<enzyme_id>.domtbl   --tblout    output_dir/<taxon_id>__<enzyme_id>.tbl   -E 1e-5   --incE 1e-3   path/to/HMMmodel/<enzyme_id>.hmm   path/to/protein_sequences/<taxon_id>.faa
```

---

### 2. Abundance and metadata files

The following CSV files are required and must be mutually consistent.

#### Abundance table

- File example: `Abbundance_simulated_low_variance.csv`
- **First column:** `Taxa` (taxon identifier)  
- **Remaining columns:** sample IDs (`STD_rep1`, `STD_rep2`, `PE_rep1`, …)  
- **Values:** normalized or absolute abundances (numeric)

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
  - `taxon_id` → taxon identifier (must match `Taxa`)
  - `Phylum` → phylum assignment

Example:

```
taxon_id,Phylum
Bifidobacterium_pseudolongum,Actinobacteria
Cellulosimicrobium_sp,Actinobacteria
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

This experimental design follows the conditions described in:

> **Zrimec et al., 2023**  
> *Plastic-degrading potential across the global microbiome*  
> Microbiome  
> https://doi.org/10.1186/s40168-023-01649-0

⚠️ **Important:** The example datasets are **not intended to represent full experimental results** and should not be used for biological inference. They are included solely to illustrate the workflow and usage of the PlasticScore pipeline.

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
hits_filt          <- df$hits_filt
plastic_score      <- df$plastic_score
enzyme_load        <- df$enzyme_load
taxon_enzyme_mat   <- df$taxon_enzyme_mat
taxon_stats        <- df$taxon_stats
abundance          <- df$abundance
metadata_groups    <- df$metadata_groups
metadata_taxa      <- df$metadata_taxa
enzyme_by_sample   <- df$enzyme_by_sample
contrib_phylum_mat <- df$contrib_phylum_mat
phylum_enzyme_mat  <- df$phylum_enzyme_mat
df_phylum_enzyme   <- df$df_phylum_enzyme
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

## 📖 Citation

If you use PlasticScore, please cite:

Zrimec et al., 2023. *Plastic-degrading potential across the global microbiome*.  
Microbiome. https://doi.org/10.1186/s40168-023-01649-0

---

## 📜 License

Specify license here (e.g. MIT, GPL-3.0).




