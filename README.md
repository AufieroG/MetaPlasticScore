# MetaPlasticScore

<p align="center">
<img width="400" height="400" alt="MetaPlasticScore" src="https://github.com/user-attachments/assets/355ba6ac-e519-4b00-bbd7-43c846d5f111" />
<p align="center">

**PlasticScore** is an R-based computational pipeline designed to assess and quantify the plastic-degrading potential of microbial taxa. By integrating HMMER `hmmsearch` outputs with normalized taxon abundance data and metadata, PlasticScore calculates degradation scores, analyzes enzyme loads, and visualizes contributions across taxonomic levels.

**Development status:** MetaPlasticScore R package is currently under development. Default parameters, and outputs may change in future releases.

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
# install and load
devtools::install_github("AufieroG/MetaPlasticScore")
library("MetaPlasticScore")
```

---

## ⚙️ Input Preparation

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

#### Pre-generated HMM models

The package includes **pre-generated HMM models** for plastic-degrading enzymes, located in:

```
PlasticScore/inst/extdata/HMMmodel/
```

These HMM profiles were built from protein sequences identified using:

> *PlasticDB*
> [https://plasticdb.org/](https://plasticdb.org/)
> *Victor Gambarini, Olga Pantos, Joanne M Kingsbury, Louise Weaver, Kim M Handley, Gavin Lear, PlasticDB: a database of microorganisms and proteins linked to plastic biodegradation. Database*
> [https://doi.org/10.1093/database/baac008](https://doi.org/10.1093/database/baac008)

and annotations described in:

> *De Filippis, F., Bonelli, M., Bruno, D. et al. Plastics shape the black soldier fly larvae gut microbiome and select for biodegrading functions. Microbiome 11, 205 (2023).*
> [https://doi.org/10.1186/s40168-023-01649-0](https://doi.org/10.1186/s40168-023-01649-0)

Users may use these HMM models directly or replace them with custom HMM profiles if desired.

The HMM models are used as input to `hmmsearch` to perform profile-to-sequence searches against the protein complement of each taxon.

#### Example `hmmsearch` command

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
### Example command to run the pipeline
```r
# locate example data loaded with the package
extdata_path <- system.file("extdata", package = "MetaPlasticScore")

hmm_dir <- file.path(extdata_path, "hmmsearch")
abundance_csv <- file.path(extdata_path, "Abundance_simulated.csv")
metadata_groups_csv <- file.path(extdata_path, "Metadata_groups.csv")
metadata_taxa_csv <- file.path(extdata_path, "Metadata_taxa.csv")


# Wrapper that executes the pipeline
result <- run_MetaPlasticScore(
  hmmsearch_files = hmm_dir,
  abundance_csv = abundance_csv,
  metadata_groups_csv = metadata_groups_csv,
  metadata_taxa_csv = metadata_taxa_csv,
  min_hmm_coverage = 0.7,
  max_dom_bias = 0.1,
  enzymes_for_plastic = "All",
  method = "abundance",
  normalize_abundance = FALSE
)
```

---

### 🧾 Main Output Objects

```r
# Filtered HMMER hits retained after coverage and domain-bias filtering
hits_filt <- result$hits_filt

# PlasticScore computed for each sample
plastic_score <- result$plastic_score

# Enzyme load per taxon (number or abundance-weighted count of plastic-degrading enzymes)
enzyme_load <- result$enzyme_load

# Taxon × enzyme matrix (presence / counts / abundance-weighted values)
taxon_enzyme_mat <- result$taxon_enzyme_mat

# Per-taxon summary statistics and contribution metrics
taxon_stats <- result$taxon_stats

# Abundance table actually used in the analysis (raw or normalized)
abundance <- result$abundance

# Sample metadata (e.g. experimental groups such as STD and PE)
metadata_groups <- result$metadata_groups

# Taxonomic metadata used for aggregation (e.g. Phylum)
metadata_taxa <- result$metadata_taxa

# Enzyme abundance aggregated at the sample level
enzyme_by_sample <- result$enzyme_by_sample

# Contribution of each phylum to the PlasticScore or enzyme load
contrib_phylum_mat <- result$contrib_phylum_mat

# Phylum × enzyme matrix (enzyme profiles aggregated by phylum)
phylum_enzyme_mat <- result$phylum_enzyme_mat

# Long-format table of enzyme abundances grouped by phylum (tidy format)
df_phylum_enzyme <- result$df_phylum_enzyme
```

---

### 📈 Plot generation

```r
# The default plotting configuration provided by the default_plot_config() function is:
 
   # # ============================================================
    # # OUTPUT / SAVING OPTIONS (GLOBAL)
    # 
    # output_dir = "Plots",      # Directory where plots are saved if save = TRUE
    # save = TRUE,               # If TRUE, plots are written to disk automatically
    # dpi = 300,                 # Resolution for raster images (publication-ready)
    # width = 8,                 # Default plot width in inches
    # height = 5,                # Default plot height in inches
    # 
    # # ============================================================
    # # GENERAL SIZES
    # 
    # base_font_size = 10,       # Base font size used by theme_minimal()
    # axis_text_size = 9,       # Font size for axis tick labels
    # title_size = 14,           # Font size for plot titles
    # 
    # # ============================================================
    # # COLOR PALETTES
    # 
    # palette_enzyme = "cividis", # Used for enzyme-level stacked bars
    # palette_group  = "cividis", # Used for group-based plots (PCoA, boxplots)
    # palette_phylum = "Set2",    # Used for phylum-level stacked contributions
    # 
    # # ============================================================
    # # GENERAL AESTHETICS
    # 
    # alpha_points = 0.8,        # Transparency for scatter points
    # jitter_width = 0.15,       # Horizontal jitter for overplotted points
    # 
    # # ============================================================
    # # HEATMAP OPTIONS (plot_enzyme_heatmap)
    # 
    # heatmap_row_level = "Taxon",   # "Taxon" or "Phylum"
    # heatmap_log = TRUE,            # Apply log10(x + 1) for visualization only
    # heatmap_cluster_rows = TRUE,   # Hierarchical clustering of rows
    # heatmap_cluster_cols = TRUE,   # Hierarchical clustering of columns
    # heatmap_top_n_taxa = NULL,     # Show only top-N rows (by total abundance)
    # heatmap_title = "Enzyme abundance heatmap",
    # heatmap_row_label_size = 8,    # Font size of row labels
    # heatmap_col_label_size = 8,    # Font size of column labels
    # 
    # # ============================================================
    # # PCoA OPTIONS (plot_pcoa_enzyme_profiles)
    # 
    # pcoa_label = TRUE,             # Draw sample labels using ggrepel
    # pcoa_point_size = 3,           # Size of points in ordination plot
    # pcoa_transform = "none",       # "none" | "log"
    # pcoa_distance  = "bray",       # "bray" | "euclidean"
    # 
    # # ============================================================
    # # STACKED CONTRIBUTION PLOTS
    # 
    # normalize_contrib = TRUE      # Show contributions as percentages per sample

# To customize the plotting parameters, modify the configuration object returned by default_plot_config() and pass it to the cfg argument of plot_MetaPlasticScore().
# config <- default_plot_config()
# config$output_dir = "Plots_result"      # Directory where plots are saved if save = TRUE
# config$save = FALSE                     # If TRUE, plots are written to disk automatically
# config$dpi = 200                        # Resolution for raster images (publication-ready)
# config$width = 15                       # Default plot width in inches
# config$height = 10                      # Default plot height in inches

# Wrapper to generate and save all plots
plot_MetaPlasticScore(
  df_phylum_enzyme = df_phylum_enzyme,
  taxon_enzyme_mat = taxon_enzyme_mat,
  phylum_enzyme_mat = phylum_enzyme_mat,
  enzyme_by_sample = enzyme_by_sample,
  plastic_score = plastic_score,
  contrib_phylum_mat = contrib_phylum_mat,
  meta = metadata_groups,
  meta_taxa = metadata_taxa,
  sample_col = "Samples",
  condition_col = "Groups",
  cfg = NULL                                   # Set to NULL to use the default configuration, or provide a custom configuration object
)
```

<img width="600" height="375" alt="stacked_bar" src="https://github.com/user-attachments/assets/50255315-dbc9-4121-845c-d1e06cc6f98e" />


<img width="600" height="375" alt="plastic_score" src="https://github.com/user-attachments/assets/64ac8df9-dfa0-41bd-8718-645d096fc4b8" />


<img width="600" height="375" alt="pcoa" src="https://github.com/user-attachments/assets/65216ea1-56df-47b8-ba54-c38274aa1c20" />


<img width="600" height="375" alt="heatmap" src="https://github.com/user-attachments/assets/51806b64-dd1e-479a-b103-3ac1b54a2501" />


<img width="600" height="375" alt="contributions" src="https://github.com/user-attachments/assets/4e480b16-9931-4926-b4e5-b074ec663baa" />


<img width="600" height="375" alt="bubble" src="https://github.com/user-attachments/assets/f8d89cfb-bd9f-46a8-bc5e-be87f572ef8a" />


<img width="600" height="375" alt="abundance_violin" src="https://github.com/user-attachments/assets/868396a4-dcb7-445f-bf5c-788d0ccf6cd8" />




