# MONTI: A multi-omics non-negative tensor decomposition framework for gene-level integrative analysis
https://www.frontiersin.org/articles/10.3389/fgene.2021.682841/abstract

Multi-omics data is frequently measured to characterize biological mechanisms underlying phenotypes. Complex relationships in multi-omics data, if mined, can lead to more accurate classification of patient samples according to the phenotypes.

MONTI (Multi-Omics Non-negative Tensor decomposition for Integrative analysis) is a tool that can be used to integrate and analyze large sets of multi-omics data. MONTI identifies gene regulatory multi-omics features specific to a group of samples that share a common biological trait. This repository extends the original MONTI code base with a reproducible pipeline that downloads The Cancer Genome Atlas (TCGA) data, converts every omics layer into the MONTI gene-centric representation, and exports MOPA-ready tensors.

### Pipeline overview

```
TCGA data → MONTI preprocessing → gene-level matrices → tensor.npy → MOPA
```

The repository is organised as follows:

```
src/                # Original MONTI python package
pipelines/          # Orchestrated data engineering scripts
data/               # Raw downloads (git-ignored)
output/             # Gene-level matrices and tensors (git-ignored)
config/config.yml   # Pipeline configuration
```

### Environment setup

1. **R packages** – required for downloading TCGA data and for quantile normalisation.

   ```bash
   Rscript install_dependencies.R
   ```

   The script installs `TCGAbiolinks`, `minfi`, `limma`, `edgeR`, `preprocessCore`, `reticulate`, and `yaml` through Bioconductor/CRAN.

2. **Python packages** – required for tensor construction.

   ```bash
   python3 -m pip install numpy pandas pyyaml
   ```

   (Optional) configure `reticulate` to use the same Python environment by setting `RETICULATE_PYTHON` before running the R scripts.

3. **Metadata resources** – the gene-centric conversion relies on ancillary mapping tables distributed with MONTI:

   - `metadata/ensembl_gene_transcript_map.txt`
   - `metadata/promoter_probes_illumina450.txt`
   - `metadata/mirna_target_mirdb.txt`

   Download these files from the original [MONTI tutorial package](http://cobi.knu.ac.kr/tools.php) and place them inside a local `metadata/` directory (or update the paths in `config/config.yml`).

### Running the pipeline

1. **Download TCGA BRCA data** (customise the project and query parameters in `config/config.yml`).

   ```bash
   Rscript pipelines/download_TCGA.R config/config.yml
   ```

   The script uses `TCGAbiolinks::GDCquery`, `GDCdownload`, and `GDCprepare` to fetch RNA-seq (HTSeq-FPKM), DNA methylation (Illumina 450K β-values), and miRNA-seq counts. Raw matrices are written to `data/` as CSV files alongside the raw `SummarizedExperiment` objects (`.rds`).

2. **Generate gene-level matrices** via the MONTI utilities.

   ```bash
   Rscript pipelines/generate_gene_level.R config/config.yml
   ```

   This step calls `make_methylation_gcentric` and `make_mir_gcentric` from the Python `src/monti.py` module (through `reticulate`). Expression, methylation, and miRNA matrices are log2-transformed, quantile-normalised, and min–max scaled to [0, 1]. Harmonised matrices covering the intersecting genes and samples are stored under `output/`:

   - `gene_expression_matrix.tsv`
   - `gene_methylation_matrix.tsv`
   - `gene_miRNA_matrix.tsv`

3. **Build the MOPA tensor**.

   ```bash
   python pipelines/build_tensor.py --config config/config.yml
   ```

   The script stacks the three matrices into `output/tensor.npy`, writes `output/genelist.txt`, and records the sample IDs in `output/samplelist.txt`.

### Outputs

The processed artefacts in `output/` are ready for downstream tensor decomposition or pathway analysis with MOPA. Every run is controlled through `config/config.yml`, which exposes the TCGA project identifier, download parameters, raw input paths, and output filenames.

Below is an illustration of the analysis workflow of MONTI.
![workflow](./images/monti_workflow.jpg)

The output of MONTI is a simple gene list with information of their associated subtypes, which can be used for further downstream analysis. For example, the Venn diagram below shows the genes that are found to be associated to colorectal cancer subtypes CMS1, CMS2, CMS3 and CMS4. These genes showed to be informative in separating the four subtypes as shown in the t-SNE plot.
<!--![example output](./images/monti_outputexample.png =250x)-->
<img src="./images/monti_outputexample.png" alt="example output" width="600"/>

---

<!-- ## Download MONTI
```bash
git clone https://github.com/inukj/MONTI.git
``` -->

## Install MONTI
MONTI is developed in python3 and can be installed as below
```bash
python3 -m pip install monti
```
## Documentation
The functions and objects used by MONTI are documented [here](https://github.com/inukj/MONTI/blob/main/documentation/documentation.md).

## Tutorial using colon cancer data (TCGA-COAD)
A brief tutorial for using MONTI can be found under the ['tutorial'](https://github.com/inukj/MONTI/tree/main/tutorial) directory.
The associated multi-omics data are included.

If the above link does not work, the tutorial data is also available [here](http://cobi.knu.ac.kr/tools.php).

Before starting the tutorial, the dataset should be downloaded.
After download decompress data by
```bash
cd <download_path>
tar -xzvf tutorial_data_coad.tar.gz
```

The *<download_path>* should also be used as the tutorial directory, or you can simply move the data to another directory to be used for the tutorial.

The data includes three omics data, 1) gene expression (mRNA), 2) methylation level and 3) miRNA expression.
They are raw data directly collected from the TCGA portal.

In the [jupyter notebook](https://github.com/inukj/MONTI/blob/main/tutorial/tutorial_coad.ipynb) shows an example of how to integrate multi-omics data in a gene-level manner and extract features that can classify the molecular subtypes of COAD.

The tutorial includes the below analysis procedures:
* gene-level transformation
* normalization
* feature selection
* classification accuracy measurement and
* plotting of the results










