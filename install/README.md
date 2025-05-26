# Installation


## Spatial-GT-seq RNA data pipeline

| Item         | Link                                                |
|--------------|-----------------------------------------------------|
| RNA_pipeline | https://github.com/MingyuYang-Yale/DBiT-seq         |



## Spatial-GT-seq DNA and single-cell DNA-seq pipeline

### Python and some key Python Packages

| Item            | Purpose      | Link                                                          |
|-----------------|--------------|---------------------------------------------------------------|
| Python (3.11.4) | programming  | https://www.anaconda.com/download                             |


* **Python 3.11** is required for compatibility with all included scripts.

Install the required Python dependencies using `pip`:

```bash
pip install pysam natsort numpy
```

### External Software Tools

The following tools are required and should be installed separately. These tools are used via absolute paths defined in the scripts.

```bash
conda install -c bioconda bwa samtools sambamba
```

You must modify the absolute paths in the scripts (`preprocess.py`, `preprocess1.py`) to reflect the location of these tools on your system:

```python
BWA_PATH = "/path/to/bwa"
SAMTOOLS_PATH = "/path/to/samtools"
SAMBAMBA_PATH = "/path/to/sambamba"
```


## R and some key R Packages

| Item               | Purpose       | Link                                                                     |
|--------------------|---------------|--------------------------------------------------------------------------|
| R language (4.3.2) | programming   | https://www.r-project.org/                                               |
| Seurat             | RNA           | https://satijalab.org/seurat/articles/install                            |
| Monocle3           | RNA           | https://cole-trapnell-lab.github.io/monocle3/                            |
| CellChat           | RNA           | https://github.com/sqjin/CellChat                                        |
| Copykit            | CNV           | https://github.com/navinlabcode/copykit                                  |
| ComplexHeatmap     | visualization | https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html  |
| ggplot2            | visualization | https://ggplot2.tidyverse.org/                                           |

---
