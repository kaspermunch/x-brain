# x-brain


```
conda create -n x-brain -c conda-forge scipy scanpy anndata jupyterlab

conda activate x-brain
conda install -c conda-forge -c bioconda -c kaspermunch  jupyterlab biopython wget goatools geneinfo matplotlib-venn openpyxl pytables ipympl ipython nodejs plotly
conda env export > x-brain.yml
```
