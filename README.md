# Single-cell FDR benchmark

Code for the main statistical methods in *Benchmarking Test Statistics and False Discovery Rate Procedures for Single-Cell Differential Expression*.

- `main.py`: run a statistic and an FDR method.
- `methods.py`: four statistics and BH, working e-BH, DS and TD.
- `official_methods.R`: DESeq2 and MAST functions.

The Python statistic names are `t`, `wilcoxon`, `deseq2` and `mast`. The last two are simplified scores, not the R packages. The default level is 0.05 and the working e-value uses c=2.

## Run

```sh
python -m pip install -r requirements.txt
python main.py --input counts.npz --statistic wilcoxon --rule BH --output selected.npz
```

`counts.npz` contains `XA` and `XB`: rows are cells and columns are genes, in the same order. Save them with `np.savez_compressed("counts.npz", XA=XA, XB=XB)`. The output is a Boolean rejection mask.

DS and TD also need `--method-seed`. The `stable_seed` function uses base seed 2025 and SHA-256 of experiment identifiers. For example, `stable_seed("ms-method", 0.3, 0, "wilcoxon", "DS")` returns `1731682018`.

For R, source `official_methods.R`. Use integer gene-by-cell counts, group levels `c("B", "A")` and aligned donor labels. Apply BH to the returned p-values.

Kang18 data: [GEO GSE96583](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583).

## Software

Python 3.12.7, NumPy 1.26.4 and SciPy 1.13.1.

R 4.5.2, [DESeq2](https://bioconductor.org/packages/DESeq2/) 1.50.2, [MAST](https://bioconductor.org/packages/MAST/) 1.36.0 and [SingleCellExperiment](https://bioconductor.org/packages/SingleCellExperiment/) 1.32.0.
