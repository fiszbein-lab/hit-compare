# hit-compare
Pairwise HIT index output comparison.

### Dependencies
- Python 3.9.13
- NumPy 1.23.2

### Usage
First, import the needed functions and parse the `.EXON` files:

```python
from hit_compare import read_exon, rep_merge, get_matrices

locus_to_exon1 = read_exon("path/to/exon1")
locus_to_exon2 = read_exon("path/to/exon2")
```
Alternatively, if a replicate is available, `rep_merge` can be used to require 
that exons maintain their positional classes over both replicates; exons that 
pass this requirement will be returned with their HIT indices averaged:
```python
locus_to_exon1 = rep_merge("path/to/exon1-rep1", "path/to/exon1-rep2")
locus_to_exon2 = rep_merge("path/to/exon2-rep1", "path/to/exon2-rep2")
```
Next, the parsed data are passed to the matrix-building function:
```python
index_matrix, id_matrix = get_matrices(locus_to_exon1, locus_to_exon2)
```

`index_matrix` and `id_matrix` can then be visualized with modules
such as [matplotlib.pyplot.imshow](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.imshow.html) or 
[seaborn.heatmap](https://seaborn.pydata.org/generated/seaborn.heatmap.html). 