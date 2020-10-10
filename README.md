# maxquant.tools
This is an R package for reading MaxQuant protein and peptide data, currently from label-free experiments only.

Installation
------------

```r
# Get the the development version from GitHub:
# install.packages("devtools")
devtools::install_github("slurmclassic/maxquant.tools")
```

Usage
-----

```r
# Read a proteinGroups.txt file.
protein_group_data = readMaxQuantProteinMeasurements('examples/proteinGroups.txt')

# Read a peptides.txt file.
peptide_data = readMaxQuantPeptideMeasurements('examples/peptides.txt')

# Merge two protein group datasets.
a = readMaxQuantProteinMeasurements('examples/proteinGroups1.txt')
b = readMaxQuantProteinMeasurements('examples/proteinGroups2.txt')
merged_dataset = a + b

# Or use the underlying merge function with a list of datasets.
merged_dataset = mergeMaxQuantProteinGroupsResults(list(a, b), 'Merged Dataset Name')

# Merging two peptide datasets is very similar.

a = readMaxQuantPeptideMeasurements('examples/peptides1.txt')
b = readMaxQuantPeptideMeasurements('examples/peptides2.txt')
merged_dataset = a + b

# Or use the underlying merge function with a list of datasets.
merged_dataset = mergeMaxQuantPeptideResults(list(a, b), 'Merged Dataset Name')
```




