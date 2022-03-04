# PDBClean F.V.
PDBClean offers curation tools for structural ensemble deposited in the Protein Data Bank.

*For installation instructions, please see [below](#installation).*

## Tutorials

### Downloading and curating a structural ensemble

The overall protocol is broken down in elementary sequential steps described in the following notebooks

[0. Download a structural ensemble from the RCSB PDB](https://github.com/csblab/PDBClean/blob/master/notebooks/0.%20Download%20a%20structural%20ensemble%20from%20RCSB%20PDB.ipynb)

![](figures/fig_download_pdb.png)

[1. Cleaning the CIF files just downloaded](https://github.com/csblab/PDBClean/blob/master/notebooks/1.%20Cleaning%20the%20CIF%20files%20just%20downloaded.ipynb)

[2. Assign MolID to the entities found in the CIF files](https://github.com/csblab/PDBClean/blob/master/notebooks/2.%20Assign%20MolID%20to%20the%20entities%20found%20in%20the%20CIF%20files.ipynb)

![](figures/fig_curate_MolID.png)

[3. Standardize chain IDs](https://github.com/csblab/PDBClean/blob/master/notebooks/3.%20Chain%20ID%20standardization.ipynb)

![](figures/fig_curate_ChainID.png)

[4. Standardize residue IDs](https://github.com/csblab/PDBClean/blob/master/notebooks/4.%20Residue%20ID%20standardization.ipynb)

![](figures/fig_curate_ResID.png)

[5. Finalize curation](https://github.com/csblab/PDBClean/blob/master/notebooks/5.%20Finalize%20curation.ipynb)


### Sharing curated dataset

We provide some ways to upload datasets to and download datasets from [OSF](osf.io), together with our examples.

[Pulling and pushing datasets on OSF](https://github.com/csblab/PDBClean/blob/master/notebooks/Datasets%20in%20the%20cloud%20-%20how%20to%20pull%20and%20push.ipynb)

[List of datasets curated by the Levitt Lab](https://github.com/csblab/PDBClean/blob/master/notebooks/List%20of%20datasets%20curated%20by%20the%20Levitt%20Lab.ipynb)

### Chain and atom selection

[Chain and atom selection](https://github.com/csblab/PDBClean/blob/master/notebooks/Chain%20and%20atom%20selection.ipynb)

### Extracting a homogeneous dataset

For many types of analysis, one would need to be able to load the dataset as a feature-by-sample array that requires all samples to exhibit the same features. This homogeneization step is not unique

[Extracting a homogeneous dataset](https://github.com/csblab/PDBClean/blob/master/notebooks/Extracting%20a%20homogeneous%20dataset.ipynb)

### Analysis of the resulting dataset

#### Conformational heterogeneity

## Installation

### Download from Pypi

For now we only uploaded the package to [TestPypi](https://test.pypi.org/project/PDBClean/), so you also need to install the required tools listed below:

`pip install --index-url  https://test.pypi.org/simple/ --no-deps PDBClean`

### Download from Github

Assuming you have the required tools and libraries listed below, just type:

`git clone https://github.com/csblab/PDBClean.git`

`python setup.py install`

#### Requirements
- [muscle](http://www.drive5.com/muscle/downloads.htm)
- [biopython](https://biopython.org/wiki/Download)
