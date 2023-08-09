# Identification of endogenous viral elements (EVE) in host genome 
## Table of Contents
- [Purpose](#purpose)
- [Installation](#installation)
- [Dependencies](#dependencies)
- [Usage](#usage)
- [Directory contents](#Directorycontents)
- [Contributing](#contributing)
- [License](#license)

### Purpose
The EVE_annotation scripts is developed to identify Endogenous Viral Elements (EVEs) within genomic sequences. It takes advantage of several established bioinformatics tools and custom scripts to accurately annotate and analyze potential EVEs.

## Installation
Clone this repository to your local machine and provide necessary permissions to the scripts using the following commands:

```bash
git clone https://github.com/carolebelliardo/EVE_annotation.git
chmod +x EVE_annotation
```

Then, download 
  1. Viral proteins on the NCBI using web page downloading service <path_to_db1>
  2. The complet NR database <path_to_db2> using the following command:

```bash
wget -Nc -o wgetNRfasta.log 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz'
gzip -d nr.gz
chmod 775 nr
```
### Dependencies ###
The EVE Detection Tool relies on the following external tools and libraries:
- [Diamond](https://github.com/bbuchfink/diamond): A sequence aligner for protein and translated DNA searches [^1].
- [Bedtools](https://bedtools.readthedocs.io/en/latest/): A powerful toolset for genome arithmetic[^2].
- [R](https://www.r-project.org/): A programming language for statistical computing and graphics[^3].
- Python 3: A programming language required for executing the provided Python scripts[^4].

The following R librairies are required:
```
library(data.table)
library(taxonomizr)
```

Ensure that these softwares are installed and accessible in your system's PATH before running the EVE annotation scripts.

## Usage
```bash
cd EVE_annotation
./EVE_annotation.sh <path_to_db1> <path_to_db2> <path_to_host_fasta>
```
Replace <path_to_db1>, <path_to_db2>, and <path_to_host_fasta> with the paths to the required database and input files.

### Directory contents ###
The EVE annotation directory consists of the following scripts:

    * EVE_annotation.sh: The main script that coordinates the entire EVE annotation process.
    * bestHitsToFasta.py: Python script to parse BLAST results for identifying best hits and returning a FASTA file.
    * absolutPosi.R: R script for absolute coordinate calculations.
    * Get_EVE_annotation_summary.r: R script to generate a summary of EVE lineages.
    * addFamily.py: Python script to add family information to EVE annotations.




### Contributing ###
Contributions to the EVE Detection Tool are welcome! If you find any issues or have suggestions for improvements, please open an issue or submit a pull request in the GitHub repository.

### Authors ###
* Carole Belliardo: carole.belliardo@inrae.fr
* Clément Gilbert: clement.gilbert1@universite-paris-saclay.fr



### Funding ###
This work was supported by the plant heath and environment departement of the INRAe and MYCOPHYTO sas.  


### Acknowlegements ###
We acknowledge the people involved in generating the transcriptomic and genome assembly data for their publications and for making the data publicly available. We also would like to acknowledge members of the ITN Insect Doctors consortium and the current members of IRBI who provided some feedback throughout the study. 

### License ###

This project is licensed under the MIT License.



[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6554302.svg)](https://doi.org/10.5281/zenodo.6554302)


## References

[^1]: Buchfink, B., Xie, C., & Huson, D. H. (2015). Fast and sensitive protein alignment using DIAMOND. Nature Methods, 12(1), 59-60. [DOI: 10.1038/nmeth.3176](https://doi.org/10.1038/nmeth.3176)

[^2]: Quinlan, A. R., & Hall, I. M. (2010). BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics, 26(6), 841-842. [DOI: 10.1093/bioinformatics/btq033](https://doi.org/10.1093/bioinformatics/btq033)

[^3]: R Core Team (2021). R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria. [URL: https://www.R-project.org/](https://www.R-project.org/)

[^4]: Python Software Foundation. (2021). Python Language Reference, version 3.9.6. [URL: https://www.python.org/](https://www.python.org/)


