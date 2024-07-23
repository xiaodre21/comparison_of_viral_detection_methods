# Viral detection in wastewaters: comparison of bioinformatic tools for virome analysis
## Read me

### Dependencies

This project requires the following Python packages:

- `os` (standard library)
- `random` (standard library)
- `networkx`
- `requests`
- `pandas`
- `ete3`
- `warnings` (standard library)
- `gzip` (standard library)
- `numpy` (standard library)


### Installation

You can install the required packages using `pip`. It is recommended to use a virtual environment to manage your dependencies.

### Using pip

```bash
pip install networkx requests pandas ete3 biopython numpy
```
<br />

# Description

This script envolves all the steps from start 





## Workflow of the script:
**1. Raw Reads:**<br />
- Read and count the raw reads from all fastq.gz files.

**2. Taxonomy Data:**<br />
- Update and download necessary taxonomy databases.

**3. Reports Processing:**<br />
- Read, store, and ensure completeness of all reports from various software tools.<br />
- Filter out invalid data.

**4. Taxonomy Enrichment:**<br />
- Add detailed taxonomy information to each sample.

**5. Pathogen Identification:**<br />
- Identify and flag potentially pathogenic genera in the samples.

**6. Data Merging:**<br />
- Combine data from all samples into a format suitable for R analysis.

**7. Host Source and Classification:**<br />
- Enrich the data with host source information and classification levels.

**8. Summary and OTU Tables:**<br />
- Create summary tables and OTU tables for further analysis and reporting.

**9. Output Writing:**<br />
- Save the processed data and tables to Excel files.
<br />

## It was created in a way such that you should follow the following steps in order for it to work correctly.

#### 1 - Download the template analysis folder "metagenomic_analysis" from this repository.

#### 2 - Open the "metadata_table.xlsx" and fill it accordingly. The file contains instructions.

#### 3 - Open the fastq.gz folder:
- For pooled sequencing files, place them onto the "pools" folders
- For replicate sequencing files, place them onto the "replicates" folders

#### 4.1 - Open the pool_results folder and navegate the different software folders:
- CZ.ID reports for pool samples should be put inside the folder "reports_czid"
- GenomeDetective reports for pool samples should be put inside the folder "reports_gd"
- INSaFLU reports for pool samples should be put inside the folder "reports_insaflu"
- Kraken2 reports for pool samples should be put inside the folder "reports_kraken2"

#### 4.2 - Open the rep_results folder and navegate the different software folders:
- CZ.ID reports for replicate samples should be put inside the folder "reports_czid"
- GenomeDetective reports for replicate samples should be put inside the folder "reports_gd"
- INSaFLU reports for replicate samples should be put inside the folder "reports_insaflu"
- Kraken2 reports for replicate samples should be put inside the folder "reports_kraken2"


#### 5 - Run the python script in commandline:
Example
```bash
python script_name_to_decide.py --update_ncbi True --filter_warnings True --tax_levels ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'] --phage_hosts ['bacteria', 'archaea']
```

### 6 - If you followed all the steps so far, your analysis folder should look like:
- fastq.gz
- main
- pool_results
- rep_results
- utils
- metadata_table.xlsx

### 7 - After you've ran the script the analysis folder should look like:
- fastq.gz
- main
- pool_results
- rep_results
- script_results (new)
- utils
- latest_virus_taxonomy.xlsx (new)
- metadata_table.xlsx


