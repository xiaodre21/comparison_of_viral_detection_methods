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
- `re` (standard library)
- `warnings` (standard library)
- `gzip` (standard library)
- `biopython` (for `Bio.SeqIO`)
- `numpy`


### Installation

You can install the required packages using `pip`. It is recommended to use a virtual environment to manage your dependencies.

### Using pip

```bash
pip install networkx requests pandas ete3 biopython numpy
```
<br />

# Description

This script envolves all the steps from start (getting the original reports from each software) to finish
(producing the heatmaps).





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

#### 1 - Create a folder for the analysis (e.g. metagenomic_analysis_01_01_2023)

#### 2.1 - Place the script (place_to_decide.py) into that folder.

#### 2.2 - Download the folder "filtering_phages" present in the GitHub rep and place it inside the folder of the analysis.

#### 3 - Create an text file called "pool_samples_reads.txt" and a file called "rep_samples_reads.txt", each file should contain the information for the raw number of reads of each sample with the following sctructure:

<sample filename> <number of raw reads>

**Example: (for pools)**<br />
A 230493 <br />
B 1234932 <br />
H 2839302 <br />
It will extract the letter (the first character).<br /> <br />
**Example: (for replicates)**<br />
A1_S1 3209246 <br />
D2_S11 1404672 <br />
F3_S18 1615053 <br />
It will extract the sample code (the first two characters).<br /> <br />

***Note:*** **If for some reason the sequencer outputted fastq/fastq.gz files with the wrong labels, (if the sample names are wrong) and you wish to swap them according to a specific correspondence, please fill a file called "correspondence.txt" with the following format**:
original_sample_code desired_sample_code<br />
**Example:**<br />
A1_S1 C2_S3<br />
D2_S11 F1_S9<br />
F3_S18 H9_S5<br />

### 4.1 Create 2 folders for each type of samples (pools and replicates):
#### &ensp; 4.1.1 Folder called "pool_results"
#### &ensp; 4.1.2 Folder called "rep_results"
Place the "pool_samples_reads.txt" and "rep_samples_reads.txt" files their correspondent folders.
    
### 4.2 Create 4 subfolders inside each of the previous folders:
#### &ensp; 4.2.1 Folder called "reports_czid"
#### &ensp; 4.2.2 Folder called "reports_insaflu"
#### &ensp; 4.2.3 Folder called "reports_gd"
#### &ensp; 4.2.4 Folder called "reports_kraken2"<br /><br />
**Note:** You should now have 2 folders, "pool_results" and "rep_results", each with 4 folders inside, one for each software, as
described in the previous point (4.2) and the "x_samples_reads.txt" file.<br /><br />
    
    
If you followed all the steps so far, your analysis folder should look like this:
- script.py
- pool_results
- rep_results
- filtering_phages

### 6. Run the script<br />
### 7. After you've ran the script, 2 files will be created:
- all_but_phage_reps_and_pools_for_R.xlsx
- viral_percent_across_workflow_reps_and_pools.xlsx

So after you run the script, the analysis folder should be the following:
- from_start_to_finish_all_reports_to_ready_for_heatmap.py
- all_but_phage_reps_and_pools_for_R.xlsx
- viral_percent_across_workflow_reps_and_pools.xlsx
- pool_results
- rep_results
- filtering_phages


### 8. Go to the GitHub repository and download the R script file with all the code to produce heatmaps (creating_heatmaps_for_phage_and_all_but_phage.R) You will have to install all the packages needed.

After you run the R script, a file called "all_but_phage_families_pools_and_reps.jpeg" will be created.

In the end, the analysis folder should contain all of the following files:
- from_start_to_finish_all_reports_to_ready_for_heatmap.py
- all_but_phage_families_pools_and_reps.jpeg
- all_but_phage_reps_and_pools_for_R.xlsx
- viral_percent_across_workflow_reps_and_pools.xlsx
- pool_results
- rep_results
- filtering_phages
