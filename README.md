# Read me

## Viral detection in wastewaters: comparison of bioinformatic tools for virome analysis

This script envolves all the steps from start (getting the original reports from each software) to finish
(producing the heatmaps).
We only care about non-phage viruses.

## Workflow of the script:
1. Get viral percent for excel sheet; each sheet will be for a software.
2. Get the viral details (number of reads, viral percent, etc).
3. Get full taxonomy for each indentification.
4. Grouping everything in one and formatting it accordingly for an excel file for for R plotting.


It was created in a way such that you should follow the following rules 
in order for it to work correctly.

## Steps:
### 1 - Create a folder for the analysis (e.g. metagenomic_analysis_01_01_2023)

### 2.1 - Place the script (from_start_to_finish_all_reports_to_ready_for_heatmap.py) into that folder.

### 2.2 - Download the folder "filtering_phages" present in the GitHub rep and place it inside the folder of the analysis.

### 3 - Create an text file called "pool_samples_reads.txt" and a file called "rep_samples_reads.txt", each file should contain the information for the raw number of reads of each sample with the following sctructure:

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
