


import os
import random
import networkx as nx
import requests
import pandas as pd
from ete3 import NCBITaxa
import warnings
import gzip



# 1. Read in all the raw number of reads
# 1.1 Function to read from a fastq.gz file
def getting_raw_reads(filename):
    """
    Decompress and read the information on fastq.gz file to count the amount of reads.

    :param filename: name of the file to read
    :return: int, number of reads
    """
    reads = 0
    with gzip.open(filename, 'rb') as read:
        for _ in read:
            seq = next(read)
            reads += 1
    return reads


# 1.2 Function to read all fastq.gz files
def read_in_all_raw_reads():
    """
    Iterate over each fastq.gz file and count the amount of reads.
    Stores the results in a dictionary where keys are the sample codes, and values are the amount of reads.
    The sample code is the first two characters of each filename, e.g. A5_S1_L001_R1_001.fastq.gz -> A5.

    :return: dic
    """
    raw_reads_for_all = {}

    # Go to
    os.chdir("./fastq.gz")

    for folder in os.listdir(os.getcwd()):
        os.chdir(f"./{folder}")
        filenames = [filename for filename in os.listdir(os.getcwd()) if "R2" not in filename]

        for filename in filenames:
            print("Processing: ", filename, "\n")
            raw_reads_for_all[filename[:2]] = [getting_raw_reads(filename)]

        os.chdir("../")

    os.chdir("../")

    raw_reads_for_all = pd.DataFrame(raw_reads_for_all)
    raw_reads_for_all = raw_reads_for_all.T
    raw_reads_for_all.columns = ["Raw_reads"]
    return raw_reads_for_all


# 2. Get the latest ICTV Taxonomy
# 2.1. Download the latest version
def download_latest_virus_taxonomy():
    """
    Makes a HTTP request to the specified URL and writes the content of the response in an excel file.

    :return: None
    """
    # Url link we wish to request
    file_url = "https://ictv.global/vmr/current"
    # Make the request
    r = requests.get(file_url)
    # Save the restul in an excel file
    with open("latest_virus_taxonomy.xlsx", 'wb') as f:
        f.write(r.content)


# 2.2. Turn it into a pandas dataframe
def get_ictv_taxonomy_into_dataframe(taxonomic_levels):
    """
    Reads the ICTV taxonomy excel file and creates a pandas dataframe of it.

    This was obtained via ---> Virus Metadata Resource (VMR), https://ictv.global/vmr

    # The VMR: Exemplar Viruses for Species
    The ICTV chooses an exemplar virus for each species and the VMR provides a list of these exemplars.
    An exemplar virus serves as an example of a well-characterized virus isolate of that species and includes
    the GenBank accession number for the genomic sequence of the isolate as well as the virus name, isolate designation,
    suggested abbreviation, genome composition, and host source.

    :param taxonomic_levels: list, the taxonomic levels which we wish to retain
    :return: pd.Dataframe (ICTV Taxonomy)
    """
    # Read the excel file
    ictv_taxonomy = pd.read_excel("latest_virus_taxonomy.xlsx")
    # Remove unwanted columns
    ictv_taxonomy.drop([col for col in ictv_taxonomy.columns if "Sub" in col], axis=1)
    # Make everything lowercase to match our data
    ictv_taxonomy.columns = [x.lower() for x in ictv_taxonomy.columns]
    # Fill missing values with "undefined"
    ictv_taxonomy[taxonomic_levels] = ictv_taxonomy[taxonomic_levels].fillna('undefined')

    return ictv_taxonomy


# 2.3. Create a graph to transverse the taxonomy and obtain the host source
def create_taxonomy_graph(ictv_taxonomy, identification_data, taxonomic_levels, long_format=False):
    """
    Create a taxonomy graph and add host source and graph path to the identification data.

    :param ictv_taxonomy: pd.Dataframe, (pandas dataframe containing the ICTV taxonomy data)
    :param identification_data: pd.Dataframe, (pandas dataframe containing the identification data)
    :param taxonomic_levels: list, (list of taxonomic levels to consider)
    :param long_format: bool, (True for identification_data in long format, False otherwise)
    :return: pd.Dataframe, (Updated identification_data pandas dataframe with host source and graph path columns)
    """

    # If in long format, remove prefixes from taxonomic level columns
    if long_format:
        identification_data.columns = [x.lower() for x in identification_data.columns]

    # Create a directed graph
    G = nx.DiGraph()

    # Add nodes and edges to the graph from ICTV taxonomy, including "undefined" nodes
    for _, row in ictv_taxonomy.iterrows():
        for i in range(len(taxonomic_levels) - 1):
            parent = row[taxonomic_levels[i]]
            child = row[taxonomic_levels[i + 1]]
            G.add_edge(parent, child)
            if 'host_source' not in G.nodes[child]:
                G.nodes[child]['host_source'] = row['host source']
            if 'undefined' not in G:
                G.add_node('undefined')
            G.add_edge(parent, 'undefined')
            G.add_edge('undefined', child)

    def search_host_source(record, G):
        for level in range(len(taxonomic_levels), 0, -1):
            for i in range(level):
                node = record[taxonomic_levels[i]][3:] if long_format else record[taxonomic_levels[i]]
                if node in G and 'host_source' in G.nodes[node]:
                    return G.nodes[node]['host_source']
        return 'unknown'

    def trace_graph_path(record, G):
        path = []
        host_source = 'unknown'
        for level in range(len(taxonomic_levels)):
            node = record[taxonomic_levels[level]][3:] if long_format else record[taxonomic_levels[level]]
            if node != 'undefined' and node in G:
                path.append(node)
                if 'host_source' in G.nodes[node]:
                    host_source = G.nodes[node]['host_source']
        path.append(f"Host Source: {host_source}")
        return " -> ".join(path) if path else "Path not found"

    host_sources = []
    graph_paths = []
    # Create a list of results for host sources and the graph path to append to the pandas dataframe
    for _, row in identification_data.iterrows():
        host_sources.append(search_host_source(row, G))
        graph_paths.append(trace_graph_path(row, G))
    # Add the host source and graph path to the pandas dataframe
    identification_data['host_source'] = host_sources
    identification_data['graph_path'] = graph_paths

    # Extract the host source
    def extract_final_result(graph_path):
        if "Host Source:" in graph_path:
            return graph_path.split("Host Source:")[-1].strip()
        return "unknown"

    identification_data = identification_data.drop(columns=['host_source'])
    identification_data['host_source'] = identification_data['graph_path'].apply(extract_final_result)

    return identification_data


# 3. Determine the classification level based on the amount of missing taxonomical levels
def determine_classification_level(row, tax_levels_for_classification):
    """
    Returns a string according to the level of classification:
    - full if every level is defined
    - partial if there's at least one level as undefined
    - undefined if all levels are undefined

    :param row: pd.Dataframe row (for usage in lambda functions)
    :param tax_levels_for_classification: list, (list of taxonomic levels to consider)
    :return: str
    """
    if all(row[tax_levels_for_classification] != 'undefined'):
        return 'full'
    elif all(row[tax_levels_for_classification] == 'undefined'):
        return 'undefined'
    else:
        return 'partial'


# 4. Read all the reports and store them into pandas dataframes
# 4.1. Read the czid reports
def read_cz_report_into_pd():
    """
    Scans through both folders (replicates and pools) for each report for CZ.ID.
    Reads each report, reformats column names and retains only specific information:
    - TaxID (NCBI taxid for the reference it mapped to)
    - Assignment (the designation of the organism)
    - Reads (amount of reads mapping to that organism)
    - Known_pathogen (whether it's potencially pathogenic for humans or not)

    :return: dic, (a dictionary where the keys are sample names (str) and values are pd.Dataframe)
    """
    res = {}

    # Save og dir
    cwd = os.getcwd()

    # Jump to CZ.ID's replicate reports folder
    os.chdir("rep_results/reports_czid")

    # Get all the sample filenames to loop
    filenames = os.listdir(os.getcwd())

    # Iterate over each report and read it as pandas dataframe
    for file in filenames:

        # Read as pandas dataframe
        data = pd.read_csv(file)

        data.drop(data.tail(1).index, inplace=True)

        # Rename the pertinent columns to a standard
        data = data.rename(columns={"taxId": "TaxID",
                                    "nt.count": "Reads",
                                    "name": "Assignment"})

        # Keep only the desired columns
        data = data.loc[:, ['TaxID', 'Assignment', 'Reads', "known_pathogen"]]

        # Reformat column dtypes to int, otherwise NaN
        data["Reads"] = pd.to_numeric(data["Reads"], errors='coerce')
        data['TaxID'] = pd.to_numeric(data['TaxID'], errors='coerce')

        # Append the pandas dataframe to the dic
        res[file[:2]] = data

    # Move back to original directory
    os.chdir(cwd)

    # Jump to CZ.ID's Pool reports folder
    os.chdir("pool_results/reports_czid")

    # Get all the sample filenames to loop
    filenames = os.listdir(os.getcwd())

    # Iterate over each report and read it as pandas dataframe
    for file in filenames:

        # Read as pandas dataframe
        data = pd.read_csv(file)

        # Rename the pertinent columns to a standard (apply for all softwares)
        data = data.rename(columns={"tax_id": "TaxID",
                                    "nt_count": "Reads",
                                    "name": "Assignment"})

        # Keep only the desired columns
        data = data.loc[:, ['TaxID', 'Assignment', 'Reads', "known_pathogen"]]

        # Reformat column dtypes to int, otherwise NaN
        data["Reads"] = pd.to_numeric(data["Reads"], errors='coerce')
        data['TaxID'] = pd.to_numeric(data['TaxID'], errors='coerce')

        # Append the pandas dataframe to the dic
        res[file[:2]] = data

    # Change back to original directory
    os.chdir(cwd)

    return res


# 4.2. Read the INSaFLU-TELEVIR reports
def read_insa_report_into_pd():
    """
    Scans through both folders (replicates and pools) for each report for INSaFLU.
    It reads them, reformats column names and retains only specific information:
    - TaxID (NCBI taxid for the reference it mapped to)
    - Assignment (the designation of the organism)
    - Reads (amount of reads mapping to that organism)

    :return: dic, (a dictionary where the keys are sample names (str) and values are pd.Dataframe)
    """
    res = {}

    # Save working directory
    cwd = os.getcwd()

    # Jump to INSaFLU's replicate reports folder
    os.chdir("rep_results/reports_insaflu")

    # Get all the sample filenames to loop
    filenames = os.listdir(os.getcwd())

    # Iterate over each replicate report and read it as pandas dataframe
    for file in filenames:

        # This file is sometimes present in MAC computers, had to filter out
        if file != '.DS_Store':

            # Read as pandas dataframe
            data = pd.read_csv(file, sep="\t")

            # Group by 'taxid' and get the index of the maximum 'counts' for each group
            idx = data.groupby('taxid')['counts'].idxmax()

            # Use these indices to filter the original dataframe
            data = data.loc[idx].reset_index(drop=True)

            if "Taxid" in data.columns:
            # Rename the pertinent columns to a standard
                data = data.rename(columns={"Taxid": "TaxID",
                                        "Mapped reads": "Reads",
                                        "Description": 'Assignment'})
            else:
                data = data.rename(columns={"taxid": "TaxID",
                                        "counts": "Reads",
                                        "description": 'Assignment'})

            # Keep only the desired columns
            data = data.loc[:, ['TaxID', 'Assignment', 'Reads']]

            # Reformat column dtypes into int, otherwise NaN
            data["Reads"] = pd.to_numeric(data["Reads"], errors='coerce')
            data['TaxID'] = pd.to_numeric(data['TaxID'], errors='coerce')

            # Append the pandas dataframe to the dic
            res[file[:2]] = data

    os.chdir(cwd)

    # Jump to INSaFLU's pool reports folder
    os.chdir("pool_results/reports_insaflu")

    # Get all the sample filenames to loop
    filenames = os.listdir(os.getcwd())

    # Iterate over each pool report and read it as pandas dataframe
    for file in filenames:
        if file != '.DS_Store':
            data = pd.read_csv(file, sep="\t")

            # Group by 'taxid' and get the index of the maximum 'counts' for each group
            idx = data.groupby('taxid')['counts'].idxmax()

            # Use these indices to filter the original dataframe
            data = data.loc[idx].reset_index(drop=True)

            if not data.empty:
                # Rename the pertinent columns to a standard (apply for all softwares)
                if "Taxid" in data.columns:
                    # Rename the pertinent columns to a standard (apply for all softwares)
                    data = data.rename(columns={"Taxid": "TaxID",
                                                "Mapped reads": "Reads",
                                                "Description": 'Assignment'})
                else:
                    data = data.rename(columns={"taxid": "TaxID",
                                                "counts": "Reads",
                                                "description": 'Assignment'})

                # Keep only pertinent columns
                data = data.loc[:, ['TaxID', 'Assignment', 'Reads']]

                # Append the pandas dataframe to the dic
                res[file[:2]] = data

    os.chdir(cwd)

    return res


# 4.3. Read the Kraken2 reports
def read_kr_report_into_pd():
    """
    Scans through both folders (replicates and pools) for each report for Kraken2.
    It reads them, reformats column names and retains only specific information:
    - TaxID (NCBI taxid for the reference it mapped to)
    - Reads (amount of reads mapping to that organism)

    :return: dic, (a dictionary where the keys are sample names (str) and values are pd.Dataframe)
    """
    res = {}

    # Save working directory
    cwd = os.getcwd()

    # Jump to Kraken2's replicates reports folder
    os.chdir("rep_results/reports_kraken2")

    # Read excel and have each sample as a separate pd
    excel_all_samples_per_sheet = pd.read_excel("results_kraken2_galaxy.xlsx", sheet_name=None)

    # Get the sample codes to iterate
    samples = excel_all_samples_per_sheet.keys()

    # Iterate over each excel sheet (sample data)
    for sample_code in samples:
        current_sample_pd = excel_all_samples_per_sheet[sample_code]

        # Rename the pertinent columns to a standard
        current_sample_pd.columns = ["Percent", "Number of Reads in Clade", "Reads", "Rank", "TaxID", "Taxonomy"]

        # Keep only pertinent columns
        current_sample_pd = current_sample_pd.loc[:, ['TaxID', 'Reads']]

        # Append the pandas dataframe to the dic
        res[sample_code[:2]] = current_sample_pd

    os.chdir(cwd)

    # Jump to Kraken2's pool reports folder
    os.chdir("pool_results/reports_kraken2")

    # Read excel and have each sample as a separate pd
    excel_all_samples_per_sheet = pd.read_excel("Results_kraken2_reports_trimmed.xlsx", sheet_name=None)

    # Get the sample codes to iterate
    samples = excel_all_samples_per_sheet.keys()

    # Iterate over each excel sheet (sample data)
    for sample_code in samples:
        current_sample_pd = excel_all_samples_per_sheet[sample_code]

        # Drop the last 2 useless columns
        current_sample_pd = current_sample_pd.drop(current_sample_pd.iloc[:, 6:], axis=1)

        # Rename the pertinent columns to a standard (apply for all softwares)
        current_sample_pd.columns = ["Percent", "Number of Reads in Clade", "Reads", "Rank", "TaxID",
                                     "Taxonomy"]  # Add column names to pd dataframe

        # Keep only pertinent columns
        current_sample_pd = current_sample_pd.loc[:, ['TaxID', 'Reads']]

        # Append the pandas dataframe to the dic
        res[sample_code] = current_sample_pd

    os.chdir(cwd)

    return res


# 4.4. Read the Genome Detective reports
def read_gd_report_into_pd():
    """
    Scans through both folders (replicates and pools) for each report for Genome Detective.
    It reads them, reformats column names and retains only specific information:
    - TaxID (NCBI taxid for the reference it mapped to)
    - Assignment (the designation of the organism)
    - Reads (amount of reads mapping to that organism)

    :return: dic, (a dictionary where the keys are sample names (str) and values are pd.Dataframe)
    """
    res = {}

    # Save working directory
    cwd = os.getcwd()

    # Jump to GenomeDetective's replicates reports folder
    os.chdir("rep_results/reports_gd")

    # Get all the sample filenames to loop
    filenames = os.listdir(os.getcwd())

    # Iterate over each replicate report
    for file in filenames:

        # Read as a pandas dataframe
        data = pd.read_csv(file)

        # Rename the pertinent columns to a standard
        data = data.rename(columns={"TaxonomyId": "TaxID",
                                    "NbReads": "Reads",
                                    "StrainName": "Assignment"})

        # Keep only pertinent columns
        data = data.loc[:, ['TaxID', 'Assignment', 'Reads']]

        # Append the pandas dataframe to the dic
        res[file[:2]] = data

    os.chdir(cwd)

    # Jump to GenomeDetective's pool reports folder
    os.chdir("pool_results/reports_gd")

    # Get all the sample filenames to loop
    filenames = os.listdir(os.getcwd())

    # Iterate over each pool report
    for file in filenames:

        # Read as a pandas dataframe
        data = pd.read_csv(file)

        # Rename the pertinent columns to a standard (apply for all softwares)
        data = data.rename(columns={"TaxonomyId": "TaxID",
                                    "NbReads": "Reads",
                                    "StrainName": "Assignment"})

        # Select only pertinent columns
        data = data.loc[:, ['TaxID', 'Assignment', 'Reads']]

        # Append the pandas dataframe to the dic
        res[file[:2]] = data

    os.chdir(cwd)

    return res

# 4.5. Function to read all reports from each software and store them in a dic
def read_all_reports_and_store_into_df():
    """
    The function calls specific report-reading functions for each software tool:
    - read_cz_report_into_pd(): Reads the report from CZ.ID and stores the resulting DataFrames into the dictionary
                                with the key "czid".
    - read_insa_report_into_pd(): Reads the report from INSaFLU, parameter, and stores the resulting DataFrames
                                into the dictionary with the key "insa".
    - read_kr_report_into_pd(): Reads the report from Kraken2 and stores the resulting DataFrames into the dictionary
                                with the key "kr".
    - read_gd_report_into_pd(): Reads the report from GenomeDetective and stores the resulting DataFrames into the
                                dictionary with the key "gd".

    :return: dic, (dictionary all the data, where keys are software names (str), values are dictionaries with keys as
                    sample names (str) and values as pd.Dataframe)
    """
    all_samples = {}

    # CZ.ID
    all_samples["czid"] = read_cz_report_into_pd()

    # INSaFLU
    all_samples["insa"] = read_insa_report_into_pd()

    # Kraken2
    all_samples["kr"] = read_kr_report_into_pd()

    # GenomeDetective
    all_samples["gd"] = read_gd_report_into_pd()

    return all_samples


# 5. Add the full taxonomy to all reports
# 5.1 Function to obtain a dictionary with all the taxonomy
def get_desired_ranks(taxid, desired_ranks, ncbi):
    """
    Retrieves the taxonomic ranks for a given taxonomic ID (taxid) using the NCBI taxonomy database.

    :param taxid: int, (the taxonomic ID for which the lineage and ranks are to be retrieved)
    :param desired_ranks: list, (list of taxonomic ranks that you want to retrieve, e.g., 'superkingdom', 'kingdom', 'phylum', etc)
    :param ncbi: object, (provides methods to interact with the NCBI taxonomy database)
    :return: dic, (dictionary where the keys are the desired ranks suffixed with '_id',
                    and the values are the corresponding taxids. If a desired rank is
                    not present in the lineage, it assigns '<not present>' as the value)
    """
    try:
        lineage = ncbi.get_lineage(taxid)
        names = ncbi.get_taxid_translator(lineage)
    except (ValueError, TypeError) as e:
        return {'superkingdom_id': 1, 'kingdom_id': 1, 'phylum_id': 1, 'class_id': 1, 'order_id': 1, 'family_id': 1,
                'genus_id': 1, 'species_id': 1}

    lineage2ranks = ncbi.get_rank(names)
    ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
    return {'{}_id'.format(rank): ranks2lineage.get(rank, '<not present>') for rank in desired_ranks}


# 5.2. Function to create an empty template for the taxonomy dictionary
def create_taxo_dic(desired_ranks):
    """
    Creates and returns a dictionary with keys representing various taxonomic ranks and values initialized
    as empty lists. This structure can be used to store taxonomic information corresponding to each rank.

    :return: dic, (returns the taxo_dic dictionary containing the taxonomic ranks as keys and empty lists as values)
    """
    taxo_dic = {}
    for rank in desired_ranks:
        taxo_dic[rank] = []
    return taxo_dic


# 5.3. Create a taxonomy dataframe for all identifications in the sample to then merge with the sample dataframe
def characterize_sample(tax_ids, ncbi):
    """
    Takes a list of taxonomic IDs and an NCBI object to retrieve and translate taxonomic information
    (their respective taxonomic ranks), populating a pandas DataFrame with the taxonomic classifications for each
    identification.
    This pd.Dataframe provides an organized representation of the taxonomic hierarchy for the input IDs.


    :param tax_ids: list, (list of taxonomic IDs to be characterized)
    :param ncbi: object, (provides methods to interact with the NCBI taxonomy database)
    :return: pd.Dataframe, (pd.Dataframe with the taxonomy for each of the desired ranks)
    """
    # Select the desired ranks
    desired_ranks = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

    # Create the empty taxonomy dictionary
    taxonomy_dic = create_taxo_dic(desired_ranks)
    results = []

    # Iterate over all the taxids identified
    for tax_id in tax_ids:
        results.append(list())
        results[-1].append(str(tax_id))
        # Get the ranks (number) for the taxid
        ranks = get_desired_ranks(tax_id, desired_ranks, ncbi)
        for key, rank in ranks.items():
            if rank != '<not present>':
                # Append the taxonomical name for that rank
                results[-1].append(list(ncbi.get_taxid_translator([rank]).values())[0])
            else:
                results[-1].append(rank)

    # Populate the taxonomy dictionary
    for result in results:
        taxonomy_dic["superkingdom"].append(result[1])
        taxonomy_dic["kingdom"].append(result[2])
        taxonomy_dic["phylum"].append(result[3])
        taxonomy_dic["class"].append(result[4])
        taxonomy_dic["order"].append(result[5])
        taxonomy_dic["family"].append(result[6])
        taxonomy_dic["genus"].append(result[7])
        taxonomy_dic["species"].append(result[8])

    # Transform the dictionary into a pandas dataframe
    taxonomy_df = pd.DataFrame(taxonomy_dic)

    return taxonomy_df

# 5.4. Add the full taxonomical classification to each sample
def get_full_taxonomy(all_samples, ncbi, software_key):
    """
    Get the full taxonomy for the specified software and also extract a list of dataframes with the organisms flagged
    as potentially pathogenic for humans (in the case of CZ.ID).

    :param all_samples: dic, (dictionary containing all samples for different softwares)
    :param ncbi: object, (provides methods to interact with the NCBI taxonomy database)
    :param software_key: str, (key for the software to process, e.g., "czid", "insa", "kr", "gd")
    :return: Updated all_samples dictionary and known_pathogens list if applicable.
    """
    # Initialize the known_pathogens list if the software_key is "czid"
    known_pathogens = [] if software_key == "czid" else None

    # Iterate over samples
    for sample in all_samples[software_key].keys():

        # Read in the sample dataframe
        data = all_samples[software_key][sample]

        # Get all the taxids into a list
        data_taxids = data["TaxID"].tolist()

        # Characterize the sample
        taxonomy_df = characterize_sample(data_taxids, ncbi)

        # Need to reset indexes so concat doesn't fail
        data.reset_index(drop=True, inplace=True)
        taxonomy_df.reset_index(drop=True, inplace=True)

        # Concatenate the dataframes horizontally
        full_df = pd.concat([data, taxonomy_df], axis=1)

        # Replace <not present> and root values for "undefined"
        for col in full_df.columns:
            full_df[col] = full_df[col].replace("<not present>", "undefined")
            full_df[col] = full_df[col].replace("root", "undefined")

        # Reformat into all numbers
        full_df['Reads'] = pd.to_numeric(data['Reads'], errors='coerce').astype('Int64')
        full_df['TaxID'] = pd.to_numeric(data['TaxID'], errors='coerce').astype('Int64')

        if software_key == "czid":
            full_df['known_pathogen'] = pd.to_numeric(data['known_pathogen'], errors='coerce')
            full_pathogen_df = full_df[["Assignment", "Reads", 'known_pathogen', 'superkingdom', 'kingdom', 'phylum',
                                        'class', 'order', 'family', 'genus', 'species']].copy()
            known_pathogens.append(full_pathogen_df[full_pathogen_df["known_pathogen"] != 0])

        if software_key == "kr":
            full_df["Assignment"] = full_df['species']

        # Convert columns from float64 to string
        full_df[['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']] = full_df[
            ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']].astype(str)

        # Update the pandas dataframe on the 'all_samples' variable
        all_samples[software_key][sample] = full_df

    if software_key == "czid":
        return all_samples, known_pathogens
    else:
        return all_samples

# 5.5. Add the full taxonomical classification to all samples for each software
def add_taxonomy_to_all_workflows(all_samples, ncbi):
    """
    Iterates through different workflows (e.g., 'czid', 'insa', 'kr', 'gd') and calls the get_full_taxonomy function
    to add taxonomic information to the sample data. It returns the updated sample data and known pathogens.

    :param all_samples: dic, (dictionary containing all samples for different softwares)
    :param ncbi: object, (provides methods to interact with the NCBI taxonomy database)
    :return: dic, list, (updated all_samples dictionary and known_pathogens list)
    """
    # CZ.ID
    all_samples, known_pathogens = get_full_taxonomy(all_samples, ncbi, "czid")
    # INSaFLU
    all_samples = get_full_taxonomy(all_samples, ncbi, "insa")
    # Kraken2
    all_samples = get_full_taxonomy(all_samples, ncbi, "kr")
    # GenomeDetective
    all_samples = get_full_taxonomy(all_samples, ncbi, "gd")

    return all_samples, known_pathogens


# 6. Add the pathogenic_or_not flag
# 6.1. Extract all the potentially pathogenic genera
def get_all_pathogenic_genera(known_pathogens):
    """
    Processes a list of DataFrames containing known pathogen information, removes duplicate and undefined entries,
    and returns a set of unique pathogenic genera.

    CZ.ID has a list of pathogens with known pathogenicity in immunocompetent human hosts as well as a few
    common causes of disease in immunocompromised human hosts (https://czid.org/pathogen_list). This list is used to
    flag each identification present in the reports. By extracting the taxid of flagged identifications,
    we can then obtain a list of all potentially pathogenic genera.

    :param known_pathogens: list, (list of pd.Dataframe)
    :return: list, (list of potentially pathogenic genera)
    """
    # Group all dataframe containing only flagged identifications (as potentially pathogenic)
    all_known_pathogens = pd.concat(known_pathogens)

    # Remove duplicate entries
    all_known_pathogens.drop_duplicates(
        subset=['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'], inplace=True)

    # Drop the "undefined" genera
    only_genus = all_known_pathogens[(all_known_pathogens["genus"] != "undefined")].copy()

    return set(only_genus["genus"].to_list())


# 6.2. Add the flag as a new column to all samples
def add_pathogen_flag(all_samples, all_known_pathogenic_genera):
    """
    Adds a new column to each sample DataFrame to flag whether the genus is known to be pathogenic based on the
    provided set of pathogenic genera.
    This helps in quickly identifying samples containing known pathogenic genera.

    :param all_samples: dic, (dictionary containing dictionaries of pd.Dataframes)
    :param all_known_pathogenic_genera: list, (list of genera that are known to be pathogenic)
    :return: all_samples: dic, (dictionary containing dictionaries of pd.Dataframes)
    """
    for software in all_samples.keys():
        for sample in all_samples[software].keys():
            all_samples[software][sample]["pathogenic_or_not"] = all_samples[software][sample]['genus'].apply(
                lambda x: 1 if x in all_known_pathogenic_genera else 0)

    return all_samples

# 7. Data cleanup
# 7.1. For the cases of empty sample reports, a dataframe still needs to exist for the loops not to conflitc
def reset_and_fill_dataframe(df):
    """
    Creates a new DataFrame with the same columns as the input DataFrame, but it resets the content to a single row
    filled with default values.

    :param df: pd.Dataframe
    :return: pd.Dataframe
    """
    # Create an empty dataframe with the same columns
    empty_df = pd.DataFrame(columns=df.columns)

    # Create a dictionary with default values based on the column types
    default_values = {}
    for col in df.columns:
        # If column is numeric, put 0
        if pd.api.types.is_numeric_dtype(df[col]):
            default_values[col] = 0
        else:
            default_values[col] = 'missing_info'

    # Create a DataFrame with the default row
    default_row_df = pd.DataFrame([default_values])

    # Concatenate the empty dataframe with the default row dataframe
    new_df = pd.concat([empty_df, default_row_df], ignore_index=True)

    return new_df


# 7.2. Check missing reports and fill in the gaps
def check_missing_data_and_fill(all_samples):
    """
    Ensures that all software datasets in all_samples have entries for all unique sample names by filling in any
    missing samples with a DataFrame that has default values.
    This ensures uniformity across the datasets.

    :param all_samples: dic, (dictionary containing dictionaries of pd.Dataframes)
    :return: all_samples: dic, (dictionary containing dictionaries of pd.Dataframes)
    """
    # Identify all unique sample names
    all_sample_names = set(sample for software in all_samples.values() for sample in software.keys())

    for software in all_samples.keys():
        for sample_name in all_sample_names:
            if sample_name not in all_samples[software].keys():
                random_sample = random.choice(list(all_samples[software].keys()))
                template = all_samples[software][random_sample].copy()
                all_samples[software][sample_name] = reset_and_fill_dataframe(template)
    return all_samples


# 7.3. Filter any identifications with 0 reads
def filter_out_ids_with_zero_reads(all_samples):
    """
    Processes each pd.Dataframe in all_samples to remove rows where the "Reads" column has a value of zero or
    is not numeric.
    This ensures that each sample DataFrame only contains rows with positive read counts.

    :param all_samples: dic, (dictionary containing dictionaries of pd.Dataframes)
    :return: all_samples: dic, (dictionary containing dictionaries of pd.Dataframes)
    """
    for software in all_samples.keys():
        for sample in all_samples[software].keys():
            if all_samples[software][sample] is not None:
                data = all_samples[software][sample].copy()
                data["Reads"] = pd.to_numeric(data["Reads"], errors='coerce')
                data = data[data['Reads'] > 0]
                all_samples[software][sample] = data

    return all_samples


# 8. Functions to extract data from the samples in lambda functions
# 8.1. Extract the replicate number
def get_replicate_number(row):
    """
    Returns the second element of a string.

    :param row: str
    :return: str
    """
    return row[1]


# 8.2. Extract the software used
def get_software(row):
    """
    Extracts the software name from a string by splitting the string at underscores and returning the second part.
    If an exception occurs, it returns 'Information missing'.

    :param row: str
    :return: str
    """
    try:
        return row.split("_")[1]
    except KeyError:
        return 'Information missing'


# 8.3. Extract the software used
def get_raw_reads(row, raw_reads_for_all):
    """
    Extracts raw number of reads information from pd.Dataframe based the row (sample name).

    :param row: str
    :param metadata: pd.Dataframe
    :return: str
    """
    try:
        return raw_reads_for_all.loc[row[:2], 'Raw_reads']
    except KeyError:
        return 'Raw_reads missing'



# 9. Function to extract data from the metadata dataframe, in lambda functions
def get_metadata_info(row, metadata, column_name):
    """
    Extracts specified information from the metadata DataFrame based on the sample name and column name.

    :param row: str, sample name
    :param metadata: pd.DataFrame, the metadata DataFrame
    :param column_name: str, the name of the column to retrieve information from
    :return: str, the extracted information or a default message if not found
    """
    try:
        if column_name == "<collection_step>":
            return "Step " + str(metadata.loc[row[:2], column_name])
        else:
            return metadata.loc[row[:2], column_name]
    except KeyError:
        return f'{column_name} missing'


# 10. Function to calculate sums of specific columns
# 10.1. Generic function which works for any column
def get_info(row, all_samples_consensus_result, from_, to_):
    """
    Calculates the sum of a specific column (from_) for each sample in the form of pd.Dataframe.
    It retrieves this information for a specific sample indicated by row

    :param row: str
    :param all_samples_consensus_result: pd.Dataframe
    :param from_: str, (column to perform sum)
    :param to_: str, (column to store the result)
    :return: int
    """
    result = all_samples_consensus_result.groupby('SampleID')[from_].sum().reset_index()
    result.columns = ['SampleID', to_]
    result = result.set_index("SampleID")
    try:
        return result.loc[row, to_]
    except KeyError:
        return 0


# 10.2. Function to calculate the sum of undefined reads for each sample
def get_undefined(row, all_samples_consensus_result):
    """
    Calculates the sum of undefined reads for each sample in the form of pd.Dataframe.
    It retrieves this information for a specific sample indicated by row

    :param row: str
    :param all_samples_consensus_result: pd.Dataframe
    :return: int
    """
    filtered_samples = all_samples_consensus_result[all_samples_consensus_result['host_source'].isin(['unknown'])]
    undefined_reads = filtered_samples.groupby('SampleID')['Reads'].sum().reset_index()
    undefined_reads.columns = ['SampleID', 'Number_of_undefined_reads']
    undefined_reads = undefined_reads.set_index("SampleID")
    try:
        return undefined_reads.loc[row, 'Number_of_undefined_reads']
    except KeyError:
        return 0


# 10.3. Function to calculate the sum of phage reads for each sample
def get_all_phage_reads(row, all_samples_consensus_result):
    """
    Calculates the sum of phage reads for each sample in the form of pd.Dataframe.
    It retrieves this information for a specific sample indicated by row

    :param row: str
    :param all_samples_consensus_result: pd.Dataframe
    :return: int
    """
    all_samples_consensus_result = all_samples_consensus_result[
        ~all_samples_consensus_result["host_source"].isin(["unknown"])].copy()
    phage_reads_per_sample = all_samples_consensus_result.groupby('SampleID')['Reads'].sum().reset_index()
    phage_reads_per_sample.columns = ['SampleID', 'Number_of_only_phage_reads']
    phage_reads_per_sample = phage_reads_per_sample.set_index("SampleID")
    try:
        return phage_reads_per_sample.loc[row, 'Number_of_only_phage_reads']
    except KeyError:
        return 0


# 10.4. Function to calculate the sum of fully identified non phage reads for each sample
def get_full_id_non_phage(row, all_samples_consensus_result):
    """
    Calculates the sum of fully identified non-phage reads for each sample in the form of pd.Dataframe.
    It retrieves this information for a specific sample indicated by row

    :param row: str
    :param all_samples_consensus_result: pd.Dataframe
    :return: int
    """
    # all_samples_consensus_result = all_samples_consensus_result.query(f'host_source not in {phage_hosts} & classification_level == "full"')
    mask_classification = all_samples_consensus_result['classification_level'] == 'full'
    filtered_samples = all_samples_consensus_result[mask_classification]
    full_id_non_phage_reads = filtered_samples.groupby('SampleID')['Reads'].sum().reset_index()
    full_id_non_phage_reads.columns = ['SampleID', 'Number_of_full_id_non_phage_reads']
    full_id_non_phage_reads = full_id_non_phage_reads.set_index("SampleID")
    try:
        return full_id_non_phage_reads.loc[row, 'Number_of_full_id_non_phage_reads']
    except KeyError:
        return 0


# 10.5. Function to calculate the sum of partially identified non phage reads for each sample
def get_partial_id_non_phage(row, all_samples_consensus_result):
    """
    Calculates the sum of partially identified non-phage reads for each sample in the form of pd.Dataframe.
    It retrieves this information for a specific sample indicated by row

    :param row: str
    :param all_samples_consensus_result: pd.Dataframe
    :return: int
    """
    mask_classification = all_samples_consensus_result['classification_level'] == 'partial'
    filtered_samples = all_samples_consensus_result[mask_classification]
    partial_id_non_phage_reads = filtered_samples.groupby('SampleID')['Reads'].sum().reset_index()
    partial_id_non_phage_reads.columns = ['SampleID', 'Number_of_partially_identified_non_phage_reads']
    partial_id_non_phage_reads = partial_id_non_phage_reads.set_index("SampleID")
    try:
        return partial_id_non_phage_reads.loc[row, 'Number_of_partially_identified_non_phage_reads']
    except KeyError:
        return 0


# 10.6. Function to calculate the sum of fully identified phage reads for each sample
def get_full_id_phage(row, all_samples_consensus_result):
    """
    Calculates the sum of fully identified phage reads for each sample in the form of pd.Dataframe.
    It retrieves this information for a specific sample indicated by row

    :param row: str
    :param all_samples_consensus_result: pd.Dataframe
    :return: int
    """
    mask_unknown_hosts = all_samples_consensus_result['host_source'].isin(['unknown'])
    mask_classification = all_samples_consensus_result['classification_level'] == 'full'
    filtered_samples = all_samples_consensus_result[mask_classification & ~mask_unknown_hosts]
    full_id_phage_reads = filtered_samples.groupby('SampleID')['Reads'].sum().reset_index()
    full_id_phage_reads.columns = ['SampleID', 'Number_of_full_id_phage_reads']
    full_id_phage_reads = full_id_phage_reads.set_index("SampleID")
    try:
        return full_id_phage_reads.loc[row, 'Number_of_full_id_phage_reads']
    except KeyError:
        return 0

# 10.7. Function to calculate the sum of partially identified phage reads for each sample
def get_partial_id_phage(row, all_samples_consensus_result):
    """
    Calculates the sum of partially identified phage reads for each sample in the form of pd.Dataframe.
    It retrieves this information for a specific sample indicated by row

    :param row: str
    :param all_samples_consensus_result: pd.Dataframe
    :return: int
    """
    mask_unknown_hosts = all_samples_consensus_result['host_source'].isin(['unknown'])
    mask_classification = all_samples_consensus_result['classification_level'] == 'partial'
    filtered_samples = all_samples_consensus_result[mask_classification & ~mask_unknown_hosts]
    partial_id_non_phage_reads = filtered_samples.groupby('SampleID')['Reads'].sum().reset_index()
    partial_id_non_phage_reads.columns = ['SampleID', 'Number_of_partially_identified_phage_reads']
    partial_id_non_phage_reads = partial_id_non_phage_reads.set_index("SampleID")
    try:
        return partial_id_non_phage_reads.loc[row, 'Number_of_partially_identified_phage_reads']
    except KeyError:
        return 0


# 10.8. Function to calculate the sum of non phage reads for each sample
def get_all_non_viral_reads(row, all_samples_consensus_result):
    """
    Calculates the sum of non-viral reads for each sample in the form of pd.Dataframe.
    It retrieves this information for a specific sample indicated by row

    :param row: str
    :param all_samples_consensus_result: pd.Dataframe
    :return: int
    """
    all_samples_consensus_result['SampleID'] = all_samples_consensus_result['Sample'] + '_' + \
                                               all_samples_consensus_result['Software']
    viral_reads_per_sample = all_samples_consensus_result.groupby('SampleID')['Reads'].sum().reset_index()
    viral_reads_per_sample.columns = ['SampleID', 'Number_of_all_non_viral_reads']
    viral_reads_per_sample = viral_reads_per_sample.set_index("SampleID")
    try:
        return viral_reads_per_sample.loc[row, 'Number_of_all_non_viral_reads']
    except KeyError:
        return 0


# 10.9. Function to calculate the sum of crAssphage reads for each sample
def get_crassphages(row, all_samples_consensus_result):
    """
    Calculates the sum of non-viral reads for each sample in the form of pd.Dataframe.
    It retrieves this information for a specific sample indicated by row

    :param row: str
    :param all_samples_consensus_result: pd.Dataframe
    :return: int
    """
    all_samples_consensus_result['SampleID'] = all_samples_consensus_result['Sample'] + '_' + \
                                               all_samples_consensus_result['Software']
    viral_reads_per_sample = all_samples_consensus_result.groupby('SampleID')['Reads'].sum().reset_index()
    viral_reads_per_sample.columns = ['SampleID', 'Number_of_crassphage_reads']
    viral_reads_per_sample = viral_reads_per_sample.set_index("SampleID")
    try:
        return viral_reads_per_sample.loc[row, 'Number_of_crassphage_reads']
    except KeyError:
        return 0

# 11. Function to create the "table_samples" dataframe,
#     a summary table with metrics and metadata for each sample
def populate_baseplate_df(all_samples_consensus_result, metadata, raw_reads_for_all, phage_hosts):
    """
    Creates and populates a DataFrame (baseplate_df) with various statistics derived from the pd.Dataframe
    all_samples_consensus_result and additional metadata.
    The function returns the populated pd.Dataframe.

    :param all_samples_consensus_result: pd.Dataframe
    :param metadata: pd.Dataframe
    :param raw_reads_for_all: pd.Dataframe
    :param phage_hosts: list, (list of phage host sources)
    :return: pd.Dataframe
    """
    # Initialize an empty pandas dataframe
    baseplate_df = pd.DataFrame()

    # Create a new column just for the sake of the functions based on row
    all_samples_consensus_result['SampleID'] = all_samples_consensus_result['Sample'] + '_' + \
                                               all_samples_consensus_result['Software']

    # Making sure its still an int
    all_samples_consensus_result["pathogenic_or_not"] = all_samples_consensus_result["pathogenic_or_not"].astype(int)

    # Filter out everything into sub-dataframes
    all_viruses = all_samples_consensus_result[all_samples_consensus_result["superkingdom"].isin(["Viruses"])].copy()
    all_but_viruses = all_samples_consensus_result[
        ~all_samples_consensus_result["superkingdom"].isin(["Viruses", "undefined"])].copy()
    undefined_stuff = all_samples_consensus_result[
        all_samples_consensus_result["superkingdom"].isin(["undefined"])].copy()
    all_non_phage = all_viruses[~all_viruses["host_source"].isin(phage_hosts)].copy()
    all_phage = all_viruses[all_viruses["host_source"].isin(phage_hosts)].copy()

    all_crassphage = all_phage[all_phage["order"] == "Crassvirales"].copy()
    all_poly = all_non_phage[all_non_phage["family"].isin(["Polyomaviridae"])].copy()
    all_adeno = all_non_phage[all_non_phage["family"].isin(["Adenoviridae"])].copy()
    all_boca = all_non_phage[all_non_phage["genus"].isin(["Bocaparvovirus"])].copy()


    baseplate_df["SampleID"] = sorted(set(list(all_samples_consensus_result['SampleID'])))
    baseplate_df["Sample"] = baseplate_df.SampleID.str[:2]
    baseplate_df["Replicate"] = baseplate_df["SampleID"].apply(get_replicate_number)
    baseplate_df["Software"] = baseplate_df["SampleID"].apply(get_software)

    # Constant columns
    constant_columns = [
        '<original sample code>',
        '<desired sample code>',
        '<sample_number>',
        '<sample_type>'
    ]

    # Identify non-constant columns in metadata
    non_constant_columns = [col for col in metadata.columns if col not in constant_columns]

    # Apply the get_metadata_info function to each non-constant column
    for column in non_constant_columns:
        baseplate_df[column[1:-1]] = baseplate_df['SampleID'].apply(lambda row: get_metadata_info(row, metadata, column))

    # Get the number of raw reads
    baseplate_df["Raw_reads"] = baseplate_df["SampleID"].apply(lambda row: get_raw_reads(row, raw_reads_for_all))

    baseplate_df["Number_of_all_reads"] = baseplate_df["SampleID"].apply(
        lambda row: get_info(row, all_samples_consensus_result, 'Reads', 'Number_of_all_reads'))
    baseplate_df["Number_of_all_viral_reads"] = baseplate_df["SampleID"].apply(
        lambda row: get_info(row, all_viruses, 'Reads', 'Number_of_all_viral_reads'))
    baseplate_df["Number_of_all_non_viral_reads"] = baseplate_df["SampleID"].apply(
        lambda row: get_all_non_viral_reads(row, all_but_viruses))
    baseplate_df["Number_of_undefined_reads"] = baseplate_df["SampleID"].apply(
        lambda row: get_undefined(row, undefined_stuff))


    baseplate_df["Number_of_non_phage_reads"] = baseplate_df["SampleID"].apply(
        lambda row: get_info(row, all_non_phage, 'Reads', 'Number_of_non_phage_reads'))
    baseplate_df["Number_of_only_phage_reads"] = baseplate_df["SampleID"].apply(
        lambda row: get_all_phage_reads(row, all_phage))


    baseplate_df["Number_of_full_id_non_phage_reads"] = baseplate_df["SampleID"].apply(
        lambda row: get_full_id_non_phage(row, all_non_phage))
    baseplate_df["Number_of_partially_identified_non_phage_reads"] = baseplate_df["SampleID"].apply(
        lambda row: get_partial_id_non_phage(row, all_non_phage))


    baseplate_df["Number_of_full_id_phage_reads"] = baseplate_df["SampleID"].apply(
        lambda row: get_full_id_phage(row, all_phage))
    baseplate_df["Number_of_partially_identified_phage_reads"] = baseplate_df["SampleID"].apply(
        lambda row: get_partial_id_phage(row, all_phage))


    baseplate_df["Number_of_total_pathogenic_IDs"] = baseplate_df["SampleID"].apply(
        lambda row: get_info(row, all_samples_consensus_result, 'pathogenic_or_not', 'Number_of_total_pathogenic_IDs'))
    baseplate_df["Number_of_viral_pathogenic_IDs"] = baseplate_df["SampleID"].apply(
        lambda row: get_info(row, all_viruses, 'pathogenic_or_not', 'Number_of_viral_pathogenic_IDs'))
    baseplate_df["Number_of_non_viral_pathogenic_IDs"] = baseplate_df["SampleID"].apply(
        lambda row: get_info(row, all_but_viruses, 'pathogenic_or_not', 'Number_of_non_viral_pathogenic_IDs'))


    baseplate_df["Number_of_non_phage_pathogenic_IDs"] = baseplate_df["SampleID"].apply(
        lambda row: get_info(row, all_non_phage, 'pathogenic_or_not', 'Number_of_non_phage_pathogenic_IDs'))
    baseplate_df["Number_of_phage_pathogenic_IDs"] = baseplate_df["SampleID"].apply(
        lambda row: get_info(row, all_phage, 'pathogenic_or_not', 'Number_of_phage_pathogenic_IDs'))


    baseplate_df["Number_of_crassphages_reads"] = baseplate_df["SampleID"].apply(
        lambda row: get_crassphages(row, all_crassphage))

    baseplate_df["Number_of_poly_reads"] = baseplate_df["SampleID"].apply(
        lambda row: get_crassphages(row, all_poly))

    baseplate_df["Number_of_adeno_reads"] = baseplate_df["SampleID"].apply(
        lambda row: get_crassphages(row, all_adeno))

    baseplate_df["Number_of_boca_reads"] = baseplate_df["SampleID"].apply(
        lambda row: get_crassphages(row, all_boca))


    return baseplate_df


# 12. Function to create the "table_samples" dataframe for the paper (simplified version)
def format_for_paper_table(tables_samples_df, metadata):
    """
    Formats a pd.Dataframe (tables_samples_df) into a table suitable for inclusion in a paper (as a table).

    :param tables_samples_df: pd.Dataframe
    :param metadata: pd.Dataframe
    :return: pd.Dataframe
    """
    def get_replicate_or_pool(row):
        if int(row) in (1, 2, 3):
            return "Replicate"
        else:
            return "Pool"

    def get_volume(row, metadata):
        return metadata.loc[row, '<volume_collected>']

    paper_table = pd.DataFrame()

    paper_table["Treament_Plant"] = tables_samples_df["treatment_plant"].copy()
    paper_table["Step"] = tables_samples_df["collection_step"]
    paper_table["SampleID"] = tables_samples_df["SampleID"].str[:2].copy()
    paper_table = paper_table.drop_duplicates(subset=['SampleID'])
    paper_table["Sample"] = tables_samples_df["SampleID"].str[:2].copy()
    paper_table["SampleID"] = tables_samples_df["SampleID"].str[0].copy()
    paper_table["Replicate/Pool"] = tables_samples_df["Replicate"].apply(lambda row: get_replicate_or_pool(row))
    paper_table["Sample Number"] = tables_samples_df["Replicate"].copy()
    paper_table["Volume Collected (L)"] = tables_samples_df["SampleID"].str[:2].apply(
        lambda row: get_volume(row, metadata))
    paper_table["Concentration"] = tables_samples_df["library_concentration"].copy()
    paper_table["Number of raw reads"] = tables_samples_df["Raw_reads"].copy()

    paper_table['Viral reads proportion CZID'] = [0] * len(paper_table)
    paper_table['Viral reads proportion GD'] = [0] * len(paper_table)
    paper_table['Viral reads proportion INSA'] = [0] * len(paper_table)
    paper_table['Viral reads proportion KR'] = [0] * len(paper_table)

    # Loop over each sample and software combination
    for sample in list(tables_samples_df["Sample"]):
        for software in ['CZID', 'GD', 'INSA', 'KR']:
            # Extract the relevant values from the original dataframe
            total_reads = tables_samples_df.loc[(tables_samples_df['Sample'] == sample) & (tables_samples_df['Software'] == software), 'Raw_reads']
            viral_reads = tables_samples_df.loc[(tables_samples_df['Sample'] == sample) & (tables_samples_df['Software'] == software), 'Number_of_all_viral_reads']

            # Calculate the viral proportion
            if not total_reads.empty and not viral_reads.empty:
                viral_proportion = viral_reads.values[0] / total_reads.values[0]
            else:
                viral_proportion = 0

            # Place the viral proportion in the correct place in df2
            paper_table.loc[paper_table['Sample'] == sample, f'Viral reads proportion {software}'] = viral_proportion

    # Sort them alphabetically to make it pretty
    paper_table = paper_table.sort_values(by=["Treament_Plant", "Step", "SampleID"], ascending=True)

    return paper_table


# 13. Retrieve the metadata file to extract info for each sample
# 13.1. Function to retrieve the metadata file and read it as an pandas dataframe
def get_metadata():
    """
    Reads the metadata table in order to get additional information about samples.

    :return: pd.Dataframe
    """
    metadata = pd.read_excel("metadata_table.xlsx", skiprows=10)
    metadata = metadata.set_index('<desired sample code>')
    return metadata


# 13.2. Function to retrieve the treatment plant for a sample, from the metadata
def get_treatment_plant(metadata, sample_code):
    """
    Extracts treatment plant information from pd.Dataframe (metadata) based the sample name.

    :param metadata: pd.Dataframe
    :param sample_code: str
    :return: str
    """
    return metadata.loc[sample_code, '<treatment plant>']


# 13.3. Function to retrieve the collection step for a sample, from the metadata
def get_collection_step(metadata, sample_code):
    """
    Extracts collection step information from pd.Dataframe (metadata) based the sample name..

    :param metadata: pd.Dataframe
    :param sample_code: str
    :return: str
    """
    return str(metadata.loc[sample_code, '<collection_step>'])


# 13.4. Add the extra details obtained from the metadata to the identification data
def extract_essencials_from_pd(virus_data, sample, sw):
    """
    Extracts essential metadata from a pd.Dataframe containing virus data, enriches it with additional information,
    and returns the modified pd.Dataframe along with the metadata

    :param virus_data: pd.Dataframe
    :param sample: str, (sample name)
    :param sw: string, (software)
    :return: pd.Dataframe, pd.Dataframe (modified pd.Dataframe and metadata)
    """
    # Constant columns
    constant_columns = [
        '<original sample code>',
        '<desired sample code>',
        '<sample_number>',
        '<sample_type>'
    ]

    metadata = get_metadata()

    virus_pd = virus_data.copy()

    software_analysis_list = [sw.upper() for _ in range(len(virus_pd))]
    sample_list = [sample for _ in range(len(virus_pd))]

    i = 0
    virus_pd.insert(i, "Location", ["1" for _ in range(len(virus_pd))])
    i+=1
    virus_pd.insert(i, "Sample", sample_list)

    # Identify non-constant columns in metadata
    non_constant_columns = [col for col in metadata.columns if col not in constant_columns]

    for idx, column in enumerate(non_constant_columns):
        i+=1
        virus_pd.insert(i, column[1:-1],
                        [metadata.loc[metadata['<original sample code>'] == row['Sample'], f'{column}'].values[0] for
                         _, row in virus_pd.iterrows()])

    i+=1
    virus_pd.insert(i, "Software", software_analysis_list)


    return virus_pd, metadata


# 14. Extract all the pertinent information and enrich it with extra details, for each dataframe
def extract_info_for_R(all_pds, all_samples, software_key):
    """
    Processes sample data for a specific software key, enriches it with metadata, and appends the processed
    pd.Dataframe's to a list.
    It returns the list of processed DataFrames and the metadata.

    :param all_pds: list, (list of pd.Dataframe)
    :param all_samples: dic, (dictionary containing all samples for different softwares)
    :param software_key:
    :return: list, pd.Dataframe
    """
    sample_names = all_samples[software_key].keys()

    for sample in sample_names:
        pd_virus, metadata = extract_essencials_from_pd(all_samples[software_key][sample], sample, software_key)
        pd_virus.reset_index(drop=True, inplace=True)
        all_pds.append(pd_virus)

    return all_pds, metadata


# 15. Merge all dataframes into one
def merge_everything_for_r(non_phage_reads):
    """
    A fucntion to extract all the information from each sample's dataframe and merge everything into one

    :param non_phage_reads: dic
    :return: pandas.Dataframe, pandas.Dataframe
    """
    all_non_phage_pds = []
    all_non_phage_pds, _ = extract_info_for_R(all_non_phage_pds, non_phage_reads, software_key="czid")
    all_non_phage_pds, _ = extract_info_for_R(all_non_phage_pds, non_phage_reads, software_key="insa")
    all_non_phage_pds, _ = extract_info_for_R(all_non_phage_pds, non_phage_reads, software_key="kr")
    all_non_phage_pds, metadata = extract_info_for_R(all_non_phage_pds, non_phage_reads, software_key="gd")

    # Merge everything
    consensus_result = pd.concat(all_non_phage_pds, ignore_index=True)
    return consensus_result, metadata


# 16. A function to concatenate taxonomic ranks into a single string
def concatenate_taxonomic_ranks(row):
    """
    Join all taxonomic ranks and separate them by a comma to simulate the OTU Table output of QIIME.
    :param row:
    :return:
    """
    return ';'.join([
        f'd__{row["superkingdom"]}',  # This is actually domain
        f'p__{row["phylum"]}',
        f'c__{row["class"]}',
        f'o__{row["order"]}',
        f'f__{row["family"]}',
        f'g__{row["genus"]}',
        f's__{row["species"]}'
    ])


# 17. A function to transform the data into long format
def pivot_dataframe(dataframe):
    """
    Transform dataframe into long format

    :param dataframe: pd.Dataframe
    :return: pd.Dataframe
    """
    # Pivot the dataframe to the desired format with taxonomic assignments as rows and sample_software as columns
    pivoted_df = dataframe.pivot_table(index='taxonomic_assignment', columns='SampleID', values='Reads',
                                       aggfunc='sum', fill_value=0)

    # Reset index to make 'taxonomic_assignment' a column again
    pivoted_df.reset_index(inplace=True)

    # Rename columns to match the target format
    pivoted_df.columns.name = None  # Remove the name of the index column
    pivoted_df.rename(columns={'taxonomic_assignment': 'Taxonomy_ID'}, inplace=True)

    return pivoted_df


# 18. A function to transform the data into OTU table, the common output of QIIME.
def get_table_otu_genus(data, viral_or_not):
    """
    Transform the data into the OTU table format by reordering the columns and only keeping pertinent information.

    :param data: pd.Dataframe
    :param viral_or_not: bool
    :return: pd.Dataframe
    """
    if viral_or_not:
        table_otu_genus_df = data[data["host_source"] != "unknown"].copy()
    else:
        table_otu_genus_df = data.copy()

    # Apply the function to each row to create a new 'taxonomic_assignment' column
    table_otu_genus_df['taxonomic_assignment'] = table_otu_genus_df.apply(concatenate_taxonomic_ranks, axis=1)

    # Pivot into desired format
    pivoted_df = pivot_dataframe(table_otu_genus_df)

    pivoted_df["Taxonomy_ID"] = pivoted_df["Taxonomy_ID"].astype("string")

    # Column that should not be ordered
    keep_col = "Taxonomy_ID"

    # Sort the rest of the column names excluding the first column
    sorted_remaining_columns = sorted(pivoted_df.columns.drop(keep_col))

    # Concatenate the first column name with the sorted list of the remaining column names
    new_column_order = [keep_col] + sorted_remaining_columns

    # Reindex the DataFrame with the new column order
    pivoted_df = pivoted_df.reindex(new_column_order, axis=1)

    # Using DataFrame.insert() to add a column
    pivoted_df.insert(0, "OTU_ID", [num for num in range(len(pivoted_df))], True)

    # Split the 'Taxonomy' column into multiple columns
    taxonomy_split_df = pivoted_df['Taxonomy_ID'].str.split(';', expand=True)

    # Define the names of the new columns. Adjust as needed.
    new_columns = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

    # Find the index of the 'Taxonomy_ID' column
    taxonomy_id_idx = pivoted_df.columns.get_loc('Taxonomy_ID')

    # Insert each new column into the DataFrame right after 'Taxonomy_ID'
    for i, col_name in enumerate(new_columns):
        pivoted_df.insert(taxonomy_id_idx + 1 + i, col_name, taxonomy_split_df[i])

    return pivoted_df


# 19. A function to transform OTU table shaped data into a long format
def get_table_otu_genus_long(dataframe, all_known_pathogenic_genera, ictv_taxonomy, taxonomic_levels):
    # Transforming the dataframe using pandas to achieve the same result as the R code
    df_long = pd.melt(dataframe,
                      id_vars=['OTU_ID', 'Taxonomy_ID', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus',
                               'Species'],
                      var_name='Sample', value_name='Abundance')

    # convert column "Abundance" to numeric
    df_long["Abundance"] = pd.to_numeric(df_long["Abundance"])

    # Calculating total reads per sample
    total_reads_per_sample = df_long.groupby('Sample')['Abundance'].sum().reset_index()
    total_reads_per_sample.columns = ['Sample', 'Total_Reads_Per_Sample']

    # Merging the total reads back to the long dataframe
    df_long = pd.merge(df_long, total_reads_per_sample, on='Sample', how='left')

    # Create new columns by extracting parts of the "Sample" column
    df_long['Sample_Code'] = df_long['Sample'].str[0]  # First character
    df_long['Sample_Number'] = df_long['Sample'].str[1]  # Second character
    df_long['Software'] = df_long['Sample'].str.split('_').str[1]  # Part after the underscore
    df_long['Sample_Code_Number'] = df_long['Sample'].str[0] + df_long['Sample'].str[1]

    # This adds the information on whether or not a virus is pathogentic (check the list for additional info)
    df_long["pathogenic_or_not"] = df_long['Genus'].apply(lambda x: 1 if x[3:] in all_known_pathogenic_genera else 0)

    df_long = create_taxonomy_graph(ictv_taxonomy, df_long, taxonomic_levels, long_format=True)

    return df_long


# 20. The main function, where the magic happens
def main(update_ncbi_taxa, tax_levels_for_classification, phage_hosts):
    # Save the original working directory (script file)
    cwd = os.getcwd()

    # Move back to the main folder where everything is
    os.chdir("../")

    print("Getting raw_reads...")
    # Get the initial raw number of reads for each file
    raw_reads_for_all = read_in_all_raw_reads()
    print("Raw_reads obtained\n")

    # Moving back a directory to store the latest phage taxonomy
    os.chdir(cwd)
    os.chdir("../")

    # To suppress the specific warning for taxid translation
    warnings.filterwarnings("ignore", category=UserWarning, message="taxid \d+ was translated into \d+")
    # To supress the specific warning for data filtering into non_phage, phage, etc
    warnings.filterwarnings("ignore", 'This pattern has match groups')

    # Create a "Results" folder
    if not os.path.exists("script_results"):
        os.mkdir("script_results")

    print("Reading all reports...")
    # Initialize a dictionary that stores the pandas dataframes per sample, per workflow
    # Read all reports and store information in pandas dataframes
    all_samples = read_all_reports_and_store_into_df()
    print("All reports have been read and stored\n")

    print("Making sure no report is missing...")
    # Make sure every software has a df, even if there's no report, it should be the same but empty
    all_samples = check_missing_data_and_fill(all_samples)


    print("\nFiltering out any identifications with 0 reads...")
    # Filter out any identifications with 0 reads
    all_samples = filter_out_ids_with_zero_reads(all_samples)
    print("Filtering finished\n")


    print("Initializing NCBI object to retrieve the taxonomy...")
    # Initialize the ncbi object to retrieve the taxonomy
    ncbi = NCBITaxa()
    # Update NCBI or not, for taxonomy
    if update_ncbi_taxa:
        print("Updating NCBI Taxa...")
        ncbi.update_taxonomy_database()
        print("NCBI Taxa updated")
    print()


    print("Downloading lastest ICTV Taxonomy...")
    # Get the lattest ICTV Virus Taxonomy
    if not os.path.exists("latest_virus_taxonomy.xlsx"):
        download_latest_virus_taxonomy()
    taxonomic_levels = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus']
    ictv_taxonomy = get_ictv_taxonomy_into_dataframe(taxonomic_levels)
    print("Lastest ICTV Taxonomy aquired\n")


    print("Adding the full taxonomy to each sample... (this takes 2-3 minutes)")
    # Add taxonomy full taxonomy to all the samples (from here onward, the colnames are standardized = ["Assignment", "Reads", 'superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    all_samples, known_pathogens = add_taxonomy_to_all_workflows(all_samples, ncbi)
    print("Taxonomy added\n")


    print("Getting all pathogenic genera...")
    all_known_pathogenic_genera = get_all_pathogenic_genera(known_pathogens)
    print("All pathogenic genera obtained\n")


    print("Adding pathogen flag to all identifications...")
    # Add the pathogen flag to all identifications
    all_samples = add_pathogen_flag(all_samples, all_known_pathogenic_genera)
    print("Pathogen flag added to all identifications\n")


    print("Merging everything for R...")
    # Merge everything for R
    all_samples_consensus_result, metadata = merge_everything_for_r(all_samples)
    print("Everything merged for R\n")


    print("Adding host source...")
    # Creating a graph with the taxonomy (kingdom to species -> host source)
    all_samples_consensus_result = create_taxonomy_graph(ictv_taxonomy, all_samples_consensus_result, taxonomic_levels, long_format=False)
    print("Host source added\n")


    print("Adding classification level...")
    all_samples_consensus_result['classification_level'] = all_samples_consensus_result.apply(
        determine_classification_level, args=(tax_levels_for_classification,), axis=1)
    print("Classification level added\n")

    # Remove the "known_pathogen"
    all_samples_consensus_result = all_samples_consensus_result.drop(["known_pathogen"], axis=1)


    print("Creating tables_samples...")
    # Gathering data for the "table_samples" excel file
    tables_samples_df = populate_baseplate_df(all_samples_consensus_result, metadata, raw_reads_for_all, phage_hosts)

    tables_samples_df_for_paper = format_for_paper_table(tables_samples_df, metadata)
    print("tables_samples has been created\n")


    # Filter out everything into sub-dataframes for OTU_genes
    all_viruses = all_samples_consensus_result[all_samples_consensus_result["superkingdom"].isin(["Viruses"])].copy()
    all_but_viruses = all_samples_consensus_result[
        ~all_samples_consensus_result["superkingdom"].isin(["Viruses", "undefined"])].copy()
    undefined_data = all_samples_consensus_result[
        all_samples_consensus_result["superkingdom"].isin(["undefined"])].copy()
    all_non_phage = all_viruses[~all_viruses["host_source"].isin(phage_hosts)].copy()
    all_phage = all_viruses[all_viruses["host_source"].isin(phage_hosts)].copy()

    print("Creating table_otu_genus...")
    all_viruses_otu_table = get_table_otu_genus(all_viruses, viral_or_not=True)
    all_but_viruses_otu_table = get_table_otu_genus(all_but_viruses, viral_or_not=False)
    non_phage_otu_table = get_table_otu_genus(all_non_phage, viral_or_not=True)
    phage_otu_table = get_table_otu_genus(all_phage, viral_or_not=True)
    print("table_otu_genus has been created\n")

    print("Creating table_otu_genus_long...")
    all_viruses_otu_table_long = get_table_otu_genus_long(all_viruses_otu_table, all_known_pathogenic_genera,
                                                          ictv_taxonomy, taxonomic_levels)
    all_but_viruses_otu_table_long = get_table_otu_genus_long(all_but_viruses_otu_table, all_known_pathogenic_genera,
                                                              ictv_taxonomy, taxonomic_levels)
    non_phage_otu_table_long = get_table_otu_genus_long(non_phage_otu_table, all_known_pathogenic_genera, ictv_taxonomy,
                                                        taxonomic_levels)
    non_phage_otu_table_long = non_phage_otu_table_long[non_phage_otu_table_long["host_source"] != "unknown"].copy()
    phage_otu_table_long = get_table_otu_genus_long(phage_otu_table, all_known_pathogenic_genera, ictv_taxonomy,
                                                    taxonomic_levels)
    phage_otu_table_long = phage_otu_table_long[phage_otu_table_long["host_source"] != "unknown"].copy()
    print("table_otu_genus_long has been created\n")

    print("Writting everything:")

    os.chdir(os.path.join(os.getcwd(), "script_results"))
    all_samples_consensus_result.to_excel("taxid_df.xlsx")  # all_samples_for_R
    tables_samples_df.to_excel("table_samples.xlsx")   # fica igual
    tables_samples_df_for_paper.to_excel("table_samples_for_paper.xlsx", index=False)  # ww_samples_df (complete metadata)
    undefined_data.to_excel("taxid_df_undefined.xlsx")  # undefined_data_for_R

    print("Writting all otu_table's...")
    all_viruses_otu_table.to_excel("taxid_df_viruses.xlsx")  # all_viruses_otu_table
    all_but_viruses_otu_table.to_excel("taxid_df_non_viruses.xlsx")  # all_but_viruses_otu_table
    non_phage_otu_table.to_excel("otu_table_non_phage.xlsx")  # non_phage_otu_table
    phage_otu_table.to_excel("otu_table_phage.xlsx")  # phage_otu_table
    print("All otu_table's written.\n")

    print("Writting all otu_table_long's...")
    all_viruses_otu_table_long.to_excel("taxid_df_viruses_long.xlsx")  # all_viruses_otu_table_long
    all_but_viruses_otu_table_long.to_excel("taxid_df_non_viruses_long.xlsx")  # all_but_viruses_otu_table_long
    non_phage_otu_table_long.to_excel("otu_table_non_phage_long.xlsx")  # non_phage_otu_table_long
    phage_otu_table_long.to_excel("otu_table_phage_long.xlsx")  # otu_table_phage_long phage_otu_table_long
    print("All otu_table_long's written.")


if __name__ == "__main__":
    update_ncbi_taxa = False

    tax_levels_for_classification = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus',
                                     'species']  # 'superkingdom'
    phage_hosts = ["bacteria", "archaea"]

    main(update_ncbi_taxa, tax_levels_for_classification, phage_hosts)