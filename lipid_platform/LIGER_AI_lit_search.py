import os
import pandas as pd
import requests
import xml.etree.ElementTree as ET
from langchain.chat_models import ChatOpenAI
from langchain.chains import LLMChain
from langchain.prompts import PromptTemplate
from langchain.tools import Tool

file_path = "/scratch/gilbreth/cbeveri/CLAW/lipid_platform/Variable_Storage/folder_path.txt"

# Read the file and store the value in path_variable
with open(file_path, "r") as file:
    path_variable = file.read().strip()

# ================== SETUP API KEYS ==================
os.environ["OPENAI_API_KEY"] = "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"

# ================== PUBMED SEARCH FUNCTIONS ==================
def search_pubmed(query, max_results=5):
    """Searches PubMed for articles matching the query."""
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "pubmed",
        "term": query,
        "retmode": "json",
        "retmax": max_results  # Limit to top results
    }
    response = requests.get(base_url, params=params)
    if response.status_code == 200:
        return response.json()["esearchresult"].get("idlist", [])
    return []

def get_paper_details(paper_ids):
    """
    Fetches details (Title, DOI, and Abstract summary) for given PubMed article IDs using EFetch.
    """
    if not paper_ids:
        return "No papers found."
    
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "pubmed",
        "id": ",".join(paper_ids),
        "retmode": "xml"
    }
    response = requests.get(base_url, params=params)
    if response.status_code != 200:
        return "No details found."
    
    summaries = []
    try:
        root = ET.fromstring(response.content)
        for article in root.findall(".//PubmedArticle"):
            # Get title
            article_title_elem = article.find(".//ArticleTitle")
            title = article_title_elem.text if article_title_elem is not None else "No title available"
            # Get DOI from ELocationID element with attribute EIdType="doi"
            doi = "Not available"
            for eloc in article.findall(".//ELocationID"):
                if eloc.attrib.get("EIdType", "").lower() == "doi":
                    doi = eloc.text
                    break
            # Get abstract (first AbstractText)
            abstract_text = "No abstract available"
            abstract_elem = article.find(".//Abstract/AbstractText")
            if abstract_elem is not None and abstract_elem.text:
                abstract_text = abstract_elem.text
            summaries.append(f"Title: {title}\nDOI: {doi}\nSummary: {abstract_text}\n")
    except Exception as e:
        return f"Error parsing paper details: {e}"
    
    return "\n".join(summaries)

def fetch_papers(query):
    """Fetches PubMed papers based on the query."""
    paper_ids = search_pubmed(query)
    return get_paper_details(paper_ids)

# ================== LLM SETUP ==================
llm = ChatOpenAI(temperature=0.7, model_name="gpt-4")

# LLMChain for filtering PubMed paper details.
filter_prompt = PromptTemplate(
    input_variables=["question", "papers"],
    template=("Below are details for papers retrieved from PubMed for the query: {question}.\n"
              "Paper Details:\n{papers}\n\n"
              "Filter out any non-relevant results and provide for each relevant paper its title, DOI, "
              "and a brief summary. If no papers are relevant, state 'No relevant papers found.'")
)
llm_chain = LLMChain(llm=llm, prompt=filter_prompt)

def process_query(query):
    """For a given query, search PubMed, get paper details, filter with LLM, and return the result."""
    print(f"Searching PubMed for: {query}\n")
    papers = fetch_papers(query)
    if papers.strip() == "No papers found.":
        return f"Query: {query}\nNo papers found.\n{'-'*40}\n"
    else:
        filtered_summary = llm_chain.run({"question": query, "papers": papers})
        result = f"Query: {query}\nOriginal Paper Details:\n{papers}\nFiltered Summary:\n{filtered_summary}\n{'-'*40}\n"
        return result

# ================== INFER DISEASE FROM FILENAME ==================
def infer_disease(filename):
    """
    Uses an LLM to infer the disease context from the given file name based on biomedical domain knowledge.
    The LLM is prompted to return only the disease name. If no clear association is found, it should return 'unknown disease'.
    """
    prompt = PromptTemplate(
        input_variables=["filename"],
        template=("Based on the file name '{filename}', using your biomedical domain knowledge, "
                  "determine the most likely disease or condition associated with this file. "
                  "Return only the name of the disease. If none is clearly associated, return 'unknown disease'.")
    )
    inference_chain = LLMChain(llm=llm, prompt=prompt)
    disease = inference_chain.run({"filename": filename}).strip()
    return disease

# ================== FILE PROCESSING PIPELINE ==================
def process_file(filepath, output_folder):
    """
    Processes a given CSV file:
      - Reads the file and finds significant rows (P-Value < 0.05; or the row with the lowest p-value if none)
      - Uses the 'Gene' column (note the capital 'G') for forward reaction queries.
      - Uses T-Statistic to determine direction:
            * T-Statistic > 0: forward genes are upregulated and the inverse reaction genes should be searched as downregulated.
            * T-Statistic < 0: forward genes are downregulated and the inverse reaction genes should be searched as upregulated.
      - For each row, also checks the lipid pair in "Lipid 1" and "Lipid 2".
            * The inverse reaction is determined by switching the lipid columns.
      - The disease context is determined by inferring from the filename using an LLM.
      - For every query, it searches PubMed and saves the details (title, DOI, and summary) into a text file.
    """
    try:
        df = pd.read_csv(filepath)
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return
    
    # Check required columns
    required_cols = ["P-Value", "Gene", "Lipid 1", "Lipid 2", "T-Statistic"]
    for col in required_cols:
        if col not in df.columns:
            print(f"File {filepath} is missing required column: {col}")
            return
    
    # Filter rows where P-Value is below 0.05; if none, choose the row with the smallest p-value
    sig_df = df[df["P-Value"] < 0.05]
    if sig_df.empty:
        sig_df = df.loc[[df["P-Value"].idxmin()]]
    
    # Use the file name and LLM-based domain knowledge to infer the disease context
    filename = os.path.basename(filepath)
    disease = infer_disease(filename)
    print(f"Inferred disease for {filename}: {disease}")
    
    query_results = ""
    # Process each significant row individually
    for index, row in sig_df.iterrows():
        t_stat = row["T-Statistic"]
        if t_stat >= 0:
            forward_direction = "increased"
            inverse_direction = "decreased"
        else:
            forward_direction = "decreased"
            inverse_direction = "increased"
        
        # Process forward reaction genes from the "Gene" column.
        forward_genes = set()
        for gene_entry in str(row["Gene"]).split(","):
            gene = gene_entry.strip()
            if gene:
                forward_genes.add(gene)
        
        for gene in forward_genes:
            query = f"{forward_direction} {gene} in {disease}"
            query_results += process_query(query) + "\n"
        
        # Process inverse reaction by switching Lipid 1 and Lipid 2.
        lipid1 = str(row["Lipid 1"]).strip().upper()
        lipid2 = str(row["Lipid 2"]).strip().upper()
        inverse_lipid1 = lipid2  # switched
        inverse_lipid2 = lipid1  # switched
        # Look up rows with the reversed lipid pair in the entire dataframe
        inverse_rows = df[
            (df["Lipid 1"].astype(str).str.strip().str.upper() == inverse_lipid1) &
            (df["Lipid 2"].astype(str).str.strip().str.upper() == inverse_lipid2)
        ]
        if not inverse_rows.empty:
            inverse_genes = set()
            for _, inv_row in inverse_rows.iterrows():
                for gene_entry in str(inv_row["Gene"]).split(","):
                    gene = gene_entry.strip()
                    if gene:
                        inverse_genes.add(gene)
            for gene in inverse_genes:
                query = f"{inverse_direction} {gene} in {disease}"
                query_results += process_query(query) + "\n"
    
    # Save the queries and their results to an output text file
    output_filename = os.path.splitext(filename)[0] + "_query_results.txt"
    output_path = os.path.join(output_folder, output_filename)
    with open(output_path, "w") as f:
        f.write(query_results)
    print(f"Saved query results for {filepath} to {output_path}")

def main():
    # Define input and output folders
    input_folder = path_variable + "sum"
    output_folder = path_variable + "Query_Results"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    # Process all CSV files in the input folder
    for file in os.listdir(input_folder):
        if file.endswith(".csv"):
            file_path = os.path.join(input_folder, file)
            print(f"Processing file: {file_path}")
            process_file(file_path, output_folder)

if __name__ == "__main__":
    main()



# import os
# import pandas as pd
# import requests
# import xml.etree.ElementTree as ET
# from langchain.chat_models import ChatOpenAI
# from langchain.chains import LLMChain
# from langchain.prompts import PromptTemplate
# from langchain.tools import Tool


# file_path = "/scratch/gilbreth/cbeveri/CLAW/lipid_platform/Variable_Storage/folder_path.txt"

# # Read the file and store the value in path_variable
# with open(file_path, "r") as file:
#     path_variable = file.read().strip()

# # ================== SETUP API KEYS ==================
# os.environ["OPENAI_API_KEY"] = "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"

# # ================== PUBMED SEARCH FUNCTIONS ==================
# def search_pubmed(query, max_results=5):
#     """Searches PubMed for articles matching the query."""
#     base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
#     params = {
#         "db": "pubmed",
#         "term": query,
#         "retmode": "json",
#         "retmax": max_results  # Limit to top results
#     }
#     response = requests.get(base_url, params=params)
#     if response.status_code == 200:
#         return response.json()["esearchresult"].get("idlist", [])
#     return []

# def get_paper_details(paper_ids):
#     """
#     Fetches details (Title, DOI, and Abstract summary) for given PubMed article IDs using EFetch.
#     """
#     if not paper_ids:
#         return "No papers found."
    
#     base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
#     params = {
#         "db": "pubmed",
#         "id": ",".join(paper_ids),
#         "retmode": "xml"
#     }
#     response = requests.get(base_url, params=params)
#     if response.status_code != 200:
#         return "No details found."
    
#     summaries = []
#     try:
#         root = ET.fromstring(response.content)
#         for article in root.findall(".//PubmedArticle"):
#             # Get title
#             article_title_elem = article.find(".//ArticleTitle")
#             title = article_title_elem.text if article_title_elem is not None else "No title available"
#             # Get DOI from ELocationID element with attribute EIdType="doi"
#             doi = "Not available"
#             for eloc in article.findall(".//ELocationID"):
#                 if eloc.attrib.get("EIdType", "").lower() == "doi":
#                     doi = eloc.text
#                     break
#             # Get abstract (first AbstractText)
#             abstract_text = "No abstract available"
#             abstract_elem = article.find(".//Abstract/AbstractText")
#             if abstract_elem is not None and abstract_elem.text:
#                 abstract_text = abstract_elem.text
#             summaries.append(f"Title: {title}\nDOI: {doi}\nSummary: {abstract_text}\n")
#     except Exception as e:
#         return f"Error parsing paper details: {e}"
    
#     return "\n".join(summaries)

# def fetch_papers(query):
#     """Fetches PubMed papers based on the query."""
#     paper_ids = search_pubmed(query)
#     return get_paper_details(paper_ids)

# # ================== LLM SETUP ==================
# llm = ChatOpenAI(temperature=0.7, model_name="gpt-4")

# # LLMChain for filtering PubMed paper details.
# filter_prompt = PromptTemplate(
#     input_variables=["question", "papers"],
#     template=("Below are details for papers retrieved from PubMed for the query: {question}.\n"
#               "Paper Details:\n{papers}\n\n"
#               "Filter out any non-relevant results and provide for each relevant paper its title, DOI, "
#               "and a brief summary. If no papers are relevant, state 'No relevant papers found.'")
# )
# llm_chain = LLMChain(llm=llm, prompt=filter_prompt)

# def process_query(query):
#     """For a given query, search PubMed, get paper details, filter with LLM, and return the result."""
#     print(f"Searching PubMed for: {query}\n")
#     papers = fetch_papers(query)
#     if papers.strip() == "No papers found.":
#         return f"Query: {query}\nNo papers found.\n{'-'*40}\n"
#     else:
#         filtered_summary = llm_chain.run({"question": query, "papers": papers})
#         result = f"Query: {query}\nOriginal Paper Details:\n{papers}\nFiltered Summary:\n{filtered_summary}\n{'-'*40}\n"
#         return result

# # ================== INFER DISEASE FROM FILENAME ==================
# def infer_disease(filename):
#     """
#     Uses an LLM to infer the disease context from the given file name based on biomedical domain knowledge.
#     The LLM is prompted to return only the disease name. If no clear association is found, it should return 'unknown disease'.
#     """
#     prompt = PromptTemplate(
#         input_variables=["filename"],
#         template=("Based on the file name '{filename}', using your biomedical domain knowledge, "
#                   "determine the most likely disease or condition associated with this file. "
#                   "Return only the name of the disease. If none is clearly associated, return 'unknown disease'.")
#     )
#     inference_chain = LLMChain(llm=llm, prompt=prompt)
#     disease = inference_chain.run({"filename": filename}).strip()
#     return disease

# # ================== FILE PROCESSING PIPELINE ==================
# def process_file(filepath, output_folder):
#     """
#     Processes a given CSV file:
#       - Reads the file and finds significant rows (P-Value < 0.05; or the row with the lowest p-value if none)
#       - Uses the 'Gene' column (note the capital 'G') for forward reaction queries.
#       - Uses T-Statistic to determine direction:
#             * T-Statistic > 0: forward genes are upregulated and the inverse reaction genes should be searched as downregulated.
#             * T-Statistic < 0: forward genes are downregulated and the inverse reaction genes should be searched as upregulated.
#       - For each row, also checks the lipid pair in "Lipid 1" and "Lipid 2".
#             * If they are 'LPI' and 'PI' respectively (case insensitive), then the inverse reaction is assumed to be PI → LPI,
#               and the script searches separately for the genes PLA2G4D and PLA2G4E.
#       - The disease context is determined by inferring from the filename using an LLM.
#       - For every query, it searches PubMed and saves the details (title, DOI, and summary) into a text file.
#     """
#     try:
#         df = pd.read_csv(filepath)
#     except Exception as e:
#         print(f"Error reading {filepath}: {e}")
#         return
    
#     # Check required columns
#     required_cols = ["P-Value", "Gene", "Lipid 1", "Lipid 2", "T-Statistic"]
#     for col in required_cols:
#         if col not in df.columns:
#             print(f"File {filepath} is missing required column: {col}")
#             return
    
#     # Filter rows where P-Value is below 0.05; if none, choose the row with the smallest p-value
#     sig_df = df[df["P-Value"] < 0.05]
#     if sig_df.empty:
#         sig_df = df.loc[[df["P-Value"].idxmin()]]
    
#     # Use the file name and LLM-based domain knowledge to infer the disease context
#     filename = os.path.basename(filepath)
#     disease = infer_disease(filename)
#     print(f"Inferred disease for {filename}: {disease}")
    
#     query_results = ""
#     # Process each significant row individually
#     for index, row in sig_df.iterrows():
#         t_stat = row["T-Statistic"]
#         if t_stat >= 0:
#             forward_direction = "increased"
#             inverse_direction = "decreased"
#         else:
#             forward_direction = "decreased"
#             inverse_direction = "increased"
        
#         # Process forward reaction genes from the "Gene" column.
#         forward_genes = set()
#         for gene_entry in str(row["Gene"]).split(","):
#             gene = gene_entry.strip()
#             if gene:
#                 forward_genes.add(gene)
        
#         for gene in forward_genes:
#             query = f"{forward_direction} {gene} in {disease}"
#             query_results += process_query(query) + "\n"
        
#         # Process inverse reaction based on lipid pair
#         lipid1 = str(row["Lipid 1"]).strip().upper()
#         lipid2 = str(row["Lipid 2"]).strip().upper()
#         # Example: if forward reaction is LPI -> PI then inverse reaction is PI -> LPI
#         if lipid1 == "LPI" and lipid2 == "PI":
#             inverse_genes = ["PLA2G4D", "PLA2G4E"]
#             for gene in inverse_genes:
#                 query = f"{inverse_direction} {gene} in {disease}"
#                 query_results += process_query(query) + "\n"
    
#     # Save the queries and their results to an output text file
#     output_filename = os.path.splitext(filename)[0] + "_query_results.txt"
#     output_path = os.path.join(output_folder, output_filename)
#     with open(output_path, "w") as f:
#         f.write(query_results)
#     print(f"Saved query results for {filepath} to {output_path}")

# def main():
#     # Define input and output folders
#     input_folder = path_variable+"sum"
#     output_folder = path_variable+"Query_Results"
#     if not os.path.exists(output_folder):
#         os.makedirs(output_folder)
    
#     # Process all CSV files in the input folder
#     for file in os.listdir(input_folder):
#         if file.endswith(".csv"):
#             file_path = os.path.join(input_folder, file)
#             print(f"Processing file: {file_path}")
#             process_file(file_path, output_folder)

# if __name__ == "__main__":
#     main()
