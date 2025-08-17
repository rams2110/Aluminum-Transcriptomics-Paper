import subprocess

# Use pip to install openpyxl
subprocess.check_call(["pip", "install", "requests"])
import requests

def get_annotations(accession):
    base_url = "https://alphafold.ebi.ac.uk/api"
    endpoint = f"/uniprot/summary/{accession}.json?key=AIzaSyCeurAJz7ZGjPQUtEaerUkBZ3TaBkXrY94"

    url = base_url + endpoint

    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception if the request was unsuccessful

        data = response.json()
        # Process the data as needed to extract annotations
        annotations = data["structures"]
        annotations = annotations[0]
        annotations = annotations["summary"]
        annotations = annotations["entities"]
        final = str(annotations[0]["description"])
        
        #annotations = data.get('structures', [])
        
        

        return annotations, final

    except requests.exceptions.HTTPError as http_err:
        print(f"HTTP error occurred: {http_err}")
    except Exception as err:
        print(f"Error occurred: {err}")

    return "None"

# Example usage:

accession_number = "I1JGM9"  # Replace this with the protein accession number you want to query
annotation = get_annotations(accession_number)[1]
if annotations:
    for annotation in annotations:
        print(annotation)
else:
    print("Annotations not found for the given accession number.")

# Load the Excel workbook
import openpyxl
workbook = openpyxl.load_workbook('Final.xlsx')

# Select the specific worksheet you want to work with (e.g., the first sheet)
worksheet = workbook.active
structure_description ={}

# Iterate through rows and columns to replace accession numbers with descriptions
for row in worksheet.iter_rows(min_row=2, values_only=True):
    clean_row = [value for value in row if value is not None] 
    description_list= []
    # Assuming accession numbers are in columns B, C, D, etc. (adjust as needed)
    protein_ID =clean_row[0]
    for accession_number in clean_row[1:]:
        description = get_annotations(accession_number)[1]
        description_list.append(description)
        print(accession_number)
        
    structure_description[protein_ID] = description_list

# Close the workbook when done
workbook.close()
import pandas as pd
df = pd.DataFrame(structure_description)
df.to_excel('structure_uniprot_output.xlsx', index=False)

