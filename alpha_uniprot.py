import pandas as pd
df = pd.read_excel('unique.xlsx')
fold_ids = df['Model names'].tolist()
import requests
annotationsl ={}
for i in fold_ids:
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
          
            return annotations

        except requests.exceptions.HTTPError as http_err:
            print(f"HTTP error occurred: {http_err}")
        except Exception as err:
            print(f"Error occurred: {err}")

        return None

    accession_number = i  # Replace this with the protein accession number you want to query
    annotations = get_annotations(accession_number)
    annotationsl[i]=annotations

import pandas as pd

dfa = pd.DataFrame(annotationsl,index=[1])
dfa = dfa.T
dfa
#convert into excel
dfa.to_excel("Annotations_uniprot_structure_based.xlsx")
#print("Dictionary converted into excel...")
