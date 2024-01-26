import requests
import xml.etree.ElementTree as ET


def get_smiles_from_chebi(chebi_id):
    # ChEBI URL for the SOAP service to get the complete entity information using the provided endpoint
    url = f'https://www.ebi.ac.uk/webservices/chebi/2.0/test/getCompleteEntity?chebiId={chebi_id}'
    
    # Send a GET request to the ChEBI API
    response = requests.get(url)
    
    # Check if the request was successful
    if response.status_code == 200:
        # Parse the XML response
        root = ET.fromstring(response.content)
        
        # Define the namespace to find the SMILES element
        # Namespaces are usually found in the XML header, but here we extract it from the provided XML structure
        ns = {
            'default': 'http://schemas.xmlsoap.org/soap/envelope/',
            'chebi': 'https://www.ebi.ac.uk/webservices/chebi'
        }
        
        # Extract the SMILES string using XPath
        # Adjust the XPath query to navigate the XML structure properly
        smiles = root.find('.//chebi:smiles', ns)
        
        if smiles is not None:
            return smiles.text
        else:
            print(f"SMILES string not found for ChEBI ID {chebi_id}")
            return None
    
    else:
        print(f"Error: Unable to retrieve data for ChEBI ID {chebi_id}. Status code: {response.status_code}")
        return None



# Example usage
chebi_id      = '57856'
smiles_string = get_smiles_from_chebi(chebi_id)

if smiles_string:
    print(f"The SMILES string for ChEBI ID {chebi_id} is: {smiles_string}")
else:
    print("Could not retrieve the SMILES string.")














