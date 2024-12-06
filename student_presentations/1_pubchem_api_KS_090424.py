**Exercise 3a:** 
# Write your code in this cell: (The solution code below will be removed later)

names = [ 'methane', 'ethane', 'propane', 'butane', 'pentane', \
          'hexane','heptane', 'octane', 'nonane', 'decane', \
          'undecane', 'dodecane']
import time

pugrest = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
pugoper = "property/XLogP"
pugout  = "txt"

for i in range(len(names)):    
    
    pugin = "compound/name/" + names[i]    
    
    url = "/".join( [pugrest, pugin, pugoper, pugout] )
    res = requests.get(url)
    print(names[i], ":", res.text)

    if  i % 5 == 4:  
        time.sleep(1)

**Exercise 3b** 
# Write your code in this cell (The solution code below will be removed later)

names = [ 'glycine', 'L-alanine', 'L-serine', 'L-threonine', 'L-cysteine', \
          'L-valine','L-leucine', 'L-isoleucine', 'L-methionine', 'L-proline', \
          'L-phenylalanine', 'L-tyrosine', 'L-tryptophan', 'L-aspartic acid', 'L-glutamic acid', 'L-asparagine', 'L-glutamine', 'L-histidine', 'L-lysine', 'L-arginine']
import time

pugrest = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
pugoper = "property/CanonicalSMILES"
pugout  = "txt"

for i in range(len(names)):    
    
    pugin = "compound/name/" + names[i]   
    
    url = "/".join( [pugrest, pugin, pugoper, pugout] )
    res = requests.get(url)
    print(names[i], ":", res.text)

    if  i % 5 == 4:  
        time.sleep(1)

## 4. Getting multiple molecular properties

pugrest = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
pugin   = "compound/cid/4485,4499,5026,5734,8082"
pugoper = "property/HBondDonorCount,HBondDonorCount,XLogP,TPSA"
pugout  = "csv"

url = "/".join([pugrest, pugin, pugoper, pugout])   # Construct the URL
print(url)
print("-" * 30)   # Print "-" 30 times (to print a line for readability)

res = requests.get(url)
print(res.text)
res.text.rstrip()
PubChem has a standard time limit of **30 seconds per request**.  When you try to retrieve too many properties for too many compounds with a single request, it can take longer than the 30-second limit and a time-out error will be returned.  Therefore, you may need to split the compound list into smaller chunks and process one chunk at a time.
cids = [ 443422,  72301,   8082,    4485,    5353740, 5282230, 5282138, 1547484, 941361, 5734,  \
         5494,    5422,    5417,    5290,    5245,    5026,    4746,    4507,    4499,   4497,  \
         4494,    4474,    4418,    4386,    4009,    4008,    3949,    3926,    3878,   3784,  \
         3698,    3547,    3546,    3336,    3333,    3236,    3076,    2585,    2520,   2351,  \
         2312,    2162,    1236,    1234,    292331,  275182,  235244,  108144,  104972, 77157, \
         5942250, 5311217, 4564402, 4715169, 5311501]
chunk_size = 10

if ( len(cids) % chunk_size == 0 ) : # check if total number of cids is divisible by 10 with no remainder
    num_chunks = len(cids) // chunk_size # sets number of chunks
else : # if divide by 10 results in remainder
    num_chunks = len(cids) // chunk_size + 1 # add one more chunk

print("# Number of CIDs:", len(cids) )
print("# Number of chunks:", num_chunks )
pugrest = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
pugoper = "property/HBondDonorCount,HBondAcceptorCount,XLogP,TPSA"
pugout  = "csv"

csv = ""   #sets a variable called csv to save the comma separated output

for i in range(num_chunks) : # sets number of requests to number of data chunks as determined above
    
    idx1 = chunk_size * i        # sets a variable for a moving window of cids to start in a data chunk
    idx2 = chunk_size * (i + 1)  # sets a variable for a moving window of cids to end ina data chunk

    pugin = "compound/cid/" + ",".join([ str(x) for x in cids[idx1:idx2] ]) # build pug input for chunks of data
    url = "/".join( [pugrest, pugin, pugoper, pugout] )   # Construct the URL
    
    res = requests.get(url)

    if i == 0: # if this is the first request, store result in empty csv variable
        csv = res.text 
    else :          # if this is a subsequent request, add the request to the csv variable adding a new line between chunks
        csv = csv + "\n".join(res.text.split()[1:]) + "\n" 
    
    if i % 5 == 4:  
        time.sleep(1)

print(csv)
print(type(csv))
**Exercise 4a:** Below is the list of CIDs of known antiinflmatory agents (obtained from PubChem via the URL: https://www.ncbi.nlm.nih.gov/pccompound?LinkName=mesh_pccompound&from_uid=68000893).  Use GPT to write a script that downloads the following properties of those compounds in a comma-separated format: Heavy atom count, rotatable bond count, molecular weight, XLogP, hydrogen bond donor count, hydrogen bond acceptor count, TPSA, and isomeric SMILES.

- Split the input CID list into small chunks (with a chunk size of 100 CIDs).
- Process one chunk at a time using a for loop.
- Do not forget to add sleep() to comply the usage policy.

What problems do you encounter in using GPT and how can you solve these?
# Write your code in this cell.

import requests
import time
import csv

def fetch_pubchem_data(cids):
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/property/HeavyAtomCount,RotatableBondCount,MolecularWeight,XLogP,HBondDonorCount,HBondAcceptorCount,TopologicalPolarSurfaceArea,IsomericSMILES/CSV"
    response = requests.get(url.format(','.join(map(str, cids))))
    response.raise_for_status()  # Raise an error if the request failed
    return response.text

def save_to_csv(data, filename="compound_properties.csv"):
    with open(filename, 'w', newline='') as file:
        file.write(data)

def main(cid_list, chunk_size=100):
    all_data = []

    # Split the list of CIDs into chunks
    for i in range(0, len(cid_list), chunk_size):
        chunk = cid_list[i:i + chunk_size]
        
        # Fetch data for the current chunk
        try:
            csv_data = fetch_pubchem_data(chunk)
            all_data.append(csv_data)
        except requests.HTTPError as e:
            print(f"An error occurred: {e}")
            continue
        
        # Sleep to comply with PubChem's rate limits
        time.sleep(1)
    
    # Combine all chunks and save to CSV
    combined_data = "".join(all_data)
    save_to_csv(combined_data)

# Example usage
if __name__ == "__main__":
    # Replace this with your actual list of CIDs
    cids = [ 170924268, 168431238, 166639011, 165360157, 163359777 ]
        
        # Add more CIDs as needed 
    main(cids)


#in above cell, I tested with 5 CIDs prior to running a longer list and received the above error. I input GPT the above error message and received this replacement script:

import requests
import time
import csv

def fetch_pubchem_data(cids):
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/property/HeavyAtomCount,RotatableBondCount,MolecularWeight,XLogP,HBondDonorCount,HBondAcceptorCount,TopologicalPolarSurfaceArea,IsomericSMILES/CSV"
    response = requests.get(url.format(','.join(map(str, cids))))
    response.raise_for_status()  # Raise an error if the request failed
    return response.text

def save_to_csv(data, filename="compound_properties.csv"):
    with open(filename, 'w', newline='') as file:
        file.write(data)

def main(cid_list, chunk_size=50):  # Reduced chunk size to 50
    all_data = []

    # Split the list of CIDs into chunks
    for i in range(0, len(cid_list), chunk_size):
        chunk = cid_list[i:i + chunk_size]
        
        try:
            csv_data = fetch_pubchem_data(chunk)
            all_data.append(csv_data)
        except requests.HTTPError as e:
            print(f"An error occurred for CIDs {chunk}: {e}")
            continue
        
        # Sleep to comply with PubChem's rate limits
        time.sleep(1)
    
    # Combine all chunks and save to CSV
    combined_data = "".join(all_data)
    save_to_csv(combined_data)

# Example usage
if __name__ == "__main__":
    cids = [
        5280931, 3672, 5280450, 445639, 643651, 5280343, 2244, 5280613, 3301, 
        1673, 5280723, 2241, 3671, 5280794, 5280720, 12365139, 10479646, 12365782,
        # Add more CIDs as needed
        
    ]
    main(cids)
import requests
from bs4 import BeautifulSoup
import time
import csv

def scrape_cids(url):
    # Fetch the HTML content from the PubChem URL
    response = requests.get(url)
    response.raise_for_status()
    
    # Parse the HTML content using BeautifulSoup
    soup = BeautifulSoup(response.text, 'html.parser')
    
    # Extract all CIDs from the page (assuming CIDs are in the <dd> tags with class 'item')
    cids = []
    for dd in soup.find_all('dd', class_='item'):
        cid = dd.get_text(strip=True)
        if cid.isdigit():
            cids.append(cid)
    
    return cids

def fetch_pubchem_data(cids):
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/property/HeavyAtomCount,RotatableBondCount,MolecularWeight,XLogP,HBondDonorCount,HBondAcceptorCount,TopologicalPolarSurfaceArea,IsomericSMILES/CSV"
    response = requests.get(url.format(','.join(map(str, cids))))
    response.raise_for_status()  # Raise an error if the request failed
    return response.text

def save_to_csv(data, filename="compound_properties.csv"):
    with open(filename, 'w', newline='') as file:
        file.write(data)

def main():
    # Step 1: Scrape CIDs from the provided URL
    cid_url = "https://www.ncbi.nlm.nih.gov/pccompound?LinkName=mesh_pccompound&from_uid=68000893"
    cids = scrape_cids(cid_url)
    
    print(f"Found {len(cids)} CIDs.")
    
    all_data = []
    errors = []
    chunk_size = 100

    # Step 2: Process the CIDs in chunks
    for i in range(0, len(cids), chunk_size):
        chunk = cids[i:i + chunk_size]
        
        try:
            csv_data = fetch_pubchem_data(chunk)
            all_data.append(csv_data)
        except requests.HTTPError as e:
            print(f"An error occurred for CIDs {chunk}: {e}")
            errors.append(chunk)
        
        # Sleep to comply with PubChem's rate limits
        time.sleep(1)
    
    # Step 3: Combine all chunks and save to CSV
    combined_data = "".join(all_data)
    save_to_csv(combined_data)

    if errors:
        print(f"These CID chunks caused errors and were skipped: {errors}")

if __name__ == "__main__":
    main()

# missing modules for "bs4" and "BeautifulSoup" prevent the script from scraping CIDs from pubchem.
# GPT uses different functions in some areas from functions learned in parts 2 and 3 of this API practice.