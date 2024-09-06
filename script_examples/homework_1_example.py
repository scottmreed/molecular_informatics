import requests
import time

def fetch_properties(cids):
    """Fetches properties for a list of CIDs from PUGREST."""
    base_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/property/HeavyAtomCount,RotatableBondCount,MolecularWeight,XLogP,HBondDonorCount,HBondAcceptorCount,TPSA,IsomericSMILES/JSON'
    url = base_url.format(','.join(map(str, cids)))
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    else:
        print(f"Error fetching data: {response.status_code}")
        return None

def process_cid_list(cids, chunk_size=100):
    """Processes CIDs in chunks and retrieves properties."""
    results = []
    for i in range(0, len(cids), chunk_size):
        chunk = cids[i:i + chunk_size]
        print(f"Processing chunk: {chunk}")
        properties = fetch_properties(chunk)
        if properties:
            results.extend(properties['PropertyTable']['Properties'])
        time.sleep(1)  # Sleep for 1 second to comply with usage policy
    return results

def main():
    cids = [
        443422, 72301, 8082, 4485, 5353740, 5282230, 5282138, 1547484, 941361, 5734, 5494, 5422, 5417, 5290, 
        5245, 5026, 4746, 4507, 4499, 4497, 4494, 4474, 4418, 4386, 4009, 4008, 3949, 3926, 3878, 3784, 3698, 
        3547, 3546, 3336, 3333, 3236, 3076, 2585, 2520, 2351, 2312, 2162, 1236, 1234, 292331, 275182, 235244, 
        108144, 104972, 77157, 5942250, 5311217, 4564402, 4715169, 5311501
    ]
    chunk_size = 100
    results = process_cid_list(cids, chunk_size)
    
    # Print or save results
    print("List of properties:")
    for result in results:
        print(result)

if __name__ == "__main__":
    main()
