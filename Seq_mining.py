import requests
import time


def fetch_interpro_sequences(interpro_id="IPR033966", output_file="interpro_sequences_all.fasta", page_size=100):
    """
    Fetch sequences from an InterPro entry using their API
    """
    base_url = "https://www.ebi.ac.uk/interpro/api/protein/UniProt/entry/InterPro"
    url = f"{base_url}/{interpro_id}/?page_size={page_size}"
    
    headers = {
        "Accept": "application/json",
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36"
    }
    
    print(f"Fetching sequences for InterPro entry {interpro_id}...")
    
    all_sequences = []
    current_url = url
    page = 1
    max_sequences = 300000  # Limit to first 1000 sequences for practical purposes
    
    while current_url and len(all_sequences) < max_sequences:
        print(f"Fetching page {page}...")
        response = requests.get(current_url, headers=headers)
        
        if response.status_code == 200:
            data = response.json()
            total_count = data.get('count', 0)
            if page == 1:
                print(f"Total entries found: {total_count}")
                print(f"Will fetch first {max_sequences} sequences")
            
            # Process results from this page
            for result in data.get('results', []):
                if len(all_sequences) >= max_sequences:
                    break
                    
                metadata = result.get('metadata', {})
                if metadata:
                    accession = metadata.get('accession', 'Unknown')
                    name = metadata.get('name', '')
                    organism = metadata.get('source_organism', {}).get('scientificName', '')
                    
                    # Get sequence using UniProt API
                    seq_url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
                    seq_response = requests.get(seq_url)
                    if seq_response.status_code == 200:
                        sequence = seq_response.text
                        all_sequences.append(sequence)
                        print(f"Retrieved sequence {len(all_sequences)}/{max_sequences}: {accession} - {name} ({organism})")
            
            # Get URL for next page
            current_url = data.get('next')
            page += 1
            
            # Be nice to the server
            time.sleep(0.5)  # Reduced sleep time but still being respectful
        else:
            print(f"Failed to fetch page. Status code: {response.status_code}")
            print(f"URL attempted: {current_url}")
            break
    
    if all_sequences:
        # Write all sequences to file
        with open(output_file, 'w') as f:
            f.writelines(all_sequences)
        
        print(f"\nTotal sequences saved: {len(all_sequences)}")
        print(f"Sequences saved to: {output_file}")
        return output_file
    else:
        print("No sequences were retrieved")
        return None

if __name__ == "__main__":
    interpro_id = "IPR033966"
    fetch_interpro_sequences(interpro_id, page_size=100)