import requests
import time
from Bio import SeqIO  # Add this import
from io import StringIO


def fetch_interpro_sequences(interpro_id="IPR033966", output_file="interpro_sequences_small_unit.fasta", 
                           page_size=100, query_sequence=None, min_distance=None, max_distance=None):
    """
    Fetch sequences from an InterPro entry using their API and 
    """
    base_url = "https://www.ebi.ac.uk/interpro/api/protein/UniProt/entry/InterPro"
    url = f"{base_url}/{interpro_id}/?page_size={page_size}"
    
    headers = {
        "Accept": "application/json",
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36"
    }
    
    print(f"Fetching sequences for InterPro entry {interpro_id}...")
    
    sequence_count = 0
    filtered_count = 0
    current_url = url
    page = 1
    
    # Open the file in write mode at the start
    with open(output_file, 'w') as fasta_out:
        while current_url:
            print(f"Fetching page {page}...")
            response = requests.get(current_url, headers=headers)
            
            if response.status_code == 200:
                data = response.json()
                if page == 1:
                    total_count = data.get('count', 0)
                    print(f"Total entries found: {total_count}")
                
                # Process results from this page
                for result in data.get('results', []):
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
                            sequence_count += 1
                            
                            # Parse the FASTA sequence to get just the sequence string
                            fasta = list(SeqIO.parse(StringIO(sequence), "fasta"))[0]
                            seq_str = str(fasta.seq)
                            
                            # If we have a query sequence, filter based on Levenshtein distance
                            if query_sequence and max_distance is not None:
                                distance = levenshtein_distance(query_sequence, seq_str)
                                # Check if distance falls within the desired range
                                if (min_distance is None or distance >= min_distance) and distance <= max_distance:
                                    fasta_out.write(sequence)
                                    fasta_out.flush()
                                    filtered_count += 1
                                    print(f"Retrieved and kept sequence {filtered_count} (distance={distance}): {accession}")
                            else:
                                fasta_out.write(sequence)
                                fasta_out.flush()
                                print(f"Retrieved sequence {sequence_count}: {accession}")
                
                # Get URL for next page
                current_url = data.get('next')
                page += 1
                
                # Be nice to the server
                time.sleep(0.5)
            else:
                print(f"Failed to fetch page. Status code: {response.status_code}")
                print(f"URL attempted: {current_url}")
                break
    
    if sequence_count > 0:
        print(f"\nTotal sequences processed: {sequence_count}")
        if query_sequence:
            print(f"Sequences within distance range: {filtered_count}")
        print(f"Sequences saved to: {output_file}")
        return output_file
    else:
        print("No sequences were retrieved")
        return None

if __name__ == "__main__":
    interpro_id = "IPR036385"
    fetch_interpro_sequences(interpro_id, page_size=100)