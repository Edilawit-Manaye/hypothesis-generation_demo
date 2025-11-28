#!/usr/bin/env python3
import requests
import sys
import time
import gzip

def get_rsid_from_ucsc(chrom, pos, genome="hg19"):
    """Get rsID from UCSC API for a given chromosome and position"""
    url = f"https://api.genome.ucsc.edu/getData/track?genome={genome};track=dbSnp155;chrom=chr{chrom};start={pos-1};end={pos}"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            if "dbSnp155" in data and data["dbSnp155"]:
                return data["dbSnp155"][0]["name"]
    except Exception as e:
        print(f"Error fetching rsID for {chrom}:{pos}: {e}", file=sys.stderr)
    return None

def main():
    if len(sys.argv) != 3:
        print("Usage: python map_snps_to_rsid.py <input_ldscore_file> <output_ldscore_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Count total SNPs for progress tracking
    print("Counting SNPs...")
    with gzip.open(input_file, 'rt') as f:
        total_snps = sum(1 for line in f) - 1  # Subtract header
    
    print(f"Processing {total_snps} SNPs...")
    
    # Process the file
    with gzip.open(input_file, 'rt') as fin, gzip.open(output_file, 'wt') as fout:
        # Write header
        header = fin.readline().strip()
        fout.write(header + "\n")
        
        processed = 0
        found_rsid = 0
        not_found = 0
        
        for line in fin:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                chrom = parts[0]
                snp_id = parts[1]
                bp = parts[2]
                l2_score = parts[3] if len(parts) > 3 else "0"
                
                # Extract position from chr:pos:ref:alt format
                if ":" in snp_id:
                    pos_parts = snp_id.split(":")
                    if len(pos_parts) >= 2:
                        try:
                            pos = int(pos_parts[1])
                            rsid = get_rsid_from_ucsc(chrom, pos)
                            
                            if rsid:
                                fout.write(f"{chrom}\t{rsid}\t{bp}\t{l2_score}\n")
                                found_rsid += 1
                            else:
                                # If no rsID found, keep original but note it
                                fout.write(f"{chrom}\t{snp_id}\t{bp}\t{l2_score}\n")
                                not_found += 1
                                print(f"No rsID found for {snp_id}")
                            
                            processed += 1
                            
                            # Progress update every 100 SNPs
                            if processed % 100 == 0:
                                print(f"Processed {processed}/{total_snps} SNPs. Found: {found_rsid}, Not found: {not_found}")
                            
                            # Be nice to the API - small delay
                            time.sleep(0.1)
                            
                        except ValueError:
                            fout.write(line)
                            processed += 1
                    else:
                        fout.write(line)
                        processed += 1
                else:
                    fout.write(line)
                    processed += 1
            else:
                fout.write(line)
                processed += 1
    
    print(f"Finished! Processed {processed} SNPs.")
    print(f"Found rsIDs for {found_rsid} SNPs.")
    print(f"No rsID found for {not_found} SNPs.")

if __name__ == "__main__":
    main()