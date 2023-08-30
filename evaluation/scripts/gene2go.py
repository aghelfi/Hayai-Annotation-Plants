import concurrent.futures

def search_in_file(search_term):
    results = []
    with open("goa_arabidopsis.gaf", 'r') as f:
        for record in f:
            fields = record.split('\t')
            if len(fields) > 10 and search_term.lower() in fields[10].lower():
                results.append((search_term, fields[4]))
    return results

if __name__ == "__main__":
    all_results = []

    # Read search terms from araport11_list.txt
    with open("araport11_list.txt", 'r') as list_file:
        search_terms = [line.strip() for line in list_file]

    # Use a ProcessPoolExecutor for parallel processing
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for matches in executor.map(search_in_file, search_terms):
            all_results.extend(matches)

    # Write all the results to gold_standards.tsv
    with open("gold_standards.tsv", 'a') as out_file:  
        for result in all_results:
            out_file.write("\t".join(result) + "\n")
