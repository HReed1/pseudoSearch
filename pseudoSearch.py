# Written by Harrison Reed

import sys
import re
import argparse
import subprocess as sub
import os
from os.path import exists


def parse_args(parser, argv):
    io = parser.add_argument_group("Input/Output")

    io.add_argument("--query",
                    help="Query sequence. To search for pseudogenes, query sequence should be",
                    required=True)

    io.add_argument("--reference",
                    help="Reference sequences.",
                    required=True)

    io.add_argument("--kmer-length-pct",
                    help="The minimum k-mer length defined as a percentage of the query sequence length.",
                    required=False,
                    default=30)

    parser = parser.parse_args(argv)

    return parser


def supermatcher(query, seq):
    with open('query', 'w') as f:
        f.write('>query\n' + query)

    with open('seq', 'w') as f:
        f.write('>seq\n' + seq)

    sub.call(
        ['supermatcher', '-asequence', 'query', '-bsequence', 'seq', '-gapopen', '10', '-gapextend', '0.5', '-aformat',
         'pair', '-outfile', 'supermatcher-output'])

    first = ''
    second = ''
    previous = ''
    last = ''
    count = 0

    with open('supermatcher-output', 'r') as f:
        for line in f:
            line = line.rstrip()

            if line.startswith("#"):
                line = line.replace("# ", '')
                if line.startswith("Length"):
                    length = re.match(r'^Length..(.*)', line).group(1)

                elif line.startswith("Identity"):
                    identity = re.match(r'^.*\((.*)%\)', line).group(1)

                elif line.startswith("Similarity"):
                    similarity = re.match(r'^.*\((.*)%\)', line).group(1)

                elif line.startswith("Gaps"):
                    gaps = re.match(r'^.*(\d.*)\/.*', line).group(1)

                elif line.startswith("Score"):
                    score = re.match(r'^.*\s(\d+\.\d+)', line).group(1)

            elif not line.startswith("#") and re.match("^.*\d", line):
                if count == 0:
                    first = line
                elif count == 1:
                    second = line
                else:
                    previous = last
                    last = line
                count += 1

    print("Query: " + query)
    print("Seq: " + seq)
    print("First: " + first)
    print("Second: " + second)

    # Ensure supermatcher was able to make a match
    try:
        query_start = re.match(r'^.*\s(\d+)\s', first).group(1)
        seq_start = re.match(r'^.*\s(\d+)\s', second).group(1)
    except AttributeError:
        return None


    # In case the match is only one line
    try:
        query_end = re.match(r'^.*\s(\d+)?', previous).group(1)
        seq_end = re.match(r'^.*\s(\d+)?', last).group(1)
    except AttributeError:
        query_end = re.match(r'^.*\s(\d+)?', first).group(1)
        seq_end = re.match(r'^.*\s(\d+)?', second).group(1)

    return [query_start, query_end, seq_start, seq_end, length, identity, similarity, gaps, score]


def coordinate(match, query, seq, seq_head, cushion):
    # Get coordinates
    g1start = match[0]
    g2start = match[1] - cushion

    # Cast into lists
    # query, seq = list(query), list(seq)
    g1end = g1start + match[2]
    g2end = g2start + match[2] + cushion

    sub_query = query[g1start:g1end]

    # set cushion on chromosome
    sub_seq = seq[abs(g2start):g2end]

    diff = ""

    return [sub_query, sub_seq, g1start, g2start, diff, g1end, g2end, seq_head]


def load_kmers(seq, kmer_length, reference=False, alt_dict=None):
    # Saving each kmer in hash.
    kmer_dict = {}

    # Parsing through every position of the sequence
    for i in range(0, len(seq)):

        # Retrieve kmer from dictionary
        # Sliding window of length kmer_length

        # Prevent sliding window index error
        if i + kmer_length < len(seq):
            kmer = seq[i:i + kmer_length]
        else:
            kmer = seq[i:len(seq)]

        # If loading kmers for reference sequence...
        # ...only keep kmers that match query sequence
        if reference == True:
            if kmer in alt_dict.keys():

                # Check for kmer in dictionary
                if kmer not in kmer_dict:

                    # Add kmer to dictionary as a key
                    # Add the start position of the kmer within the sequence as the value
                    kmer_dict[kmer] = [i]

                # Add subsequent observations to the kmer count
                else:
                    kmer_dict[kmer].append(i)

        # Check for kmer in dictionary
        elif kmer not in kmer_dict:

            # Add kmer to dictionary as a key
            # Add the start position of the kmer within the sequence as the value
            kmer_dict[kmer] = [i]

        # Append multiple start positions of the same kmer
        else:
            kmer_dict[kmer].append(i)

    return kmer_dict


def match_kmers(kmer_dict, kmer_length, seq):
    # Initialize and Append Matches
    with open('matches.tsv', 'a') as f:
        if not os.path.exists('matches.tsv'):
            f.write("Query_Start\tMatch_Start\tKmer_Length\n")

        matches = []

        # Random high number
        highest_match = 999999

        # kmer match counter
        hit_count = 0

        # Match Status Initialization
        m_status = False

        # Parse Chromosome using a Sliding Window with length of kmer_length
        for i in range(0, len(seq)):
            sub = seq[i:i + kmer_length]

            # Capture Match
            if sub in kmer_dict:

                hit_count += 1
                highest_match = kmer_length
                m_status = True  # Boolean to check for matches

                # Printing update
                if hit_count % 10000 == 0:
                    print("HITS: " + str(hit_count))

                # Capture hits
                for hit in kmer_dict[sub]:

                    # hit: match locus for query
                    # i: match query for reference
                    matches.append([hit, i, kmer_length])

                    try:
                        f.write(str(hit) + "\t" + str(i) + "\t" + str(kmer_length) + "\n")
                    except ValueError:
                        print("Value Error - Exiting")
                        print(hit)
                        print(i)
                        print(kmer_length)
                        exit()

    return matches


def collapse_matches(matches):
    matches.sort()

    dele = False

    ## TODO: Replace this with a while Loop
    for i in range(len(matches) - 1):
        try:
            # Find intersecting matches
            while matches[i + 1][1] <= matches[i][1] + matches[i][2]:

                # Combine intersecting matches
                matches[i][2] = matches[i][2] + matches[i + 1][2]

                # After combining recrod 1, delete record 2
                del matches[i + 1]

        # Prevents deletions from ruining indexes
        except IndexError:
            break

    return matches


def pseudoSearch(query, fasta_list_fn, kmer_length):
    """
    Create dictionary of all matching kmers 
    between both genes with the length of 'n'" 
    """

    cushion = 500

    coord_list = []
    chrom = ""
    seq = ""
    current = {}

    print("Retrieve all k-mers of length " + str(kmer_length) + " from query")
    # Retrieve all query k-mers of the length of kmer_length)
    kmer_dict = load_kmers(query, kmer_length)

    print("Length of Query Dictionary: " + str(len(kmer_dict)))

    ref_fasta = {}

    # Stream Reference Sequence
    with open(fasta_list_fn, 'r') as f:
        for line in f:
            line = line.rstrip()

            # Skipping empty lines
            if line == "":
                continue

            # Handle headers
            elif line.startswith(">"):
                current = line

            # Handle sequences
            else:
                ref_fasta[current] = line

                # Skip chunk on the first iteration to initialize current
                if current != {}:

                    # Retrieve completed nucleotide sequence
                    ref = current.get(chrom)[0]

                    if (ref is None):
                        print("Ref is NONE. Printing Line. Exiting")
                        print(line)
                        print(ref)
                        print(current)
                        exit()

                    print("Retrieve all k-mers of length " + str(kmer_length) + " from " + chrom)

                    # Retrieve all reference k-mers of the length of kmer_length)
                    ref_dict = load_kmers(ref, kmer_length, reference=True, alt_dict=kmer_dict)

                    print("Length of Reference Dictionary: " + str(len(ref_dict)))

                    print("Finding k-mer matches between query and reference: " + chrom)
                    matches = match_kmers(kmer_dict, kmer_length, ref)

                    matches = collapse_matches(matches)

                    # Create a header for the file on the first go around
                    if not exists("collapsed_matches.tsv"):
                        with open("collapsed_matches.tsv", 'w') as fn:
                            fn.write("Query Start\tRef Start\tMatch Length\tChromosome\n")

                    with open("collapsed_matches.tsv", 'a') as fn:

                        for match in matches:
                            coord_list.append(coordinate(match, query, ref, chrom, cushion))
                            fn.write(str(match[0]) + "\t" + str(match[1]) + "\t" + str(match[2]) + "\t" + chrom + "\n")

                    current = {chrom: []}

                # Capture Chromosome
                chrom = line.replace('>', '')

                # Reset Current with Latest chromosome
                current = {chrom: []}

            # Handle sequence
            else:
                # Loading in current fasta record
                current[chrom].append(line.upper())


    return coord_list


def main(argv):
    # Parsing Args and Creating Help Dialog
    parser = argparse.ArgumentParser(
        description="Find every k-mer aligned between a query and reference sequence - starting with a k-mer length of 'n' to the length of the query sequence.")
    parser = parse_args(parser, argv)

    query_fn = parser.query
    fasta_list_fn = parser.reference

    query = ""

    # Read gene1 into string variable
    with open(query_fn, 'r') as f:
        for line in f:
            if not line.startswith(">"):
                line = line.rstrip()
                query += line

    kmer_length = int((float(parser.kmer_length_pct) / 100) * len(query))
    match_list = []

    # Get all matching k-mers between reference and query
    match_list = pseudoSearch(query, fasta_list_fn, kmer_length)

    print("All K-mers Generated. Checking for Alignment...")

    with open(query_fn + '_out.tsv', 'w') as fn:
        fn.write(
            'Chromosome\tQuery_Start\tQuery_End\tMatch_Start\tMatch_End\tAlign_Length\tIdentity\tSimilarity\tGaps\tScore\n')
        for record in match_list:
            print("PRINTING MATCH LIST RECORD")
            print(record)
            exit()
            query = record[0]
            gene2 = record[1]
            query_start = record[2]
            g2start = record[3]
            query_end = record[5]
            g2end = record[6]
            seq_head = record[7]

            # H = matrix(query, gene2)
            # ans = traceback(H, gene2)

            print("QUERY LENGTH: " + str(len(query)))
            print("SEQ LENGTH " + str(len(gene2)))

            # query_start, query_end, g2start, g2end, length, identity, similarity, gaps, score = e(query, gene2)
            match_stats = supermatcher(query, gene2)
            # print(match_stats)

            # Skip for empty alignment
            if match_stats is None:
                continue

            # Adjust for relative seq start position
            match_stats[2] = str(g2start + int(match_stats[2]))

            # Adjust for relative seq end position
            match_stats[3] = str(int(match_stats[2]) + int(match_stats[4]))

            # Add chromosome information
            match_stats = [seq_head] + match_stats
            fn.write('\t'.join(match_stats) + '\n')

            # fn.write(seq_head + "\t" + str(g2start) + "\t" + str(end + g2start) + "\t" + str(g2end-g2start) + "\t" + str(score)+ '\n')

    print("Processing Complete!")

if __name__ == '__main__':
    main(sys.argv[1:])
