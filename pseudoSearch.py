# Written by Harrison Reed

import sys
import re
import argparse
import itertools
import numpy as np
import subprocess as sub
import os

def parse_args(parser, argv):

    io = parser.add_argument_group("Input/Output")

    io.add_argument("--query",
                    help = "Gene 1",
                    required = True)
    io.add_argument("--fasta_list",
                    help = "Gene 2",
                    required = True)

    io.add_argument("--nmer",
                    help = "Percent Minimum Perfect Match",
                    required = False,
                    default = 30)

    parser = parser.parse_args(argv)

    return parser


def supermatcher(query, seq):

    with open('query', 'w') as f:
        f.write('>query\n' + query)

    with open('seq', 'w') as f:
        f.write('>seq\n' + seq)

    sub.call(['supermatcher', '-asequence', 'query', '-bsequence', 'seq', '-gapopen',  '10', '-gapextend',  '0.5', '-aformat',  'pair',  '-outfile',  'test'])
    
    first = ''
    second = ''
    previous = ''
    last = ''
    count = 0

    with open('test', 'r') as f:
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
            
    
    query_start = re.match(r'^.*\s(\d+)\s', first).group(1)
    seq_start = re.match(r'^.*\s(\d+)\s', second).group(1)

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
    #query, seq = list(query), list(seq)
    g1end = g1start + match[2]
    g2end = g2start + match[2] + cushion

    sub_query = query[g1start:g1end]

    # set cushion on chromosome
    sub_seq = seq[abs(g2start):g2end]
    
    diff = ""

    return [sub_query, sub_seq, g1start, g2start, diff, g1end, g2end, seq_head]


def load_kmers(query, lkmer, seq):
    
    matches = []
    s_nmer = 0

    # matchMer = {}

    with open('matches.tsv', 'a') as f:
        if not os.path.exists('matches.tsv'):
            f.write("Query_Start\tMatch_Start\tKmer_Length\n")
        # Saving each kmer in hash.
        # for i in range(0, len(query), nmer):
        # for nmer in range(20, len(query)):
        for nmer in range(lkmer, len(query)):
            kmer_dict = {}
            
            if nmer % 100 == 0:
                print("Loading nmer of " + str(nmer))

            for i in range(0, len(query)):

                # Left over nmers
                if i + nmer < len(query): 
                    kmer = query[i:i+nmer]

                else:
                    s_nmer = len(query) - i
                    kmer = query[i:i+s_nmer]
                if kmer not in kmer_dict:
                    kmer_dict[kmer] = [i]
                else:
                    kmer_dict[kmer].append(i)
                    
            highest_match = 999999
            hit_count = 0
            # Match Status Initialization
            m_status = False

                
            # Parse Chromosome for each kmer length of nmer
            for i in range(0, len(seq), nmer):
                sub = seq[i:i+nmer]

                 # Capture Match
                if sub in kmer_dict:

                    hit_count += 1
                    highest_match = nmer
                    m_status = True

                    if hit_count % 10000 == 0:
                        print("HITS: " + str(hit_count))

                    for hit in kmer_dict[sub]:
                        matches.append([hit, i, nmer])
                        try:
                            f.write(str(hit) + "\t" + str(i) + "\t" + str(nmer) + "\n")
                        except ValueError:
                            print(hit)
                            print(i)
                            print(nmer)
                            exit()

            # If no matches are found, no need to extend
            if m_status is False:
                break

    return matches


def collapse_matches(matches):
    
    matches.sort()

    dele = False
    for i in range(len(matches) -1 ):
        try:
            # Find intersecting matches
            while matches[i+1][1] <= matches[i][1] + matches[i][2]:

                # Combine intersecting matches
                matches[i][2] = matches[i][2] + matches[i+1][2]
            
                # After combining recrod 1, delete record 2
                del matches[i+1]

        # Prevents deletions from ruining indexes
        except IndexError:
            break

    return matches


def pseudoSearch(query, fasta_list_fn, nmer):

    
    """
    Create dictionary of all matching kmers 
    between both genes with the length of 'n'" 
    """
    
    cushion = 500

    coord_list = []
    chrom = ""
    seq = ""

    # Stream Reference Sequence
    with open(fasta_list_fn, 'r') as f:
        for line in f:
            line = line.rstrip()

            if line.startswith(">"):

                chrom = line.replace('>','')
                                            
            else:
                seq = line.upper()

                matches = load_kmers(query, nmer, seq)
    
                print("ALL LOADED " + str(nmer) + " on " + chrom)

                matches = collapse_matches(matches)
                
                with open("collapsed_matches.tsv", 'a') as fn:

                    for match in matches:
                        coord_list.append(coordinate(match, query, seq, chrom, cushion))
                        fn.write(str(match[0]) + "\t" + str(match[1]) + "\t" + str(match[2]) + "\n")

    return coord_list

def main(argv):


    # Parsing Args and Creating Help Dialog
    parser = argparse.ArgumentParser(description="Finding the Longest kmer")
    parser = parse_args(parser, argv)

    query_fn = parser.query
    fasta_list_fn = parser.fasta_list

    query = ""

    # Read gene1 into string variable
    with open(query_fn, 'r') as f:
        for line in f:
            if not line.startswith(">"):
                line = line.rstrip()
                query += line

    nmer = int((float(parser.nmer)/100) * len(query))
    match_list = []

    match_list = pseudoSearch(query, fasta_list_fn, nmer)

        
    with open(query_fn + '_out.tsv', 'w') as fn:
        fn.write('Chromosome\tQuery_Start\tQuery_End\tMatch_Start\tMatch_End\tAlign_Length\tIdentity\tSimilarity\tGaps\tScore\n')
        for record in match_list:
            #print(match_list)
            query = record[0]
            gene2 = record[1]
            query_start = record[2]
            g2start = record[3]
            query_end = record[5]
            g2end = record[6]
            seq_head = record[7]


            #H = matrix(query, gene2)
            #ans = traceback(H, gene2)

            print("ALIRNGMENT")
            print("QUERY LENGTH: " + str(len(query)))
            print("SEQ LENGTH " + str(len(gene2)))

            # query_start, query_end, g2start, g2end, length, identity, similarity, gaps, score = supermatcher(query, gene2)
            match_stats = supermatcher(query, gene2)
            #print(match_stats)

            # Adjust for relative seq start position
            match_stats[2] = str(g2start + int(match_stats[2]))

            # Adjust for relative seq end position
            match_stats[3] = str(int(match_stats[2])+ int(match_stats[4]))

            # Add chromosome information
            match_stats = [seq_head] + match_stats
            fn.write('\t'.join(match_stats) + '\n')

            #fn.write(seq_head + "\t" + str(g2start) + "\t" + str(end + g2start) + "\t" + str(g2end-g2start) + "\t" + str(score)+ '\n')


if __name__ == '__main__':
    main(sys.argv[1:])
