# pseudoSearch.py

##### System Requirements 
- Python 3.5+
- MacOS or Linux OS
- EMBOSS: http://emboss.sourceforge.net/download/

##### Python Packages
- numpy: https://scipy.org/install.html

##### Program Prequisites
- reference (fasta-list) must be converted to single-line mulit-fasta format: https://www.biostars.org/p/9262/

### Usage
```
usage: pseudoSearch.py [-h] --query QUERY --fasta_list FASTA_LIST
                       [--nmer NMER]

Finding the Longest kmer

optional arguments:
  -h, --help            show this help message and exit

Input/Output:
  --query QUERY         Gene 1
  --fasta_list FASTA_LIST
                        Gene 2
  --nmer NMER           Percent Minimum Perfect Match
```

### Example Usage
```
python pseudoSearch.py --query /Data/Homo_sapiens_NCF1B_sequence.fa --fasta_list /Reference/hg19.fa --nmer 5
```
