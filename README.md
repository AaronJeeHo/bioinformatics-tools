# bioinformatics-tools
A collection of useful algorithms and scripts written for bioinformatics and genetics research.


&nbsp;
## Requirements
- Python 3.7 +


&nbsp;
## Table of Contents
- [Alignment](#Alignment)
    - [global_align.py](#global_alignpy)
    - [global_align_aa.py](#global_align_aapy)
    - [local_align.py](#local_alignpy)
    - [local_align_aa.py](#local_align_aapy)
- [File Tools](#File-Tools)
    - [parse_fasta.py](#parse_fastapy)
    - [parse_fastq.py](#parse_fastqpy)
- [Hidden Markov Model](#Hidden-Markov-Model)
    - [hmm.py](#hmmpy)
- [Pattern Matching](#Pattern-Matching)
    - [aho-corasick.py](#aho-corasickpy)
- [Seq Tools](#Seq-Tools)
    - [dna_map.py](#dna_mappy)
    - [orf_finder.py](#orf_finderpy)  
    - [gene_finder.py](#gene_finderpy)
    

<br/><br/>
## Alignment
Programs that find optimal local or global alignments for nucleotide and amino acid sequences.

&nbsp;
### global_align.py
Needleman-Wunsch algorithm for aligning two nucleotide sequences.

#### Usage
```shell script
global_align.py [fasta] -m [match] -s [mis] -d [indel]
```
Required Arguments

- fasta: File containing sequence data for alignment, denoted by ">"
- -m / --match: Alignment score per match
- -s / --mismatch: Penalty per mismatch
- -d / --indel: Penalty per insertion or deletion

Optional Arguments

- -a: Output Alignment


&nbsp;
### global_align_aa.py
Needleman-Wunsch algorithm for aligning two amino acid sequences.  
Scored using BLOSUM62 scoring matrix.

#### Usage
```shell script
global_align_aa.py [fasta] -d [indel]
```
Required Arguments

- fasta: File containing sequence data for alignment, denoted by ">"
- -d / --indel: Penalty per insertion or deletion

Optional Arguments

- -a: Output Alignment


&nbsp;
### local_align.py
Smith-Waterman algorithm for aligning two nucleotide sequences.

#### Usage
```shell script
local_align.py [fasta] -m [match] -s [mis] -d [indel]
```
Required Arguments

- fasta: File containing sequence data for alignment, denoted by ">"
- -m / --match: Alignment score per match
- -s / --mismatch: Penalty per mismatch
- -d / --indel: Penalty per insertion or deletion

Optional Arguments

- -a: Output Alignment


&nbsp;
### local_align_aa.py
Smith-Waterman algorithm for aligning two amino acid sequences.  
Scored using PAM250 scoring matrix.

#### Usage
```shell script
local_align_aa.py [fasta] -d [indel]
```
Required Arguments

- fasta: File containing sequence data for alignment, denoted by ">"
- -m / --match: Alignment score per match
- -s / --mismatch: Penalty per mismatch
- -d / --indel: Penalty per insertion or deletion

Optional Arguments

- -a: Output Alignment


<br/><br/>
## File Tools
Programs for parsing sequence data or converting file types

&nbsp;
### parse_fasta.py
Program for parsing sequence data from FASTA files

#### Usage
```shell script
parse_fasta.py [fasta]
```
Required Arguments

- fasta: File containing sequence data for alignment, denoted by ">"

Optional Arguments

- -m / --multi: File contains multiple sequences

#### Callable Functions

- `parse_fasta(fasta)`
    - fasta (str): String denoting file name
    - return (str): Sequence data
    
- `parse_multiseq(fasta)`
    - fasta (str): String denoting file name
    - return (list): List containing sequence data


&nbsp;
### parse_fastq.py
Program for parsing sequence data from FASTQ files, and converting FASTQ to FASTA

#### Usage
```shell script
parse_fastq.py [fastq]
```
Required Arguments

- fastq: FASTQ file storing sequence data and quality scores. 

Optional Arguments

- -f / --fasta: Output sequences in FASTA format

#### Callable Functions

- `parse_fastq(fastq)`
    - fastq (str): String denoting file name
    - return (dict): Dictionary containing sequence data using id's as keys
        - Sequence data denoted as tuple containing (Sequence, Quality Scores)


<br/><br/>
## Hidden Markov Model
Programs for predicting hidden states and probabilities in a sequence using Hidden Markov Models

&nbsp;
### hmm.py
Program that uses Viterbi Algorithm to create HMM's and predict optimal hidden paths and probabilities

#### Usage
```shell script
hmm.py [hmm_file] [parse_order] -p [path] -s [sequence] [action]
```
Required Arguments

- hmm_file: HMM file containing *states, initial state probabilities, transition probabilities,
 symbols emitted, and symbol emission probability* separated by lines starting with "-"
- parse_order: Five character string denoting order in which data is presented in file, eg. "qitse"
    - q: States
    - i: Initial state probabilities
    - t: State transition probabilities
    - s: Symbols emitted
    - e: Emission probabilities

Input Data
- -p / --path: Hidden path HMM will follow
    - Required For: (-d / --dprob) and (-o / --oprob)
    
- -s / --seq: Sequence of symbols HMM will emit
    - Required For: (-v / --viterbi), (-e / --eprob) and (-o / --oprob)
    
Actions
- -v / --viterbi: Viterbi algorithm for finding optimal path given an emitted sequence
- -e / --eprob: Find probability an HMM outputs a given sequence
- -d / --dprob: Find probability an HMM outputs a given hidden path
- -o / --oprob: Find probability an HMM outputs a given sequence following a given path

#### HMM File Format
*Example HMM File*  

Parse Order: "qitse"
```text
A   B
--------
A   B
0.5 0.5
--------
    A   B
A   0.641   0.359
B   0.729   0.271
--------
x   y   z
--------
    x   y   z
A   0.117   0.691   0.192   
B   0.097   0.42    0.483
```

*Input Format*  

- States: Possible Hidden States
    - Format: Single Line Tab-Delimited

- Initial State Probabilities: Probability of starting in state
    - Format: Double Line Tab-Delimited.
        - Top Row: States
        - Bottom Row: Probabilities

- Transition Probabilities: Probability of transitioning between states
    - Format: Tab-Delimited Matrix
        - Row: Current State
        - Column: Transition State

- Symbols: Possible Symbols Emitted 
    - Format: Single Line Tab-Delimited  
    
- Transition Probabilities: Probability of state emitting a given symbol
    - Format: Tab-Delimited Matrix
        - Row: Current State
        - Column: Symbol Emitted



<br/><br/>
## Pattern Matching
Programs for querying reads and fragments against sequences and databases

&nbsp;
### aho-corasick.py
Program that implements an Aho-Corasick Trie to rapidly query a set of reads 
against every position in a sequence database

#### Usage
```shell script
aho-corasick.py -d [database] -q [query] -o [output]
```

Required Arguments
- -d / --database: Single or Multi-Line Text or FASTA file that contains our database sequence
- -q / --query: File containing one read per line to be queried against our database sequence
- -o / --output: Name for output files
    - "output".tsv: Tab-separated file containing matched reads, start index, and end index
    - "output"_stats.tsv: Tab-separated file containing Expected vs Actual matches in total and per read


<br/><br/>
## Seq Tools
Programs that manipulate or extract data from a sequence

&nbsp;
### dna_map.py
Script that returns data or manipulates a nucleotide sequence

#### Usage
```shell script
dna_map.py [sequence] [actions]
```

Required Arguments
- sequence: Sequence to run script on

Actions
- -l / --length: Output sequence length
- -n / --nuc: Output nucleotide counts
- -r / --rna: Convert DNA sequence to RNA
- -c / --comp: Output reverse-complementary strand
- -p / --protein: Convert nucleotide sequence to amino acids*
    - *Requires sequence to be able split into codons (length divisible by 3)


&nbsp;
### orf_finder.py
Program that finds Open Reading Frames in a DNA sequence.  

#### Usage
```shell script
orf_finder.py -f [fasta] -m [min_size]
```

Required Arguments
- -f / --file: FASTA file containing DNA sequence
- -m / --minbp: Minimum base-pair length for ORF's

Optional Arguments
- -n / --nested: Include nested ORF's

Output
- Tab Separated Values containing Start Index, End Index, and Frame


&nbsp;
### gene_finder.py
Program that locates ORF's then predicts Genes in a DNA sequence

#### Usage
```shell script
gene_finder.py -f [fasta] -m [min_size]
```

Required Arguments
- -f / --file: FASTA file containing DNA sequence
- -m / --minbp: Minimum base-pair length for ORF's

Optional Arguments
- -n / --nested: Include nested Genes

Output
- Tab Separated Values containing Gene Label, Forward/Reverse Strand, Frame, Start Index, End Index,
 and Amino Acid Sequence
