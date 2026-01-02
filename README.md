# Bioinformatics Algorithms in Python 🧬

This repository contains basic Python implementations of common
bioinformatics algorithms used in genome analysis.

These functions are implemented **from scratch** for learning purposes,
without using external bioinformatics libraries or AI.

## 📌 Algorithms Included

### 1. Pattern Count
Counts how many times a specific DNA pattern appears in a sequence.

```python
PatternCounts("CGCGATACGTTACATACATGATAGACCGCGCGCGATCATATCGCGATTATC","CGCG")
```
5

### 2. Frequent k-mers
Finds the most frequent k-mer in a sequence

```python
Pattern_k_mer("TAAACGTGAGAGAAACGTGCTGATTACACTTGTTCGTGTGGTAT",3)
```
(['GTG'], 4)

### 3. Starting Index Finder
Finds all the starting indexes of a specific pattern inside a sequence 
```python
starting_index_finder("ATGACTTCGCTGTTACGCGC","CGC")
```
7 15 17 

### 4. Clump Finder (Replication of Origin Problem)
Finds all the k-mers that appears atleast t times within any window of length L in a sequence.
This algorithm is aimed to be used to identify potential candidate regions specially for oriC in bacterial genomes.
Note that L is usually 500 and k is usually 9 for bacterial genomes.
```python
replicatio("CGATGATGATGACGATGATG",3,12,2)
```
{'TGA', 'GAT', 'ATG'}

### 5. GC skew Identifier
This algorithm starts from 0 and increase the value by 1 each time it hits a Guanine
and decrease the value by 1 each time it hits a Cytosine resulting in a series of numbers
showing the skew value at each genome location. The index of the base having the minimum
value of GC skew is also presented which helps understand the direction of the replication
and whether we are near the location of Ori or Ter.
```python
GCskew("ATAATCGATGCTGCCGCTATCGTAAGTATTAGCTAGccc")
```
0 0 0 0 0 0 -1 0 0 0 1 0 0 1 0 -1 0 -1 -1 -1 -1 -2 -1 -1 -1 -1 0 0 0 0 0 0 1 0 0 0 1 0 -1 -2

Min GC skew indexes are: 21 39

### 6. Hamming Distance Calculator
This algorithm calculates the hamming distance between two given DNA sequences and returns
a number corresponding to the amount of nucleotides that are different (i.e. mismatch)
```python
Hamming_Distance("ATGCGC","ATGAGT")
```
2

### 7. Approximate pattern
Identify the indexes of a given pattern in a sequence having a number of "mismatch" (variation/mutated nucleotide) with a value of d. For example if seq = ATGC, pattern = AG, and d = 1. The return of the function whould be TG with index of 1 (zero based indexing); telling us that there exist a pattern in the sequence with index = 1, which is different from the given pattern by 1 nucleotide.
```python
ApproxPattern("CTC","CTGTCCCTGTAGCCTCAGGACTTTTGAAATTCTGAGGTCCTTCATAAAATTTGAA",1)
````
indexes are: 
0 2 4 6 13 20 29 31 36 39 40

The total number of fittig nucleotide base is: 11

### 8. Neighbors (mismatches) of a given pattern
Generates a neighborhood ,for a given pattern, consisting of all the neighbors that differ from the original nucleotide by
a "d" number of mismatches.
```python
Neighbors("ACG",1)
```
{'AAG', 'ACA', 'ACC', 'ACG', 'ACT', 'AGG', 'ATG', 'CCG', 'GCG', 'TCG'}

### 9. Reverse DNA strand
A simple function to produce the reverse DNA strand of a given sequence 
````python
ReverseDNA("ATCGCGATGCTAGCTGTATAGCCGGCATTCA")
````
TGAATGCCGGCTATACAGCTAGCATCGCGAT

### 10. THE DNAa box finder
This is one of the main functions in this project which aims at solving the problem of "Frequent Words with Mismatches and Reverse Complements" problem. The algorithm identifies the most frequent k-mer in a sequence with 2 main parametres:

1. Mutations: It asks for a "d" value representing the allowed mismatch number of a k-mer
2. Double-Strandedness: Makes sure the nucleotide pattern that is found appears on both the reverse and forward strands

This algortihm aims to find the potential DNAa boxes within the OriC, considering the fact that DnaA protein can bind to these boxes even if with slight mutations and can bind to both the target sequence or its reverse complement.
````python
kmer_with_d("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, 1)
````
ATGT ACAT



