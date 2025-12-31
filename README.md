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
