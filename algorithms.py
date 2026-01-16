# NOTE: Please enter your integres using double quotes annotation ("ATAATCGGCTAGCGT", "ACGC" etc.)

# Given a specific pattern like "GTC" Count how many of the pattern exists in the sequence 
def PatternCounts(seq, pattern):
    seq = str(seq)
    pattern = str(pattern)
    count = 0
    it = len(pattern)

    for i in range(0, len(seq) - it + 1):
        if seq[i:i+it] == pattern:
            count += 1
    return count

# Identify and return the maximum count of the most frequent k-mers in a list 
def max_frequent(slices):
        max_count = 0
        frequent = []

        for kmer in set(slices):
            count = slices.count(kmer)

            if count > max_count:
                max_count = count
                frequent = [kmer]
            elif count == max_count:
                frequent.append(kmer)

        return frequent, max_count

# Given a k value find the most frequent k-mer in a sequence
def Pattern_k_mer(seq, k):
    slices = [seq[i:i+k] for i in range(0, len(seq) - k + 1)]

    def max_frequent(slices):
        max_count = 0
        frequent = []

        for kmer in set(slices):
            count = slices.count(kmer)

            if count > max_count:
                max_count = count
                frequent = [kmer]
            elif count == max_count:
                frequent.append(kmer)

        return frequent, max_count

    return max_frequent(slices)

# Find all the starting indexes of a spesific pattern inside a sequence 
index_list = []
def starting_index_finder(seq, pattern):
    splices = [seq[i:i+len(pattern)] for i in range(0,len(seq)-len(pattern)+1)]
    index_list = ([i for i, n in enumerate(splices) if n == pattern])
    for i in index_list:
        print(i, sep=" ",end = " ")

# Find all the k-mers that appears atleast t times within any window of length L in the sequence
def replicatio(seq, k, L, t):
    count = {}
    frequency = set()
    
    firstwindow = seq[0:L]
    for i in range(len(firstwindow)-k+1):
        pattern = firstwindow[i:i+k]

        if pattern in count:
            count[pattern] += 1
        else:
            count[pattern] = 1
        if count[pattern] >= t:
            frequency.add(pattern)

    for i in range(1,len(seq) - L + 1):
        left_leaving = seq[i-1 : i-1+k]
        count[left_leaving] -= 1

        right_entering = seq[i+L-k : i+L]
        if right_entering in count:
            count[right_entering] += 1
        else:
            count[right_entering] = 1
        if count[right_entering] >= t:
            frequency.add(right_entering)
    return(frequency)
    
# Calcualte GC skew and find the index of the minimum value for the GC skew     
def GCskew(seq):
    skew = 0
    GC_skew_list = [0]

    for base in seq:
        if base == "G":
            skew += 1
        elif base == "C":
            skew -= 1
        GC_skew_list.append(skew)
    print("GC skew is: ")
    print(*GC_skew_list)
    min_skew = min(GC_skew_list)
    indexes = [i for i, value in enumerate(GC_skew_list) if value == min_skew]
    print("Min GC skew indexes are: ")
    print(*indexes)

# Calculate the Hamming Distance between two sequences
def Hamming_Distance(A, B):
  count = 0
  for i in range(len(A)):
    if A[i] != B[i]:
        count += 1

# Identify the indexes of a given pattern in a sequence for each "mismatch" having a variation/mutated nucleotide with a d parametre
# for example if seq = ATGC, pattern = AG, and d = 1
# the return of the function whould be TG with index of 1 (zero based indexing); telling us that there exist a pattern in the sequence 
# which is different from the given pattern by 1 nucleotide 
def ApproxPattern(pattern, seq, d):
  index = []
  k = len(pattern)
  for i in range(len(seq) - k + 1):
    splices = seq[i:i+k]
    if Hamming_Distance(pattern, splices) <= d:
      index.append(i)
  print("indexes are: ")
  print(*index)
  occurances = len(index)
  print("The total number of fittig nucleotide base is: {}".format(occurances))
  return

# Identify the "Neighbors" for each given "pattern" and a "d" value 
# construct a neighborhood having all the neighbors 
def Neighbors(pattern, d):
  if d == 0:
    return {pattern}
  if len(pattern) == 1:
    return {"A", "C", "G", "T"}
  neighborhood = set()
  suffix_neighbors = Neighbors(pattern[1:], d)
  for suffix in suffix_neighbors:
    if Hamming_Distance(suffix, pattern[1:]) < d:
      for Nucleotide in {"A", "T", "G", "C"}:
        neighborhood.add(Nucleotide + suffix)
    else:
      neighborhood.add(pattern[0] + suffix)
  print(*neighborhood)
  return neighborhood

# Generate the Reverse strand of a given DNA sequence 
def ReverseDNA(seq):
  pattern = seq[::-1]
  pattern = pattern.replace("A", "t")
  pattern = pattern.replace("T", "a")
  pattern = pattern.replace("G", "c")
  pattern = pattern.replace("C", "g")
  pattern = pattern.upper()
  return pattern


# ONE OF THE MAIN FUNCTIONS IN THIS PROJECT 
# given a sequence (seq) a k-mer value (k which is any number specified) and a nucleotide mutation value, can be set to zero if wanted,(d)
# The function returns the pattern which satisfies these conditions, has the maximum amount of appearance in the sequence and has its reverse
# This is crucial for DNAA box identification which triggers the initiation of the replication.
def kmer_with_d(seq, k, d):
  pattern_counts = {}

  for i in range(len(seq) - k + 1):
    kmer = seq[i:i+k]
    neighbors = Neighbors(kmer, d)
    for neighbor in neighbors:
      if neighbor not in pattern_counts:
        pattern_counts[neighbor] = 1
      else:
        pattern_counts[neighbor] += 1

  count_final = {}

  for pattern in pattern_counts:
    rev_seq = ReverseDNA(pattern)

    count_normal = pattern_counts[pattern]
    count_reverse = pattern_counts.get(rev_seq,0) # set at zero if not existing

    total_count = count_normal + count_reverse
    count_final[pattern] = total_count 
    count_final[rev_seq] = total_count

  max_count = 0
  results = []

  for pattern, count in count_final.items():
    if count > max_count:
      max_count = count
      results = [pattern]
    elif count == max_count:
      results.append(pattern)

  print(*results)
  return 

def MotifEnumeration(seq, k, d):
  patterns = set() # set to handle the duplciates

  def ApproxPattern2(pattern, seq, d):
    count = 0
    k = len(pattern)

    for i in range(len(seq) - k + 1):
      splices = seq[i:i+k]
      if Hamming_Distance(pattern, splices) <= d:
        count += 1
    return count

  first_seq = seq.partition(" ")[0]
  rest_seq = seq.partition(" ")[2]

  for i in range(len(first_seq) - k + 1):
    kmer = first_seq[i:i+k]
    neighbors_first_Seq = Neighbors(kmer, d)

    for potential in neighbors_first_Seq:
      motif = True
      for seqs in rest_seq.split():
        if not ApproxPattern2(potential, seqs, d) > 0:
          motif = False
          break
      if motif:
        patterns.add(potential)

  return " ".join(patterns)

filename = input("Enter the name of your file: ")

with open(filename, "r") as f:
  lines = f.readlines()

seq = lines[0].strip()
k = int(lines[1].strip())

profile_matrix = []
for i in range(2, 6):
  row = lines[i].strip().split()
  row_float = [float(x) for x in row]
  profile_matrix.append(row_float)

print("\nThe profile matrix is: {}".format(profile_matrix))
print("\nSequence is: {}".format(seq))
print("\nk-mer is: {}".format(k))

def ProbFunction(k_mer, p_matrix):
  probability = 1

  map = {"A": 0, "C": 1, "G": 2, "T": 3}
  for a in range(len(k_mer)):
    letter = k_mer[a]
    row = map[letter]
    probability *= profile_matrix[row][a]
  return probability

max_prob = -1
best_kmer = seq[0:k] # start at the first kmer

for i in range(len(seq) - k + 1):
  kmer = seq[i:i+k]

  probability = ProbFunction(kmer, profile_matrix)

  if probability > max_prob:
    max_prob = probability
    best_kmer = kmer
print("\nBest matching k-mer is : {}".format(best_kmer))
