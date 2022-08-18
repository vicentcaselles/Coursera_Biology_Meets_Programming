def SymbolArray(Genome, Symbol):
    output = {}
    ExtendedGenome = Genome + Genome[0:len(Genome)//2]
    for i in range(len(ExtendedGenome) - (len(Genome)//2)):
        output.update({i: ExtendedGenome[i:i+len(Genome)//2].count(Symbol)})
    return output

def SymbolArrayLite(Genome, Symbol):
    output = {}
    ExtendedGenome = Genome + Genome[0:len(Genome)//2]
    output[0] = ExtendedGenome[0:len(Genome)//2].count(Symbol)
    for i in range(1, len(Genome)):
        output[i] = output[i - 1]
        if ExtendedGenome[i-1] == Symbol:
            output.update({i: output[i] - 1})
        if ExtendedGenome[i + len(Genome)//2 - 1] == Symbol:
            output.update({i: output[i] + 1})
    return output

f = open("/Users/vicentcaselles/PycharmProjects/Biology_Meets_Programming_Coursera/Week_2/E_coli.txt", "r")
E_Coli = f.read()

# array = SymbolArrayLite(E_Coli, Symbol = "C")

## To plot the results
# import matplotlib.pyplot as plt
# plt.plot(*zip(*sorted(array.items())))
# plt.show()

def GC_counter(Genome):
    G_minus_C = {}
    ExtendedGenome = Genome + Genome[0:len(Genome)//2]
    G_minus_C[0] = ExtendedGenome[0:len(Genome)//2].count("G") - ExtendedGenome[0:len(Genome)//2].count("C")
    for i in range(1, len(Genome)):
        G_minus_C[i] = G_minus_C[i-1]
        if ExtendedGenome[i-1] == "G":
            G_minus_C.update({i: G_minus_C[i] - 1})
        if ExtendedGenome[i-1] == "C":
            G_minus_C.update({i: G_minus_C[i] + 1})
        if ExtendedGenome[i + len(Genome)//2 - 1] == "G":
            G_minus_C.update({i: G_minus_C[i] + 1})
        if ExtendedGenome[i + len(Genome)//2 - 1] == "C":
            G_minus_C.update({i: G_minus_C[i] - 1})
    return G_minus_C

#results_1 = GC_counter(E_Coli)
#import matplotlib.pyplot as plt
#plt.plot(*zip(*sorted(results_1.items())))
#plt.show()

def SkewArray(Genome):
    skewarray = {}
    skewarray[0] = 0
    for i in range(1, len(Genome) + 1):
        if Genome[i-1] == "A":
            skewarray[i] = skewarray[i-1]
        if Genome[i-1] == "C":
            skewarray[i] = skewarray[i-1] - 1
        if Genome[i-1] == "T":
            skewarray[i] = skewarray[i-1]
        if Genome[i-1] == "G":
            skewarray[i] = skewarray[i-1] + 1
    return skewarray

import matplotlib.pyplot as plt
#results_2 = SkewArray(E_Coli)
#plt.plot(*zip(*sorted(results_2.items())))
#plt.show()

def MinSkewArray(Genome):
    array = SkewArray(Genome)
    minGen = []
    for i in range(len(Genome)):
        if array[i] == min(array.values()):
            minGen.append(i)
    return minGen

print(f'MinSkewArray = {MinSkewArray("CATTCCAGTACTTCGATGATGGCGTGAAGA")}')

def HammingDistance(p, q):
    HamDist = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            HamDist += 1
    return HamDist

print(HammingDistance("GGGCCGTTGGT", "GGACCGTTGAC"))

def ApproximatePatternMatching(Text, Pattern, d):
    positions = []
    for i in range(len(Text) - len(Pattern) + 1):
        if HammingDistance(Text[i:i+len(Pattern)], Pattern) <= d:
            positions.append(i)
    return positions

print(ApproximatePatternMatching("CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT", "ATTCTGGA", 3))

def ApproximatePatternCount(Pattern, Text, d):
    count = 0
    for i in range(len(Text) - len(Pattern) + 1):
        if HammingDistance(Text[i:i+len(Pattern)], Pattern) <= d:
            count += 1
    return count

print(ApproximatePatternCount("GAGG", "TTTAGAGCCTTCAGAGG", 2))

## Identify the most frequent k-mers with 1 mismatch
def count_kmers(seq, k):
    counts = {}
    for i in range(len(seq) - k + 1):
        if seq[i:i + k] in counts:
            counts[seq[i:i + k]] += 1
        else:
            counts[seq[i:i + k]] = 1
    return counts

def approx_count_kmers(seq, k, d):
    counts = count_kmers(seq, k)
    for i in range(len(seq)-k+1):
        for key in counts:
            if HammingDistance(key, seq[i:i+k]) <= d:
                counts.update({f'{key}': seq.count(key)})
    return counts

def reverse_complementary(seq):
    rev_com = ""
    for i in range(len(seq)):
        if seq[i] == "A":
            rev_com += "T"
        elif seq[i] == "T":
            rev_com += "A"
        elif seq[i] == "C":
            rev_com += "G"
        elif seq[i] == "G":
            rev_com += "C"
    return rev_com[::-1]