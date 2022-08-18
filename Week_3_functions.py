# Entropy calculation

import math

import numpy as np

P = np.array([[0.2, 0.2, 0, 0, 0, 0, 0.9, 0.1, 0.1, 0.1, 0.3, 0], [0.1,0.6,0,0,0,0,0,0.4,0.1,0.2,0.4,0.6],
              [0,0,1,1,0.9,0.9,0.1,0,0,0,0,0], [0.7,0.2,0,0,0.1,0.1,0,0.5,0.8,0.7,0.3,0.4]])

H = 0
for i in range(12):
    for j in range(4):
        if P[j, i] == 0:
            H = H
        if P[j, i] != 0:
            H += -(P[j, i] * math.log(P[j,i], 2))
print(H)

# Motif count
def MotifCount(input):
    counts = {}
    for symbol in "ACTG":
        counts[symbol] = []
        for i in range(len(input[0])):
            counts[symbol].append(0)
    for i in range(len(input)):
        for j in range(len(input[0])):
            if input[i][j] == "A":
                counts['A'][j] +=1
            if input[i][j] == "G":
                counts['G'][j] += 1
            if input[i][j] == "C":
                counts['C'][j] += 1
            if input[i][j] == "T":
                counts['T'][j] += 1
    return counts


#Motif Profile
def MotifProfile(input):
    counts = {}
    for symbol in "ACTG":
        counts[symbol] = []
        for i in range(len(input[0])):
            counts[symbol].append(0)
    for i in range(len(input)):
        for j in range(len(input[0])):
            if input[i][j] == "A":
                counts['A'][j] +=1
            if input[i][j] == "G":
                counts['G'][j] += 1
            if input[i][j] == "C":
                counts['C'][j] += 1
            if input[i][j] == "T":
                counts['T'][j] += 1
    for i in "ACTG":
        for j in range(len(input[0])):
            counts[i][j] = counts[i][j]/len(input)
    return counts

def ConsensusMotif(Motifs):
    counts = MotifCount(Motifs)
    consensus = ""
    for j in range(len(Motifs[0])):
        m = 0
        frequentsymbol = ""
        for symbol in "ACTG":
            if counts[symbol][j] > m:
                m = counts[symbol][j]
                frequentsymbol = symbol
        consensus += frequentsymbol
    return consensus

def Score(Motifs):
    Consensus = ConsensusMotif(Motifs)
    Score = 0
    for i in range(len(Motifs)):
        for j in range(len(Motifs[0])):
            if Motifs[i][j] != Consensus[j]:
                Score += 1
    return Score


# function to calculate the probability of finding a sequence, given a profile of nucleotide probabilities
# for a consensus region


def computeprobability(Text, Profile):
    prob = 1
    for i in range(len(Text)):
        prob = prob * Profile[Text[i]][i]
    return prob



def ProfileMostProbablePattern(Text, k, Profile):
    MostProbablePattern = ""
    prob = {}
    for i in range(len(Text) - k + 1):
        prob.update({Text[i:i+k]: computeprobability(Text[i:i+k], Profile)})
    for key in prob:
        if prob[key] == max(prob.values()):
            MostProbablePattern = key
            return(MostProbablePattern) #Weird indentation to get the first substring with max value



def GreedyMotifSearch(Dna, k):
    BestMotifs = []
    for i in range(len(Dna)):
        BestMotifs.append(Dna[i][0:k])
    for i in range(len(Dna[0]) - k + 1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, len(Dna)):
            P = MotifProfile(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

Dna = ['GAGGCGCACATCATTATCGATAACGATTCGCCGCATTGCC', 'TCATCGAATCCGATAACTGACACCTGCTCTGGCACCGCTC', 'TCGGCGGTATAGCCAGAAAGCGTAGTGCCAATAATTTCCT',
       'GAGTCGTGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTG', 'GACGGCAACTACGGTTACAACGCAGCAACCGAAGAATATT', 'TCTGTTGTTGCTAACACCGTTAAAGGCGGCGACGGCAACT',
       'AAGCGGCCAACGTAGGCGCGGCTTGGCATCTCGGTGTGTG', 'AATTGAAAGGCGCATCTTACTCTTTTCGCTTTCAAAAAAA']


print(GreedyMotifSearch(Dna, k=5))

Miau = {'A': [0.4, 0.3, 0.0, 0.1, 0.0, 0.9],
           'C': [0.2, 0.3, 0.0, 0.4, 0.0, 0.1],
           'G': [0.1, 0.3, 1.0, 0.1, 0.5, 0.0],
           'T': [0.3, 0.1, 0.0, 0.4, 0.5, 0.0]}
print(computeprobability("CAGTGA", Miau))

DNA = ['TTACCTTAAC', 'GATGTCTGTC', 'ACGGCGTTAG', 'CCCTAACGAG', 'CGTCAGAGGT']
print(GreedyMotifSearch(Dna, 5))