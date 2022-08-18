def CountWithPseudocounts(Motifs):
    counts = {}
    for symbol in "ACTG":
        counts[symbol] = []
        for i in range(len(Motifs[0])):
            counts[symbol].append(1)
    for i in range(len(Motifs)):
        for j in range(len(Motifs[0])):
            counts[Motifs[i][j]][j] += 1
    return counts

def ProfileWithPseudocounts(Motifs):
    profile = CountWithPseudocounts(Motifs)
    for i in "ACTG":
        for j in range(len(Motifs[0])):
            profile[i][j] = profile[i][j]/(len(Motifs)+4)
    return profile

def computeprobability(Text, Profile):
    prob = 1
    for i in range(len(Text)):
        prob = prob * Profile[Text[i]][i]
    return prob

def ProfileMostProbablePattern(Text, k, Profile):
    prob = {}
    for i in range(len(Text) - k + 1):
        prob[Text[i:i+k]] = computeprobability(Text[i:i+k], Profile)
    for key in prob:
        if prob[key] == max(prob.values()):
            return key

def ConsensusMotif(Motifs):
    counts = CountWithPseudocounts(Motifs)
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
    consensus = ConsensusMotif(Motifs)
    score = 0
    for i in range(len(Motifs)):
        for j in range(len(Motifs[0])):
            if Motifs[i][j] != consensus[j]:
                score += 1
    return score

def GreedyMotifSearchWithPseudocounts(Dna, k):
    bestmotifs = []
    for i in range(len(Dna)):
        bestmotifs.append(Dna[i][0:k])
    for i in range(len(Dna[0]) - k + 1):
        motifs = []
        motifs.append(Dna[0][i:i+k])
        for j in range(1, len(Dna)):
            P = ProfileWithPseudocounts(motifs[0:j])
            motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        if Score(motifs) < Score(bestmotifs):
            bestmotifs = motifs
    return bestmotifs


def Motifs(Dna, Profile, k):
    output = []
    for i in range(len(Dna)):
        for j in range(len(Dna[0]) - k + 1):
            if Dna[i][j:j+k] == ProfileMostProbablePattern(Dna[i], k, Profile):
                output.append(Dna[i][j:j+k])
    return output



import random as rnd

def RandomMotifs(Dna, k):
    output = []
    for i in range(len(Dna)):
        a = rnd.randint(0, len(Dna[0]) - k)
        output.append(Dna[i][a:a+k])
    return output

def RandomizedMotifSearch(Dna, k):
    m = RandomMotifs(Dna, k)
    bestmotifs = m
    while True:
        profile = ProfileWithPseudocounts(m) # Generate a Profile matrix from the motifs of the random set of strings
        m = Motifs(Dna, profile, k) # Obtain the best motifs according to the odds from the Profile matrix
        if Score(m) < Score(bestmotifs):
            bestmotifs = m # If the motifs have a better score than the randomly selected, substitute them
        else:
            return bestmotifs # When the motifs no longer improve, you now have the BEST motifs in the set of strings

def callRandomizedMotifSearch(Dna, k, N):
    bestmotifs = RandomizedMotifSearch(Dna, k)
    for i in range(N):
        newmotifs = RandomizedMotifSearch(Dna, k)
        if Score(newmotifs) < Score(bestmotifs):
            bestmotifs = newmotifs
    return bestmotifs, Score(bestmotifs)

dna = ['TTACCTTAAC', 'GATGTCTGTC', 'CCGGCGTTAG', 'CACTAACGAG', 'CGTCAGAGGT']

#print(callRandomizedMotifSearch(dna, k = 4, N = 100))

def Motifs_2(Dna, Profile, k):
    D = []
    for i in range(len(Dna)):
        km = []
        sc = []
        for j in range(len(Dna[0]) - k + 1):
            km.append(Dna[i][j:j+k])
        for item in km:
            sc.append(computeprobability(item, Profile))
        D.append(km[sc.index(max(sc))])
    return D

Dna = ["GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC", "CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG", "ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC", "GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC", "GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG", "CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA", "GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA", "GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG", "GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG", "TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC"]

# print(callRandomizedMotifSearch(Dna = Dna, k = 15, N = 100000))

# Chapter 1.4 Gibbs Sampling
def Normalize(Probabilities):
    normalized_probabilities = {}
    for i in Probabilities:
        normalized_probabilities.update({i: Probabilities[i]/sum(Probabilities.values())})
    return normalized_probabilities

def WeightedDie(Probabilities):
    n = rnd.uniform(0, 1)
    for i in range(len(Probabilities)):
        values = list(Probabilities.values())
        if n <= sum(values[:i+1]):
            return list(Probabilities.keys())[i]


Probabilities = {
'AA' : 0.2,
'TT' : 0.2,
'CC' : 0.1,
'GG' : 0.1,
'AT' : 0.4
}
# print(WeightedDie(Probabilities))

def ProfileGeneratedString(Text, Profile, k):
    probabilities = {}
    for i in range(len(Text) - k + 1):
        probabilities.update({Text[i:i+k]: computeprobability(Text[i:i+k], Profile)})
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)


Text = "AAACCCAAACCC"
Profile = {'A': [0.5, 0.1], 'C': [0.3, 0.2], 'G': [0.2, 0.4], 'T': [0.0, 0.3]}

# print(ProfileGeneratedString(Text, Profile, 2))

def GibbsSampler(Dna, k, N):
    motifs = RandomMotifs(Dna, k)
    bestmotifs = motifs
    for i in range(N):
        n = rnd.randint(0, len(Dna)-1)
        motifs_new = motifs
        motifs_new.pop(n)
        profile = ProfileWithPseudocounts(motifs_new)
        motifs_new.insert((n), ProfileGeneratedString(Dna[n], profile, k))
        if Score(motifs_new) < Score(bestmotifs):
            bestmotifs = motifs_new
    return bestmotifs

# Another way of deleting random k-mer
def GibbsSampler_2(Dna, k, t, N):
    bestmotifs = []
    motifs = RandomMotifs(Dna, k)
    bestmotifs = motifs
    for i in range(N):
        n = rnd.randint(0, len(Dna)-1)
        motifs_new = []
        for j in range(len(Dna)):
            if j != n:
                motifs_new.append(motifs[j])
        profile = ProfileWithPseudocounts(motifs_new)
        motifs_new.insert(n, ProfileGeneratedString(Dna[n], profile, k))
        if Score(motifs_new) < Score(bestmotifs):
            bestmotifs = motifs_new
    return bestmotifs

Dna_final =["GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC", "CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG", "ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC", "GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC", "GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG", "CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA", "GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA", "GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG", "GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG", "TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC"]

def GibbsSamplerCaller(Dna, k, N, N_2):
    BestMotifs = GibbsSampler(Dna, k, N)
    for i in range(N_2):
        NewBestMotifs = GibbsSampler(Dna, k, N)
        if Score(NewBestMotifs) < Score(BestMotifs):
            BestMotifs = NewBestMotifs
    return BestMotifs

#result = GibbsSamplerCaller(Dna_final, k = 15, N = 100, N_2 = 200)
#print(result)
#print(Score(result))

dna_quiz = ["TGACGTTC", "TAAGAGTT", "GGACGAAA", "CTGTTCGC"]
motifs_quiz = ["TGA", "GTT", "GAA", "TGT"]
profile = ProfileWithPseudocounts(motifs_quiz)
m = Motifs(dna_quiz, profile, 3)
print(m)

y=rnd.randint(1,10)
print(y)

probabilities_quiz = {'1': 0.22, '2': 0.54, '3': 0.58, '4': 0.36, '5': 0.3}
print(Normalize(probabilities_quiz))