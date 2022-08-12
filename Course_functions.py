# Function to count how many times a pattern appears in a sequence
def PatternCount(Pattern, Text):
    count = 0
    for i in range(0, len(Text) - len(Pattern) + 1):
        if Text[i:i + len(Pattern)] == Pattern:
            count += 1
    return count


# Frequency map function (not mine):
def FrequencyMap(Text, k):
    freq = {}
    for i in range(len(Text) - k + 1):
        Pattern = Text[i:i + k]
        if Pattern in freq:
            freq[Pattern] += 1
        else:
            freq[Pattern] = 1
    return freq


# Function to count most frequent kmers (does same as FrequencyMap)
def count_kmers(seq, k):
    counts = {}
    for i in range(len(seq) - k + 1):
        if seq[i:i + k] in counts:
            counts[seq[i:i + k]] += 1
        else:
            counts[seq[i:i + k]] = 1
    return counts


# Function to obtain the most frequent k-mers
def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
    return words


# Function to return kmers with max count (same as FrequentWords)

def max_kmers(seq, k):
    max_kmers = []
    for key in count_kmers(seq, k):
        if count_kmers(seq, k)[key] == max(count_kmers(seq, k).values()):
            max_kmers.append(key)
    return max_kmers


# function to calculate the probability of a given k-mer in a n lenght seq
def calc_probability(seq, k):
    prob_k_mer = (1 / 4 ** k) * len(seq) / k
    return (prob_k_mer)


# Reverse complementary (one step)
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


# Reverse complementary (Two steps)
def Reverse(Pattern):
    Rev = Pattern[::-1]
    return Rev


def Complement(Pattern):
    Comp = ""
    for i in Pattern:
        if i == "A":
            Comp += "T"
        if i == "T":
            Comp += "A"
        if i == "C":
            Comp += "G"
        if i == "G":
            Comp += "C"
    return Comp


def ReverseComplement(Pattern):
    Pattern = Reverse(Pattern)
    Pattern = Complement(Pattern)
    return Pattern


# Another way of Reversing DNA strand
def Rev_3(seq):
    rev = ""
    for i in range(len(seq)):
        rev += seq[len(seq) - i - 1]
    return rev


# Function to obtain starting positions where a given pattern is located
def PatternMatching(Pattern, Genome):
    Positions = []
    for i in range(len(Genome) - len(Pattern) + 1):
        if Genome[i:i + len(Pattern)] == Pattern:
            Positions.append(i)
    return Positions


#Shows you the most frequent k-mers with > min_k frequency
def Frequent_words_2(seq, k, min_k):
    max_kmers_2 = {}
    for key, value in count_kmers(seq, k).items():
        if count_kmers(seq, k)[key] >= min_k:
            max_kmers_2.update({f'{key}': value})
    return max_kmers_2

#Shows the counts of k-mers + their reverse complementary
def ForwardReverseCount(seq, k, min_k):
    pair = {}
    freqwords = Frequent_words_2(seq, k, min_k)
    for key, value in freqwords.items():
        if reverse_complementary(key) in freqwords.keys():
            pair.update({f'{key}/{reverse_complementary(key)}': value + freqwords.get(reverse_complementary(key))})
            if f'{reverse_complementary(key)}/{key}' in pair.keys():
                pair.pop(f'{reverse_complementary(key)}/{key}')
    return pair