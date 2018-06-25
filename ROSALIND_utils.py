from collections import OrderedDict

def read_fasta(fp):
    seqs = OrderedDict()
    current = None
    buff = []
    for line in fp:
        if line.startswith('>'):
            if current:
                seqs[current] = ''.join(buff).rstrip()
                buff = []
            current = line.rstrip()
        else:
            buff.append(line.rstrip())
    seqs[current] = ''.join(buff).rstrip()
    return seqs

def gc_content(seq):
    """GC content of the sequence"""
    return sum([seq.count(x) for x in 'CG'])/len(seq)

def RNA_codon(seq):
    """Translating individual RNA codons into amino acids"""
    RNA_codon_table = {'UUU':'F','UUC':'F','UUA':'L','UUG':'L','UCU':'S','UCC':'S','UCA':'S',
                      'UCG':'S','UAU':'Y','UAC':'Y','UAA':'Stop','UAG':'Stop','UGU':'C','UGC':'C',
                      'UGA':'Stop','UGG':'W','CUU':'L','CUC':'L','CUA':'L','CUG':'L','CCU':'P',
                      'CCC':'P','CCA':'P','CCG':'P','CAU':'H','CAC':'H','CAA':'Q','CAG':'Q',
                      'CGU':'R','CGC':'R','CGA':'R','CGG':'R','AUU':'I','AUC':'I','AUA':'I',
                      'AUG':'M','ACU':'T','ACC':'T','ACA':'T','ACG':'T','AAU':'N','AAC':'N',
                      'AAA':'K','AAG':'K','AGU':'S','AGC':'S','AGA':'R','AGG':'R','GUU':'V',
                      'GUC':'V','GUA':'V','GUG':'V','GCU':'A','GCC':'A','GCA':'A','GCG':'A',
                      'GCG':'A','GAU':'D','GAC':'D','GAA':'E','GAG':'E','GGU':'G','GGC':'G',
                      'GGA':'G','GGG':'G'}
    return RNA_codon_table[seq]