from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

'''
#text = open('cp_teste.txt', 'r', encoding='UTF8')
#text.close()
#seq = 'GATCGATGGGCCTATATAGGATCGAAAATCGC'#text.read()

seq = Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG', IUPAC.unambiguous_dna)
print('DNA:  ' + seq)

#seq1 = seq.transcribe()
seq1 = seq.reverse_complement().transcribe()

print('mRNA: ' + seq1)'''


def slicesequence(x):
    print(x)
    aux = x[0:3]
    final = x[3:]
    print('\n' + aux)
    for i in x:
        while aux:
            aux = final[:3]
            final = final[3:]
            print(aux)

coding_dna = Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG', IUPAC.unambiguous_dna)  # Sequência DNA
print('Seq. Original: DNA = {}' .format(coding_dna))
coding_rna = coding_dna.transcribe()  # Trenscrever DNA para RNAm
print('Seq. Original: RNA = {}' .format(coding_rna))
print('Caracteres Seq. Orig: {}'.format(len(coding_dna)))

print(slicesequence(coding_rna))  # Chama a função slicesequence() passando como parâmetro a sequência RNAm
