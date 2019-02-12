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


def seq(x):
    print('i: {} \nf: {}'.format(x[0:3], x[3:]))
    for i in range(len(x)):
        slicef = x[3:]
        slicei = slicef[0:3]
        print(slicef)
    return slicei


coding_dna = Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG', IUPAC.unambiguous_dna)
print('Seq. Original: ' + coding_dna)
print('Caracteres Seq. Orig: {}'.format(len(coding_dna)))
'''
#print(coding_dna.translate())

#print(coding_dna.translate(to_stop=True))#table=#2, stop_symbol='@'))#"Vertebrate Mitochondrial"))

slice1 = coding_dna[:3]
slice2 = coding_dna[3:]

print('Slice 1 - {}'.format(slice1))
print('Slice 2 = {}'.format(slice2))
print('Caracteres slice 2 = {}'.format(len(slice2)))

slice3 = slice2[:3]
slice4 = slice2[3:]

print('Slice 3 = {}'.format(slice3))
print('Slice 4 = {}'.format(slice4))
print('Caracteres slice 4 = {}'.format(len(slice4)))
'''
print(seq(coding_dna))
