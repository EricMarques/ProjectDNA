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
    '''
    print('1')
    aux = x[0:3]
    print(aux)
    final = x[3:]
    print(final)
    print('2')
    aux = final[:3]
    print(aux)
    final = final[3:]
    print(final)
    ...
    '''


coding_dna = Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG', IUPAC.unambiguous_dna)  # Sequência DNA
print('Seq. Original: DNA = {}' .format(coding_dna))
coding_rna = coding_dna.transcribe()  # Trenscrever DNA para RNAm
print('Seq. Original: RNA = {}' .format(coding_rna))
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
print(slicesequence(coding_rna))  # Chama a função slicesequence() passando como parâmetro a sequência RNAm
