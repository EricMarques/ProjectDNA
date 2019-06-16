# -*- coding: utf-8 -*-

import os
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO

path = "Samples/"

dirs = os.listdir(path)

query = input("Arquivo com sequência para alinhamento: ")

arq_file_q = path + query

queryT = arq_file_q + ".fasta"

for seq_record in SeqIO.parse(queryT, 'fasta'):
    seqDna = seq_record.seq


def slicesequence(x):
    # print(x)
    aux = x[0:3]
    final = x[3:]
    print('\n' + aux)
    for i in x:
        while aux:
            aux = final[:3]
            final = final[3:]
            print(aux)

    return ''


'''
text = open('cp_teste.txt', 'r', encoding='UTF8')
text.close()

coding_dna = Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG', IUPAC.unambiguous_dna)  # Sequência DNA
print('Seq. Original: DNA =  {}' .format(coding_dna))
messenger_rna = coding_dna.transcribe()  # Trenscrever DNA para RNAm
print('Seq. Original: RNAm = {}' .format(messenger_rna))
'''

seqDNA = str(seqDna)
codeSeqDNA = Seq(seqDNA, IUPAC.unambiguous_dna)
messengerRNA = codeSeqDNA.transcribe()

print('DNA: {}\n'.format(seqDna))
print('Caracteres DNA: {}'.format(len(seqDNA)))
print('\nRNAm: {}\n' .format(messengerRNA))
print('Caracteres RNAm: {}'.format(len(messengerRNA)))
print(messengerRNA.translate())

# print(slicesequence(seqDNA))  # Chama a função slicesequence() passando como parâmetro a sequência RNAm
