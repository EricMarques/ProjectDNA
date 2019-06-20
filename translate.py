import os

from Bio import SeqIO


path = "./Samples/CALR/CALRs/Nucleotides/"
dirs = os.listdir(path)

query = input("Arquivo com sequência para alinhamento: ")
seq3 = open(path + query + ".fasta", "r")

for seq_record in SeqIO.parse(seq3, 'fasta'):
	seqDna = seq_record.seq

se = str(seqDna)

print("Sequência Original: {}" .format(seqDna))
print("Sequência Translate: {}" .format(seqDna.translate()))
print("Len: {}".format(len(seqDna.translate())))
