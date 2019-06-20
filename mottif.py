import os

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


def long_mottif(data):
	substr = ''
	if len(data) > 1 and len(data[0]) > 0:
		for i in range(len(data[0])):
			for j in range(len(data[0])-i + 1):
				if j > len(substr) and all(data[0][i:i+j] in x for x in data):
					substr = data[0][i:i+j]
	return substr


path = "./Samples/CALR/"
dirs = os.listdir(path)

seq3 = open(path + "seqTransB.fasta", "r")

content = seq3.read()

lines = []

for data in content.split('>'):
	line = ''.join(data.split('\n')[1:]).replace('\n', '')
	if line:
		lines.append(line)

tr = Seq(long_mottif(lines), IUPAC.unambiguous_dna)
print(lines)
print("Seq Comum: {} \nTamanho: {} caracteres" .format(long_mottif(lines), len(long_mottif(lines))))
# print("Prot: {}\nTamanho: {} caracteres" .format(tr.translate(), len(tr.translate())))

seq3.close()
