import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from Bio.pairwise2 import *
from Bio.Blast import NCBIWWW
from Bio.Blast.Applications import *

from Bio.SeqUtils import GC
import  pylab

path = "./Samples/"

dirs = os.listdir(path)

for file in dirs:
    if os.path.exists(path + "oi.txt"):
        print(file)

print("\n")


'''
# Buscar proteínas = blastp
# Buscar nucleotídeos = blastn
print("Buscando...")
blast_result = NCBIWWW.qblast("blastn", "nr", "CM000663.2")
blast_out = open("calrN_blastn.out", "w")
blast_out.write(blast_result.read())
blast_out.close()
blast_result.close()
print("Fim da busca.")


# converter arquivos
SeqIO.convert("nome_arquivo.formato_antigo", "formato_antigo" / "nome_arquivo.formato_novo", "formato_novo")
'''

'''
Alinhamento
align = align.globalxx("ACGTACGTAGCA", "ACG")
for a in align:
    print(format_alignment(*a))
'''

'''
arq = SeqIO.read("exemplo_fasta", format="fasta")
print("Buscando...")

result = NCBIWWW.qblast("blastn", "nt" / arq.seq, format_type="Text")
print(result.re)
print("Fim")
'''

'''
# Visualização de dados
seq = Seq("GTGTGTGTGTCGCGCGCGCGGCGGGCGCACTCCGTGTGTGTGTGTTGTACATC")
gc = GC(seq)
at = 100 - gc

pylab.poly([gc, at])
pylab.title("GC Content")
pylab.xlabel("GC: %0.1f\nAT: %0.1f" % (gc, at))
pylab.show()
'''

'''
# BLAST = Biopython
query = input("Arquivo com sequência para alinhamento: ")
subject = input("Arquivo com sequência original: ")  # arquivo com o qual será alinhada a 'query'
arq_out = input("Nome arquivo saída: ")

comand = NcbiblastnCommandline(
	query=query,
	subject=subject,
	outfmt=0,
	out=arq_out+".txt"
)

stdout, stderr = comand()
result = open(arq_out+".txt", "r")

print(result.read())

# BLAST = linha de comando
# MODELO # os.system("blastp -query a.fasta -subject b.fasta -outfmt 6 > resultado.txt")  # alinhamento de proteínas
# MODELO # os.system("blastn -query a.fasta -subject b.fasta > test.txt")  # alinhamento de nucleotídeos

# os.system("blastn -query " + query + " -subject " + subject + " > resultado.txt")
# os.system("blastn -query " + query + " -subject " + subject + " -outfmt 6 > resumo.txt")

for seq_record in SeqIO.parse(query, 'fasta'):
	seqDna = seq_record.seq

gene = Seq(str(seqDna), IUPAC.unambiguous_dna)

print(len(gene))
print(len(gene.translate()))
'''

'''
linhas = open("resultado.txt").readlines()

for linha in linhas:
	print(linha)
'''