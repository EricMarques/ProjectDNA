from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Blast.Applications import *


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

linhas = open("resultado.txt").readlines()

for linha in linhas:
    print(linha)
