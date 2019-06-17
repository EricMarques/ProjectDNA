# -*- coding: utf-8 -*-

import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Blast.Applications import *

path = "Samples/"
dirs = os.listdir(path)

# BLAST = Biopython

query = input("Arquivo com sequência para alinhamento: ")
subject = input("Arquivo com sequência original: ")  # arquivo com o qual será alinhada a 'query'
arq_out = input("Nome arquivo saída: ")

arq_file_q = path + query
arq_file_s = path + subject
arq_out_f = path + arq_out

queryT = arq_file_q + ".fasta"
subjectT = arq_file_s + ".fasta"
arq_outT = arq_out_f + ".txt"

comand = NcbiblastpCommandline(
    query=queryT,
    subject=subjectT,
    outfmt=0,
    out=arq_outT
)

stdout, stderr = comand()
result = open(arq_outT, "r")

print(result.read())

result.close()

# BLAST = linha de comando
# MODELO # os.system("blastp -query a.fasta -subject b.fasta -outfmt 6 > resultado.txt")  # alinhamento de proteínas
# MODELO # os.system("blastn -query a.fasta -subject b.fasta > test.txt")  # alinhamento de nucleotídeos

# os.system("blastn -query " + query + " -subject " + subject + " > resultado.txt")
# os.system("blastn -query " + query + " -subject " + subject + " -outfmt 6 > resumo.txt")

for seq_record in SeqIO.parse(queryT, 'fasta'):
    seqDna = seq_record.seq


gene = Seq(str(seqDna), IUPAC.unambiguous_dna)

print(len(gene))
# print(len(gene.translate()))

linhas = open(arq_outT).readlines()

for linha in linhas:
    print(linha)

import menu
menu.choice


'''
    Forçar a mutação da sequencia
    "gerar" novamente aminoácidos
    Verificar as proteínas gerada ,se permaneceu a mesma
    comparar com "original"
    
    
    Identifica algumas proteinas pela sequencia de bases
    nas diferencas, tenta identificar possiveis mutações que fariam a sequencia deixar de ser uma proteina e virar outra
'''