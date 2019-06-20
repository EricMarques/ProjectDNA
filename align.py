# -*- coding: utf-8 -*-

import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Blast.Applications import *

pathI = "./Samples/CALR/CALRs/Proteins/"
pathO = "./Samples/CALR/CALRs/Results/"
dirs = os.listdir(pathI)

# BLAST = Biopython
query = input("Arquivo com sequência para alinhamento: ")
subject = input("Arquivo com sequência original: ")
arq_out = input("Nome arquivo saída: ")

arq_file_q = pathI + query
arq_file_s = pathI + subject
arq_out_f = pathO + arq_out

queryT = arq_file_q + ".fasta"
subjectT = arq_file_s + ".fasta"
arq_outT = arq_out_f + ".txt"
arq_outT_Resume = arq_out_f + "_resume.txt"

comand = NcbiblastpCommandline(
    query=queryT,
    subject=subjectT,
    outfmt=0,
    out=arq_outT
)

comand2 = NcbiblastpCommandline(
    query=queryT,
    subject=subjectT,
    outfmt=6,
    out=arq_outT_Resume
)

stdout, stderr = comand()
stdout, stderr = comand2()

result = open(arq_outT, "r")
result2 = open(arq_outT_Resume, "r")

print(result.read())

result.close()
result2.close()

for seq_record in SeqIO.parse(queryT, 'fasta'):
    seqDna = seq_record.seq


gene = Seq(str(seqDna), IUPAC.unambiguous_dna)

print(len(gene))

linhas = open(arq_outT).readlines()

for linha in linhas:
    print(linha)
