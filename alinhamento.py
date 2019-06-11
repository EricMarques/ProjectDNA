import os
import xmltodict

import Exceptions as e
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from Bio.pairwise2 import *
from Bio.Blast.Applications import *

from Bio.SeqUtils import GC
import  pylab

path = "Samples/"

dirs = os.listdir(path)
'''
query = input("Arquivo com sequência para alinhamento: ")

arq_file = path + query


if os.path.isfile(arq_file):
    try:
        print("Arquivo '" + query + "' encontrado!")
        arq = open(path + query, "r")
        if query not in arq_file:
            raise NameError(path+query)

    except NameError:
        print("OPA")
        raise

    finally:
        
        arq = open(path + query, "r")
        print(arq.read())


print(os.path.isfile(arq_file))
print(arq_file)
'''

with open(path + "calr_blastp.xml") as fd:
    doc = xmltodict.parse(fd.read())

# print(doc["BlastOutput"]["BlastOutput_iterations"]["Iteration"]["Iteration_query-def"])
# print(doc["BlastOutput"]["BlastOutput_iterations"]["Iteration"]["Iteration_hits"]["Hit"])
print(doc["BlastOutput"]["BlastOutput_iterations"]["Iteration"]["Iteration_hits"]["Hit"][0]["Hit_num"])


# calr_blastp.xml

'''
# converter arquivos
SeqIO.convert(path + "carl_blastp.xml", "xml", "carl_blastp.fasta", "fasta")



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
