import os
from Bio.Blast import NCBIWWW
import time

path = "Samples/"
dirs = os.listdir(path)

query = input("Digite o nome da Proteína a ser buscada: ")

arq_file = path + query
queryT = arq_file + ".xml"


def buscaNcbi(query):
    # Buscar proteínas = blastp
    # Buscar nucleotídeos = blastn
    try:
        print("Buscando arquivo...")
        blast_result = NCBIWWW.qblast("blastn", "nr", query)
        blast_out = open(arq_file + ".xml", "w")
        blast_out.write(blast_result.read())
        blast_out.close()
        blast_result.close()
        print("Fim da busca.\nArquivo " + query + ".xml encontra-se disponível no diretório '" +
              path + "' para análise")

    except ValueError:
        print("\nProteína inexistente ou inválida\n\n\n")
        time.sleep(5)


if os.path.exists(queryT):
    print("Arquivo '" + query + ".xml' já existente no diretório '" + path.strip("/") + "'!")

else:
    buscaNcbi(query)

import menu
menu.choice