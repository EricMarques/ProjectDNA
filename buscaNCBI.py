from Bio.Blast import NCBIWWW

path = "Samples/"

query = input("Digite o nome da Proteína a ser buscada: ")

'''
arq_file = path + query

if os.path.isfile(arq_file):
    try:
        print("Arquivo '" + query + "' encontrado!")
        arq = open(path + query, "r")
        if query not in arq_file:
            raise NameError(path + query)

    except NameError:
        print("OPA")
        raise

    finally:

        arq = open(path + query, "r")
        print(arq.read())

print(os.path.isfile(arq_file))
print(arq_file)
'''

# Buscar proteínas = blastp
# Buscar nucleotídeos = blastn
print("Buscando...")
blast_result = NCBIWWW.qblast("blastp", "nr", query)
blast_out = open(path + query + ".xml", "w")
blast_out.write(blast_result.read())
blast_out.close()
blast_result.close()
print("Fim da busca.")
