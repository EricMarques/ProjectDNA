from Bio.Blast import NCBIWWW

path = "Samples/"

file_request = input("Digite o nome da Proteína a ser buscada: ")

# Buscar proteínas = blastp
# Buscar nucleotídeos = blastn
print("Buscando...")
blast_result = NCBIWWW.qblast("blastp", "nr", file_request)
blast_out = open(path + file_request + ".xml", "w")
blast_out.write(blast_result.read())
blast_out.close()
blast_result.close()
print("Fim da busca.")
