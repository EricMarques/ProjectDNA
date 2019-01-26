fileName = str(input("Arquivo: "))

arq = open(fileName + '.txt', 'r')
arqCopy = open('cp_' + arq + '.txt', 'w')

arqCopy.writelines(arq)

print(arq.read())

arqCopy.close()
arq.close()

