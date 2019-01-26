fileName = str(input("Arquivo: "))

arq = open(fileName + '.txt', 'r')
arqCopy = open('cp_' + arq + '.txt', 'w')

arq.writelines(arqCopy)

print(arq.read())

arqCopy.close()
arq.close()

