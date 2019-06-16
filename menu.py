# -*- coding: utf-8 -*-

import time

print("********Software para busca e alinhamento de sequências de Aminoácidos********\n")
print("1) Download de arquivo de uma Proteína específica.\n"
      "2) Alinhamento de Aminoácidos.\n"
      "3) Hello World.\n"
      "9) Sair")
choice = int(input("Escolha uma opção: "))
choices = (1, 2, 3, 9)

if choice == 1:
    import searchNCBI as bNcbi

    bNcbi.query

elif choice == 2:
    import align as a

    a.query

elif choice == 3:
    import slice
    slice.query

elif choice == 9:
    print("\n\nObrigado pela utilização!")
    time.sleep(5)
    exit()

elif choice != choices:
    print("Opção inválida!")
    time.sleep(3)
