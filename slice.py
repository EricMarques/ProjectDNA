def slicesequence(x):
    print(x)
    aux = x[0:3]
    final = x[3:]
    print('\n' + aux)
    for i in x:
        aux = final[:3]
        final = final[3:]
        print(aux)
    '''
    print('1')
    aux = x[0:3]
    print(aux)
    final = x[3:]
    print(final)
    print('2')
    aux = final[:3]
    print(aux)
    final = final[3:]
    print(final)
    ...
    '''

seq = "ABCDEFGHIJKLMNOPQR"

slicesequence(seq)
