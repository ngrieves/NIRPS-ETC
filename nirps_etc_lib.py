# Library for NIRPS ETC
# Last change Dez/2017
#


# TODO: Documentar o que e cada funcao
#

def reading_table(file_name):
    TABLE_ARCHIVE=open(file_name,'rt')
    lines=TABLE_ARCHIVE.readlines()
    M=[]

    for line in lines:
        if not line.strip()[0] == '#':
            v=line.strip().split()
            M.append(list(map(float,v)))

    TABLE_ARCHIVE.close()
    return M

def linear_interpolation(M,wc):

    asc = True
    if (M[0][0]>M[1][0]):
        asc = False

    f = -1

    idx1 = 0
    idx2 = len(M)-1
    while not idx2 == idx1+1:
        idx = int((idx1+idx2)/2)
        w = M[idx][0]
        if wc>w:
            idx1 = idx if asc else idx1
            idx2 = idx if not asc else idx2
        else:
            idx2 = idx if asc else idx2
            idx1 = idx if not asc else idx1

    if not asc:
        tmp = idx1
        idx1 = idx2
        idx2 = tmp

    w_a = M[idx1][0]
    w_p = M[idx2][0]
    f_a = M[idx1][1]
    f_p = M[idx2][1]
    f = f_a + (wc-w_a)*(f_p-f_a)/(w_p-w_a)

    return f

# def linear_interpolation(M,wc):
#     w_anterior=M[0][0] # comprimento de onda: elemento 0 da linha 0
#     f_anterior=M[0][1] # fluxo: elemento 1 da linha zero
#
#     f1=-1
#
#     for m in M:
#         w_posterior=m[0]
#         f_posterior=m[1]
#         if w_posterior > w_anterior:
#             if w_posterior >= wc and w_anterior <= wc:
#                 w1=w_anterior
#                 f1=f_anterior
#                 w2=w_posterior
#                 f2=f_posterior
#         else:
#             if w_posterior <= wc and w_anterior >= wc:
#                 w1=w_posterior
#                 f1=f_posterior
#                 w2=w_anterior
#                 f2=f_anterior
#         w_anterior=w_posterior
#         f_anterior=f_posterior
#
#     if f1 < 0:
#         print("########### ERRO ###############\n")
#         print(M)
#         print(wc)
#
#     f=f1+(wc-w1)/(w2-w1)*(f2-f1)
#
#     return f

def atmosphere_efficiencies(airmass):

    file_name='tapas_wcentral_interpol.txt'
    lines=reading_table(file_name)
    file_name2='tapas_wcentral10_interpol.dat'
    lines2=reading_table(file_name2)
#    file_name3='tapas_wcentral10_interpol.dat'
#    lines3=reading_table(file_name3)

    ATM_EFF=[]
    ATM_EFF10=[]
#    ATM_EFF20=[]

    for line in lines:

        if airmass == 1.0:
            eff=line[1]
            ATM_EFF.append(eff)

        elif airmass > 1.0 and airmass < 1.02:
            f=line[1]+(line[2]-line[1])*(airmass-1.00)/0.02
            eff=f
            ATM_EFF.append(eff)

        elif airmass == 1.02:
            eff=line[2]
            ATM_EFF.append(eff)

        elif airmass > 1.02 and airmass < 1.06:
            f=line[2]+(line[3]-line[2])*(airmass-1.02)/0.04
            eff=f
            ATM_EFF.append(eff)

        elif airmass == 1.06:
            eff=line[3]
            ATM_EFF.append(eff)

        elif airmass > 1.06 and airmass < 1.15:
            f=line[3]+(line[4]-line[3])*(airmass-1.06)/0.09
            eff=f
            ATM_EFF.append(eff)

        elif airmass == 1.15:
            eff=line[4]
            ATM_EFF.append(eff)

        elif airmass > 1.15 and airmass < 1.31:
            f=line[4]+(line[5]-line[4])*(airmass-1.15)/0.16
            eff=f
            ATM_EFF.append(eff)

        elif airmass == 1.31:
            eff=line[5]
            ATM_EFF.append(eff)

        elif airmass > 1.31 and airmass < 1.56:
            f=line[5]+(line[6]-line[5])*(airmass-1.31)/0.25
            eff=f
            ATM_EFF.append(eff)

        elif airmass == 1.56:
            eff=line[6]
            ATM_EFF.append(eff)

        elif airmass > 1.56 and airmass < 2.00:
            f=line[6]+(line[7]-line[6])*(airmass-1.56)/0.44
            eff=f
            ATM_EFF.append(eff)

        elif airmass == 2.00:
            eff=line[7]
            ATM_EFF.append(eff)

#    print(ATM_EFF)

    for line in lines2:

        if airmass == 1.0:
            eff10=line[1]
            ATM_EFF10.append(eff10)

        elif airmass > 1.0 and airmass < 1.02:
            f=line[1]+(line[2]-line[1])*(airmass-1.00)/0.02
            eff10=f
            ATM_EFF10.append(eff10)

        elif airmass == 1.02:
            eff10=line[2]
            ATM_EFF10.append(eff10)

        elif airmass > 1.02 and airmass < 1.06:
            f=line[2]+(line[3]-line[2])*(airmass-1.02)/0.04
            eff10=f
            ATM_EFF10.append(eff10)

        elif airmass == 1.06:
            eff10=line[3]
            ATM_EFF10.append(eff10)

        elif airmass > 1.06 and airmass < 1.15:
            f=line[3]+(line[4]-line[3])*(airmass-1.06)/0.09
            eff10=f
            ATM_EFF10.append(eff10)

        elif airmass == 1.15:
            eff10=line[4]
            ATM_EFF10.append(eff10)

        elif airmass > 1.15 and airmass < 1.31:
            f=line[4]+(line[5]-line[4])*(airmass-1.15)/0.16
            eff10=f
            ATM_EFF10.append(eff10)

        elif airmass == 1.31:
            eff10=line[5]
            ATM_EFF10.append(eff10)

        elif airmass > 1.31 and airmass < 1.56:
            f=line[5]+(line[6]-line[5])*(airmass-1.31)/0.25
            eff10=f
            ATM_EFF10.append(eff10)

        elif airmass == 1.56:
            eff10=line[6]
            ATM_EFF10.append(eff10)

        elif airmass > 1.56 and airmass < 2.00:
            f=line[6]+(line[7]-line[6])*(airmass-1.56)/0.44
            eff10=f
            ATM_EFF10.append(eff10)

        elif airmass == 2.00:
            eff10=line[7]
            ATM_EFF10.append(eff10)

#    print(ATM_EFF10)

    return ATM_EFF,ATM_EFF10
