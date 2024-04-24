# Martin AYMÉ
# ENS de Lyon
# Python, UTF-8
# January 2023

## Librairies

import matplotlib.pyplot as plt
import numpy as np
import sklearn.metrics as sk
# Si,Ca ou Ti, Fe, Rb ou Sr


## Courbes

def trace(voltage,calib):
    resultsPB = loadData("avg_results_calib_PB.csv")
    resultsCC = loadData(f"avg_results_{voltage}_{calib}.csv")

    resultsLab = loadData("INCONNUS-CONNUS_NoAt.csv")

    N = 25 # default nb of element studied (LV & HV)

    if voltage == "LV":
        resultsLab = [resultsLab[i][:13] for i in range (len(resultsLab))]
        resultsPB = [resultsPB[i][:13] for i in range (len(resultsPB))]
        N = 12 # CC LV calib considers 12 elements (Na-Fe)
    else :
        resultsLab = [[resultsLab[i][0]]+resultsLab[i][12:] for i in range (len(resultsLab))]
        resultsPB = [[resultsPB[i][0]]+resultsPB[i][12:] for i in range (len(resultsPB))]
        N = 14 # CC HV calib considers 14 elements (Fe-Ba)
    # sample the list of elements matching with either the LV or HV list of elements
    # of the CC calib

    plt.figure()

    for i in range(1,N+1):
        plt.subplot(3,5,i)
        plt.title(resultsLab[0][i])
        CC = [resultsCC[j][i] for j in range(1,len(resultsCC))]
        PB = [resultsPB[j][i] for j in range(1,len(resultsPB))]
        Lab = [resultsLab[j][i] for j in range(1,len(resultsLab))]
        plt.scatter(Lab,CC, label="CC", marker='+', color='b')
        plt.scatter(Lab,PB, label="PB", marker='x', color='orange')
        mini = min(Lab+CC+PB)
        maxi = max(Lab+CC+PB)
        marge = (maxi - mini)*0.05
        plt.plot([mini,maxi],[mini,maxi],linestyle='--', color='grey',label='Droite 1:1')
        plt.xlim(mini-marge,maxi+marge)
        plt.ylim(mini-marge,maxi+marge)
        #plt.legend()


    plt.show()

def traceTwoCalib(voltage1,calib1,voltage2,calib2,element):
    resultsPB = loadData("avg_results_calib_PB.csv")
    resultsCC1 = loadData(f"avg_results_{voltage1}_{calib1}.csv")
    resultsCC2 = loadData(f"avg_results_{voltage2}_{calib2}.csv")
    resultsM = loadData("avg_results_mining.csv")

    resultsLab = loadData("INCONNUS-CONNUS_NoAt.csv")


    if voltage1 == "LV":
        resultsLab = [resultsLab[i][:13] for i in range (len(resultsLab))]
        resultsPB = [resultsPB[i][:13] for i in range (len(resultsPB))]
    else :
        resultsLab = [[resultsLab[i][0]]+resultsLab[i][12:] for i in range (len(resultsLab))]
        resultsPB = [[resultsPB[i][0]]+resultsPB[i][12:] for i in range (len(resultsPB))]


    resultsCC1,resultsCC2,resultsPB,resultsM,resultsLab = supprZerosDbl(resultsCC1,resultsCC2,resultsPB,resultsM,resultsLab)

    allElements = resultsLab[0][1:]
    elementIndex = allElements.index(element)

    plt.title(element)
    plt.subplot(1,2,1)
    CC1 = [resultsCC1[j][elementIndex+1] for j in range(1,len(resultsCC1))]
    CC2 = [resultsCC2[j][elementIndex+1] for j in range(1,len(resultsCC2))]
    PB = [resultsPB[j][elementIndex+1] for j in range(1,len(resultsPB))]
    Lab = [resultsLab[j][elementIndex+1] for j in range(1,len(resultsLab))]
    plt.scatter(Lab,CC1, label=f"CC_{voltage1}_{calib1}", marker='+', color='b')
    plt.scatter(Lab,PB, label="PB", marker='x', color='orange')
    mini = min(Lab+CC1+CC2+PB)
    maxi = max(Lab+CC1+CC2+PB)
    marge = (maxi - mini)*0.05
    plt.plot([mini,maxi],[mini,maxi],linestyle='--', color='grey',label='Droite 1:1')
    plt.xlim(mini-marge,maxi+marge)
    plt.ylim(mini-marge,maxi+marge)
    plt.legend()

    plt.subplot(1,2,2)
    plt.scatter(Lab,CC2, label=f"CC_{voltage2}_{calib2}", marker='+', color='b')
    plt.scatter(Lab,PB, label="PB", marker='x', color='orange')
    plt.plot([mini,maxi],[mini,maxi],linestyle='--', color='grey',label='Droite 1:1')
    plt.xlim(mini-marge,maxi+marge)
    plt.ylim(mini-marge,maxi+marge)
    plt.legend()


    plt.show()


def traceOneElem(voltage,calib,element):
    resultsPB = loadData("avg_results_calib_PB.csv")
    resultsCC = loadData(f"avg_results_{voltage}_{calib}.csv")
    resultsM = loadData("avg_results_mining.csv")

    resultsLab = loadData("INCONNUS-CONNUS_NoAt.csv")


    if voltage == "LV":
        resultsLab = [resultsLab[i][:13] for i in range (len(resultsLab))]
        resultsPB = [resultsPB[i][:13] for i in range (len(resultsPB))]
        resultsM = [resultsM[i][:13] for i in range (len(resultsM))]
    else :
        resultsLab = [[resultsLab[i][0]]+resultsLab[i][12:] for i in range (len(resultsLab))]
        resultsPB = [[resultsPB[i][0]]+resultsPB[i][12:] for i in range (len(resultsPB))]
        resultsM = [[resultsM[i][0]]+resultsM[i][12:] for i in range (len(resultsM))]


    allElements = resultsLab[0][1:]
    elementIndex = allElements.index(element)

    CC = [resultsCC[j][elementIndex+1] for j in range(1,len(resultsCC))]
    PB = [resultsPB[j][elementIndex+1] for j in range(1,len(resultsPB))]
    M = [resultsM[j][elementIndex+1] for j in range(1,len(resultsM))]
    Lab = [resultsLab[j][elementIndex+1] for j in range(1,len(resultsLab))]

    CC,PB,M,Lab = supprZeros(CC,PB,M,Lab)

    r2 = [getR2(Lab,CC),getR2(Lab,PB),getR2(Lab,M)]

    plt.figure()
    plt.subplot(111)
    plt.scatter(Lab,CC, label=f"CloudCal_{voltage}_{calib}\nR2={r2[0]}", marker='+', color='b')
    plt.scatter(Lab,PB, label=f"Profile Builder\nR2={r2[1]}", marker='x', color='orange')
    plt.scatter(Lab,M, label=f"Mining\nR2={r2[2]}", marker='.', color='green')
    mini = min(Lab+CC+PB+M)
    maxi = max(Lab+CC+PB+M)
    marge = (maxi - mini)*0.05
    plt.plot([mini,maxi],[mini,maxi],linestyle='--', color='grey',label='Droite 1:1')
    plt.xlim(mini-marge,maxi+marge)
    plt.ylim(mini-marge,maxi+marge)
    plt.legend()
    plt.title(element, fontsize="30")
    plt.xlabel("Concentration de référence (wt %)", fontsize="20")
    plt.ylabel("Concentration mesurée (wt %)", fontsize="20")
    plt.subplots_adjust(left=0.16)

    #plt.savefig(f"{element}_validation2.svg")
    plt.show()


def traceMining(voltage,calib):
    resultsPB = loadData("avg_results_calib_PB.csv")
    resultsCC = loadData(f"avg_results_{voltage}_{calib}.csv")
    resultsM = loadData("avg_results_mining.csv")

    resultsLab = loadData("INCONNUS-CONNUS_NoAt.csv")

    N = 25 # default nb of element studied (LV & HV)

    if voltage == "LV":
        resultsLab = [resultsLab[i][:13] for i in range (len(resultsLab))]
        resultsPB = [resultsPB[i][:13] for i in range (len(resultsPB))]
        resultsM = [resultsM[i][:13] for i in range (len(resultsM))]
        N = 12 # CC LV calib considers 12 elements (Na-Fe)
    else :
        resultsLab = [[resultsLab[i][0]]+resultsLab[i][12:] for i in range (len(resultsLab))]
        resultsPB = [[resultsPB[i][0]]+resultsPB[i][12:] for i in range (len(resultsPB))]
        resultsM = [[resultsM[i][0]]+resultsM[i][12:] for i in range (len(resultsM))]
        N = 14 # CC HV calib considers 14 elements (Fe-Ba)
    # sample the list of elements matching with either the LV or HV list of elements
    # of the CC calib

    plt.figure()

    for i in range(1,N+1):
        plt.subplot(5,3,i)
        plt.title(resultsLab[0][i])
        CC = [resultsCC[j][i] for j in range(1,len(resultsCC))]
        PB = [resultsPB[j][i] for j in range(1,len(resultsPB))]
        M = [resultsM[j][i] for j in range(1,len(resultsM))]
        Lab = [resultsLab[j][i] for j in range(1,len(resultsLab))]
        plt.scatter(Lab,CC, label="CC", marker='+', color='b')
        plt.scatter(Lab,PB, label="PB", marker='x', color='orange')
        plt.scatter(Lab,M, label="M", marker='.', color='green')
        mini = min(Lab+CC+PB+M)
        maxi = max(Lab+CC+PB+M)
        marge = (maxi - mini)*0.05
        plt.plot([mini,maxi],[mini,maxi],linestyle='--', color='grey',label='Droite 1:1')
        plt.xlim(mini-marge,maxi+marge)
        plt.ylim(mini-marge,maxi+marge)
        #plt.legend()


    plt.show()


def supprZeros(resultsCC,resultsPB,resultsM,resultsLab):
    notZerosIndex = [-1]*len(resultsLab)
    Lab = [-1]*len(resultsLab)
    CC = [-1]*len(resultsLab)
    PB = [-1]*len(resultsLab)
    M = [-1]*len(resultsLab)

    notZerosIndex = [j for j in range(len(resultsLab)) if resultsLab[j] != 0]

    Lab = [resultsLab[j] if j != -1 else -1 for j in notZerosIndex]
    CC = [resultsCC[j] for j in notZerosIndex if j != -1]
    PB = [resultsPB[j] for j in notZerosIndex if j != -1]
    M = [resultsM[j] for j in notZerosIndex if j != -1]

    return CC, PB, M, Lab


def supprZerosDbl(resultsCC1,resultsCC2,resultsPB,resultsM,resultsLab):
    notZerosIndex = [[-1]*len(resultsLab[0])]*len(resultsLab)
    Lab = [resultsLab[0]] + [[-1]*len(resultsLab[0])]*(len(resultsLab)-1)
    CC1 = [resultsCC1[0]] + [[-1]*len(resultsLab[0])]*(len(resultsLab)-1)
    CC2 = [resultsCC2[0]] + [[-1]*len(resultsLab[0])]*(len(resultsLab)-1)
    PB = [resultsPB[0]] + [[-1]*len(resultsLab[0])]*(len(resultsLab)-1)
    M = [resultsM[0]] + [[-1]*len(resultsLab[0])]*(len(resultsLab)-1)

    for i in range(1,len(resultsLab)):

        notZerosIndex[i] = [j for j in range(len(resultsLab[i])) if resultsLab[i][j] != 0]

        Lab[i] = [resultsLab[i][j] for j in notZerosIndex[i] if j != -1]
        CC1[i] = [resultsCC1[i][j] for j in notZerosIndex[i] if j != -1]
        CC2[i] = [resultsCC2[i][j] for j in notZerosIndex[i] if j != -1]
        PB[i] = [resultsPB[i][j] for j in notZerosIndex[i] if j != -1]
        M[i] = [resultsM[i][j] for j in notZerosIndex[i] if j != -1]

    return CC1, CC2, PB, M, Lab


def rechIndexEleCom(resultsPB,resultsCC,resultsLab,eleCommun):
    n = len(eleCommun)
    indexPB, indexCC, indexLab = [0]*n,[0]*n,[0]*n

    for i in range(n):
        indexPB[i] = resultsPB.index(eleCommun[i])
        indexCC[i] = resultsCC.index(eleCommun[i])
        indexLab[i] = resultsLab.index(eleCommun[i])

    return indexPB, indexCC, indexLab


def loadData(file):
    f = open(file)
    tout = f.readlines()

    for i in range(len(tout)):
        tout[i] = tout[i].strip()
        tout[i] = tout[i].split(",")
        if i>0:
            tout[i] = [tout[i][0]] + [float(data) if data != '' else 0.0 for data in tout[i][1:]]
    f.close()
    return tout


def rechEleCom(resultsPB,resultsCC,resultsLab):
    eleCommun = []

    eleCC = resultsCC[0][1:]
    elePB = resultsPB[0][1:]

    for i in range(len(1,resultsLab)):
        if resultsLab[0][i] in eleCC and resultsLab[0][i] in elePB:
            eleCommun.append(resultsLab[0][i])

    return eleCommun


## Quantification de l'erreur avec coef de détermination

def getR2(resultsExp,resultsLab):
    return round(sk.r2_score(resultsLab,resultsExp),5)

def getAllR2(voltage,calib):
    resultsPB = loadData("avg_results_calib_PB.csv")
    resultsCC = loadData(f"avg_results_{voltage}_{calib}.csv")
    resultsM = loadData("avg_results_mining.csv")

    resultsLab = loadData("INCONNUS-CONNUS_NoAt.csv")

    if voltage == "LV":
        resultsLab = [resultsLab[i][:13] for i in range (len(resultsLab))]
        resultsPB = [resultsPB[i][:13] for i in range (len(resultsPB))]
        resultsM = [resultsM[i][:13] for i in range (len(resultsM))]
    else :
        resultsLab = [[resultsLab[i][0]]+resultsLab[i][12:] for i in range (len(resultsLab))]
        resultsPB = [[resultsPB[i][0]]+resultsPB[i][12:] for i in range (len(resultsPB))]
        resultsM = [[resultsM[i][0]]+resultsM[i][12:] for i in range (len(resultsM))]

    allElements = resultsLab[0][1:]
    R2 = [0]*len(allElements)

    f = open(f"all_R2_{voltage}.csv","w")

    for i in range(len(allElements)):

        elementIndex = allElements.index(allElements[i])

        CC = [resultsCC[j][elementIndex+1] for j in range(1,len(resultsCC))]
        PB = [resultsPB[j][elementIndex+1] for j in range(1,len(resultsPB))]
        M = [resultsM[j][elementIndex+1] for j in range(1,len(resultsM))]
        Lab = [resultsLab[j][elementIndex+1] for j in range(1,len(resultsLab))]

        CC,PB,M,Lab = supprZeros(CC,PB,M,Lab)

        try :
            R2[i] = [getR2(Lab,CC),getR2(Lab,PB),getR2(Lab,M)]
            R2[i] = [round(R2Cal,5) for R2Cal in R2[i]]
        except :
            R2[i] = None

        f.write(allElements[i]+",")
        try :
            f.write(f"{R2[i][0]},{R2[i][1]},{R2[i][2]}\n")
        except :
            f.write("None\n")

    f.close()

    return R2


## Réarrengement CC et PB (No de spectre → Nom du kn-unkn) + moyenne des 3 spectres


def RAZresultsCC(file):

    name = ['GR2', 'GR3', 'GR4', 'GR5', 'GR6', 'CAVA2018-1', 'CAVA2018-4', 'LAV2018-1', 'LAV2018-4', 'KH04-14', 'SAN2019-1', 'KG03B']
    ranges = [[139, 140, 141], [142, 143, 144], [145, 146, 147], [172, 173, 174], [148, 149, 150], [151, 152, 153], [154, 155, 156], [157, 158, 159], [160, 161, 162], [163, 164, 165], [166, 167, 168], [169, 170, 171]]

    tout = loadData(file)

    f = open(f"modif_{file}","w")

    k = 0
    f.write(tout[0][0])
    for i in range(1,len(tout[0])):
        f.write(","+tout[0][i])
    f.write("\n")
    for i in range(1,len(tout),3):
        f.write(name[k]+",")

        for j in range(1,len(tout[0])):
            f.write(str((tout[i][j]+tout[i+1][j]+tout[i+2][j])/3))

            if j != len(tout[0])-1:
                f.write(",")
            else :
                f.write("\n")

        k += 1

    f.close()


def RAZresults_spec(file):

    name = ['GR2', 'GR3', 'GR4', 'GR5', 'GR6', 'CAVA2018-1', 'CAVA2018-4', 'LAV2018-1', 'LAV2018-4', 'KH04-14', 'SAN2019-1', 'KG03B']
    spectrumNb = [[139, 140, 141], [142, 143, 144], [145, 146, 147], [172, 173, 174], [148, 149, 150], [151, 152, 153], [154, 155, 156], [157, 158, 159], [160, 161, 162], [163, 164, 165], [166, 167, 168], [169, 170, 171]]

    tout = loadData(file)

    f = open(f"avg_{file}","w")

    tout[0][0] = "Spectrum_#"
    if tout[1][0][:4] == 'Test':
        for i in range(1,len(tout)):
            tout[i][0] = int(tout[i][0].split("_")[1][1:]) # take the spectrum nb in the
            # complete name ('Test_#159_B1-10.0keV...' → 159)

    allSpectrumNb = [0]+[int(tout[i][0]) for i in range(1,len(tout))] # [0] is
    # there to replace the 'Spectrum_#'

    f.write(tout[0][0])
    for i in range(1,len(tout[0])):
        f.write(","+tout[0][i])
    f.write("\n")

    k = 0
    for i in range(len(spectrumNb)):
        f.write(name[k]+",")
        index = [0,0,0] # position of the spectrum corresponding to the 'name[k]'
        # sample

        for j in range(3):
            index[j] = allSpectrumNb.index(spectrumNb[i][j])

        for j in range(1,len(tout[0])):
            f.write(str((tout[index[0]][j]+tout[index[1]][j]+tout[index[2]][j])/3))

            if j != len(tout[0])-1:
                f.write(",")
            else :
                f.write("\n")

        k += 1

    f.close()


def corresp_nb_id(file):

    name = ['GR2', 'GR3', 'GR4', 'GR5', 'GR6', 'CAVA2018-1', 'CAVA2018-4', 'LAV2018-1', 'LAV2018-4', 'KH04-14', 'SAN2019-1', 'KG03B']
    spectrumNb = [[139, 140, 141], [142, 143, 144], [145, 146, 147], [172, 173, 174], [148, 149, 150], [151, 152, 153], [154, 155, 156], [157, 158, 159], [160, 161, 162], [163, 164, 165], [166, 167, 168], [169, 170, 171]]

    tout = loadData(file)

    f = open(f"ID_{file}","w")

    tout[0][0] = "Spectrum_#"
    if tout[1][0][:4] == 'Test':
        for i in range(1,len(tout)):
            tout[i][0] = int(tout[i][0].split("_")[1][1:]) # take the spectrum nb in the
            # complete name ('Test_#159_B1-10.0keV...' → 159)

    allSpectrumNb = [0]+[int(tout[i][0]) for i in range(1,len(tout))] # [0] is
    # there to replace the 'Spectrum_#'

    f.write(tout[0][0])
    for i in range(1,len(tout[0])):
        f.write(","+tout[0][i])
    f.write("\n")

    k = 0
    for i in range(len(spectrumNb)):

        for j in range(3):
            f.write(name[k]+",")
            index = allSpectrumNb.index(spectrumNb[i][j])

            for l in range(1,len(tout[0])):
                f.write(str(tout[index[0]][l]))


                if l != len(tout[0])-1:
                    f.write(",")
                else :
                    f.write("\n")

        k += 1

    f.close()


def RAZresults_mining(file):

    name = ['GR2', 'GR3', 'GR4', 'GR5', 'GR6', 'CAVA2018-1', 'CAVA2018-4', 'LAV2018-1', 'LAV2018-4', 'KH04-14', 'SAN2019-1', 'KG03B']
    spectrumNb = [[231, 232, 233], [234, 235, 236], [237, 238, 239], [240, 241, 242], [243, 244, 245], [246, 247, 248], [249, 250, 251], [252, 253, 254], [255, 256, 257], [258, 259, 260], [261, 262, 263], [264, 265, 266]]

    tout = loadData(file)

    f = open(f"avg_{file}","w")

    tout[0][0] = "Spectrum_#"
    if tout[1][0][:4] == 'Test':
        for i in range(1,len(tout)):
            tout[i][0] = int(tout[i][0].split("_")[1][1:]) # take the spectrum nb in the
            # complete name ('Test_#159_B1-10.0keV...' → 159)

    allSpectrumNb = [int(tout[i][0]) for i in range(1,len(tout))] # all spectrum nb
    # present in the results file

    correctUnit = np.zeros([len(tout)-1,len(tout[0])-1]) # to convert the ppm in %
    # correctUnit[0,:] = tout[0] # not possible to put string in a np array
    for i in range(1,len(tout)):
        for j in range(1,len(tout[0])):
            if int(tout[i][j]) - float(tout[i][j]) != 0: # true if tout[i][j] in %
            # i.e. is a float (ppm results are integers)
                correctUnit [i-1,j-1] = tout[i][j]
            else:
                correctUnit [i-1,j-1] = tout[i][j]/10000

    f.write(tout[0][0])
    for i in range(1,len(tout[0])):
        f.write(","+tout[0][i])
    f.write("\n")


    for i in range(len(spectrumNb)):
        f.write(name[i]+",")
        index = [0,0,0] # positions of spectra corresponding to the 'name[i]'
        # sample

        for j in range(3):
            index[j] = allSpectrumNb.index(spectrumNb[i][j])

        for j in range(len(tout[0])-1):
            f.write(str((correctUnit[index[0],j]+correctUnit[index[1],j]+correctUnit[index[2],j])/3))

            if j != len(tout[0])-2:
                f.write(",")
            else :
                f.write("\n")


    f.close()



