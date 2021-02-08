import os
import sys
import matplotlib
import pylab as plt
import numpy as np


def PlotEnergyData(annihilator):
    grabbedenergydict=annihilator.masterdict['energy']
    patharray=[]
    solvHarray=[]
    solvHerrarray=[]
    solvSarray=[]
    solvSerrarray=[]
    solvGarray=[]
    solvGerrarray=[]
    compHarray=[]
    compHerrarray=[]
    compSarray=[]
    compSerrarray=[]
    compGarray=[]
    compGerrarray=[]
    for path in grabbedenergydict.keys():
        head,tail=os.path.split(path)
        patharray.append(tail)
        pathdict=grabbedenergydict[path]
        if u'ΔHˢᵒˡᵛ' in pathdict.keys():
            solvHarray.append(pathdict[u'ΔHˢᵒˡᵛ'])
        else:
            solvHarray.append(0)
        if u'ΔHˢᵒˡᵛᵉʳʳ' in pathdict.keys():
            solvHerrarray.append(pathdict[u'ΔHˢᵒˡᵛᵉʳʳ'])
        else:
            solvHerrarray.append(0)
        if u'ΔSˢᵒˡᵛ' in pathdict.keys():
            solvSarray.append(pathdict[u'ΔSˢᵒˡᵛ'])
        else:
            solvSarray.append(0)
        if u'ΔSˢᵒˡᵛᵉʳʳ' in pathdict.keys():
            solvSerrarray.append(pathdict[u'ΔSˢᵒˡᵛᵉʳʳ'])
        else:
            solvSerrarray.append(0)
        if u'ΔGˢᵒˡᵛ' in pathdict.keys():
            solvGarray.append(pathdict[u'ΔGˢᵒˡᵛ'])
        else:
            solvGarray.append(0)
        if u'ΔGˢᵒˡᵛᵉʳʳ' in pathdict.keys():
            solvGerrarray.append(pathdict[u'ΔGˢᵒˡᵛᵉʳʳ'])
        else:
            solvGerr.append(0)
        if u'ΔHᶜᵒᵐᵖ' in pathdict.keys():
            compHarray.append(pathdict[u'ΔHᶜᵒᵐᵖ'])
        else:
            compHarray.append(0)
        if u'ΔHᶜᵒᵐᵖᵉʳʳ' in pathdict.keys():
            compHerrarray.append(pathdict[u'ΔHᶜᵒᵐᵖᵉʳʳ'])
        else:
            compHerrarray.append(0)
        if u'ΔSᶜᵒᵐᵖ' in pathdict.keys():
            compSarray.append(pathdict[u'ΔSᶜᵒᵐᵖ'])
        else:
            compSarray.append(0)
        if u'ΔSᶜᵒᵐᵖᵉʳʳ' in pathdict.keys():
            compSerrarray.append(pathdict[u'ΔSᶜᵒᵐᵖᵉʳʳ'])
        else:
            compSerrarray.append(0)
        if u'ΔGᶜᵒᵐᵖᶜᵒʳʳ' in pathdict.keys():
            compGarray.append(pathdict[u'ΔGᶜᵒᵐᵖᶜᵒʳʳ'])
        else:
            compGarray.append(0)
        if u'ΔGᶜᵒᵐᵖᶜᵒʳʳᵉʳʳ' in pathdict.keys():
            compGerrarray.append(pathdict[u'ΔGᶜᵒᵐᵖᶜᵒʳʳᵉʳʳ'])
        else:
            compGerrarray.append(0)

    x = 3       # the width of the bars: can also be len(x) sequence
    width=1
    fig1=plt.figure(figsize=(8, 10))
    
    solvSarray=[0 if v is None or v is '' else v for v in solvSarray]
    solvSarray=[float(i) for i in solvSarray]
    solvSerrarray=[0 if v is None or v is '' else v for v in solvSerrarray]
    solvSerrarray=[float(i) for i in solvSerrarray]

    solvGerrarray=[0 if v is None or v is '' else v for v in solvGerrarray] 
    solvGerrarray=[float(i) for i in solvGerrarray]

    solvGarray=[0 if v is None or v is '' else v for v in solvGarray] 
    solvGarray=[float(i) for i in solvGarray]

    solvSarray=[0 if v is None or v is '' else v for v in solvSarray]
    solvSarray=[float(i) for i in solvSarray]

    solvHerrarray=[0 if v is None or v is '' else v for v in solvHerrarray]
    solvHerrarray=[float(i) for i in solvHerrarray]

    solvHarray=[0 if v is None or v is '' else v for v in solvHarray]
    solvHarray=[float(i) for i in solvHarray]

    compHerrarray=[0 if v is None or v is '' else v for v in compHerrarray]
    compHerrarray=[float(i) for i in compHerrarray]

    compSerrarray=[0 if v is None or v is '' else v for v in compSerrarray]
    compSerrarray=[float(i) for i in compSerrarray]

    compGerrarray=[0 if v is None or v is '' else v for v in compGerrarray]
    compGerrarray=[float(i) for i in compGerrarray]

    compHarray=[0 if v is None or v is '' else v for v in compHarray]
    compHarray=[float(i) for i in compHarray]

    compSarray=[0 if v is None or v is '' else v for v in compSarray]
    compSarray=[float(i) for i in compSarray]

    compGarray=[0 if v is None or v is '' else v for v in compGarray]
    compGarray=[float(i) for i in compGarray]

    solvSarray=np.array(solvSarray)*(int(annihilator.roomtemp))
    compSarray=np.array(compSarray)*(int(annihilator.roomtemp))
    Hind = np.arange(0,3*len(patharray)-2,x)    # the x locations for the groups
    Sind = [i+1 for i in Hind]
    Gind = [i+2 for i in Hind] # plot _TDeltaS instead of DeltaS

    if annihilator.solvation==True or annihilator.binding==True:    
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True,figsize=(8,10))
        Hp = ax2.bar(Hind, solvHarray, yerr=solvHerrarray,color='red')
        Sp = ax2.bar(Sind, solvSarray, yerr=solvSerrarray,color='blue')
        Gp = ax1.bar(Gind, solvGarray, yerr=solvGerrarray,color='green')
        ax1.set_ylabel('Solvation Absolute Energies (kcal/mol)')
        ax1.set_title('Absolute Solvation Energies')
        plt.xticks(Sind, tuple(patharray))
        plt.legend((Hp, Sp,Gp), (u'ΔHˢᵒˡᵛ', u'TΔSˢᵒˡᵛ',u'ΔGˢᵒˡᵛ'))
        fig.savefig(annihilator.outputpath+'SolvAbsBarPlot.png')

    
    if annihilator.complexation==True or annihilator.binding==True:    

        fig2, (ax1, ax2) = plt.subplots(2, 1,sharex=True,figsize=(8,10))
        
        Hp = ax2.bar(Hind, compHarray, yerr=compHerrarray,color='red')
        Sp = ax2.bar(Sind, compSarray, yerr=compSerrarray,color='blue')
        Gp = ax1.bar(Gind, compGarray, yerr=compGerrarray,color='green')
        ax1.set_ylabel('Complexation Absolute Energies (kcal/mol)')
        ax1.set_title('Absolute Complexation Energies')
        plt.xticks(Sind, tuple(patharray))
        plt.legend((Hp,Sp,Gp), (u'ΔHᶜᵒᵐᵖ', u'TΔSᶜᵒᵐᵖ',u'ΔGᶜᵒᵐᵖ'))
        fig2.savefig(annihilator.outputpath+'CompAbsBarPlot.png')
    
       
    if annihilator.binding==True:
        relHarray=-np.array(solvHarray)+np.array(compHarray)
        relHerrarray=np.sqrt(np.square(np.array(solvHerrarray))+np.square(np.array(compHerrarray)))
        relSarray=-np.array(solvSarray)+np.array(compSarray)
        relSerrarray=np.sqrt(np.square(np.array(solvSerrarray))+np.square(np.array(compSerrarray)))
        relGarray=-np.array(solvGarray)+np.array(compGarray)
        relGerrarray=np.sqrt(np.square(np.array(solvGerrarray))+np.square(np.array(compGerrarray)))

        Hind = np.arange(0,3*len(patharray)-2,x)    # the x locations for the groups
        Sind = [i+1 for i in Hind]
        Gind = [i+2 for i in Hind] # plot _TDeltaS instead of DeltaS
        fig3, (ax1, ax2) = plt.subplots(2, 1, sharex=True,figsize=(8,10))
        Hp = ax2.bar(Hind, relHarray, yerr=relHerrarray,color='red')
        Sp = ax2.bar(Sind, relSarray, yerr=relSerrarray,color='blue')
        Gp = ax1.bar(Gind, relGarray, yerr=relGerrarray,color='green')
        ax1.set_ylabel('Binding Absolute Energies (kcal/mol)')
        ax1.set_title('Binding Energies')
        plt.xticks(Gind, tuple(patharray))
        plt.legend((Hp, Sp,Gp), (u'ΔHᵇᶦⁿᵈ', u'TΔSᵇᶦⁿᵈ',u'ΔGᵇᶦⁿᵈ'))
        fig3.savefig(annihilator.outputpath+'BindBarPlot.png')

def PlotHeatmap(annihilator,matrix,xaxis,yaxis):
    shp=matrix.shape
    if len(shp)==2:
        fig, ax = plt.subplots()
        im = ax.imshow(matrix)
        figfname='FreeEnergyMatrix.png'
        cbar = ax.figure.colorbar(im, ax=ax)
        cbar.ax.set_ylabel('', rotation=-90, va="bottom")
        ax.set_xticks(numpy.arange(len(xaxis)))
        ax.set_yticks(numpy.arange(len(yaxis)))
        ax.set_xticklabels(xangles)
        ax.set_yticklabels(yangles)
        ax.set_ylim(len(matrix),-0.5, -0.5) 
        # Loop over data dimensions and create text annotations.
        for i in range(len(yaxis)):
            for j in range(len(xaxis)):
                energyvalue=str(round(matrix[i,j]))
                text = ax.text(j, i, energyvalue,ha="center", va="center", color="w")
    


        fig.tight_layout()
        plt.show() 
        fig.savefig(figfname)


