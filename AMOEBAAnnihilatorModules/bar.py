import submitjobs
import tables
import terminate as term
import os
import sys
import time
import numpy as np
import csv
import restraints

def ExecuteBARSecondOption(annihilator,joblistname):
    jobtolog={}
    jobtojobpath={}
    annihilator.WriteToLog('Computing free energy from .bar files ',prin=True)
    os.chdir(annihilator.outputpath+annihilator.simfoldname)
    for i in range(len(annihilator.thermooutputfilepath)): 
        outputfilepath=annihilator.thermooutputfilepath[i]
        barfilepath=annihilator.barfilepath[i]
        path=os.path.split(outputfilepath)[0]
        maxframenum=int(os.path.getsize(path+'/'+annihilator.proddynarcwaterboxfilename)/os.path.getsize(path+'/'+annihilator.proddynwaterboxfilename))
        firstframenum=1
        cmdstr=BARSecondOptionCommand(annihilator,barfilepath,firstframenum,maxframenum,outputfilepath)
        terminate=term.CheckFileTermination(annihilator,outputfilepath)
        if terminate==False:
            jobtolog[cmdstr]=annihilator.outputpath+joblistname
            jobtojobpath[cmdstr]=path
    os.chdir('..')
    submitjobs.SubmitJobs(annihilator,jobtolog,jobtojobpath,annihilator.outputpath+joblistname)

def BARSecondOptionCommand(annihilator,barfilepath,firstframenum,maxframenum,outputname):
    cmdstr=annihilator.truebarpath+' '+'2'+' '+barfilepath+' '+str(firstframenum)+' '+str(maxframenum)+' '+'1'+' '+str(firstframenum)+' '+str(maxframenum)+' '+'1'+' > '+outputname
    return cmdstr

def ExecuteBAR(annihilator):
    jobtolog={}
    jobtojobpath={}
    annihilator.WriteToLog('Running BAR for '+annihilator.simfoldname,prin=True)
    os.chdir(annihilator.outputpath+annihilator.simfoldname)
    for i in range(len(annihilator.baroutputfilepath)): 
        outputfilepath=annihilator.baroutputfilepath[i]
        path=os.path.split(outputfilepath)[0]
        secondarcpath=annihilator.secondarcpaths[i]
        firstarcpath=annihilator.firstarcpaths[i]
        cmdstr=BARCommand(annihilator,secondarcpath,firstarcpath,outputfilepath)
        terminate=term.CheckFileTermination(annihilator,outputfilepath)
        if terminate==False:
            jobtolog[cmdstr]=annihilator.outputpath+annihilator.barjobsfilename
            jobtojobpath[cmdstr]=path
    os.chdir('..')
    submitjobs.SubmitJobs(annihilator,jobtolog,jobtojobpath,annihilator.outputpath+annihilator.barjobsfilename)

def BARCommand(annihilator,secondarcpath,firstarcpath,outputfilepath):
    cmdstr=annihilator.truebarpath+' '+'1'+' '+secondarcpath+' '+str(annihilator.equilibriatescheme[-1])+' '+firstarcpath+' '+str(annihilator.equilibriatescheme[-1])+' > '+outputfilepath
    return cmdstr

def SumTheFreeEnergyStepsFromBAR(annihilator):
    firststate=False
    laststate=False
    annihilator.freeenergy=0
    annihilator.freeenergylist=[]
    annihilator.freeenergyerrorlist=[]
    annihilator.freeenergyfwd=0
    annihilator.freeenergylistfwd=[]
    annihilator.freeenergyerrorlistfwd=[]
    annihilator.freeenergybwd=0
    annihilator.freeenergylistbwd=[]
    annihilator.freeenergyerrorlistbwd=[]

    annihilator.freeenergyviabariter=0
    annihilator.freeenergylistviabariter=[]
    annihilator.freeenergyerrorlistviabariter=[]
    annihilator.freeenergyviabootstrap=0
    annihilator.freeenergylistviabootstrap=[]
    annihilator.freeenergyerrorlistviabootstrap=[]

    annihilator.enthalpy=0
    annihilator.enthalpylist=[]
    annihilator.enthalpyerrorlist=[]
    annihilator.enthalpyerrorlisttotal=[]
    annihilator.entropy=0
    annihilator.entropylist=[]
    annihilator.entropyerrorlist=[]
    annihilator.entropyerrorlisttotal=[]
    os.chdir(annihilator.outputpath+annihilator.simfoldname)
    for i in range(len(annihilator.thermooutputfilepath)):
        outputfilepath=annihilator.thermooutputfilepath[i]
        temp=open(outputfilepath,'r')
        foundfreeenergyline=False
        for line in temp.readlines():
            if 'Free Energy via BAR Iteration' in line:
                foundfreeenergyline=True
                linesplit=line.split()
                freeenergy=float(linesplit[5])
                index=line.find('+/-')
                newstring=line[index+3:]
                newstringlinesplit=newstring.split()
                freenergyerror=float(newstringlinesplit[0])
                annihilator.freeenergy+=freeenergy
                annihilator.freeenergylist.append(freeenergy)
                annihilator.freeenergyerrorlist.append(freenergyerror)
                annihilator.freeenergyviabariter+=freeenergy
                annihilator.freeenergylistviabariter.append(freeenergy)
                annihilator.freeenergyerrorlistviabariter.append(freenergyerror)

            elif 'Free Energy via BAR Bootstrap' in line: 
                foundfreeenergyline=True
                linesplit=line.split()
                freeenergy=float(linesplit[5])
                index=line.find('+/-')
                newstring=line[index+3:]
                newstringlinesplit=newstring.split()
                freenergyerror=float(newstringlinesplit[0])
                annihilator.freeenergyviabootstrap+=freeenergy
                annihilator.freeenergylistviabootstrap.append(freeenergy)
                annihilator.freeenergyerrorlistviabootstrap.append(freenergyerror)

            elif 'Enthalpy via BAR Estimate' in line:
                linesplit=line.split()
                enthalpy=float(linesplit[4])
                index=line.find('+/-')
                newstring=line[index+3:]
                newstringlinesplit=newstring.split()
                enthalpyerror=float(newstringlinesplit[0])
                annihilator.enthalpy+=enthalpy
                annihilator.enthalpylist.append(enthalpy)
                annihilator.enthalpyerrorlisttotal.append(enthalpyerror)
            elif 'Entropy via BAR Estimate' in line:
                linesplit=line.split()
                entropy=float(linesplit[4])
                entropyerror=np.sqrt((enthalpyerror)**2+(freenergyerror)**2)
                annihilator.entropy+=entropy
                annihilator.entropylist.append(entropy)
                annihilator.entropyerrorlisttotal.append(entropyerror)
            elif 'Free Energy via Forward FEP' in line:
                linesplit=line.split()
                freeenergy=float(linesplit[5])
                index=line.find('+/-')
                newstring=line[index+3:]
                newstringlinesplit=newstring.split()
                freenergyerror=float(newstringlinesplit[0])
                annihilator.freeenergyfwd+=freeenergy
                annihilator.freeenergylistfwd.append(freeenergy)
                annihilator.freeenergyerrorlistfwd.append(freenergyerror)
            elif 'Free Energy via Backward FEP' in line:
                linesplit=line.split()
                freeenergy=float(linesplit[5])
                index=line.find('+/-')
                newstring=line[index+3:]
                newstringlinesplit=newstring.split()
                freenergyerror=float(newstringlinesplit[0])
                annihilator.freeenergybwd+=freeenergy
                annihilator.freeenergylistbwd.append(freeenergy)
                annihilator.freeenergyerrorlistbwd.append(freenergyerror)


        temp.close()


        annihilator.freeenergyerror=np.sqrt(np.sum(np.square(annihilator.freeenergyerrorlist)))
        annihilator.enthalpyerror=np.sqrt(np.sum(np.square(annihilator.enthalpyerrorlist)))
        annihilator.entropyerror=np.sqrt(np.sum(np.square(annihilator.entropyerrorlist)))
        if annihilator.solvation==True:
            annihilator.tabledict[u'ΔGˢᵒˡᵛ']=annihilator.freeenergy
            annihilator.tabledict[u'ΔGˢᵒˡᵛᵉʳʳ']=annihilator.freeenergyerror
            annihilator.tabledict[u'ΔHˢᵒˡᵛ']=annihilator.enthalpy
            annihilator.tabledict[u'ΔHˢᵒˡᵛᵉʳʳ']=annihilator.enthalpyerror
            annihilator.tabledict[u'ΔSˢᵒˡᵛ']=annihilator.entropy
            annihilator.tabledict[u'ΔSˢᵒˡᵛᵉʳʳ']=annihilator.entropyerror
        elif annihilator.complexation==True:
            annihilator.tabledict[u'ΔGᶜᵒᵐᵖᵘⁿᶜᵒʳʳ']=annihilator.freeenergy
            annihilator.tabledict[u'ΔGᶜᵒᵐᵖᶜᵒʳʳᵉʳʳ']=annihilator.freeenergyerror
            annihilator.tabledict[u'ΔHᶜᵒᵐᵖ']=annihilator.enthalpy
            annihilator.tabledict[u'ΔHᶜᵒᵐᵖᵉʳʳ']=annihilator.enthalpyerror
            annihilator.tabledict[u'ΔSᶜᵒᵐᵖ']=annihilator.entropy
            annihilator.tabledict[u'ΔSᶜᵒᵐᵖᵉʳʳ']=annihilator.entropyerror
            tables.WriteTableUpdateToLog(annihilator)

       
        tempname='BARResults.csv'
        with open(annihilator.outputpath+tempname, mode='w') as energy_file:
            energy_writer = csv.writer(energy_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            energy_writer.writerow(['Ele-Lambda']+annihilator.estatlambdascheme)
            energy_writer.writerow(['Vdw-Lambda']+annihilator.vdwlambdascheme)
            energy_writer.writerow(['Rest-Lambda']+annihilator.restlambdascheme)
            if annihilator.solvation==True:
                energy_writer.writerow([u'ΔGˢᵒˡᵛ']+annihilator.freeenergylist)
                energy_writer.writerow([u'ΔGˢᵒˡᵛᵉʳʳ']+annihilator.freeenergyerrorlist)
                energy_writer.writerow([u'ΔGˢᵒˡᵛᵇᵃʳᶦᵗᵉʳ']+annihilator.freeenergylistviabariter)
                energy_writer.writerow([u'ΔGˢᵒˡᵛᵉʳʳᵇᵃʳᶦᵗᵉʳ']+annihilator.freeenergyerrorlistviabariter)
                energy_writer.writerow([u'ΔGˢᵒˡᵛᵇᵒᵒᵗˢᵗʳᵃᵖ']+annihilator.freeenergylistviabootstrap)
                energy_writer.writerow([u'ΔGˢᵒˡᵛᵉʳʳᵇᵒᵒᵗˢᵗʳᵃᵖ']+annihilator.freeenergyerrorlistviabootstrap)
                energy_writer.writerow([u'ΔGˢᵒˡᵛᶠʷᵈ']+annihilator.freeenergylistfwd)
                energy_writer.writerow([u'ΔGˢᵒˡᵛᵉʳʳᶠʷᵈ']+annihilator.freeenergyerrorlistfwd)
                energy_writer.writerow([u'ΔGˢᵒˡᵛᵇʷᵈ']+annihilator.freeenergylistbwd)
                energy_writer.writerow([u'ΔGˢᵒˡᵛᵉʳʳᵇʷᵈ']+annihilator.freeenergyerrorlistbwd)
                energy_writer.writerow([u'ΔHˢᵒˡᵛ']+annihilator.enthalpylist)
                energy_writer.writerow([u'ΔHˢᵒˡᵛᵉʳʳ']+annihilator.enthalpyerrorlisttotal)
                energy_writer.writerow([u'ΔSˢᵒˡᵛ']+annihilator.entropylist)
                energy_writer.writerow([u'ΔSˢᵒˡᵛᵉʳʳ']+annihilator.entropyerrorlisttotal)
            elif annihilator.complexation==True:
                energy_writer.writerow([u'ΔGᶜᵒᵐᵖ']+annihilator.freeenergylist)
                energy_writer.writerow([u'ΔGᶜᵒᵐᵖᵉʳʳ']+annihilator.freeenergyerrorlist)
                energy_writer.writerow([u'ΔGᶜᵒᵐᵖᵇᵃʳᶦᵗᵉʳ']+annihilator.freeenergylistviabariter)
                energy_writer.writerow([u'ΔGᶜᵒᵐᵖᵉʳʳᵇᵃʳᶦᵗᵉʳ']+annihilator.freeenergyerrorlistviabariter)
                energy_writer.writerow([u'ΔGᶜᵒᵐᵖᵇᵒᵒᵗˢᵗʳᵃᵖ']+annihilator.freeenergylistviabootstrap)
                energy_writer.writerow([u'ΔGᶜᵒᵐᵖᵉʳʳᵇᵒᵒᵗˢᵗʳᵃᵖ']+annihilator.freeenergyerrorlistviabootstrap)
                energy_writer.writerow([u'ΔGᶜᵒᵐᵖᶠʷᵈ']+annihilator.freeenergylistfwd)
                energy_writer.writerow([u'ΔGᶜᵒᵐᵖᵉʳʳᶠʷᵈ']+annihilator.freeenergyerrorlistfwd)
                energy_writer.writerow([u'ΔGᶜᵒᵐᵖᵇʷᵈ']+annihilator.freeenergylistbwd)
                energy_writer.writerow([u'ΔGᶜᵒᵐᵖᵉʳʳᵇʷᵈ']+annihilator.freeenergyerrorlistbwd)
                energy_writer.writerow([u'ΔHᶜᵒᵐᵖ']+annihilator.enthalpylist)
                energy_writer.writerow([u'ΔHᶜᵒᵐᵖᵉʳʳ']+annihilator.enthalpyerrorlisttotal)
                energy_writer.writerow([u'ΔSᶜᵒᵐᵖ']+annihilator.entropylist)
                energy_writer.writerow([u'ΔSᶜᵒᵐᵖᵉʳʳ']+annihilator.entropyerrorlisttotal)

def BARProtocol(annihilator):
    if term.CheckFilesTermination(annihilator,annihilator.baroutputfilepath)==False:         
        ExecuteBAR(annihilator)
    
    annihilator.WriteToLog('BAR is running',prin=True)
    messages=[]
    while term.CheckFilesTermination(annihilator,annihilator.baroutputfilepath)==False:
        msg='BAR is not complete '
        if msg not in messages:
            annihilator.WriteToLog(msg,prin=True)
            messages.append(msg)
        time.sleep(annihilator.waitingtime)

    if term.CheckFilesTermination(annihilator,annihilator.thermooutputfilepath)==False:
        ExecuteBARSecondOption(annihilator,annihilator.freeenergyjobsfilename)
    messages=[]
    while term.CheckFilesTermination(annihilator,annihilator.thermooutputfilepath)==False:
        msg='System BAR option 2 is not complete '
        if msg not in messages:
            annihilator.WriteToLog(msg,prin=True)
            messages.append(msg)
        time.sleep(annihilator.waitingtime)
    SumTheFreeEnergyStepsFromBAR(annihilator)
    if annihilator.complexation==True: 
        annihilator.freeenergy=annihilator.freeenergy+annihilator.rescorrection
        annihilator.tabledict[u'ΔGᵃⁿᵃᶜᵒᵐᵖᶜᵒʳʳ']=annihilator.rescorrection
        annihilator.tabledict[u'ΔGᶜᵒᵐᵖᶜᵒʳʳ']=annihilator.freeenergy

