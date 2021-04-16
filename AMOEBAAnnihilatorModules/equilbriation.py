import productiondynamics as prod 
import submitjobs
import tables
import terminate as term
import time
import numpy as np
import shutil
import keyfilemodifications as keymods
import os
import restraints
import sys


def AverageBoxSizeFromNPTArc(annihilator,arcpath,firstframe,lastframe,framestep,totalnumberframes):
    firstline=True
    framecount=0
    framestoextract=np.arange(firstframe,lastframe+1,framestep)
    framearray=[]
    extractingframe=False
    aaxisarray=[]
    baxisarray=[]
    caxisarray=[]
    with open(arcpath) as infile:
        for line in infile:
            if firstline==True:
                firstlinesplit=line.split()
                framestring=firstlinesplit[0].lstrip().rstrip()
                firstline=False
            linesplit=line.split()
            if '90.000000' in line and framecount in framestoextract:
                aaxis=float(linesplit[0])
                baxis=float(linesplit[1])
                caxis=float(linesplit[2])
                aaxisarray.append(aaxis)
                baxisarray.append(baxis)
                caxisarray.append(caxis)
            if framestring in line and (len(linesplit)==1):
                extractingframe=False
                framecount+=1
    aaxisaverage=np.mean(np.array(aaxisarray))
    baxisaverage=np.mean(np.array(baxisarray))
    caxisaverage=np.mean(np.array(caxisarray))
  
    annihilator.tabledict['Average Box Size']=aaxisaverage
    tables.WriteTableUpdateToLog(annihilator)
    return aaxisaverage,baxisaverage,caxisaverage



def EquilbriateDynamicCommand(annihilator,steps,ensemble,temp,outputfilename,NPT=False):
    if NPT==False:
        cmdstr=annihilator.truedynamicpath+' '+ annihilator.outputpath+annihilator.equilwaterboxfilename+' '+ '-k'+' '+ annihilator.outputpath+annihilator.configkeyfilename+' '+str(steps)+' '+ str(annihilator.equiltimestep)+' '+ str(annihilator.equilwritefreq)+' '+str(ensemble)+' '+str(temp)+' '+ 'N > '+outputfilename  
    else:
        cmdstr=annihilator.truedynamicpath+' '+ annihilator.outputpath+annihilator.equilwaterboxfilename+' '+ '-k'+' '+ annihilator.outputpath+annihilator.configkeyfilename+' '+str(steps)+' '+ str(annihilator.equiltimestep)+' '+ str(annihilator.equilwritefreq)+' '+str(ensemble)+' '+str(temp)+' '+str(annihilator.pressure)+' '+ 'N > '+outputfilename  

    return cmdstr         

def ExecuteEquilibriation(annihilator):
    jobtolog={}
    jobtoidx={}
    jobtojobpath={}
    jobtorestconstant={}
    for idx in range(len(annihilator.equilibriatescheme)):
        temp=annihilator.equilibriatescheme[idx]
        restrainpositionconstant=annihilator.equilibriaterestscheme[idx]
        cmdstr=EquilbriateDynamicCommand(annihilator,annihilator.equilstepsNVT,annihilator.NVTensem,temp,annihilator.equiloutputarray[idx])
        jobtolog[cmdstr]=annihilator.equiloutputarray[idx]
        jobtoidx[cmdstr]=idx
        jobtojobpath[cmdstr]=annihilator.outputpath
        jobtorestconstant[cmdstr]=restrainpositionconstant
    temp=annihilator.equilibriatescheme[-1]
    cmdstr=EquilbriateDynamicCommand(annihilator,annihilator.equilstepsNPT,annihilator.NPTensem,temp,annihilator.equiloutputarray[-1],NPT=True)
    jobtolog[cmdstr]=annihilator.equiloutputarray[-1]
    jobtoidx[cmdstr]=-1
    jobtojobpath[cmdstr]=annihilator.outputpath
    jobtorestconstant[cmdstr]=0

    for job in jobtolog.keys():
        newjobtolog={}
        newjobtojobpath={}
        newjobtolog[job]=jobtolog[job] # only submit one equilibriate job at a time (increasing temp...)
        newjobtojobpath[job]=jobtojobpath[job]
        restrainpositionconstant=jobtorestconstant[job]
        outarray=[annihilator.equiloutputarray[jobtoidx[job]]]
        finished=term.CheckFilesTermination(annihilator,outarray)
        if finished==False:
            if annihilator.complexation==True:
                keymods.RemoveKeyWord(annihilator,annihilator.configkeyfilename,'restrain-position')
                if restrainpositionconstant!=0 and annihilator.dontrestrainreceptorligand==False: # dont add for NPT only NVT
                    resposstring='restrain-position -'+str(1)+' '+str(annihilator.totalatomnumberxyzfilename-len(annihilator.ligandindices))+' '+str(restrainpositionconstant)+' '+str(annihilator.equilrestrainsphereradius)+'\n'
                    keymods.AddKeyWord(annihilator,annihilator.outputpath+annihilator.configkeyfilename,resposstring)

            submitjobs.SubmitJobs(annihilator,newjobtolog,newjobtojobpath,annihilator.outputpath+annihilator.equiljobsfilename)
            messages=[]
            while finished==False:
                msg='System equilibriation is not complete '
                if msg not in messages:
                    annihilator.WriteToLog(msg,prin=True)
                    messages.append(msg)
                time.sleep(annihilator.waitingtime)
                finished=term.CheckFilesTermination(annihilator,outarray)
            annihilator.WriteToLog('Single system equilibriation job is complete ',prin=True)

    annihilator.WriteToLog('System equilibriation is complete ',prin=True)


def ExtractLastTinkerFrame(annihilator,arcfilename):
    annihilator.WriteToLog('Extracting the last frame of '+arcfilename,prin=True)
    framenum=int(os.path.getsize(annihilator.outputpath+arcfilename)/os.path.getsize(annihilator.outputpath+arcfilename.replace('.arc','.xyz')))
    ExtractTinkerFrames(annihilator,annihilator.outputpath+arcfilename,framenum,framenum,1,framenum)

def ExtractTinkerFrames(annihilator,arcpath,firstframe,lastframe,framestep,totalnumberframes):
    firstline=True
    framecount=0
    framestoextract=np.arange(firstframe,lastframe+1,framestep)
    framearray=[]
    extractingframe=False
    with open(arcpath) as infile:
        for line in infile:
            if firstline==True:
                firstlinesplit=line.split()
                framestring=firstlinesplit[0].lstrip().rstrip()
                firstline=False
           
            linesplit=line.split()
            if framestring in line and (len(linesplit)==1):
                extractingframe=False
                framecount+=1
                if len(framearray)!=0:
                    numberofzeroes=len(str(totalnumberframes))-len(str(framecount))
                    zerostring=''
                    for i in range(numberofzeroes):
                        zerostring+='0'
                    framename=arcpath.replace('.arc','.'+zerostring+str(framecount))
                    temp=open(framename,'w')
                    for saveline in framearray:
                        temp.write(saveline)
                    temp.close()
                    framearray=[]
                if framecount in framestoextract:
                    extractingframe=True
            if(extractingframe):
                framearray.append(line)
            if len(framearray)!=0: # for last frame
                numberofzeroes=len(str(totalnumberframes))-len(str(framecount))
                zerostring=''
                for i in range(numberofzeroes):
                    zerostring+='0'
                framename=arcpath.replace('.arc','.'+zerostring+str(framecount))
                temp=open(framename,'w')
                for saveline in framearray:
                    temp.write(saveline)
                temp.close()
    os.rename(framename,annihilator.proddynwaterboxfilename)
    return
  
def EquilibriationProtocol(annihilator):
    if annihilator.restrainatomgroup1==None and annihilator.restrainatomgroup2==None and annihilator.complexation==True and annihilator.dontrestrainreceptorligand==False:
        restraints.ComputeIdealGroupRestraints(annihilator,annihilator.minwaterboxfilename)
    if annihilator.complexation==True and annihilator.dontrestrainreceptorligand==False: 
        dist=restraints.AverageCOMGroups(annihilator,annihilator.minwaterboxfilename)
        annihilator.restraintdistance=dist
        if annihilator.dontrestrainreceptorligand==False:
            restraints.AddHarmonicRestrainGroupTermsToKeyFile(annihilator,annihilator.outputpath+annihilator.configkeyfilename,dist)

    shutil.copy(annihilator.minwaterboxfilename,annihilator.equilwaterboxfilename)
    ExecuteEquilibriation(annihilator)
    firstframe=annihilator.equilframenum-annihilator.equilframenumNPT
    lastframe=annihilator.equilframenum
    annihilator.aaxis,annihilator.baxis,annihilator.caxis=AverageBoxSizeFromNPTArc(annihilator,annihilator.outputpath+annihilator.equilarcwaterboxfilename,firstframe,lastframe,1,annihilator.equilframenumNPT)
    keymods.RemoveKeyWord(annihilator,annihilator.configkeyfilename,'axis')
    keymods.AddKeyWord(annihilator,annihilator.configkeyfilename,'aaxis'+' '+str(annihilator.aaxis)+'\n')
    if not os.path.isfile(annihilator.outputpath+annihilator.proddynwaterboxfilename):
        ExtractLastTinkerFrame(annihilator,annihilator.equilarcwaterboxfilename)
    keymods.RemoveKeyWord(annihilator,annihilator.configkeyfilename,'restrain')
    keymods.RemoveKeyWord(annihilator,annihilator.configkeyfilename,'group')
    if annihilator.complexation==True and annihilator.proddyngrprests==True:
        equildist=restraints.AverageCOMGroups(annihilator,annihilator.equilarcwaterboxfilename)
        annihilator.restraintdistance=equildist
        if annihilator.dontrestrainreceptorligand==False:
            restraints.AddHarmonicRestrainGroupTermsToKeyFile(annihilator,annihilator.outputpath+annihilator.configkeyfilename,equildist)
            restraints.GroupRestraintFreeEnergyFix(annihilator)
    elif annihilator.complexation==True and annihilator.proddyngrprests==False:
        if annihilator.dontrestrainreceptorligand==False:
            restraints.ComputeIdealRestraints(annihilator,annihilator.equilarcwaterboxfilename)

    if not os.path.isfile(annihilator.outputpath+annihilator.proddynwaterboxfilenamepymol):
        annihilator.PymolReadableFile(annihilator.proddynwaterboxfilename,annihilator.proddynwaterboxfilenamepymol)

    if annihilator.equilonly==True:
        sys.exit()
