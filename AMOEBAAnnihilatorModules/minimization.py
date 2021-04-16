import submitjobs
import os
import terminate as term
import time
import keyfilemodifications as keymods
import restraints as res
import sys

def ExecuteLooseMinimization(annihilator):

    cmd=MinimizeCommand(annihilator,annihilator.outputpath+annihilator.waterboxfilename,annihilator.outputpath+annihilator.configkeyfilename,annihilator.loosemincriteria,annihilator.looseminoutput)
    jobtolog={cmd:annihilator.outputpath+annihilator.looseminjobsfilename}
    jobtojobpath={cmd:annihilator.outputpath}
    submitjobs.SubmitJobs(annihilator,jobtolog,jobtojobpath,annihilator.outputpath+annihilator.looseminjobsfilename)

def ExecuteTightMinimization(annihilator):
    cmd=MinimizeCommand(annihilator,annihilator.outputpath+annihilator.waterboxfilename+'_2',annihilator.outputpath+annihilator.configkeyfilename,annihilator.tightmincriteria,annihilator.tightminoutput)
    jobtolog={cmd:annihilator.outputpath+annihilator.tightminjobsfilename}
    jobtojobpath={cmd:annihilator.outputpath}
    submitjobs.SubmitJobs(annihilator,jobtolog,jobtojobpath,annihilator.outputpath+annihilator.tightminjobsfilename)
           
def MinimizeCommand(annihilator,xyzfilename,keyfilename,gradrms,outputfilename):
    cmd=annihilator.minimizepath+' '+xyzfilename+' -k '+keyfilename+' '+str(gradrms)+' '+'> '+annihilator.outputpath+outputfilename
    annihilator.WriteToLog('Calling '+cmd)
    return cmd

def ExpensiveMinimizationProtocol(annihilator):

    if annihilator.restrainatomsduringminimization:
        resposstring='restrain-position -'+str(1)+' '+str(annihilator.totalatomnumberxyzfilename)+' '+str(annihilator.restrainpositionconstant)+' '+'0'+'\n'
        keymods.AddKeyWords(annihilator,annihilator.outputpath+annihilator.configkeyfilename,resposstring)

    if not os.path.isfile(annihilator.outputpath+annihilator.minwaterboxfilename): # currently, only minimize and equilibraite default prm set
        annihilator.WriteToLog('Minimizing system ',prin=True)
        string='polarizeterm none'+'\n'
        keymods.AddKeyWord(annihilator,annihilator.configkeyfilename,string)
        ExecuteLooseMinimization(annihilator)
    messages=[]
    while term.CheckFileTermination(annihilator,annihilator.outputpath+annihilator.looseminoutput)==False:
        msg='Loose minimization is not complete '
        if msg not in messages:
            annihilator.WriteToLog(msg,prin=True)
            messages.append(msg)
        time.sleep(annihilator.waitingtime)
    if not os.path.isfile(annihilator.outputpath+annihilator.minwaterboxfilename): 
        keymods.RemoveKeyWord(annihilator,annihilator.configkeyfilename,'polarizeterm none') 
        ExecuteTightMinimization(annihilator)
    messages=[]
    while term.CheckFileTermination(annihilator,annihilator.outputpath+annihilator.tightminoutput)==False:
        msg='Tight minimization is not complete '
        if msg not in messages:
            annihilator.WriteToLog(msg,prin=True)
            messages.append(msg)
        time.sleep(annihilator.waitingtime)
    if not os.path.isfile(annihilator.outputpath+annihilator.minwaterboxfilename):
        if os.path.exists(annihilator.outputpath+annihilator.waterboxfilename+'_3'):
            os.rename(annihilator.outputpath+annihilator.waterboxfilename+'_3',annihilator.outputpath+annihilator.minwaterboxfilename)
            #os.remove(annihilator.outputpath+annihilator.waterboxname+'_2')

    if annihilator.restrainatomsduringminimization:
        keymods.RemoveKeyWord(annihilator,annihilator.outputpath+annihilator.configkeyfilename,'restrain-position')

def CheapMinimizationProtocol(annihilator):
    if annihilator.restrainatomsduringminimization:
        resposstring='restrain-position -'+str(1)+' '+str(annihilator.totalatomnumberxyzfilename)+' '+str(annihilator.restrainpositionconstant)+' '+'0'+'\n'
        keymods.AddKeyWord(annihilator,annihilator.outputpath+annihilator.configkeyfilename,resposstring)

    if not os.path.isfile(annihilator.outputpath+annihilator.minwaterboxfilename): # currently, only minimize and equilibraite default prm set
        annihilator.WriteToLog('Minimizing system ',prin=True)
        ExecuteLooseMinimization(annihilator)
    messages=[]
    while term.CheckFileTermination(annihilator,annihilator.outputpath+annihilator.looseminoutput)==False:
        msg='Loose minimization is not complete '
        if msg not in messages:
            annihilator.WriteToLog(msg,prin=True)
            messages.append(msg)
        time.sleep(annihilator.waitingtime)
    if not os.path.isfile(annihilator.outputpath+annihilator.minwaterboxfilename):
        if os.path.exists(annihilator.outputpath+annihilator.waterboxfilename+'_2'):
            os.rename(annihilator.outputpath+annihilator.waterboxfilename+'_2',annihilator.outputpath+annihilator.minwaterboxfilename)

    if not os.path.isfile(annihilator.outputpath+annihilator.minwaterboxfilenamepymol):
        annihilator.PymolReadableFile(annihilator.minwaterboxfilename,annihilator.minwaterboxfilenamepymol)

    if annihilator.restrainatomsduringminimization:
        keymods.RemoveKeyWord(annihilator,annihilator.outputpath+annihilator.configkeyfilename,'restrain-position')
    if annihilator.minonly==True:
        sys.exit()
