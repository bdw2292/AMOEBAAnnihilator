import submitjobs
import keyfilemodifications as keymods
import os
import shutil
import boxsetup as box
import terminate as term
import sys
import time
import numpy as np
import mutation as mutate


def ExecuteProductionDynamics(annihilator):
   jobtolog={}
   jobtojobpath={}
   annihilator.WriteToLog('Running production dynamics',prin=True)
   os.chdir(annihilator.outputpath+annihilator.simfoldname)
   subfolders=os.listdir(os.getcwd())
   for i in range(len(annihilator.proddynoutfilepath)):
       outputfilepath=annihilator.proddynoutfilepath[i]
       path,tail=os.path.split(outputfilepath)
       os.chdir(path)
       cmdstr=ProductionDynamicsCommand(annihilator,annihilator.proddynwaterboxfilename,annihilator.proddynwaterboxkeyfilename,annihilator.proddynsteps,annihilator.proddynensem,outputfilepath)
       terminate=term.CheckFileTermination(annihilator,outputfilepath)
       if terminate==False and os.path.exists(path+annihilator.proddynarcwaterboxfilename): 
           shutil.copy(outputfilepath,path+fold+'_PreviousOutput.out')
           stepstaken=CheckLastNumberDynamicStepsCompleted(annihilator,outputfilepath) # make copy of old output file and save
           newstepstotake=int(annihilator.proddynsteps)-stepstaken
           cmdstr=ProductionDynamicsCommand(annihilator,annihilator.proddynarcwaterboxfilename,annihilator.proddynwaterboxkeyfilename,newstepstotake,annihilator.proddynensem,outputfilepath)
       if terminate==False:
           jobtolog[cmdstr]=annihilator.outputpath+annihilator.proddynjobsfilename
           jobtojobpath[cmdstr]=path

       os.chdir('..')
   os.chdir('..')
   submitjobs.SubmitJobs(annihilator,jobtolog,jobtojobpath,annihilator.outputpath+annihilator.proddynjobsfilename)


def ProductionDynamicsCommand(annihilator,inputxyzname,keyfile,steps,ensemble,outputpath):
    cmdstr=annihilator.truedynamicpath+' '+inputxyzname+' -k ' + keyfile + ' '+str(steps)+' '+str(annihilator.proddyntimestep)+' '+str(annihilator.proddynwritefreq)+' '+str(ensemble)+' '+str(annihilator.roomtemp)+' '+'N'+' > '+outputpath
    return cmdstr
   


def SetupProductionDynamics(annihilator):
   annihilator.WriteToLog('Setting up dynamics for '+annihilator.simfoldname,prin=True)
   if not os.path.isdir(annihilator.outputpath+annihilator.simfoldname):
       os.mkdir(annihilator.outputpath+annihilator.simfoldname)
   os.chdir(annihilator.outputpath+annihilator.simfoldname)
   if len(annihilator.mutlambdascheme)==0:
       mut=False
   else:
       mut=True
   for i in range(len(annihilator.lambdafolderlist)):
       fold=annihilator.lambdafolderlist[i]
       if mut==False:
           elelamb=annihilator.estatlambdascheme[i]
           vdwlamb=annihilator.vdwlambdascheme[i]
       else:
           mutlambda=annihilator.mutlambdascheme[i]

       if annihilator.complexation==True and annihilator.dontrestrainreceptorligand==False:
           reslambda=annihilator.restlambdascheme[i]
       else:
           reslambda=0
       if not os.path.isdir(fold):
           os.mkdir(fold)
       os.chdir(fold)
       newfoldpath=os.getcwd()+'/'
       newtempkeyfile=annihilator.proddynwaterboxfilename.replace('.xyz','.key')
       if mut==False:
           shutil.copyfile(annihilator.outputpath+annihilator.lambdakeyfilename,newfoldpath+newtempkeyfile)
       else:
           arrayoflinearrays=mutate.MutateAllParameters(annihilator,annihilator.bgnlinetoendline,mutlambda)
           mutate.GenerateKeyFile(annihilator,arrayoflinearrays,newfoldpath+newtempkeyfile)
           keymods.InsertKeyfileHeader(annihilator,newfoldpath+newtempkeyfile)
    
       outputboxname=newfoldpath+annihilator.proddynwaterboxfilename
       shutil.copyfile(annihilator.outputpath+annihilator.proddynwaterboxfilename,outputboxname)
       if mut==False:
           ModifyLambdaKeywords(annihilator,newfoldpath,newtempkeyfile,elelamb,vdwlamb,reslambda)

       os.chdir('..')
   os.chdir('..')


def ModifyLambdaKeywords(annihilator,newfoldpath,newtempkeyfile,elelamb,vdwlamb,reslambda):
    temp=open(newfoldpath+newtempkeyfile,'r')
    results=temp.readlines()
    temp.close()
    newkeyfile=open(newfoldpath+newtempkeyfile,'w')
    for line in results:
        linesplit=line.split()
        if "ele-lambda" in line:
            newline=line.replace('\n'," "+str(elelamb)+'\n')
            newkeyfile.write(newline)
        elif "vdw-lambda" in line:
            newline=line.replace('\n'," "+str(vdwlamb)+'\n')
            newkeyfile.write(newline)
        elif "restrain-groups 1 2" in line:
            group1index=linesplit[1]
            group2index=linesplit[2]
            constant=linesplit[3]
            teatherdist=linesplit[-1]
            constant=str(float(constant)*float(reslambda))
            if annihilator.flatbotrest==False:
                restrainstring='restrain-groups '+str(group1index)+' '+ str(group2index)+' '+str(constant)+' '+str(teatherdist)+' '+str(teatherdist)
            else:
                restrainstring='restrain-groups '+str(group1index)+' '+ str(group2index)+' '+str(constant)+' '+'0'+' '+str(teatherdist)

            
            newkeyfile.write(restrainstring+'\n')        
        elif 'restrain-distance' in line:
            index1=linesplit[1]
            index2=linesplit[2]
            constant=linesplit[3]
            constant=str(float(constant)*float(reslambda))
            lowerdist=linesplit[4]
            upperdist=linesplit[5]
            string='restrain-distance '+index1+' '+index2+' '+constant+' '+lowerdist+' '+upperdist
            newkeyfile.write(string+'\n')        

        elif 'restrain-angle' in line:
            index1=linesplit[1]
            index2=linesplit[2]
            index3=linesplit[3]
            constant=linesplit[4]
            constant=str(float(constant)*float(reslambda))
            lowerangle=linesplit[5]
            upperangle=linesplit[6]
            string='restrain-angle '+index1+' '+index2+' '+index3+' '+constant+' '+lowerangle+' '+upperangle
            newkeyfile.write(string+'\n')        

        elif 'restrain-torsion' in line:
            index1=linesplit[1]
            index2=linesplit[2]
            index3=linesplit[3]
            index4=linesplit[4]
            constant=linesplit[5]
            constant=str(float(constant)*float(reslambda))
            lowerangle=linesplit[6]
            upperangle=linesplit[7]
            string='restrain-torsion '+index1+' '+index2+' '+index3+' '+index4+' '+constant+' '+lowerangle+' '+upperangle
            newkeyfile.write(string+'\n')        

        elif 'multipole' in line:
            typeindex=linesplit[1]
            if '-' in typeindex:
                xframeindex=linesplit[2]
                yframeindex=linesplit[3]
                chg=float(linesplit[4])
                chg=chg*elelamb
                newline='multipole '+xframeindex+' '+yframeindex+' '+str(chg)+'\n'
                newkeyfile.write(newline) 
            else:
                newkeyfile.write(line)

        else:
            newkeyfile.write(line)
        newkeyfile.flush()
        os.fsync(newkeyfile.fileno())
        sys.stdout.flush()

    newkeyfile.close()

 
def CheckLastNumberDynamicStepsCompleted(annihilator,outputfilepath):
   steps=None
   if os.path.isfile(outputfilepath):
       temp=open(outputfilepath,'r')
       results=temp.readlines()
       temp.close()
       for line in results:
           if 'Instantaneous Values for Frame Saved at' in line:
               linesplit=line.split()
               if linesplit[6].isdigit():
                   steps=int(linesplit[6])
               else:
                   steps=int(linesplit[7])
       return steps
   return steps

def DetermineIonIndicesToModifyCharge(annihilator):
    chg=int(annihilator.systemcharge)
    ionindexes=[]
    charge=0
    if chg>0:
        Clions=np.absolute(chg)
        iontypenumber=annihilator.elementsymtotinktype['Cl']
        ionindexes=GrabIonIndexes(annihilator,Clions,annihilator.proddynwaterboxfilename,iontypenumber)
        charge=-1
    else:
        Kions=np.absolute(chg)
        iontypenumber=annihilator.elementsymtotinktype['K']
        ionindexes=GrabIonIndexes(annihilator,Kions,annihilator.proddynwaterboxfilename,iontypenumber)
        charge=1
    return ionindexes,charge


def GrabIonIndexes(annihilator,ionnumber,boxfilename,iontypenumber):
    ionindexes=[]
    count=0     
    temp=open(boxfilename,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if count==ionnumber:
            break
        linesplit=line.split()
        if len(linesplit)!=1 and '90.000000' not in line: # not line containing number of atoms
            index=int(linesplit[0])
            typenum=int(linesplit[5]) 
            if typenum==int(iontypenumber):
                ionindexes.append(index)
                count+=1 
    return ionindexes


def ProductionDynamicsProtocol(annihilator):
    if len(annihilator.mutlambdascheme)==0:
        if not os.path.isfile(annihilator.outputpath+annihilator.lambdakeyfilename):
            shutil.copyfile(annihilator.outputpath+annihilator.configkeyfilename,annihilator.outputpath+annihilator.lambdakeyfilename)
            string='ligand'+' '
            firstligidx=str(annihilator.ligandindices[0])
            lastligidx=str(annihilator.ligandindices[-1])
            string+='-'+firstligidx+' '+lastligidx+'\n'
            keymods.AddKeyWord(annihilator,annihilator.lambdakeyfilename,string)
            string='ele-lambda'+'\n'
            keymods.AddKeyWord(annihilator,annihilator.lambdakeyfilename,string)
            string='vdw-lambda'+'\n'
            keymods.AddKeyWord(annihilator,annihilator.lambdakeyfilename,string)

    ionindexes,charge=DetermineIonIndicesToModifyCharge(annihilator)
    keymods.AddMultipoleDefinitionsForIonIndices(annihilator,ionindexes,charge,annihilator.lambdakeyfilename)
    if term.CheckFilesTermination(annihilator,annihilator.proddynoutfilepath)==False:
        SetupProductionDynamics(annihilator)
        ExecuteProductionDynamics(annihilator)
    messages=[]
    while term.CheckFilesTermination(annihilator,annihilator.proddynoutfilepath)==False:
        msg='System dynamics is not complete '
        if msg not in messages:
            annihilator.WriteToLog(msg,prin=True)
            messages.append(msg)
        time.sleep(annihilator.waitingtime)
    annihilator.WriteToLog('System dynamics is complete ',prin=True)

