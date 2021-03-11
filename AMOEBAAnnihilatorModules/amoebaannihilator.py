import os
import sys
import annihilation as ann
import tables 
import boxsetup as box
import time
import pdbxyz
import restraints
import openbabel
import plots
import mutation as mutate

class Annihilator():

    def __init__(self,endstatekey=None,bgnstatekey=None,endstatexyz=None,bgnstatexyz=None,dontrestrainreceptorligand=False,mutlambdascheme=[],minonly=False,usegpu=True,truedynamicpath=None,truebarpath=None,equilonly=False,binding=False,proddyngrprests=False,equilrestrainsphereradius=2,restrainpositionconstant=1,ligandfilename=None,tightmincriteria=1,loosemincriteria=10,rescorrection=0,anglerestraintconstant=0.003046,pdbxyzpath='pdbxyz.x',distancerestraintconstant=10,minimizepath='minimize.x',tinkerdir=None,averageenergies=False,roomtemp=300,complexedproteinpdbname=None,uncomplexedproteinpdbname=None,addphysioions=True,equilibriatescheme=[50,100,150,200,300],equilibriaterestscheme=[5,2,1,.1,0],prmfilepath=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+ "/ParameterFiles/amoebabio18.prm",keyfilename=None,xyzfilename=None,externalapi=None,bashrcpath=None,restrainatomsduringminimization=True,restrainatomgroup1=None,restrainatomgroup2=None,ligandxyzfilename=None,receptorligandxyzfilename=None,xyzeditpath='xyzedit.x',lowerperf=7,upperperf=12,simpathlist=None,fixedboxsize=None,extendelelambdaoneside=[],estatlambdascheme=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1,1],extendelelambdazeroside=[],vdwlambdascheme=[0,.45,.52,.56,.58,.6,.62,.64,.67,.7,.75,.8,.85,.9,.95,1,1,1,1,1,1,1,1,1,1,1,1],extendvdwlambdaoneside=[],extendvdwlambdazeroside=[],restlambdascheme=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0],extendrestlambdaoneside=[],extendrestlambdazeroside=[],waitingtime=5,boxbufferlength=3,receptorcharge=0,ligandcharge=0,barpath='bar.x',dynamicpath='dynamic.x',barommpath='bar_omm.x',dynamicommpath='dynamic_omm.x',complexation=False,solvation=False,flatbotrest=True,logname='TINKER.log',equilwritefreq=100,proddynwritefreq=2,equiltimeNVT=5,equiltimeNPT=2,equiltimestep=2,proddyntimestep=2,proddyntime=5,pressure=1,NVTensem=2,NPTensem=4,vdwcutoff=12,ewaldcutoff=7,polareps=0.0001,barostatmethod='montecarlo',integrator='RESPA',thermostat='BUSSI',listofsaltcons='[KCl]=100'):
        self.endstatekey=endstatekey
        self.bgnstatekey=bgnstatekey
        self.endstatexyz=endstatexyz
        self.bgnstatexyz=bgnstatexyz
        self.dontrestrainreceptorligand=dontrestrainreceptorligand
        self.mutlambdascheme=mutlambdascheme
        self.minonly=minonly
        self.usegpu=usegpu
        self.truedynamicpath=truedynamicpath
        self.truebarpath=truebarpath
        self.equilonly=equilonly
        self.binding=binding
        self.proddyngrprests=proddyngrprests
        self.equilrestrainsphereradius=equilrestrainsphereradius
        self.equilibriaterestscheme=equilibriaterestscheme
        self.restrainpositionconstant=restrainpositionconstant
        self.ligandfilename=ligandfilename
        self.loosemincriteria=loosemincriteria
        self.tightmincriteria=tightmincriteria
        self.rescorrection=rescorrection
        self.anglerestraintconstant=anglerestraintconstant
        self.pdbxyzpath=pdbxyzpath
        self.distancerestraintconstant=distancerestraintconstant
        self.minimizepath=minimizepath
        self.tinkerdir=tinkerdir
        self.averageenergies=averageenergies
        self.roomtemp=roomtemp
        self.complexedproteinpdbname=complexedproteinpdbname
        self.uncomplexedproteinpdbname=uncomplexedproteinpdbname
        self.addphysioions=addphysioions
        self.equilibriatescheme=equilibriatescheme
        self.prmfilepath=prmfilepath
        self.keyfilename=keyfilename
        self.xyzfilename=xyzfilename
        self.bashrcpath=bashrcpath
        self.externalapi=externalapi
        self.restrainatomgroup1=restrainatomgroup1
        self.restrainatomgroup2=restrainatomgroup2
        self.restrainatomsduringminimization=restrainatomsduringminimization
        self.tabledict={}
        self.ligandxyzfilename=ligandxyzfilename
        self.receptorligandxyzfilename=receptorligandxyzfilename
        self.lowerperf=lowerperf
        self.upperperf=upperperf 
        self.simpathlist=simpathlist
        self.barommpath=barommpath
        self.dynamicommpath=dynamicommpath
        self.barpath=barpath
        self.dynamicpath=dynamicpath

        self.xyzeditpath=xyzeditpath
        self.simpath=os.getcwd()        
        self.simname=os.path.basename(self.simpath)

        self.complexation=complexation
        self.solvation=solvation
        self.flatbotrest=flatbotrest

        self.logname=logname

        self.equilwritefreq=equilwritefreq # in picoseconds
        self.proddynwritefreq=proddynwritefreq # in picoseconds
        self.equiltimeNVT=equiltimeNVT # in ns
        self.equiltimeNPT=equiltimeNPT # ns
        self.equiltimestep=equiltimestep # fs
        self.proddyntimestep=proddyntimestep # fs
        self.proddyntime=proddyntime # ns
        self.pressure=pressure 
        self.NVTensem=NVTensem
        self.NPTensem=NPTensem
        self.proddynensem=NVTensem
        self.vdwcutoff=vdwcutoff
        self.ewaldcutoff=ewaldcutoff
        self.polareps=polareps

        self.barostatmethod=barostatmethod
        self.integrator=integrator
        self.thermostat=thermostat

        self.listofsaltcons=listofsaltcons   # default physiological concentrations
        self.ligandcharge=ligandcharge
        self.receptorcharge=receptorcharge
        self.boxbufferlength=boxbufferlength # angstroms
        self.waitingtime=waitingtime # time in minutes, to wait to check if all jobs have terminated
        self.fixedboxsize=None
        self.extendelelambdaoneside=extendelelambdaoneside
        self.estatlambdascheme=estatlambdascheme
        self.extendelelambdazeroside=extendelelambdazeroside
        self.vdwlambdascheme=vdwlambdascheme
        self.extendvdwlambdaoneside=extendvdwlambdaoneside
        self.extendvdwlambdazeroside=extendvdwlambdazeroside
        self.restlambdascheme=restlambdascheme
        self.extendrestlambdaoneside=extendrestlambdaoneside
        self.extendrestlambdazeroside=extendrestlambdazeroside

        self.elementsymtocharge={'K':1,'Cl':-1,'Mg':2,'Li':1,'Na':1,'Rb':1,'Cs':1,'Be':2,'Ca':2,'Zn':2}
        self.elementsymtomass={'H':1.00794,'HN':1.00794,'HO':1.00794,'O':15.9994,'OH':15.9994,'N':14.0067,'C':12.0107,'CA':12.0107,'F':18.9984032,'Cl':35.453,'Cl-':35.453,'S':32.065,'P':30.9737,'Na':22.98976,'K':39.0983,'Ca':40.078,'Mg':24.305,'Mg+':24.305,'K+':39.0983}

        temp=open(os.getcwd()+r'/'+'AMOEBA.ini','r')
        results=temp.readlines()
        temp.close()
        for line in results:
            if '#' in line:
                continue
            if '=' in line:
                linesplit=line.split('=',1)
                a=linesplit[1].replace('\n','').rstrip().lstrip()
                commalist=a.split(',')
                commalist=[i.rstrip().lstrip() for i in commalist]
            if 'uncomplexedproteinpdbname' in line:
                self.uncomplexedproteinpdbname=a
            elif 'complexedproteinpdbname' in line:
                self.complexedproteinpdbname=a
            elif 'externalapi' in line:
                self.externalapi=a
            elif 'bashrcpath' in line:
                self.bashrcpath=a
            elif 'restrainatomgroup1' in line:
                self.restrainatomgroup1=[int(i) for i in commalist]
            elif 'restrainatomgroup2' in line:
                self.restrainatomgroup2=[int(i) for i in commalist]
            elif 'simpathlist' in line:
                self.simpathlist=commalist
                templist=[]
                for ele in self.simpathlist:
                    paths=ele.lstrip().rstrip().split()
                    temp=[]
                    for e in paths:
                        temp.append(e)
                    templist.append(temp)
                self.simpathlist=templist
            elif "flatbotrest" in line:
                self.flatbotrest=False
            elif "proddyngrprests" in line:
                self.proddyngrprests=True
            elif "dontrestrainreceptorligand" in line:
                self.dontrestrainreceptorligand=True
            elif "extendelelambdaoneside" in line:
                self.extendelelambdaoneside=commalist
            elif "extendvdwlambdaoneside" in line:
                self.extendvdwlambdaoneside=commalist
            elif "extendrestlambdaoneside" in line:
                self.extendrestlambdaoneside=commalist
            elif "extendprmsetoneside" in line:
                self.extendprmsetoneside=commalist
            elif "extendelelambdazeroside" in line:
                self.extendelelambdazeroside=commalist
            elif "extendvdwlambdazeroside" in line:
                self.extendvdwlambdazeroside=commalist
            elif "extendrestlambdazeroside" in line:
                self.extendrestlambdazeroside=commalist
            elif "proddynensem" in line:
                self.proddynensem = a
            elif ("keyfilename") in line:
                self.keyfilename = a
            elif ("ligandcharge") in line:
                self.ligandcharge = int(a)
            elif ("fixedboxsize") in line:
                self.fixedboxsize = float(a)
            elif ("tightmincriteria") in line:
                self.tightmincriteria = float(a)
            elif ("restrainpositionconstant") in line:
                self.restrainpositionconstant = float(a)
            elif ("loosemincriteria") in line:
                self.loosemincriteria = float(a)
            elif ("equilrestrainsphereradius") in line:
                self.equilrestrainsphereradius = float(a)
            elif ("barostatmethod") in line:
                self.barostatmethod = a
            elif ("receptorcharge") in line:
                self.receptorcharge = float(a)
            elif ("waitingtime") in line:
                self.waitingtime = float(a)
            elif ("listofsaltcons") in line:
                self.listofsaltcons = a
            elif ("ligandxyzfilename") in line and "receptor" not in line:
                self.ligandxyzfilename = a
            elif ("ligandfilename") in line:
                self.ligandfilename = a
            elif ("receptorligandxyzfilename") in line:
                self.receptorligandxyzfilename = a
            elif ("boxbufferlength") in line:
                self.boxbufferlength = float(a)
            elif ("integrator") in line:
                self.integrator = a
            elif ("thermostat") in line:
                self.thermostat = a
            elif ("vdwcutoff") in line:
                self.vdwcutoff = a
            elif ("ewaldcutoff") in line:
                self.ewaldcutoff = a
            elif ("polareps") in line:
                self.polareps = a
            elif ("bgnstatexyz") in line:
                self.bgnstatexyz = a
            elif ("endstatexyz") in line:
                self.endstatexyz = a
            elif ("bgnstatekey") in line:
                self.bgnstatekey = a
            elif ("endstatekey") in line:
                self.endstatekey = a
            elif ("distancerestraintconstant") in line:
                self.distancerestraintconstant= float(a)
            elif ("anglerestraintconstant") in line:
                self.anglerestraintconstant= float(a)
            elif ("equilibriatescheme") in line:
                self.equilibriatescheme=commalist
            elif ("mutlambdascheme") in line:
                self.mutlambdascheme=commalist
            elif ("restlambdascheme") in line:
                self.restlambdascheme=commalist
            elif ("vdwlambdascheme") in line:
                self.vdwlambdascheme=commalist
            elif ("estatlambdascheme") in line:
                self.estatlambdascheme=commalist
            elif ("equilwritefreq") in line:
                self.equilwritefreq= a
            elif ("equiltimestep") in line:
                self.equiltimestep= int(a)
            elif ("proddyntimestep") in line:
                self.proddyntimestep= int(a)
            elif ("equiltimeNPT") in line:
                self.equiltimeNPT= float(a)
            elif ("proddynwritefreq") in line:
                self.proddynwritefreq= a
            elif ("proddyntime") in line:
                self.proddyntime= float(a)
            elif ("complexation") in line:
                self.complexation=True
            elif ("solvation") in line:
                self.solvation=True
            elif ("equiltimeNVT") in line:
                self.equiltimeNVT= float(a)
            elif ("equilonly") in line:
                self.equilonly=True
            elif ("minonly") in line:
                self.minonly=True
            elif ("equilrestlambdascheme") in line:
                self.equilrestlambdascheme=commalist
        if self.ligandxyzfilename==None:
            self.ligandxyzfilename='None'
        if self.uncomplexedproteinpdbname==None:
            self.uncomplexedproteinpdbname='None'
        self.outputpath=os.path.join(os.getcwd(),'')
        self.logfh=open(self.outputpath+self.logname,'a+')
        self.SanitizeMMExecutables()
        if self.simpathlist!=None:
            self.binding=True
            tables.GrabSimDataFromPathList(self)
            plots.PlotEnergyData(self)
            sys.exit()

        if self.complexation==True and self.uncomplexedproteinpdbname!=None and self.complexedproteinpdbname!=None: 

            self.uncomplexedxyzname=self.uncomplexedproteinpdbname.replace('.pdb','.xyz')
            self.complexedxyzname=self.uncomplexedxyzname.replace('.xyz','_comp.xyz')
            self.receptorligandxyzfilename=self.complexedxyzname
            self.ReadReceptorCharge()
            pdbxyz.GenerateProteinTinkerXYZFile(self)
        elif self.complexation==True and self.uncomplexedproteinpdbname==None and self.complexedproteinpdbname==None and self.receptorligandxyzfilename==None: 
            raise ValueError('Missing either uncomplexed or complexed proteinpdbname')
        if self.ligandfilename!=None:
            self.ReadLigandCharge()

        
        if self.complexation==True:
            self.systemcharge=self.receptorcharge+self.ligandcharge
        elif self.solvation==True:
            self.systemcharge=self.ligandcharge
        if self.solvation==True:
            self.xyzfilename=self.ligandxyzfilename
        elif self.complexation==True:
            self.xyzfilename=self.receptorligandxyzfilename

        self.elementsymtotinktype={'K':box.GrabTypeNumber(self,'Potassium Ion K+'),'Cl':box.GrabTypeNumber(self,'Chloride Ion Cl-'),'Mg':box.GrabTypeNumber(self,'Magnesium Ion Mg+2'),'Li':box.GrabTypeNumber(self,'Lithium Ion Li+'),'Na':box.GrabTypeNumber(self,'Sodium Ion Na+'),'Rb':box.GrabTypeNumber(self,'Rubidium Ion Rb+'),'Cs':box.GrabTypeNumber(self,'Cesium Ion Cs+'),'Be':box.GrabTypeNumber(self,'Beryllium Ion Be+2'),'Ca':box.GrabTypeNumber(self,'Calcium Ion Ca+2'),'Zn':box.GrabTypeNumber(self,'Zinc Ion Zn+2'),'Mg+':box.GrabTypeNumber(self,'Magnesium Ion Mg+2')}
        self.waterOtypenum=349
        self.waterHtypenum=350 
        self.estatlambdascheme=self.extendelelambdaoneside+self.estatlambdascheme+self.extendelelambdazeroside
        self.vdwlambdascheme=self.extendvdwlambdaoneside+self.vdwlambdascheme+self.extendvdwlambdazeroside
        self.restlambdascheme=self.extendrestlambdaoneside+self.restlambdascheme+self.extendrestlambdazeroside
        self.equilstepsNVT=str(int((self.equiltimeNVT*1000000)/self.equiltimestep/len(self.equilibriatescheme))) # convert ns to fs divide by the length of temperature scheme
        self.equilstepsNPT=str(int((self.equiltimeNPT*1000000)/self.equiltimestep))
        self.proddynframenum=str(int(float(self.proddyntime)/(float(self.proddynwritefreq)*0.001)))
        self.proddynsteps=str(int((self.proddyntime*1000000)/self.proddyntimestep))
        self.equilframenumNPT=int((float(self.equiltimeNPT))/(float(self.equilwritefreq)*0.001))
        self.equilframenumNVT=int((float(self.equiltimeNVT))/(float(self.equilwritefreq)*0.001))

        self.equilframenum=self.equilframenumNPT+self.equilframenumNVT
        self.proddynframenum=int((self.proddyntime)/float(self.proddynwritefreq)*0.001)
         
        self.simfoldname='Sim'
        self.waterboxfilename='waterbox.xyz'
        if self.complexation==True:
            self.simfoldname='Comp'+self.simfoldname
            self.waterboxfilename='comp'+self.waterboxfilename
            self.WriteToLog('Complexation job')
        elif self.solvation==True:
            self.simfoldname='Solv'+self.simfoldname
            self.waterboxfilename='solv'+self.waterboxfilename
            self.WriteToLog('Solvation job')
        self.minwaterboxfilename=self.waterboxfilename.replace('.xyz','min.xyz')
        self.minwaterboxfilenamepymol=self.minwaterboxfilename.replace('.xyz','.pdb')
        self.equilwaterboxfilename=self.waterboxfilename.replace('.xyz','equil.xyz')
        self.proddynwaterboxfilename=self.waterboxfilename.replace('.xyz','proddyn.xyz')
        self.proddynwaterboxfilenamepymol=self.proddynwaterboxfilename.replace('.xyz','.pdb')
        self.equilarcwaterboxfilename=self.equilwaterboxfilename.replace('.xyz','.arc')
        self.proddynarcwaterboxfilename=self.proddynwaterboxfilename.replace('.xyz','.arc')
        self.proddynwaterboxkeyfilename=self.proddynwaterboxfilename.replace('.xyz','.key')
        self.equiljobsfilename='equiljobs.txt'
        self.proddynjobsfilename='proddynamicsjobs.txt'
        self.barjobsfilename='barjobs.txt'
        self.freeenergyjobsfilename='freeenergyjobs.txt'
        self.tightminjobsfilename='tightboxminjobs.txt'
        self.looseminjobsfilename='looseboxminjobs.txt'
        self.analyzejobsfilename='analyzejobs.txt'
        self.configkeyfilename=self.keyfilename.replace('.key','_config.key')
        self.lambdakeyfilename=self.keyfilename.replace('.key','_lambda.key')
        self.tightminoutput='tightmin.out'
        self.looseminoutput='loosemin.out'
        self.iontypetoionnumberphysio={} # key is ion type (string), value is number of computed ions 
        self.iontypetoionnumberneut={} # key is ion type (string), value is number of computed ions 
        self.equiloutputarray=[]
        for temp in self.equilibriatescheme:
            self.equiloutputarray.append(self.outputpath+self.simname+'_'+str(temp)+'_'+self.equilstepsNVT+'.out')
        self.nextfiletofinish=self.equiloutputarray[0]
        self.lambdafolderlist=[]
        self.proddynoutfilepath=[]
        mut=False
        if len(self.mutlambdascheme)==0:
            for i in range(len(self.vdwlambdascheme)):
                elelamb=self.estatlambdascheme[i]
                vdwlamb=self.vdwlambdascheme[i]
                if 'Comp' in self.simfoldname:
                    rest=self.restlambdascheme[i]
                    fold=self.simfoldname+"Ele%s_Vdw%s_Res%s"%(elelamb,vdwlamb,rest)
                else:
                    fold=self.simfoldname+"Ele%s_Vdw%s"%(elelamb,vdwlamb)
                self.lambdafolderlist.append(fold)
                outputfilepath=os.getcwd()+'/'+self.simfoldname+'/'+fold+'/'
                outputfilepath+=fold+'.out'
                self.proddynoutfilepath.append(outputfilepath)
        else:
            mut=True
            for i in range(len(self.mutlambdascheme)):
                mutlambda=self.mutlambdascheme[i]
                fold=self.simfoldname+"Mut%s"%(mutlambda)
                self.lambdafolderlist.append(fold)
                outputfilepath=os.getcwd()+'/'+self.simfoldname+'/'+fold+'/'
                outputfilepath+=fold+'.out'
                self.proddynoutfilepath.append(outputfilepath)
        if mut==True:
            mutate.SingleTopologyMutationProtocol(self)

        self.baroutputfilepath=[]
        self.barfilepath=[]
        self.thermooutputfilepath=[]
        self.secondarcpaths=[]
        self.firstarcpaths=[]
        for i in range(len(self.lambdafolderlist)-1): # stop before last one because using i+1 for grabbing next index
            firstfoldname=self.lambdafolderlist[i]
            secondfoldname=self.lambdafolderlist[i+1]
            firstarcpath=self.outputpath+self.simfoldname+r'/'+firstfoldname+'/'+self.proddynarcwaterboxfilename
            secondarcpath=self.outputpath+self.simfoldname+r'/'+secondfoldname+'/'+self.proddynarcwaterboxfilename
            self.secondarcpaths.append(secondarcpath)
            self.firstarcpaths.append(firstarcpath)
            baroutputfilepath=self.outputpath+self.simfoldname+r'/'+secondfoldname+'/'+firstfoldname+secondfoldname+'BAR1.out'
            thermooutputfilepath=baroutputfilepath.replace('BAR1.out','BAR2.out')
            barfilepath=secondarcpath.replace('.arc','.bar')
            self.barfilepath.append(barfilepath)
            self.baroutputfilepath.append(baroutputfilepath)
            self.thermooutputfilepath.append(thermooutputfilepath)


                
            self.tabledictkeysused=[]

        if not (self.which(self.dynamicommpath)) and self.externalapi==None:
            self.usegpu=False

        if self.usegpu==True:
            self.truedynamicpath=dynamicommpath
            self.truebarpath=barommpath
        else:
            self.truedynamicpath=dynamicpath
            self.truebarpath=barpath
 

        self.tabledict['Prod MD Ensemb']=self.proddynensem
        self.tabledict['Prod MD Time']=self.proddyntime
        self.tabledict['Prod MD Steps']=self.proddynsteps
        self.tabledict['Dynamic Writeout Frequency (ps)']=self.proddynwritefreq
        self.tabledict['Dynamic Time Step (fs)']=self.equiltimestep
        self.tabledict['Equil Time NPT']=self.equiltimeNPT
        self.tabledict['Equil Time NVT']=self.equiltimeNVT
        self.tabledict['Ligand Charge']=self.ligandcharge
        self.tabledict['Ligand Name']=self.ligandxyzfilename.replace('.xyz','')
        if self.complexation==True:
            self.tabledict['Receptor Charge']=self.receptorcharge
            self.tabledict['Receptor Name']=self.uncomplexedproteinpdbname.replace('.pdb','')


        tables.WriteTableUpdateToLog(self)

                     
    def TraceBack(self):
        tb = traceback.format_exc()
        self.WriteToLog(tb,prin=True)
 
    
    def WriteToLog(self,string,prin=False):
        now = time.strftime("%c",time.localtime())
        self.logfh.write(now+' '+ string + "\n")
        self.logfh.flush()
        os.fsync(self.logfh.fileno())
        if prin==True:
            print(now+' '+ string + "\n")

    def SanitizeMMExecutables(self):
        self.xyzeditpath=self.SanitizeMMExecutable(self.xyzeditpath)
        self.barpath=self.SanitizeMMExecutable(self.barpath)
        self.dynamicpath=self.SanitizeMMExecutable(self.dynamicpath)
        self.minimizepath=self.SanitizeMMExecutable(self.minimizepath)
        self.pdbxyzpath=self.SanitizeMMExecutable(self.pdbxyzpath)


    def SanitizeMMExecutable(self, executable):
        # Try to find Tinker executable with/without suffix
        if self.tinkerdir is None:
            self.tinkerdir = os.getenv("TINKERDIR", default="")
        exe = os.path.join(self.tinkerdir, executable)
        if self.which(exe) is None:
            exe = exe[:-2] if exe.endswith('.x') else exe + '.x'
            if self.which(exe) is None:
                print("ERROR: Cannot find Tinker {} executable".format(executable))
                sys.exit(2)
        return exe

    def which(self,program):
        def is_exe(fpath):
            try:
                 return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
            except:
                 return None
    
        fpath, fname = os.path.split(program)
        if fpath:
            if is_exe(program):
                return program
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                exe_file = os.path.join(path, program)
                if is_exe(exe_file):
                    return exe_file
    
        return None

    def ReadReceptorCharge(self):
        pdbmol=openbabel.OBMol()
        obConversion = openbabel.OBConversion()
        obConversion.SetInFormat('pdb')
        obConversion.ReadFile(pdbmol,self.uncomplexedproteinpdbname)
        chg=pdbmol.GetTotalCharge()
        self.receptorcharge=chg 

    def ReadLigandCharge(self):
        structfname=self.ligandfilename
        tmpconv = openbabel.OBConversion()
        inFormat = openbabel.OBConversion.FormatFromExt(structfname)
        tmpconv.SetInFormat(inFormat)
        tmpmol = openbabel.OBMol()
        tmpconv.ReadFile(tmpmol, structfname)
        chg=tmpmol.GetTotalCharge()
        self.ligandcharge=chg 

    def MakeTinkerXYZFileBabelReadable(self,tinkerxyzfilename):
        temptinkerxyzfilename=tinkerxyzfilename.replace('.xyz','_temp.xyz')
        temp=open(tinkerxyzfilename,'r')
        results=temp.readlines()
        temp.close()
        temp=open(temptinkerxyzfilename,'w')
        for lineidx in range(len(results)):
            line=results[lineidx]
            if lineidx==0:
                linesplit=line.split()
                newline=linesplit[0]+' '+'generic comment'+'\n'
            else:
                newline=line
            temp.write(newline)
        temp.close()
        os.remove(tinkerxyzfilename)
        os.rename(temptinkerxyzfilename,tinkerxyzfilename) 

    def PymolReadableFile(self,tinkerxyzfilename,outputname):
        tmpconv = openbabel.OBConversion()
        tmpconv.SetInFormat('txyz')
        mol=openbabel.OBMol()
        tmpconv.ReadFile(mol,tinkerxyzfilename)
        tmpconv.SetOutFormat('pdb')
        tmpconv.WriteFile(mol,outputname)

                
if __name__ == '__main__':
    annihilator=Annihilator()
    ann.main(annihilator)
    
