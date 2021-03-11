import numpy
import re
import time
import keyfilemodifications as keymods
import productiondynamics as prod
import terminate as term
import sys

def SingleTopologyMutationProtocol(annihilator):
   bgnstatexyzatominfo,bgnstateindextotypeindex,bgnstateatomnum=GrabXYZInfo(annihilator,annihilator.bgnstatexyz)
   endstatexyzatominfo,endstateindextotypeindex,endstateatomnum=GrabXYZInfo(annihilator,annihilator.endstatexyz)
   annihilator.bgnstatetypeindextoendstatetypeindex={}
   for bgnstateindex in bgnstateindextotypeindex.keys():
       bgnstatetypeindex=int(bgnstateindextotypeindex[bgnstateindex])
       endstatetypeindex=int(endstateindextotypeindex[bgnstateindex])
       annihilator.bgnstatetypeindextoendstatetypeindex[bgnstatetypeindex]=endstatetypeindex
   annihilator.bgnlinetoendline=GrabBgnToEndPrms(annihilator,annihilator.bgnstatetypeindextoendstatetypeindex)


def GrabXYZInfo(annihilator,xyzfile):
   temp=open(annihilator.outputpath+xyzfile,'r')
   xyzfileresults=temp.readlines()
   temp.close()
   xyzatominfo=[]
   indextotypeindex={}
   for lineidx in range(len(xyzfileresults)):
       line=xyzfileresults[lineidx]
       linesplit=line.split()
       if lineidx==0:
           xyzatomnum=int(linesplit[0])
       else:
           index=linesplit[0]
           type=linesplit[5]
           indextotypeindex[index]=type
           xyzatominfo.append(linesplit[2:])
   return xyzatominfo,indextotypeindex,xyzatomnum



def GenerateKeyFile(annihilator,arrayoflinearrays,keyfilename):
    temp=open(keyfilename,'w')
    for array in arrayoflinearrays:
        for line in array:
            temp.write(line)
        temp.write('\n')
    temp.close()

def RepresentsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

def GrabBgnToEndPrms(annihilator,bgnstatetypeindextoendstatetypeindex):
    bgnatomdefs,bgnbondprms,bgnangleprms,bgntorsionprms,bgnstrbndprms,bgnpitorsprms,bgnmpoleprms,bgnpolarizeprms,bgnopbendprms,bgnvdwprms=GrabParameters(annihilator,annihilator.bgnstatekey,bgnstatetypeindextoendstatetypeindex.keys(),bgnstatetypeindextoendstatetypeindex.keys())
    bgnlinetoendline={}
    for line in bgnatomdefs:
        if 'atom' in line:
            linesplit=re.split(r'(\s+)', line)
            bgnatomtype=int(linesplit[2])
            if bgnatomtype not in bgnstatetypeindextoendstatetypeindex.keys():
                continue
            endatomtype=bgnstatetypeindextoendstatetypeindex[bgnatomtype]
            linesplit[2]=str(endatomtype)
            linesplit[4]=str(endatomtype)
            newline=''.join(linesplit[:4+1])
            bgnlinetoendline[line]=newline


    for line in bgnbondprms:
        if 'bond' in line or 'pitors' in line:
            linesplit=re.split(r'(\s+)', line)
            bgnatomtype1=int(linesplit[2])
            bgnatomtype2=int(linesplit[4])
            if bgnatomtype1 not in bgnstatetypeindextoendstatetypeindex.keys():
                continue
            if bgnatomtype2 not in bgnstatetypeindextoendstatetypeindex.keys():
                continue

            endatomtype1=bgnstatetypeindextoendstatetypeindex[bgnatomtype1]
            endatomtype2=bgnstatetypeindextoendstatetypeindex[bgnatomtype2]
            linesplit[2]=str(endatomtype1)
            linesplit[4]=str(endatomtype2)
            newline=''.join(linesplit[:4+1])
            bgnlinetoendline[line]=newline

    for line in bgnangleprms:
        if 'angle' in line or 'strbnd' in line:
            linesplit=re.split(r'(\s+)', line)
            bgnatomtype1=int(linesplit[2])
            bgnatomtype2=int(linesplit[4])
            bgnatomtype3=int(linesplit[6])
            if bgnatomtype1 not in bgnstatetypeindextoendstatetypeindex.keys():
                continue
            if bgnatomtype2 not in bgnstatetypeindextoendstatetypeindex.keys():
                continue
            if bgnatomtype3 not in bgnstatetypeindextoendstatetypeindex.keys():
                continue

            endatomtype1=bgnstatetypeindextoendstatetypeindex[bgnatomtype1]
            endatomtype2=bgnstatetypeindextoendstatetypeindex[bgnatomtype2]
            endatomtype3=bgnstatetypeindextoendstatetypeindex[bgnatomtype3]
            linesplit[2]=str(endatomtype1)
            linesplit[4]=str(endatomtype2)
            linesplit[6]=str(endatomtype3)
            newline=''.join(linesplit[:6+1])
            bgnlinetoendline[line]=newline
               
    for line in bgntorsionprms:
        if 'torsion' in line:
            line=line.lstrip()
            linesplit=re.split(r'(\s+)', line)
            bgnatomtype1=int(linesplit[2])
            bgnatomtype2=int(linesplit[4])
            bgnatomtype3=int(linesplit[6])
            bgnatomtype4=int(linesplit[8])
            if bgnatomtype1 not in bgnstatetypeindextoendstatetypeindex.keys():
                continue
            if bgnatomtype2 not in bgnstatetypeindextoendstatetypeindex.keys():
                continue
            if bgnatomtype3 not in bgnstatetypeindextoendstatetypeindex.keys():
                continue
            if bgnatomtype4 not in bgnstatetypeindextoendstatetypeindex.keys():
                continue

            endatomtype1=bgnstatetypeindextoendstatetypeindex[bgnatomtype1]
            endatomtype2=bgnstatetypeindextoendstatetypeindex[bgnatomtype2]
            endatomtype3=bgnstatetypeindextoendstatetypeindex[bgnatomtype3]
            endatomtype4=bgnstatetypeindextoendstatetypeindex[bgnatomtype4]
            linesplit[2]=str(endatomtype1)
            linesplit[4]=str(endatomtype2)
            linesplit[6]=str(endatomtype3)
            linesplit[8]=str(endatomtype3)
            newline=''.join(linesplit[:8+1])
            bgnlinetoendline[line]=newline


    for lineidx in range(len(bgnmpoleprms)):
        line=bgnmpoleprms[lineidx]
        if 'multipole' in line:
            linesplit=re.split(r'(\s+)', line)
            framedefsplit=line.split()[1:-1]
            sametypes=True
            newframedefssplit=[]
            continbool=False
            for val in framedefsplit:
                if val[0]=='-':
                    typenum=int(val[1:])
                else:
                    typenum=int(val)
                if typenum not in bgnstatetypeindextoendstatetypeindex.keys():
                    continbool=True
                    continue
                newtype=bgnstatetypeindextoendstatetypeindex[typenum]
                if newtype!=typenum:
                    sametypes=False
                if val[0]=='-':
                    newframedefssplit.append('-'+str(newtype))
                else:
                    newframedefssplit.append(str(newtype))
            if continbool==True:
                continue
            count=0
            for j in range(len(linesplit)):
                string=linesplit[j]
                if RepresentsInt(string)==True:
                    newtype=newframedefssplit[count]
                    linesplit[j]=newtype
                    count+=1
                    continue
            newline=''.join(linesplit[:-1])
            bgnlinetoendline[line]=newline


            
    for line in bgnpolarizeprms:
        if 'polarize' in line:
            linesplit=re.split(r'(\s+)', line)
            bgnatomtype1=int(linesplit[2])
            if bgnatomtype1 not in bgnstatetypeindextoendstatetypeindex.keys():
                continue

            endatomtype1=bgnstatetypeindextoendstatetypeindex[bgnatomtype1]
            linesplit[2]=str(endatomtype1)
            
            newline=''.join(linesplit[:3])
            bgnlinetoendline[line]=newline


    for line in bgnopbendprms:
        if 'opbend' in line:
            linesplit=re.split(r'(\s+)', line)
            bgnatomtype1=int(linesplit[2])
            bgnatomtype2=int(linesplit[4])
            if bgnatomtype1 not in bgnstatetypeindextoendstatetypeindex.keys():
                continue
            if bgnatomtype2 not in bgnstatetypeindextoendstatetypeindex.keys():
                continue

            endatomtype1=bgnstatetypeindextoendstatetypeindex[bgnatomtype1]
            endatomtype2=bgnstatetypeindextoendstatetypeindex[bgnatomtype2]
            linesplit[2]=str(endatomtype1)
            linesplit[4]=str(endatomtype2)
            newline=''.join(linesplit[:4+1])
            bgnlinetoendline[line]=newline


    for line in bgnvdwprms:
        if 'vdw' in line:
            linesplit=re.split(r'(\s+)', line)
            bgnatomtype1=int(linesplit[2])
            if bgnatomtype1 not in bgnstatetypeindextoendstatetypeindex.keys():
                continue

            endatomtype1=bgnstatetypeindextoendstatetypeindex[bgnatomtype1]
            linesplit[2]=str(endatomtype1)
            newline=''.join(linesplit[:2+1])
            bgnlinetoendline[line]=newline

    finalbgnlinetoendline={}
    prmarrays=GrabParameters(annihilator,annihilator.outputpath+annihilator.endstatekey,bgnstatetypeindextoendstatetypeindex.values(),bgnstatetypeindextoendstatetypeindex.values())
    for array in prmarrays:
        for lineidx in range(len(array)):
            line=array[lineidx]
            for key in bgnlinetoendline.keys():
                delim=bgnlinetoendline[key]
                if delim in line:
                    finalbgnlinetoendline[key]=line
                    if 'multipole' in line:
                        mpolearrayidx=bgnmpoleprms.index(key)
                        finalbgnlinetoendline[bgnmpoleprms[mpolearrayidx+1]]=array[lineidx+1]
                        finalbgnlinetoendline[bgnmpoleprms[mpolearrayidx+2]]=array[lineidx+2]
                        finalbgnlinetoendline[bgnmpoleprms[mpolearrayidx+3]]=array[lineidx+3]
                        finalbgnlinetoendline[bgnmpoleprms[mpolearrayidx+4]]=array[lineidx+4]

    return finalbgnlinetoendline

def AppendParametersIfRightTypeOrClass(annihilator,listtocheck,refindexes,listtoappend,line):
    for idx in listtocheck:
        if idx in refindexes and line not in listtoappend:
            listtoappend.append(line)
    return listtoappend

def MutateParameter(annihilator,bgnprm,endprm,mutlambda):
    return float(mutlambda)*float(bgnprm)+(1-float(mutlambda))*float(endprm)


def ModifyParameterLine(linesplit):
    line=' '.join(linesplit)
    newlinesplit=line.split()
    if 'torsion' in line:
        prms=newlinesplit[5:]
        newprms=prms[0::3]
        folds=prms[2::3]
        folds=[int(i) for i in folds]
        foldtoprms=dict(zip(folds,newprms))
        for i in range(1,4):
            if i not in foldtoprms.keys():
                foldtoprms[i]=str(0)
        foldtophase={}
        for fold,prm in foldtoprms.items():
            if fold==2:
                phase=str(180)
            else:
                phase=str(0)
            foldtophase[fold]=phase
        
        torline = ' torsion %7s %4s %4s %4s   ' % (newlinesplit[1],newlinesplit[2],newlinesplit[3],newlinesplit[4])
        for nfold in sorted(foldtoprms.keys()):
            prm=foldtoprms[nfold]
            phase=foldtophase[nfold]
            torline += ' %7.3s %.1s %d' % (prm,phase, nfold)
        torline += '\n'
        linesplit=re.split(r'(\s+)', torline)
    return linesplit



def MutateParameterInLine(annihilator,bgnlinesplit,endlinesplit,lineidxs,mutlambda):
    bgnlastindex=len(bgnlinesplit)-1
    endlastindex=len(endlinesplit)-1
    if bgnlastindex<lineidxs[-1]:
        bgnlinesplit=ModifyParameterLine(bgnlinesplit)

    elif endlastindex<lineidxs[-1]:
        endlinesplit=ModifyParameterLine(endlinesplit)

    for lineidx in lineidxs:
       bgnprm=float(bgnlinesplit[lineidx])
       endprm=float(endlinesplit[lineidx])
       mutprm=MutateParameter(annihilator,bgnprm,endprm,mutlambda)
       bgnlinesplit[lineidx]=str(mutprm)
    mutline=' '.join(bgnlinesplit)
    return mutline

def MutateAllParameters(annihilator,bgnlinetoendline,mutlambda):
    arrayofprmarrays=GrabParameters(annihilator,annihilator.outputpath+annihilator.bgnstatekey,annihilator.bgnstatetypeindextoendstatetypeindex.keys(),annihilator.bgnstatetypeindextoendstatetypeindex.keys())
    mutarrayofprmarrays=[]
    for array in arrayofprmarrays:
        mutarray=MutateParameters(annihilator,array,bgnlinetoendline,mutlambda)
        mutarrayofprmarrays.append(array)
    sys.exit()
    return mutarrayofprmarrays

def MutateParameters(annihilator,bgnstateprmlines,bgnlinetoendline,mutlambda):
    mutprmlines=[]
    for i in range(len(bgnstateprmlines)):
        bgnline=bgnstateprmlines[i]
        bgnlinesplit=re.split(r'(\s+)', bgnline)
        if bgnline not in bgnlinetoendline.keys():
            continue
        endline=bgnlinetoendline[bgnline]
        endlinesplit=re.split(r'(\s+)', endline)
        if 'bond' in bgnline:
            mutline=MutateParameterInLine(annihilator,bgnlinesplit,endlinesplit,[6],mutlambda)      
            mutprmlines.append(mutline)
        elif 'angle' in bgnline:
            mutline=MutateParameterInLine(annihilator,bgnlinesplit,endlinesplit,[8],mutlambda)
            mutprmlines.append(mutline)
        elif 'strbnd' in bgnline:
            mutline=MutateParameterInLine(annihilator,bgnlinesplit,endlinesplit,[8,10],mutlambda)
            mutprmlines.append(mutline)
        elif 'opbend' in bgnline:
            mutline=MutateParameterInLine(annihilator,bgnlinesplit,endlinesplit,[10],mutlambda)
            mutprmlines.append(mutline)
        elif 'torsion' in bgnline:
            mutline=MutateParameterInLine(annihilator,bgnlinesplit,endlinesplit,[10,16,22],mutlambda)
            mutprmlines.append(mutline)
        elif 'vdw' in bgnline:
            mutline=MutateParameterInLine(annihilator,bgnlinesplit,endlinesplit,[6],mutlambda)
            mutprmlines.append(mutline)
        elif 'polarize' in bgnline:
            mutline=MutateParameterInLine(annihilator,bgnlinesplit,endlinesplit,[4],mutlambda)
            mutprmlines.append(mutline)
        elif 'multipole' in bgnline:
            mutline=MutateParameterInLine(annihilator,bgnlinesplit,endlinesplit,[8],mutlambda)
            mutprmlines.append(mutline)
            bgnlinesplit=re.split(r'(\s+)',bgnstateprmlines[i+1])
            endlinesplit=re.split(r'(\s+)',endstateprmlines[i+1])
            mutline=MutateParameterInLine(annihilator,bgnlinesplit,endlinesplit,[1,3,5],mutlambda)
            mutprmlines.append(mutline)
            bgnlinesplit=re.split(r'(\s+)',bgnstateprmlines[i+2])
            endlinesplit=re.split(r'(\s+)',endstateprmlines[i+2])
            mutline=MutateParameterInLine(annihilator,bgnlinesplit,endlinesplit,[1],mutlambda)
            mutprmlines.append(mutline)
            bgnlinesplit=re.split(r'(\s+)',bgnstateprmlines[i+3])
            endlinesplit=re.split(r'(\s+)',endstateprmlines[i+3])
            mutline=MutateParameterInLine(annihilator,bgnlinesplit,endlinesplit,[1,3],mutlambda)
            mutprmlines.append(mutline)
            bgnlinesplit=re.split(r'(\s+)',bgnstateprmlines[i+4])
            endlinesplit=re.split(r'(\s+)',endstateprmlines[i+4])
            mutline=MutateParameterInLine(annihilator,bgnlinesplit,endlinesplit,[1,3,5],mutlambda)
            mutprmlines.append(mutline)
    return mutprmlines

def GrabParameters(annihilator,filename,types,classes):
    atomdefs=[]
    bondprms=[]
    angleprms=[]
    torsionprms=[]
    strbndprms=[]
    pitorsprms=[]
    mpoleprms=[]
    polarizeprms=[]
    opbendprms=[]
    vdwprms=[]
    temp=open(filename,'r')
    results=temp.readlines()
    temp.close()
    for lineidx in range(len(results)):
        line=results[lineidx]
        linesplit=line.split()
        linesplitall=re.split(r'(\s+)', line)
        if 'cubic' in line or 'quartic' in line or 'sextic' in line or 'pentic' in line or 'type' in line or 'unit' in line or '#' in line:
            continue
        if 'atom' in line:
            atomclasslist=[int(linesplit[1]),int(linesplit[2])]
            atomdefs=AppendParametersIfRightTypeOrClass(annihilator,atomclasslist,types,atomdefs,line)     
        if 'bond' in line:
            bondclasslist=[int(linesplit[1]),int(linesplit[2])]
            bondprms=AppendParametersIfRightTypeOrClass(annihilator,bondclasslist,classes,bondprms,line)      
        elif ('angle' in line or 'anglep' in line) :
            angleclasslist=[int(linesplit[1]),int(linesplit[2]),int(linesplit[3])]
            angleprms=AppendParametersIfRightTypeOrClass(annihilator,angleclasslist,classes,angleprms,line)
        elif 'torsion' in line:
            torsionclasslist=[int(linesplit[1]),int(linesplit[2]),int(linesplit[3]),int(linesplit[4])]
            torsionprms=AppendParametersIfRightTypeOrClass(annihilator,torsionclasslist,classes,torsionprms,line)
        elif 'strbnd' in line:
            angleclasslist=[int(linesplit[1]),int(linesplit[2]),int(linesplit[3])]
            strbndprms=AppendParametersIfRightTypeOrClass(annihilator,angleclasslist,classes,strbndprms,line)
        elif 'pitors' in line:
            bondclasslist=[int(linesplit[1]),int(linesplit[2])]
            pitorsprms=AppendParametersIfRightTypeOrClass(annihilator,bondclasslist,classes,pitorsprms,line)
        elif 'opbend' in line:
            bondclasslist=[int(linesplit[1]),int(linesplit[2])]
            opbendprms=AppendParametersIfRightTypeOrClass(annihilator,bondclasslist,classes,opbendprms,line)
        elif 'polarize' in line: 
            polarizelist=[int(linesplit[1])]
            polarizeprms=AppendParametersIfRightTypeOrClass(annihilator,polarizelist,types,polarizeprms,line)
        elif 'vdw' in line:
            vdwlist=[int(linesplit[1])]
            vdwprms=AppendParametersIfRightTypeOrClass(annihilator,vdwlist,classes,vdwprms,line)
        elif 'multipole' in line:
            newlinesplit=linesplit[1:-1]
            frames=[numpy.abs(int(i)) for i in newlinesplit]
            firstidx=frames[0]
            if firstidx in types:
                mpolelist=[line,results[lineidx+1],results[lineidx+2],results[lineidx+3],results[lineidx+4]]
                for mpoleline in mpolelist:
                    mpoleprms.append(mpoleline)

    return atomdefs,bondprms,angleprms,torsionprms,strbndprms,pitorsprms,mpoleprms,polarizeprms,opbendprms,vdwprms


