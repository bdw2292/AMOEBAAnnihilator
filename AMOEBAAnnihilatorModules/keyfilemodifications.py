import os
    
def AddKeyWord(annihilator,keypath,string):
    annihilator.WriteToLog('Adding key words to '+keypath+' '+string,prin=True)
    read=open(keypath,'r')
    results=read.readlines()
    read.close()
    tempkeyname=keypath.replace('.key','-t.key')
    temp=open(tempkeyname,'w')
    found=CheckIfStringAlreadyInKeyfile(annihilator,keypath,string)
    if found==False:
        temp.write(string)
        for line in results:
            temp.write(line)
    temp.close()
    os.remove(keypath)
    os.rename(tempkeyname,keypath)

def CheckIfStringAlreadyInKeyfile(annihilator,keypath,string):
    found=False
    temp=open(keypath,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if string in line:
            found=True
    return found
    
def RemoveKeyWord(annihilator,keypath,keystring):
    annihilator.WriteToLog('Removing key word from '+keypath+' '+keystring,prin=True)
    read=open(keypath,'r')
    results=read.readlines()
    read.close()
    tempname=keypath.replace('.key','-t.key')
    temp=open(tempname,'w')
    for line in results:
        if keystring not in line:
            temp.write(line)
    temp.close()
    os.remove(keypath)
    os.rename(tempname,keypath)  
    
def InsertKeyfileHeader(annihilator,keyfilename):
    annihilator.WriteToLog('Adding header to key file '+keyfilename,prin=True)
    string='parameters '+annihilator.prmfilepath+'\n'
    if not CheckIfStringAlreadyInKeyfile(annihilator,keyfilename,string):
        AddKeyWord(annihilator,keyfilename,string)
    string='archive'+'\n'
    if not CheckIfStringAlreadyInKeyfile(annihilator,keyfilename,string):
        AddKeyWord(annihilator,keyfilename,string)
    string='integrator '+annihilator.integrator+'\n'
    if not CheckIfStringAlreadyInKeyfile(annihilator,keyfilename,string):
       AddKeyWord(annihilator,keyfilename,string)
    string='thermostat '+annihilator.thermostat+'\n'
    if not CheckIfStringAlreadyInKeyfile(annihilator,keyfilename,string):
        AddKeyWord(annihilator,keyfilename,string)
    string='ewald'+'\n'
    if not CheckIfStringAlreadyInKeyfile(annihilator,keyfilename,string):
        AddKeyWord(annihilator,keyfilename,string)
    string='vdw-cutoff '+str(annihilator.vdwcutoff)+'\n'
    if not CheckIfStringAlreadyInKeyfile(annihilator,keyfilename,string):
        AddKeyWord(annihilator,keyfilename,string)
    string='ewald-cutoff '+str(annihilator.ewaldcutoff)+'\n'
    if not CheckIfStringAlreadyInKeyfile(annihilator,keyfilename,string):
        AddKeyWord(annihilator,keyfilename,string)
    string='polar-eps '+str(annihilator.polareps)+'\n'
    if not CheckIfStringAlreadyInKeyfile(annihilator,keyfilename,string):
        AddKeyWord(annihilator,keyfilename,string)
    string='polar-predict'+'\n'
    if not CheckIfStringAlreadyInKeyfile(annihilator,keyfilename,string):
        AddKeyWord(annihilator,keyfilename,string)
    string='barostat'+' '+annihilator.barostatmethod+'\n'
    if not CheckIfStringAlreadyInKeyfile(annihilator,keyfilename,string):
        AddKeyWord(annihilator,keyfilename,string)
    string='heavy-hydrogen'+'\n'
    if annihilator.equiltimestep>=3 and annihilator.proddyntimestep>=3:
       if not CheckIfStringAlreadyInKeyfile(annihilator,keyfilename,string):
           AddKeyWord(annihilator,keyfilename,string)





def AddMultipoleDefinitionsForIonIndices(annihilator,ionindexes,charge,keyfilepath):
    for ionindex in ionindexes:
        AddIonMultipoleDefinition(annihilator,ionindex,charge,keyfilepath)

def AddIonMultipoleDefinition(annihilator,ionindex,charge,keyfilepath):
    temp=open(keyfilepath,'a')
    chgline='multipole '+'-'+str(ionindex)+'   0   0                '+str(charge)+'\n'
    temp.write(chgline)
    temp.write('                                       0.00000   0.00000  -0.00000'+'\n')
    temp.write('                                       0.00000'+'\n')
    temp.write('                                       0.00000  -0.00000'+'\n')
    temp.write('                                      -0.00000   0.00000  -0.00000'+'\n')
    temp.close()


