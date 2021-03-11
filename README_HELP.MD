## Help File


### Minimum Input Example Complexation
* Input inside AMOEBA.ini
* Make sure pdb files (complexed and uncomplexed) have no missing residues or atoms
* For complexed pdb file, make sure the ligand is appended at the end of the file.
```
uncomplexedproteinpdbname=uncomplexed.pdb
complexedproteinpdbname=anilinecomp.pdb
complexation
keyfilename=aniline.key
ligandxyzfilename=aniline.xyz
```

### Minimum Input Example Solvation
* Input inside AMOEBA.ini
* Make sure that extended conformation is used, to ensure correct box size to avoid interacting with self images.
```
solvation
keyfilename=aniline.key
ligandxyzfilename=aniline.xyz
```

### Usage
* Navigate to file with AMOEBA.ini
```
nohup python ~/AMOEBAAnnihhilator/AMOEBAAnnihilatorModules/amoebaannihilator.py &
 
```


### Binding Inputs
```
simpathlist                             List of paths to completed solvation and complexation simulation data. Useful for combining data in spreadsheets and plotting binding free energy estimates. If in input file, program will only combine data from input folders and then exit. Usage (simpathlist= path_to_comp1 path_to_solv1 , path_to_comp2 path_to_solv2 ...
                                        Default=None
```


### Complexation Inputs
```
uncomplexedproteinpdbname               Uncomplexed protein PDB file, needed for total charge and also total atom number.
                                        Default=None
complexedproteinpdbname                 Complexed protein PDB file, needed for extracting the ligand and converting to tinker XYZ file. 
                                        Default=None
receptorligandxyzfilename               This tinker XYZ file contains ligand and protein, automatically is generated when the above two inputs are given.
                                        Default=None 
```
### Box Setup Inputs
```
fixedboxsize                            Fixed box size for cubic PBC box. If not used, will be determined automatically. 
                                        Default=None
listofsaltcons                          Salt concentration for physiological systems. Units in mM
                                        Default=[KCl]=100                   
boxbufferlength                         Extra box length added in case protein unfolds a bit and becomes larger (dont want to interact with periodic images). Units in angstroms.
                                        Default=3
```
### Minimization Inputs
```
restrainpositionconstant                Restraint position force constant during minimization for freezing the ligand and protein. Units in kcal/mol/Ang^2.  
                                        Default=1 
loosemincriteria                        Convergence criteria for minimization, loose criteria. Units kcal/mol.
                                        Default=10
tightmincriteria                        Convergence criteria for minimization, tight criteria. Units kcal/mol.
                                        Default=1
```
### Restraint Inputs
```
restrainatomgroup1                      Group of atoms (protein or ligand) used for tinker group restraints. If not specified, will be automatically determine (for group restraints during equilbriation).
                                        Default=None 
restrainatomgroup2                      Group of atoms (protein or ligand) used for tinker group restraints. If not specified, will be automatically determine (for group restraints during equilbriation). 
                                        Default=None
flatbotrest                             Determines whether to use flatbottom potential with group restraints.
                                        Default=True 
distancerestraintconstant               Restraint constant for restrain-distance keyword during production dynamics. Units in kcal/mol/Ang^2. 
                                        Default=10.
anglerestraintconstant                  Restraint constant for restrain-angle and restrain-torsion during production dynamics. Units in kcal/mol/deg^2.
                                        Default=0.003046
```

### Equilibriation Inputs
```
equilrestrainsphereradius               Sphere radius for protein atoms to restrain position during equilibriation. Units angstroms.
                                        Default=2
equiltimeNVT                            Total time in NVT ensemble during equilbriation (for all temperatures combined). Units in ns.
                                        Default=5
equilrestlambdascheme                   Group restraint force constants used during equilibriation (NVT).
                                        Default=5,2,1,.1,0
equilibriatescheme                      Temperature used for each NVT step during equilibriation. Units in K. 
                                        Default=50,100,150,200,300 
equilwritefreq                          How often to write out new frame during equilbriation MD. Units in ps.
                                        Default=100                            
equiltimestep                           MD time step during equilbriation. Units in fs.
                                        Default=2
equiltimeNPT                            Total time in NPT ensemble for determining average box length. Units in ns.
                                        Default=2 
barostatmethod                          The barostat method used during NPT.
                                        Default=montecarlo
```
### Production Dynamics Inputs
```
proddyngrprests                         Determiens whether to use group restraints or not during production dynamics.If not, then use distance, angle and torsion restraints.
                                        Default=False
restlambdascheme                        Restraint lambda scheme during production dynamics.
                                        Default=1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0
vdwlambdascheme                         Van der Waals lambda scheme during production dynamics. 
                                        Default=0,.45,.52,.56,.58,.6,.62,.64,.67,.7,.75,.8,.85,.9,.95,1,1,1,1,1,1,1,1,1,1,1,1
estatlambdascheme                       Electrostatics lambda scheme during production dynamics.
                                        Default=0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1,1
proddynwritefreq                        How often to write out new frame during production dynamics MD. Units in ps.
                                        Default=2
proddyntime                             Total time to run production dynamics (in NVT). Units in ns.
                                        Default=5
proddyntimestep                         MD time step during production dynamics. Units in fs.
                                        Default=2
integrator                              MD integrator method.
                                        Default=RESPA
thermostat                              Thermostat method during NVT MD.
                                        Default=BUSSI
vdwcutoff                               Cutoff radius for van der Waals interactions. Units in angstroms.
                                        Default=12
ewaldcutoff                             Cutoff radius for electrostatic interactions. Units in angstroms.
                                        Default=7
polareps                                Grid spacing for electrostatic potential surface.
                                        Default=.001
```
### Generic Inputs 
```
complexation                            Specifies to do complexation simulations or not.
                                        Default=False
solvation                               Specifies to do solvation simulations or not.
                                        Default=False
ligandxyzfilename                       Input tinker XYZ filename for the ligand.
                                        Default=None 
keyfilename                             Input tinker key file for the ligand parameters.
                                        Default=None
ligandcharge                            Input charge for ligand. If ligandfilename is given, charge will be read from file instead.
                                        Default=None
ligandfilename                          Ligand structure file containing charge of ligand.
                                        Default=None
receptorcharge                          Input charge for protein. If uncomplexedproteinpdbname is given, receptorcharge is read from file.
                                        Default=None
```

### Misc. Inputs
```
externalapi                             External script path for submitting MD jobs. Used with submitjobs.py. 
                                        submitjobs.py needs to be edited for user specific subimssion script.
                                        Default=None

dontrestrainreceptorligand              Dont use any restraints between ligand and protein (for example if ligand                                            is covalently bonded).
                                        Default=False
```

### Mutation Inputs
```
bgnstatekey                             First state key file for mutation.
                                        Default=None

endstatekey                             End state key file for mutation.
                                        Default=None

bgnstatexyz                             First state xyz file for mutation.
                                        Default=None

endstatexyz                             End state xyz file for mutation.
                                        Default=None

mutlambdascheme                         Lambda scheme (from 0 to 1) of mutation of any parameters between begin and 
                                        end state keyfiles.
                                        Default=[]
```
