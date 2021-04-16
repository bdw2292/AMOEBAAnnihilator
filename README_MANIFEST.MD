## Manifest of Files and Modules


### Modules

* terminate.py                          Termination checks for output of Tinker programs.
* plots.py                              Plot generation from free energy data.
* submitjobs.py                         Job submisison to external program provided by the user.
* equilbriation.py                      Handles equilbriating PBC box system to room temperature. 
* pdbxyz.py                             Handles generating Tinker XYZ file from input pdb structure. 
* bar.py                                Handles free energy estimation using BAR method. 
* annihilation.py                       Handles calling other modules for free energy estimations.  
* tables.py                             Handles CSV table generation of simulation input and output parameters (UTF-8) format. 
* keyfilemodifications.py               Handles modifications made to Tinker key files. 
* productiondynamics.py                 Handles molecular dynamic setup and simulations for alchemical pertubation. 
* restraints.py                         Handles restraints between complexed ligand and protein.  
* boxsetup.py                           Handles generation of PBC box system.
* minimization.py                       Handles minimization of box sytem.  
* amoebaannihilator.py                  Handles initial simulation parameter inputs. 

### Files
* if _pymol.xyz file in name, then it is pymol readable

## PBC Box Setup
* *waterbox.xyz                         Final PBC box prior to minimization. 
* xyzedit_*.in                          Input file for Tinker xyzedit program. 
* *.key                                 Input Tinker key file (for ligand) for using Tinker programs.
* *.xyz                                 Input Tinker xyz file (for ligand) for using Tinker programs.  
## Equilibriation
* *_Temperature_StepNumber.out          Output file from MD simulation during equilbriation.  
* equiljobs.txt                         Input file for input into calling Tinker dynamics (can be used by external cluster subimssion script)
* *waterboxequil.dyn                    File for velocity of each atom at each output simulation frame.  
* *waterboxequil.arc                    File for position of each atom at each output simulation frame.
* *waterboxproddyn.xyz                  Tinker XYZ File for initial production dynamics (last frame from waterboxequil.arc ). 
* *waterboxequil.xyz                    Last frame from waterboxequil.arc.
* *_config.key                          Keyfile used during minimization and equilibriation.

## Production Dynamics
* *_lambda.key                          Keyfile with ele-lambda, vdw-lambda keywords added (added from *_config.key).
* proddynamicsjobs.txt                  Input file for calling Tinker dynamics. 
* *Sim                                  Folder for production dynamics simulations. 
* *SimEle_elelambda_Vdw_vdwlambda       Folder for production dynamics simulations for set of lambda.  
* *waterboxproddyn.xyz                  Starting frame for production dynamic simulations.  
* *waterboxproddyn.key                  Key file for production dynamic simulations (from _lambda.key). 
* *waterboxproddyn.dyn                  Velocity's of every atom for production dynamic simulations. 
* *waterboxproddyn.arc                  Positions of every aotm for production dynamic simulations.  
* *SimEle_elelambda_Vdw_vdw-lambda.out  Output file from production dynamic simulations.


## Free Energy Estimation
* barjobs.txt                                               Input file for calling Tinker BAR (compute .bar files).    
* freeenergyjobs.txt                                        Input file for calling Tinker BAR (compute free energy).
* *waterboxproddyn.bar                                      Output file from Tinker bar (forward and reverse energies).  
* *SimEle0_Vdw1SolvSimEle_elelambda_Vdw_vdwlambda_BAR1.out  Output file from Tinker BAR (computing .bar file).
* *SimEle0_Vdw1SolvSimEle_elelambda_Vdw_vdwlambda_BAR2.out  Output file from Tinker BAR (computing free energy).
* BARResults.csv                                            CSV file of free energies for each ele,vdw lambda combination (UTF-8) format. 
* *AbsBarPlot.png                                           Free energy plot for either solvation or complexation.
* BindBarPlot.png                                           Binding free energy plot (solvation and complexation).

## Program Input, Output
* AMOEBA.ini                            Input keyword file for program.                                            
* TINKER.log                            Output logfile for program. 
* SimData.csv                           Output simulation parameter and output free energy data for program (UTF-8) format. 


## Examples

