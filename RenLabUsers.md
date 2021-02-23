## Installation and Quick Start for Ren Lab Users

## Install Annihilator and Daemon
```
mkdir ~./AMOEBAAnnihilator~./AMOEBAAnnihilator
cd 
git clone https://github.com/bdw2292/AMOEBAAnnihilator.git
mkdir ~./RenLabDaemon
cd ~./RenLabDaemon
git clone https://github.com/bdw2292/Ren-Lab-Daemon.git
```

## Tinker
We do not have to install Tinker. Paths are included in the following bashrc files.
https://docs.google.com/spreadsheets/d/1EOlUwFpdNU2uBZ5XrHSvZnRYCUCOw5tCSkFisTm3big/edit#gid=1357198709


## Environment Bashrc Example For CPU Tinker
```
conda activate /home/bdw2292/miniconda3/envs/amoebamd
export PATH=/home/bdw2292/TinkerCPU/bin:$PATH
export PATH=/home/bdw2292/miniconda3/envs/p4env/bin/:$PATH
export PYTHONPATH=/home/bdw2292/:$PYTHONPATH
```

## Environment Bashrc Example For GPU Tinker
```
conda activate /home/bdw2292/miniconda3/envs/amoebamd
export CUDA_HOME=/usr/local/cuda-10.2/
export PATH=$PATH:$CUDA_HOME/bin/
export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$CUDA_HOME/lib:$LD_LIBRARY_PATH
export OPENMM_CUDA_COMPILER=$CUDA_HOME/bin/nvcc
export OPENMMHOME=/home/liuchw/OpenMM-Nov-2020/tinker/source/
export LD_LIBRARY_PATH=/home/liuchw/OpenMM-Nov-2020/bin/lib:$LD_LIBRARY_PATH
export OPENMM_PLUGIN_DIR=/home/liuchw/OpenMM-Nov-2020/bin/lib/plugins
export PATH=/home/liuchw/OpenMM-Nov-2020/tinker/source/:$PATH
alias dynamic_omm.x=$OPENMMHOME/dynamic_omm.x
alias bar_omm.x=$OPENMMHOME/bar_omm.x
alias analyze_omm.x=$OPENMMHOME/analyze_omm.x

```

To use Daemon, include a keyword like this in AMOEBA.ini
```
externalapi=/home/bdw2292/ExternalAPIRenLab/submit.py
```

* Make sure to follow readme instructions in Daemon before running Annihilator program, https://github.com/bdw2292/Ren-Lab-Daemon/blob/main/README_HELP.MD
* Daemon does not need to be called directly, the annihilator program will do this for you. However, make sure to comment out hostnames that are reserved by other people (or add CPUONLY/GPUONLY).

