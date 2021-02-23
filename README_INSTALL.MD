## Tinker CPU Installation
```
git clone git@github.com:TinkerTools/Tinker.git TinkerCPU
cd TinkerCPU
cd fftw
./configure --prefix=/path_to_TinkerCPU/fftw/ --enable-openmp --enable-threads
make -j 8
make install
cd ..
cp ./make/Makefile ./source
cd source
```
* Open Makefile and comment out default operating system (MAC, unless you have mac)
* Uncomment the default linux installation lines
* Keep the RENAME line commented out
```
make -j 8
make install
```

## Tinker GPU Installation

```
git clone git@github.com:TinkerTools/Tinker-OpenMM.git OpenMM
mkdir -p buildWithCuda10.0 bin
touch openmm_cuda_10.0.bashrc
```
* These are example paths for whatever cuda version you have
* Open openmm_cuda_10.0.bashrc

```
export CUDA_HOME=/usr/local/cuda-10.0/
export PATH=$PATH:$CUDA_HOME/bin/
export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$CUDA_HOME/lib:$LD_LIBRARY_PATH
export OPENMM_CUDA_COMPILER=$CUDA_HOME/bin/nvcc
```

```
source openmm_cuda_10.0.bashrc
cd ./buildWithCuda10.0
ccmake -i ../OpenMM
```
* Type c to configure
* Type t for advanced mode
* Set the following flags to their shown values (example path given)

```
CMAKE_INSTALL_PREFIX             /home/bdw2292/OpenMM-TinkerWithCuda10.0/bin/
CMAKE_CXX_FLAGS                  --std=c++0x
OPENMM_BUILD_PYTHON_WRAPPERS     OFF
```
* Type c to confugre
* Type g to generate and exit
```
make -j 8 && make install
```
## Compile Tinker With GPU
```
cd ..
cp -r Tinker TinkerWithCuda10.0
cd TinkerWithCuda10.0
cp openmm/* ./source
cd source
```
* Open Makefile and comment out default operating system (MAC, unless you have mac)
* Uncomment the default linux installation lines
* Keep the RENAME line commented out
* Focus on the following lines
```
TINKERDIR = $(HOME)/TinkerWithCuda10.0
FFTWDIR = /home/bdw2292/TinkerWithCuda10.0/fftw/
OPENMMDIR = /home/liuchw/OpenMM-TinkerWithCuda10.0/bin
CUDA_DIR = /usr/local/cuda/
CUDA_LIB = $(CUDA_DIR)/lib64
NVML_INCLUDE = /usr/include/nvidia/gdk
NVML_LIB = /usr/local/cuda/lib64/stubs
```
* When installing with older versions of cuda may need to delete -Ofast flag from OPTFLAGS
* Now we are ready to compile with gpu content

```
make -j 8 && make install
```
* Open openmm_cuda_10.0.bashrc
* Add the following lines

```
export LD_LIBRARY_PATH=/home/bdw2292/OpenMM-TinkerWithCuda10.0/bin/lib:$LD_LIBRARY_PATH
export OPENMM_PLUGIN_DIR=/home/bdw2292/OpenMM-TinkerWithCuda10.0/lib/plugins
export OPENMMHOME=/home/bdw2292/TinkerWithCuda10.0/source
alias dynamic_omm.x=$OPENMMHOME/dynamic_omm.x
alias bar_omm.x=$OPENMMHOME/bar_omm.x
alias analyze_omm.x=$OPENMMHOME/analyze_omm.x
```
* Close the openmm_cuda_10.0.bashrc
* We are ready to test if dynamic_omm.x is working (make sure there are no errors in calling executable)

```
source ~/openmm_cuda_10.0.bashrc
dynamic_omm.x 
cd ../../
rm -r buildWithCuda10.0
```


## CPU Tinker & Python Envioronment
```
conda create -n amoebamd python=3.6 --yes
conda activate amoebamd
conda install -c conda-forge openbabel=2.4.1 --yes
conda install -c conda-forge mdtraj --yes

```

* Environment Bashrc Example For CPU Tinker
```
conda activate amoebamd
export PATH=/home/bdw2292/TinkerCPU/bin:$PATH
export PATH=/home/bdw2292/miniconda3/envs/p4env/bin/:$PATH
export PYTHONPATH=/home/bdw2292/:$PYTHONPATH
```

* Environment Bashrc Example For GPU Tinker
```
conda activate amoebamd
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

* Example keyword for Ren-Lab-Daemon submission script (in AMOEBA.ini file) 
* Make sure to have understanding of how Daemon works before running Annihilator program, https://github.com/bdw2292/Ren-Lab-Daemon/blob/main/README_HELP.MD
* Daemon does not need to be called directly, the annihilator program will do this for you. However, make sure to comment out hostnames that are reserved by other people (or add CPUONLY/GPUONLY).
```
externalapi=/home/bdw2292/ExternalAPIRenLab/submit.py

```
