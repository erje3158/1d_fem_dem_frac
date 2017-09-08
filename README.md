# 1D Uniaxial Finite Strain, FEM-DEM Heirarchical Multiscale Modeling Code with Particle Fracture Model
* MPI/OpenMP Formulation 
* DEM - ellip3D, 600 particle assembly from Colorado Mason sand
* Simulates Split Hopkinson Pressure Bar experiments
* Particle fracture model developed by Drs. Boning Zhang and Richard Regueiro

## Boundary Conditions:
---> +----+----+----+----+----|  
&nbsp;d(t) &emsp;&emsp;&emsp;&ensp; el &emsp; <--dof &emsp;&emsp;&ensp; 0  
1. Applied strain rate from Split Hopkinson Pressure Bar Experiments (Dr. Hongbing Lu, UT Dallas).
2. Other side fixed in space.

## Initial Conditions:
d(0) = v(0) = a(0) for all degrees of freedom

## Code Overview:
./bin - contains pre-compiled programs required for ellip3D  
./include - hooks for armadillo library required for FEM  
./inputs - input files for FEM and DEM  
./lib - aramdillo library files  
./outputs - simulation outputs populated from PBS submissions  
./source - source code

## Input Files:
Four input files are required to run the code all found in ./inputs:
* dem_input - input parameters for the ellip3D  
* fem_input - input parameters for the FEM code  
* input_particle_file - provides the initial assembly of DEM particles  
* input_boundary_file - provides the initial location of the DEM assembly boundarys

### _dem_input_: Template and Variable Definitions
$Particle Overlap  
Maximum Overlap - overlap over which the force-displacement relationship is allowed to increase during particle-to-particle contact

$DEM Constitutive  
Young's Modulus - elasticity parameter for the DEM particles  
Poisson's Ratio - elasticity parameter for the DEM particles

$Time Paramters  
dt - DEM timestep

$Particle Contact  
Damping - Interparticle damping ratio

### _fem_input_: Template and Variable Definitions
$FEM Constitutive  
lambda - elasticity paramter for the FEM (not used)  
mu - elasticity parameter for the FEM (not used)  
rho - density of the FEM (not used)  

$Gravity  
gravity - acceleration due to gravity  

$FEM Geometry  
diameter - Diameter of the specimen  
L/D - Ratio of the length to the diameter of the specimen  

$DEM Geometry  
length - length of the DEM assembly  
width - width of the DEM assembly  
height - height of the DEM assembly  

$FEM Constants  
numips - number of integration points per element  
nstress - number of saved stress and strain values  
nisv - number of internal state variables  
ndof - number of global degrees of freedom  
nel - number of elements  
neldof - number of element degrees of freedom  

$Time Parameters  
t0 - start time of the simulation  
dt - global timestep  
iterations - number of DEM timesteps per global timestep  
print - number of DEM snapshots per global timestep  
time total - end time of the simlation  

$Strain Rate  
rate - strain rate of the split Hopkinson pressure bar experiment  

$Mass Damping  
alpha - mass proportional damping value  

## Output Files:
Output files will be generated in ./outputs in the directory ./outputs/\<PBS Name\>/\<Job Number\>. \<PBS Name\> is the name of the simulation specified in the PBS script, and \<Job Number\> is the job number assigned by the PBS queue.

## Code Compilation:
To compile the code, "cd" into ./source, and run ./Make_Script.sh. "Make_Script.sh" takes the platform as an argument. Current platforms are "topaz" and "excalibur":  
./source/Make_Script.sh \<platform\>

## Code Submission:
To submit the code to the PBS queue system of the specified platform, run ./Submit_Script.sh. "Submit_Script.sh" takes two arguments: the platform and the name of the PBS script. Curent platforms are "topaz" and "excalibur", and the PBS scripts are created and stored in ./source/\<platform\>/. Each PBS script needs to specify the queue, number of nodes, walltime, job code, your email, and location of each input file.  
./Submit_Script.sh \<platform> \<PBS Script\>  

### _Example PBS Script_:
#!/bin/bash  
#PBS -N \<PBS Name\>  
#PBS -l select=\<Number of Nodes\>:ncpus=32:mpiprocs=1  
#PBS -l walltime=\<walltime\>  
#PBS -l place=scatter:excl  
#PBS -q \<queue\>  
#PBS -j oe  
#PBS -V  
#PBS -A \<job code\>  
#PBS -m be  
#PBS -M \<email\>  

#--- USER INPUT ---  
export PREFIX=$PBS_JOBNAME #PBS_JOBNAME is the name of the job that's been submitted  
export LOCAL_DIR=$PBS_O_WORKDIR #the directory where the script was run  
export TMPD=$WORKDIR #my personal work directory on excalibur - data here is temp (15 days after run done)  
JOBNUM=`echo ${PBS_JOBID} | cut -d "." -f 1` #create variable out of job number PBS assigns this job  

export NSLOTS=`wc -l $PBS_NODEFILE | cut -f1 -d" "`  
echo NSLOTS = $NSLOTS  

export OMP_NUM_THREADS=32  
echo OMP_NUM_THREADS = $OMP_NUM_THREADS  

#--- HARDCODED DIRECTORIES ---  

export EXE=${LOCAL_DIR}/source  
export BIN=${LOCAL_DIR}/bin  
export LIB=${LOCAL_DIR}/lib  
export INP=${LOCAL_DIR}/inputs  
export OUT=${LOCAL_DIR}/outputs  

#--- WORKING DIRECTORY ---  

export TMP_DIR=${TMPD}/${JOBNUM} #create directory to run the job in $WORKDIR/$JOBNUM  
mkdir -p ${TMP_DIR}  
mkdir -p ${OUT}/${PREFIX}/${JOBNUM}  
cp -r ${EXE}/hu_code ${TMP_DIR} #copies everything from the place this script is run into the work dir  
cd ${TMP_DIR}  
ln -s ${TMP_DIR} ${OUT}/${PREFIX}/${JOBNUM}/${JOBNUM} #create link to the work dir  
pwd  

#--- LD_LIBRARY_PATH ---  

export LD_LIBRARY_PATH="${LIB}:$LD_LIBRARY_PATH"  
echo $LD_LIBRARY_PATH  

#--- MACHINE SPECIFIC ---  
module swap PrgEnv-intel PrgEnv-gnu  
module list  

echo Simulation started at `date`  
aprun -n $NSLOTS ./hu_code ${INP}/input_boundary_file ${INP}/input_particle_file ${BIN}/qdelaunay . ${INP}/fem_input_2el ${INP}/dem_input  
echo Simulation finished at `date`  

