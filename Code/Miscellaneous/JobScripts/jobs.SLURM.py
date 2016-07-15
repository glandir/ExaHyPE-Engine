from subprocess import call

#processes = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 1680]
#pdegrees = [1, 2, 3, 4, 5, 6, 7, 8]
#hmaxs = ["0.4", "0.2", "0.04", "0.02", "0.005", "0.002", "0.0005"]
#compilers = ["GNU", "Intel"]

processes = [1, 2, 4, 8, 16, 32]
pdegrees = [1, 2, 3, 4, 5, 6, 7, 8]
hmaxs = ["0.4", "0.2", "0.04", "0.02"]
compilers = ["GNU"]


for process in processes:
  for pdegree in pdegrees:
    for hmax in hmaxs:
      for compiler in compilers:
        name = 'ExaHyPE_Euler__process_' + `process` + '__pdegree_' + `pdegree` + '__hmax_' + hmax + '__compiler_' + compiler
        filename = name + '.slurm'
        
        file = open(filename, 'w')
        
        file.write("#!/bin/bash"                                                                                                         + "\n")
        file.write("#SBATCH -o %j." + name + ".out "                                                                                     + "\n")
        file.write("#SBATCH -D /home/hpc/pr63so/gu89tik2/scratch/logs"                                                                   + "\n")
        file.write("#SBATCH -J " + name                                                                                                  + "\n")
        file.write("#SBATCH --get-user-env "                                                                                             + "\n")
        file.write("#SBATCH --clusters=mpp2"                                                                                             + "\n")
        file.write("# alternatively, use mpp1 "                                                                                          + "\n")
        file.write("#SBATCH --ntasks=" + `process`                                                                                       + "\n")
        file.write("# multiples of 28 for mpp2"                                                                                          + "\n")
        file.write("# multiples of 16 for mpp1"                                                                                          + "\n")
        # file.write("#SBATCH --cpus-per-task=28"                                                                                         + "\n")
        file.write("#SBATCH --cpus-per-task=1"                                                                                           + "\n")
        file.write("#SBATCH --mail-type=all "                                                                                            + "\n")
        file.write("#SBATCH --mail-user=varduhn@tum.de "                                                                                 + "\n")
        file.write("#SBATCH --export=NONE "                                                                                              + "\n")
        file.write("#SBATCH --time=01:00:00 "                                                                                            + "\n")
        file.write("source /etc/profile.d/modules.sh"                                                                                    + "\n")
        
        file.write(""                                                                                                                    + "\n")

        file.write("echo \"*1**********************\""                                                                                   + "\n")
        file.write("module load gcc"                                                                                                     + "\n")
        if compiler == "GNU":
          file.write("module unload mpi.intel"                                                                                           + "\n")
          file.write("module unload intel/15.0"                                                                                          + "\n")
          file.write("module unload intel"                                                                                               + "\n")
          file.write("module load mpi.ompi/1.10/gcc"                                                                                     + "\n")
          # file.write("module load mpi.intel/5.1_gcc"                                                                                    + "\n")
        else:
          file.write("module load intel"                                                                                                 + "\n")
          file.write("module load mpi.intel"                                                                                             + "\n")
        
        # file.write(""                                                                                                                    + "\n")
        # file.write("export"                                                                                                              + "\n")

        file.write(""                                                                                                                    + "\n")
        file.write("echo \"*2**********************\""                                                                                   + "\n")
        file.write("cd $SCRATCH/jobs"                                                                                                    + "\n")

        file.write(""                                                                                                                    + "\n")
        file.write("echo \"*3**********************\""                                                                                   + "\n")
        file.write("mkdir $SLURM_JOBID." + name                                                                                          + "\n")
        file.write("cd $SLURM_JOBID." + name                                                                                             + "\n")
        file.write("cp -r ~/storage/ExaHyPE ."                                                                                           + "\n")
        file.write("cd ExaHyPE/Code"                                                                                                     + "\n")

        file.write(""                                                                                                                    + "\n")
        file.write("echo \"*4**********************\""                                                                                   + "\n")
        file.write("cd Peano"                                                                                                            + "\n")
        file.write("tar xfvz peano.tar.gz"                                                                                               + "\n")
        file.write("git checkout .gitignore"                                                                                             + "\n")
        file.write("cd .."                                                                                                               + "\n")

        file.write(""                                                                                                                    + "\n")
        file.write("echo \"*5**********************\""                                                                                   + "\n")
        file.write("cd Toolkit"                                                                                                          + "\n")
        file.write("./build.sh"                                                                                                          + "\n")
        file.write("cd .."                                                                                                               + "\n")

        file.write(""                                                                                                                    + "\n")
        file.write("echo \"*6**********************\""                                                                                   + "\n")
        file.write("cd ApplicationExamples"                                                                                              + "\n")
        file.write("echo \"\" > myUserSpec.exahype"                                                                                      + "\n")
        file.write("echo \"/**                                                                 \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"                                                                    \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \" 2D Euler Flow                                                      \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"                                                                    \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \" A simple project                                                   \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"                                                                    \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \" */                                                                 \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"exahype-project  Euler                                              \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"                                                                    \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"  peano-path                 = ./Peano/peano                        \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"  tarch-path                 = ./Peano/tarch                        \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"  multiscalelinkedcell-path  = ./Peano/multiscalelinkedcell         \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"  sharedmemoryoracles-path   = ./Peano/sharedmemoryoracles          \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"  exahype-path               = ./ExaHyPE                            \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"  output-directory           = ./ApplicationExamples/EulerFlow      \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"  architecture               = noarch                               \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"                                                                    \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"  computational-domain                                              \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"    dimension                = 2                                    \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"    width                    = 1.0, 1.0                             \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"    offset                   = 0.0, 0.0                             \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"    end-time                 = 0.4                                  \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"  end computational-domain                                          \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"                                                                    \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"/*                                                                  \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"  shared-memory                                                     \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"    identifier               = autotuning                           \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"    cores                    = 2                                    \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"    properties-file          = sharedmemory.properties       		   \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"  end shared-memory                                                 \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"*/                                                                  \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"                                                                    \" >> myUserSpec.exahype"                 + "\n")
        if process > 1:
          file.write("echo \"  distributed-memory                                                \" >> myUserSpec.exahype"                 + "\n")
          file.write("echo \"    identifier               = static_load_balancing                \" >> myUserSpec.exahype"                 + "\n")
          file.write("echo \"    configure                = {greedy,FCFS}                        \" >> myUserSpec.exahype"                 + "\n")
          file.write("echo \"    buffer-size              = 64                                   \" >> myUserSpec.exahype"                 + "\n")
          file.write("echo \"    timeout                  = 120                                  \" >> myUserSpec.exahype"                 + "\n")
          file.write("echo \"  end distributed-memory                                            \" >> myUserSpec.exahype"                 + "\n")
          file.write("echo \"                                                                    \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"  optimisation                                                      \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"    fuse-algorithmic-steps        = off                             \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"    fuse-algorithmic-steps-factor = 0.99                            \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"  end optimisation                                                  \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"                                                                    \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"  solver ADER-DG MyEulerSolver                                      \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"    variables          = 5                                          \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"    parameters         = 0                                          \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"    order              = " + `pdegree` + "                          \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"    maximum-mesh-size  = " + hmax + "                               \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"    time-stepping      = global                                     \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"    kernel             = generic::fluxes::nonlinear                 \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"    language           = C                                          \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"                                                                    \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"    plot vtk::ascii                                                 \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"      time     = 0.0                                                \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"      repeat   = 0.4                                                \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"      output   = ./solution                                         \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"      select   = {all}                                              \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"    end plot                                                        \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"  end solver                                                        \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"                                                                    \" >> myUserSpec.exahype"                 + "\n")
        file.write("echo \"end exahype-project                                                 \" >> myUserSpec.exahype"                 + "\n")
        file.write("cd .."                                                                                                               + "\n")

        file.write(""                                                                                                                    + "\n")
        file.write("echo \"*7**********************\""                                                                                   + "\n")
        file.write("java -jar Toolkit/dist/ExaHyPE.jar --not-interactive ApplicationExamples/myUserSpec.exahype"                         + "\n")
        file.write(""                                                                                                                    + "\n")
        file.write("cd ApplicationExamples/EulerFlow"                                                                                    + "\n")

        file.write(""                                                                                                                    + "\n")
        file.write("echo \"*8**********************\""                                                                                   + "\n")
        file.write("export COMPILER=" + compiler                                                                                         + "\n")
        file.write("make clean"                                                                                                          + "\n")
        file.write("make -j56"                                                                                                           + "\n")
        if process == 1:
          file.write("./ExaHyPE-Euler ../myUserSpec.exahype > " + name + ".$SLURM_JOBID.out"                                             + "\n")
        else:
          file.write("mpiexec -np " + `process` + " ./ExaHyPE-Euler ../myUserSpec.exahype > " + name + ".$SLURM_JOBID.out"               + "\n")
          
        file.write("cat " + name + ".$SLURM_JOBID.out"                                                                                   + "\n")
        file.write("cd ../.."                                                                                                            + "\n")

        file.write(""                                                                                                                    + "\n")
        file.write("echo \"*9******\"$SLURM_JOBID\"****************\""                                                                   + "\n")

      
      
      file.close()
      
      call(["sbatch", filename])
      
      #print "sbatch " + filename
      
      call(["rm", filename])
      
