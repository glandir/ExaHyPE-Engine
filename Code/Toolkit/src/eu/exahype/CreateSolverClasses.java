package eu.exahype;

import eu.exahype.analysis.DepthFirstAdapter;
import eu.exahype.node.AAderdgSolver;
import eu.exahype.node.AComputationalDomain;
import eu.exahype.node.AFiniteVolumesSolver;
import eu.exahype.node.AProfiling;
import eu.exahype.node.AProject;

public class CreateSolverClasses extends DepthFirstAdapter {
  public Boolean valid = true;

  private DirectoryAndPathChecker _directoryAndPathChecker;

  private String _projectName;

  private String _microarchitecture;

  private java.util.List<String> _supportedMicroarchitectures;

  private java.util.Set<String>  _definedSolvers;

  private String _pathToLibxsmm;

  private int _dimensions;
  
  private boolean _enableProfiler;

  public CreateSolverClasses(DirectoryAndPathChecker directoryAndPathChecker) {
    _directoryAndPathChecker = directoryAndPathChecker;
    _supportedMicroarchitectures =
        java.util.Arrays.asList("wsm", "snb", "hsw", "knc", "knl", "noarch");
    _enableProfiler = false;
  }

  @Override
  public void inAProject(AProject node) {
    _projectName     = node.getName().toString().trim();
    _definedSolvers  = new java.util.HashSet<String>();

    if (node.getSolver().size() == 0) {
      System.out.println("there are no solvers in the specification file ... nothing to be done");
    }

    _microarchitecture = node.getArchitecture().toString().trim().toLowerCase();
    if (!_supportedMicroarchitectures.contains(_microarchitecture)) {
      System.out.println("Unknown architecture specified ... fallback solution \"noarch\" taken");
      _microarchitecture = "noarch";
    }
  }

  @Override
  public void inAPaths(eu.exahype.node.APaths node) {
    if (node.getLibxsmmPath() == null) {
      // attribute 'libxsmm-path' did not occur in spec file
      _pathToLibxsmm = "";
    } else {
      _pathToLibxsmm = node.getLibxsmmPath().toString().trim();
    }
  };

  @Override
  public void inAComputationalDomain(AComputationalDomain node) {
    _dimensions = Integer.parseInt( node.getDimension().toString().trim() );
    if (_dimensions!=2 && _dimensions!=3) {
      System.err.println( "ERROR: dimension has to be either 2 or 3.");
    }
  }

  
  @Override
  public void inAProfiling(AProfiling node) {
    _enableProfiler = !node.getProfiler().toString().trim().equals("NoOpProfiler");
  };

  @Override
  public void inAAderdgSolver(AAderdgSolver node) {
    String solverName = node.getName().toString().trim();
    
    if (_definedSolvers.contains(solverName)) {
      System.err.println( "ERROR: Solver " + solverName + " multiply defined" );
      valid = false;
    }
    else {
      _definedSolvers.add(solverName);
    }

    java.io.File headerFile = new java.io.File(
        _directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/" + solverName + ".h");
    java.io.File userImplementationFile = new java.io.File(
        _directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/" + solverName + ".cpp");
    java.io.File userPDEFile = null;
    java.io.File userTypesDefFile = null;
    java.io.File generatedImplementationFile =
        new java.io.File(_directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/"
            + solverName + "_generated.cpp");


    boolean isFortran = false;
    if (node.getLanguage().getText().trim().equals("C")) {
      isFortran = false;
    } else if (node.getLanguage().getText().trim().equals("Fortran")) {
      isFortran = true;
      userPDEFile =
          new java.io.File(_directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/PDE.f90");
      userTypesDefFile = new java.io.File(
          _directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/typesDef.f90");
    } else {
      System.err.println("ERROR: unknown language for solver " + node.getName().getText()
          + ". Supported languages are C and Fortran");
      valid = false;
      return;
    }

    String  kernel             = node.getKernel().toString().trim();
    boolean isLinear           = kernel.substring(kernel.lastIndexOf("::")).equalsIgnoreCase("::linear");
    String generalKernel       = kernel.substring(0, kernel.lastIndexOf("::"));
    int     numberOfVariables  = Integer.parseInt(node.getVariables().toString().trim());
    int     numberOfParameters = Integer.parseInt(node.getParameters().toString().trim());
    int     order              = Integer.parseInt(node.getOrder().toString().trim());
    boolean hasConstants       = node.getConstants()!=null;

    if (numberOfParameters != 0) {
      System.err.println("ERROR: At the moment, parameters are not yet supported. " + 
          " Please add the parameters as additional quantities to your PDE formulation.");
      valid = false;
      return;
    }
    
    if (order < 1 || order > 9) {
      System.err.println("ERROR: Only polynomial degrees of 1..9 are supported.");
      valid = false;
      return;
    }
    
    eu.exahype.solvers.Solver solver = null;

    if (isFortran && kernel.equals( eu.exahype.solvers.UserDefinedADER_DGinFortran.Identifier )) {
      solver = new eu.exahype.solvers.UserDefinedADER_DGinFortran();
    }
    else if (!isFortran && kernel.equals( eu.exahype.solvers.UserDefinedADER_DGinC.Identifier )) {
      solver = new eu.exahype.solvers.UserDefinedADER_DGinC(numberOfVariables,
        numberOfParameters, order, hasConstants, _enableProfiler);
    }
    // TODO(Dominic): I replaced this.
//    else if (generalKernel.equals( eu.exahype.solvers.GenericFluxesADER_DG.Identifier )) {
//      solver = new eu.exahype.solvers.GenericFluxesADER_DG(_dimensions,
//        numberOfVariables, numberOfParameters, order, _enableProfiler, hasConstants, isLinear, isFortran );
//    }
    else if (generalKernel.equals( eu.exahype.solvers.GenericADERDG.Identifier )) {
      solver = new eu.exahype.solvers.GenericADERDG(_dimensions,
        numberOfVariables, numberOfParameters, order, _enableProfiler, hasConstants, isLinear, isFortran );
    }
    else if (!isFortran && kernel.equals( eu.exahype.solvers.OptimisedFluxesLinearADER_DGinC.Identifier )) {
      solver = new eu.exahype.solvers.OptimisedFluxesLinearADER_DGinC(_dimensions,
        numberOfVariables, numberOfParameters, order, _microarchitecture, _pathToLibxsmm,
        _enableProfiler, hasConstants);
    }
    else if (!isFortran && kernel.equals( eu.exahype.solvers.OptimisedFluxesNonlinearADER_DGinC.Identifier )) {
      solver = new eu.exahype.solvers.OptimisedFluxesNonlinearADER_DGinC(_dimensions,
        numberOfVariables, numberOfParameters, order, _microarchitecture, _pathToLibxsmm,
        _enableProfiler, hasConstants);
    }
    else if (!isFortran && kernel.equals( eu.exahype.solvers.KernelEuler2d.Identifier )) {
       solver = new eu.exahype.solvers.KernelEuler2d();
    }

    if (solver == null) {
      System.err.println("creation solver " + solverName + " ... failed as kernel " + kernel
          + " for language " + node.getLanguage().getText().trim() + " is not supported");
      valid = false;
      return;
    }

    try {
      // =====================
      // Write all the headers
      // =====================
      if (headerFile.exists()) {
        System.out.println("create header of solver " + solverName + " ... header "
            + headerFile.getAbsoluteFile()
            + " does exist already. Remove to allow toolkit to regenerate it (changes will be lost)");
      } else {
        java.io.BufferedWriter headerWriter =
            new java.io.BufferedWriter(new java.io.FileWriter(headerFile));
        solver.writeHeader(headerWriter, solverName, _projectName);
        System.out.println("create header of solver " + solverName + " ... ok");
        headerWriter.close();
      }

      if (userImplementationFile.exists()) {
        System.out.println("user's implementation file of solver " + solverName
            + " ... does exist already. Is not overwritten");
      } else {
        java.io.BufferedWriter userImplementationWriter =
            new java.io.BufferedWriter(new java.io.FileWriter(userImplementationFile));
        solver.writeUserImplementation(userImplementationWriter, solverName, _projectName);
        System.out.println(
            "create user implementation template of solver " + solverName + " ... please complete");
        userImplementationWriter.close();
      }

      if (generatedImplementationFile.exists()) {
        System.out.println("generated implementation file of solver " + solverName
            + " ... does exist already. Is overwritten");
      }

      java.io.BufferedWriter generatedImplementationWriter =
          new java.io.BufferedWriter(new java.io.FileWriter(generatedImplementationFile));
      solver.writeGeneratedImplementation(generatedImplementationWriter, solverName, _projectName);
      System.out.println("create generated implementation of solver " + solverName + " ... ok");
      generatedImplementationWriter.close();
    } catch (Exception exc) {
      System.err.println("ERROR: " + exc.toString());
      valid = false;
    }
  }



  @Override
  public void inAFiniteVolumesSolver(AFiniteVolumesSolver node) {
    String solverName = node.getName().toString().trim();
    
    if (_definedSolvers.contains(solverName)) {
      System.err.println( "ERROR: Solver " + solverName + " multiply defined" );
      valid = false;
    }
    else {
      _definedSolvers.add(solverName);
    }

    java.io.File headerFile = new java.io.File(
        _directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/" + solverName + ".h");
    java.io.File userImplementationFile = new java.io.File(
        _directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/" + solverName + ".cpp");
    java.io.File userPDEFile = null;
    java.io.File userTypesDefFile = null;
    java.io.File generatedImplementationFile =
        new java.io.File(_directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/"
            + solverName + "_generated.cpp");


    boolean isFortran = false;
    if (node.getLanguage().getText().trim().equals("C")) {
      isFortran = false;
    } else if (node.getLanguage().getText().trim().equals("Fortran")) {
      isFortran = true;
      userPDEFile =
          new java.io.File(_directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/PDE.f90");
      userTypesDefFile = new java.io.File(
          _directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/typesDef.f90");
    } else {
      System.err.println("ERROR: unknown language for solver " + node.getName().getText()
          + ". Supported languages are C and Fortran");
      valid = false;
      return;
    }

    String kernel = node.getKernel().toString().trim();

    int numberOfVariables  = Integer.parseInt(node.getVariables().toString().trim());
    int numberOfParameters = Integer.parseInt(node.getParameters().toString().trim());
    int patchSize          = Integer.parseInt(node.getPatchSize().toString().trim());
    boolean hasConstants   = node.getConstants()!=null;

    if (numberOfParameters != 0) {
      System.err.println("ERROR: At the moment, parameters are not yet supported. " + 
          " Please add the parameters as additional quantities to your PDE formulation.");
      valid = false;
      return;
    }
    
    eu.exahype.solvers.Solver solver = null;

    if (isFortran && kernel.equals( eu.exahype.solvers.UserDefinedFiniteVolumesinFortran.Identifier )) {
      solver = new eu.exahype.solvers.UserDefinedFiniteVolumesinFortran(_dimensions,numberOfVariables, numberOfParameters, patchSize, _enableProfiler, hasConstants);
    }
    if (!isFortran && kernel.equals( eu.exahype.solvers.UserDefinedFiniteVolumesinC.Identifier )) {
      solver = new eu.exahype.solvers.UserDefinedFiniteVolumesinC(_dimensions,numberOfVariables, numberOfParameters, patchSize, _enableProfiler, hasConstants);
    }
    if (isFortran && kernel.equals( eu.exahype.solvers.GenericFiniteVolumesMUSCLinFortran.Identifier )) {
      solver = new eu.exahype.solvers.GenericFiniteVolumesMUSCLinFortran(_dimensions,numberOfVariables, numberOfParameters, patchSize, _enableProfiler, hasConstants);
    }
    if (!isFortran && kernel.equals( eu.exahype.solvers.GenericFiniteVolumesMUSCLinC.Identifier )) {
      solver = new eu.exahype.solvers.GenericFiniteVolumesMUSCLinC(_dimensions,numberOfVariables, numberOfParameters, patchSize, _enableProfiler, hasConstants);
    }
    if (!isFortran && kernel.equals( eu.exahype.solvers.GenericFiniteVolumesGodunovInC.Identifier )) {
    	solver = new eu.exahype.solvers.GenericFiniteVolumesGodunovInC(_dimensions,numberOfVariables, numberOfParameters, patchSize, _enableProfiler, hasConstants);
    }

    if (solver == null) {
      System.err.println("creation solver " + solverName + " ... failed as kernel " + kernel
          + " for language " + node.getLanguage().getText().trim() + " is not supported");
      valid = false;
      return;
    }

    try {
      // =====================
      // Write all the headers
      // =====================
      if (headerFile.exists()) {
        System.out.println("create header of solver " + solverName + " ... header "
            + headerFile.getAbsoluteFile()
            + " does exist already. Remove to allow toolkit to regenerate it (changes will be lost)");
      } else {
        java.io.BufferedWriter headerWriter =
            new java.io.BufferedWriter(new java.io.FileWriter(headerFile));
        solver.writeHeader(headerWriter, solverName, _projectName);
        System.out.println("create header of solver " + solverName + " ... ok");
        headerWriter.close();
      }

      if (userImplementationFile.exists()) {
        System.out.println("user's implementation file of solver " + solverName
            + " ... does exist already. Is not overwritten");
      } else {
        java.io.BufferedWriter userImplementationWriter =
            new java.io.BufferedWriter(new java.io.FileWriter(userImplementationFile));
        solver.writeUserImplementation(userImplementationWriter, solverName, _projectName);
        System.out.println(
            "create user implementation template of solver " + solverName + " ... please complete");
        userImplementationWriter.close();
      }


      if (generatedImplementationFile.exists()) {
        System.out.println("generated implementation file of solver " + solverName
            + " ... does exist already. Is overwritten");
      }

      java.io.BufferedWriter generatedImplementationWriter =
          new java.io.BufferedWriter(new java.io.FileWriter(generatedImplementationFile));
      solver.writeGeneratedImplementation(generatedImplementationWriter, solverName, _projectName);
      System.out.println("create generated implementation of solver " + solverName + " ... ok");
      generatedImplementationWriter.close();
    } catch (Exception exc) {
      System.err.println("ERROR: " + exc.toString());
      valid = false;
    }
  }
}
