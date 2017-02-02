package eu.exahype;

import java.io.IOException;

import eu.exahype.analysis.DepthFirstAdapter;
import eu.exahype.node.AAderdgSolver;
import eu.exahype.node.AComputationalDomain;
import eu.exahype.node.AFiniteVolumesSolver;
import eu.exahype.node.ALimitingAderdgSolver;
import eu.exahype.node.AProfiling;
import eu.exahype.node.AProject;
import eu.exahype.node.PSolver;
import eu.exahype.solvers.Solver;
import eu.exahype.solvers.SolverFactory;
import eu.exahype.variables.Variables;

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
    _projectName     = node.getName().getText();
    _definedSolvers  = new java.util.HashSet<String>();

    if (node.getSolver().size() == 0) {
      System.out.println("there are no solvers in the specification file ... nothing to be done");
    }

    // Only one optimised solver can be used (optimised kernel would be overwritten by the latest solver otherwise)
    if (node.getSolver().size() > 1) {
      int optimisedCount = 0;
      for(PSolver psolver : node.getSolver()) {
        if(psolver instanceof AAderdgSolver) {
          AAderdgSolver asolver = (AAderdgSolver) psolver;
          if(    asolver.getKernel().getText().equals( eu.exahype.solvers.OptimisedFluxesNonlinearADER_DGinC.Identifier )
              || asolver.getKernel().getText().equals( eu.exahype.solvers.OptimisedFluxesLinearADER_DGinC.Identifier )
              ){
            optimisedCount++;
          }
        }
      }
      if(optimisedCount > 1) {
        System.err.println("ERROR: Only one optimised solver can be used at a time. Currently "+optimisedCount+" are defined.");
        valid = false;
        return;
      }
    }

    _microarchitecture = node.getArchitecture().getText().toLowerCase();
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
      _pathToLibxsmm = node.getLibxsmmPath().getText();
    }
  };

  @Override
  public void inAComputationalDomain(AComputationalDomain node) {
    _dimensions = Integer.parseInt( node.getDimension().getText() );
    if (_dimensions!=2 && _dimensions!=3) {
      System.err.println( "ERROR: dimension has to be either 2 or 3.");
    }
  }


  @Override
  public void inAProfiling(AProfiling node) {
    _enableProfiler = !node.getProfiler().getText().equals("NoOpProfiler");
  };

  private boolean validate(
      Variables variables,
      int order,String kernel,String language,
      String solverName,eu.exahype.solvers.Solver solver) {
    if (_definedSolvers.contains(solverName)) {
      System.err.println( "ERROR: Solver " + solverName + " multiple definitions." );
      return false;
    }
    
    if (!language.equals("C") && !language.equals("Fortran")) {
      System.err.println("ERROR: unknown language for solver " + solverName
          + ". Supported languages are C and Fortran");
      return false;
    }
    
    if (variables.getNumberOfParameters() != 0) {
      System.err.println("ERROR: At the moment, parameters are not supported. " + 
          " Please add the parameters as additional quantities to your PDE formulation.");
      return false;
    }
    if (order < 1 || order > 9) {
      System.err.println("ERROR: Only polynomial degrees of 1..9 are supported.");
      return false;
    }
    if (solver == null) {
      System.err.println("ERROR: creation solver " + solverName + " ... failed as kernel " + kernel
          + " for language " + language + " is not supported");
      return false;
    }
    return true;
  }
  
  @Override
  public void inAAderdgSolver(AAderdgSolver node) {
    String solverName    = node.getName().getText();
    String  kernel       = node.getKernel().getText();
    String  language     = node.getLanguage().getText();
    int     order        = Integer.parseInt(node.getOrder().getText());
    boolean hasConstants = node.getConstants()!=null;
    Variables variables  = new Variables(node, _dimensions);
    boolean isFortran    = language.equals("Fortran");
    
    SolverFactory solverFactory = new SolverFactory(_dimensions, _enableProfiler, _microarchitecture, _pathToLibxsmm);
    eu.exahype.solvers.Solver solver = solverFactory.createADERDGSolver(
        kernel, isFortran, variables.getNumberOfVariables(), variables.getNumberOfParameters(), order, hasConstants);

    valid = validate(variables,order,kernel,language,solverName,solver);
    
    if (valid) {
      _definedSolvers.add(solverName);

      // write the files
      try {
//        TODO(Dominic): Unused
//        java.io.File userPDEFile = null;
//        java.io.File userTypesDefFile = null;
//        
//        if (isFortran) {
//          userPDEFile =
//              FileSearch.relocatableFile(_directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/PDE.f90");
//          userTypesDefFile = FileSearch.relocatableFile(
//              _directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/typesDef.f90");
//        }
        
        
        tryWriteSolverHeader(solver, solverName);
        tryWriteSolverUserImplementation(solver, solverName);
        tryWriteSolverGeneratedImplementation(solver, solverName);

        if (solver.supportsVariables()) {
          tryWriteVariablesHeader(variables, solverName);
        }
      } catch (Exception exc) {
        System.err.println("ERROR: " + exc.toString());
        valid = false;
      }
    }
  }

  @Override
  public void inAFiniteVolumesSolver(AFiniteVolumesSolver node) {
    String solverName    = node.getName().getText();
    String  kernel       = node.getKernel().getText();
    String  language     = node.getLanguage().getText();
    int     patchSize    = Integer.parseInt(node.getPatchSize().getText());
    boolean hasConstants = node.getConstants()!=null;
    Variables variables  = new Variables(node, _dimensions);
    boolean isFortran    = language.equals("Fortran");
    
    SolverFactory solverFactory = new SolverFactory(_dimensions, _enableProfiler, _microarchitecture, _pathToLibxsmm);
    eu.exahype.solvers.Solver solver = solverFactory.createFiniteVolumesSolver(
        kernel, isFortran, variables.getNumberOfVariables(), variables.getNumberOfParameters(), patchSize, hasConstants);

    valid = validate(variables,1/*patchSize is always supported*/,kernel,language,solverName,solver);
    
    if (valid) {
      _definedSolvers.add(solverName);

      // write the files
      try {
// TODO(Dominic): Unused
//        java.io.File userPDEFile = null;
//        java.io.File userTypesDefFile = null;
//        
//        if (isFortran) {
//          userPDEFile =
//              FileSearch.relocatableFile(_directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/PDE.f90");
//          userTypesDefFile = FileSearch.relocatableFile(
//              _directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/typesDef.f90");
//        }
        
        
        tryWriteSolverHeader(solver, solverName);
        tryWriteSolverUserImplementation(solver, solverName);
        tryWriteSolverGeneratedImplementation(solver, solverName);

        if (solver.supportsVariables()) {
          tryWriteVariablesHeader(variables, solverName);
        }
      } catch (Exception exc) {
        System.err.println("ERROR: " + exc.toString());
        valid = false;
      }
    }
  }
  
  @Override
  public void inALimitingAderdgSolver(ALimitingAderdgSolver node) {
    String solverName    = node.getName().getText();
    String  kernel       = node.getKernel().getText();
    String  language     = node.getLanguage().getText();
    int     order        = Integer.parseInt(node.getOrder().getText());
    boolean hasConstants = node.getConstants()!=null;
    Variables variables  = new Variables(node, _dimensions);
    boolean isFortran    = language.equals("Fortran");
    
    SolverFactory solverFactory = new SolverFactory(_dimensions, _enableProfiler, _microarchitecture, _pathToLibxsmm);
    eu.exahype.solvers.Solver solver = solverFactory.createADERDGSolver(
        kernel, isFortran, variables.getNumberOfVariables(), variables.getNumberOfParameters(), order, hasConstants);

    valid = validate(variables,order,kernel,language,solverName,solver);
    
    if (valid) {
      _definedSolvers.add(solverName);

      // write the files
      try {
//        TODO(Dominic): Unused
//        java.io.File userPDEFile = null;
//        java.io.File userTypesDefFile = null;
//        
//        if (isFortran) {
//          userPDEFile =
//              FileSearch.relocatableFile(_directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/PDE.f90");
//          userTypesDefFile = FileSearch.relocatableFile(
//              _directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/typesDef.f90");
//        }
        
        
        tryWriteSolverHeader(solver, solverName);
        tryWriteSolverUserImplementation(solver, solverName);
        tryWriteSolverGeneratedImplementation(solver, solverName);

        if (solver.supportsVariables()) {
          tryWriteVariablesHeader(variables, solverName);
        }
      } catch (Exception exc) {
        System.err.println("ERROR: " + exc.toString());
        valid = false;
      }
    }
  }
  
  private void tryWriteSolverHeader(Solver solver,String solverName) throws IOException {
    java.io.File solverHeaderFile = FileSearch.relocatableFile(
        _directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/" + solverName + ".h");
    
    if (solverHeaderFile.exists()) {
      System.out.println("create header of solver " + solverName + " ... header "
          + solverHeaderFile.getAbsoluteFile()
          + " does exist already. Remove to allow toolkit to regenerate it (changes will be lost)");
    } else {
      java.io.BufferedWriter headerWriter =
          new java.io.BufferedWriter(new java.io.FileWriter(solverHeaderFile));
      solver.writeHeader(headerWriter, solverName, _projectName);
      System.out.println("create header of solver " + solverName + " ... ok");
      headerWriter.close();
    }
  }
  
  private void tryWriteSolverUserImplementation(Solver solver, String solverName) throws IOException {
    java.io.File solverUserImplementationFile = FileSearch.relocatableFile(
        _directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/" + solverName + ".cpp");
    
    if (solverUserImplementationFile.exists()) {
      System.out.println("user's implementation file of solver " + solverName
          + " ... does exist already. Is not overwritten");
    } else {
      java.io.BufferedWriter userImplementationWriter =
          new java.io.BufferedWriter(new java.io.FileWriter(solverUserImplementationFile));
      solver.writeUserImplementation(userImplementationWriter, solverName, _projectName);
      System.out.println(
          "create user implementation template of solver " + solverName + " ... please complete");
      userImplementationWriter.close();
    }
  }
  
  private void tryWriteSolverGeneratedImplementation(Solver solver, String solverName) throws IOException {
    java.io.File solverGeneratedImplementationFile = FileSearch.relocatableFile(
        _directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/" + solverName + "_generated.cpp");
    
    if (solverGeneratedImplementationFile.exists()) {
      System.out.println("generated implementation file of solver " + solverName
          + " ... does exist already. Is overwritten");
    }

    java.io.BufferedWriter generatedImplementationWriter =
        new java.io.BufferedWriter(new java.io.FileWriter(solverGeneratedImplementationFile));
    solver.writeGeneratedImplementation(generatedImplementationWriter, solverName, _projectName);
    System.out.println("create generated implementation of solver " + solverName + " ... ok");
    generatedImplementationWriter.close();
  }
  
  
  private void tryWriteVariablesHeader(Variables variables,String solverName) throws IOException {
    java.io.File solverHeaderFile = FileSearch.relocatableFile(
        _directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/" + solverName + "_Variables.h");
    
    if (solverHeaderFile.exists()) {
      java.io.BufferedWriter headerWriter =
          new java.io.BufferedWriter(new java.io.FileWriter(solverHeaderFile));
      variables.writeHeader(headerWriter, solverName, _projectName);
      System.out.println("create header of variables for solver " + solverName + " ... header "
          + solverHeaderFile.getAbsoluteFile()
          + " does exist already and will be overwritten");
      headerWriter.close();
    } else {
      java.io.BufferedWriter headerWriter =
          new java.io.BufferedWriter(new java.io.FileWriter(solverHeaderFile));
      variables.writeHeader(headerWriter, solverName, _projectName);
      System.out.println("create header of variables for solver " + solverName + " ... ok");
      headerWriter.close();
    }
  }
}
