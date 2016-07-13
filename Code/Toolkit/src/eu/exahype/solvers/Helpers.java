package eu.exahype.solvers;

import java.io.IOException;

public class Helpers {
  public static void writeMinimalADERDGSolverHeader(
      String solverName, java.io.BufferedWriter writer, String projectName) throws IOException {
    writeHeaderCopyright(writer);
    writeHeaderIncludesAndDefines(writer, solverName, projectName);
    writeHeaderMinimalADERDGClassSignature(writer, solverName, projectName);
  }

  /**
   * Creates all the public operations that are mandatory for any solver.
   */
  private static void writeHeaderMinimalADERDGClassSignature(
      java.io.BufferedWriter writer, String solverName, String projectName) throws IOException {
    writer.write(
        "class " + projectName + "::" + solverName + ": public exahype::solvers::Solver {\n");
    writer.write("  public:\n");
    writer.write("    " + solverName + "(int kernelNumber, double maximumMeshSize, exahype::solvers::Solver::TimeStepping timeStepping, std::unique_ptr<exahype::profilers::Profiler> profiler);\n");

    writer.write(
        "    void spaceTimePredictor(double* lQi, double* lFi, double* lQhi, double* lFhi, double* lQhbnd, double* lFhbnd, const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx, const double dt ) override; \n");
    writer.write(
        "    void solutionUpdate(double* luh, const double* const lduh, const double dt) override;\n");
    writer.write(
        "    void volumeIntegral(double* lduh, const double* const lFhi, const tarch::la::Vector<DIMENSIONS,double>& dx) override;\n");
    writer.write(
        "    void surfaceIntegral(double* lduh, const double* const lFhbnd, const tarch::la::Vector<DIMENSIONS,double>& dx) override;\n");
    writer.write(
        "    void riemannSolver(double* FL, double* FR, const double* const QL, const double* const QR, const double dt, const int normalNonZeroIndex) override;\n");
    writer.write(
        "    double stableTimeStepSize(const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx ) override;\n");
    writer.write(
        "    void solutionAdjustment(double *luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) override;\n");
    writer.write(
        "    bool hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t) override;\n");
    writer.write(
        "    exahype::solvers::Solver::RefinementControl refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) override;\n");
    writer.write(
        "    void faceUnknownsProlongation(double* lQhbndFine,double* lFhbndFine,const double* lQhbndCoarse,const double* lFhbndCoarse,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex) override;\n");
    writer.write(
        "    void faceUnknownsRestriction(double* lQhbndCoarse,double* lFhbndCoarse,const double* lQhbndFine,const double* lFhbndFine,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex) override;\n");
    writer.write(
        "    void volumeUnknownsProlongation(double* luhFine, const double* luhCoarse, const int coarseGridLevel, const int fineGridLevel, const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) override;\n");
    writer.write(
        "    void volumeUnknownsRestriction(double* luhCoarse, const double* luhFine, const int coarseGridLevel, const int fineGridLevel, const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) override;\n");
  }

  /**
   * Write header with ExaHyPE copyright. Should be inserted for any solver's
   * header.
   */
  private static void writeHeaderCopyright(java.io.BufferedWriter writer) throws IOException {
    writer.write("// This file is generated by the ExaHyPE toolkit.\n");
    writer.write("// Please do not modify - it will be overwritten by the next\n");
    writer.write("// ExaHyPE toolkit call.\n");
    writer.write("// \n");
    writer.write("// ========================\n");
    writer.write("//   www.exahype.eu\n");
    writer.write("// ========================\n");
  }

  /**
   * Adds all the default includes of any solver as well as the solver define.
   * Is used by all solvers.
   */
  private static void writeHeaderIncludesAndDefines(
      java.io.BufferedWriter writer, String solverName, String projectName) throws IOException {
    writer.write("\n\n");
    writer.write("#include <memory>\n\n");
    writer.write("#include \"exahype/profilers/Profiler.h\"\n");
    writer.write("#include \"exahype/solvers/Solver.h\"");
    writer.write("\n\n\n");

    writer.write("namespace " + projectName + "{\n");
    writer.write("  class " + solverName + ";\n");
    writer.write("}\n\n\n");
  }

  public static void writeMinimalADERDGSolverUserImplementation(String solverName,
      java.io.BufferedWriter writer, String projectName, int numberOfVariables, int numberOfParameters, int order)
      throws IOException {
    writer.write("#include \"" + solverName + ".h\"\n\n");
    writer.write("#include <memory>\n\n");
    writer.write(projectName + "::" + solverName + "::" + solverName + "(const std::string& identifier, int kernelNumber, int numberOfVariables, int numberOfParameters, int nodesPerCoordinateAxis, double maximumMeshSize, exahype::solvers::Solver::TimeStepping timeStepping, std::unique_ptr<exahype::profilers::Profiler> profiler):\n");
    writer.write("  exahype::solvers::Solver("
        + " "+solverName+", exahype::solvers::Solver::Type::ADER_DG, kernelNumber, "+numberOfVariables+", "+numberOfParameters+", "+(order+1)+", maximumMeshSize, timeStepping, std::move(profiler)) {\n");
    writer.write("  // @todo Please implement/augment if required\n");
    writer.write("}\n");
    writer.write("\n\n\n");
  }

  static public void invokeCodeGenerator(String solverName, int numberOfUnknowns, int numberOfParameters, int order,
      boolean isLinear, int dimensions, String microarchitecture, String pathToLibxsmm)
      throws IOException {
    String currentDirectory = System.getProperty("user.dir");
    java.nio.file.Path pathToCodeGenerator =
        java.nio.file.Paths.get(currentDirectory + "/Miscellaneous/CodeGenerator/Driver.py");
    if (java.nio.file.Files.notExists(pathToCodeGenerator)) {
      System.err.println("ERROR: Code generator not found. Can't generated optimised kernels.");
      return;
    }

    String numericsParameter = isLinear ? "linear" : "nonlinear";

    // set up the command to execute the code generator
    String args = " " + solverName + " " + numberOfUnknowns + " " + order + " "
        + Integer.toString(dimensions) + " " + numericsParameter + " " + microarchitecture + " "
        + pathToLibxsmm + " "
        + "--precision=DP"; // double precision

    String bashCommand = "python " + pathToCodeGenerator + args;

    Runtime runtime = Runtime.getRuntime();

    // execute the command line program
    Process codeGenerator = runtime.exec(bashCommand);

    // capture any output that is produced by the code generator and print it line-by-line
    java.io.InputStream stdout = codeGenerator.getInputStream();
    java.io.BufferedReader stdoutReader =
        new java.io.BufferedReader(new java.io.InputStreamReader(stdout));
    String line = "";
    while ((line = stdoutReader.readLine()) != null) {
      System.out.println("CodeGenerator: " + line);
    }
    java.io.InputStream stderr = codeGenerator.getErrorStream();
    java.io.BufferedReader stderrReader =
        new java.io.BufferedReader(new java.io.InputStreamReader(stderr));
    while ((line = stderrReader.readLine()) != null) {
      System.out.println("CodeGenerator: " + line);
    }
  }
}
