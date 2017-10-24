package eu.exahype.solvers;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Set;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

// template engine
import minitemp.Context;
import minitemp.TemplateEngine;

import eu.exahype.CodeGeneratorHelper;
import eu.exahype.kernel.ADERDGKernel;
import eu.exahype.io.IOUtils;


public class OptimisedADERDG implements Solver {
  //Internal states
  //--------------- 
  private String         solverName;
  private Context        context;
  private TemplateEngine templateEngine;
  
  private boolean        useConverterDebug;
  
  public OptimisedADERDG(String projectName, String solverName, int dimensions, int numberOfVariables, int numberOfParameters, Set<String> namingSchemeNames,
      int order,String microarchitecture, boolean enableProfiler, boolean enableDeepProfiler, boolean hasConstants, ADERDGKernel kernel) 
      throws IOException, IllegalArgumentException {    
    
    this.solverName                 = solverName;
    this.useConverterDebug          = kernel.useConverterDebug();
    
    final boolean isLinear          = kernel.isLinear();
    final boolean useFlux           = kernel.useFlux();
    final boolean useSource         = kernel.useSource();
    final boolean useNCP            = kernel.useNCP();
    final boolean usePointSource    = kernel.usePointSource();
    final boolean useMaterialParam  = kernel.useMaterialParameterMatrix();
    final boolean noTimeAveraging   = kernel.noTimeAveraging();
    final boolean patchwiseAdjust   = kernel.patchwiseAdjust();
    final int numberOfPointSources  = kernel.getNumberOfPointSources();
    
    //generate the optimised kernel, can throw IOException
    final String optKernelPath = CodeGeneratorHelper.getInstance().invokeCodeGenerator(projectName, solverName, numberOfVariables, numberOfParameters, order, isLinear, dimensions,
        microarchitecture, enableDeepProfiler, useFlux, useSource, useNCP, noTimeAveraging);
    final String optNamespace = CodeGeneratorHelper.getInstance().getNamespace(projectName, solverName);
    
    templateEngine = new TemplateEngine();
    context = new Context();
    
    //String
    context.put("project"           , projectName);
    context.put("solver"            , solverName);
    context.put("abstractSolver"    , getAbstractSolverName());
    context.put("optKernelPath"     , optKernelPath);
    context.put("optNamespace"      , optNamespace);
    
    //int
    context.put("dimensions"        , dimensions);
    context.put("order"             , order);
    context.put("numberOfVariables" , numberOfVariables);
    context.put("numberOfParameters", numberOfParameters);
    context.put("numberOfPointSources", numberOfPointSources);
    
    //boolean
    context.put("enableProfiler"    , enableProfiler);
    context.put("enableDeepProfiler", enableDeepProfiler);
    context.put("hasConstants"      , hasConstants);
    context.put("isLinear"          , isLinear);
    context.put("useFlux"           , useFlux);
    context.put("useSource"         , useSource);
    context.put("useNCP"            , useNCP);
    context.put("usePointSource"    , usePointSource);
    context.put("useMaterialParam"  , useMaterialParam);
    context.put("noTimeAveraging"   , noTimeAveraging);
    context.put("patchwiseAdjust"   , patchwiseAdjust);
    
    //boolean as String
    context.put("useFlux_s"         , boolToTemplate(useFlux));
    context.put("useSource_s"       , boolToTemplate(useSource));
    context.put("useNCP_s"          , boolToTemplate(useNCP));

    //Set<String>
    context.put("namingSchemes"     , namingSchemeNames.stream().map(s -> s.substring(0, 1).toUpperCase()+s.substring(1)).collect(Collectors.toSet())); //capitalize
    
    //List<Integer> , range used by for loops
    context.put("range_0_nDim"      , IntStream.range(0, dimensions)                          .boxed().collect(Collectors.toList()));
    context.put("range_0_nVar"      , IntStream.range(0, numberOfVariables)                   .boxed().collect(Collectors.toList()));
    context.put("range_0_nVarParam" , IntStream.range(0, numberOfVariables+numberOfParameters).boxed().collect(Collectors.toList()));
    
    //TODO JMG: linear kernels unsupported for now
    if(isLinear) 
      throw new IllegalArgumentException("Linear kernels not supported yet");
  }
    
  @Override
  public String getSolverName() {
    return solverName;
  }
  
  private String getAbstractSolverName() {
    return "Abstract"+getSolverName();
  }
  
  private String boolToTemplate(boolean b) {
    return b? "true" : "false";
  }
  
  @Override
  public void writeHeader(java.io.BufferedWriter writer) throws IOException, IllegalArgumentException {
    //reuse the generic template
	  final String template = IOUtils.convertRessourceContentToString("eu/exahype/solvers/templates/GenericADERDGSolverHeader.template");
	  writer.write(templateEngine.render(template, context));
  }
  
  @Override
  public void writeUserImplementation(java.io.BufferedWriter writer) throws java.io.IOException, IllegalArgumentException {
      //reuse the generic template
    final String template = IOUtils.convertRessourceContentToString("eu/exahype/solvers/templates/GenericADERDGSolverInCUserCode.template");
    writer.write(templateEngine.render(template, context));
  }
  
  @Override
  public void writeAbstractHeader(java.io.BufferedWriter writer) throws java.io.IOException, IllegalArgumentException {      
    final String template = IOUtils.convertRessourceContentToString("eu/exahype/solvers/templates/AbstractOptimisedADERDGSolverHeader.template");
    writer.write(templateEngine.render(template, context));
  }
  
  @Override
  public void writeAbstractImplementation(java.io.BufferedWriter writer) throws java.io.IOException, IllegalArgumentException {
    
    final String template = IOUtils.convertRessourceContentToString(
        (useConverterDebug ? 
            "eu/exahype/solvers/templates/AbstractOptimisedADERDGSolverImplementation_withConverter.template" 
          : "eu/exahype/solvers/templates/AbstractOptimisedADERDGSolverImplementation.template"
        )
      ); 
    writer.write(templateEngine.render(template, context));
  }
  
  @Override
  public boolean supportsVariables() {
    return true;
  }
}
