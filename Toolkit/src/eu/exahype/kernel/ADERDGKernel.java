package eu.exahype.kernel;

import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import eu.exahype.node.PIds;
import eu.exahype.node.AIds;
import eu.exahype.node.AIdentifierId;
import eu.exahype.node.AIdentifierWithMultId;

import eu.exahype.node.AAderdgSolver;
import eu.exahype.node.ALimitingAderdgSolver;
import eu.exahype.node.PSolver;

/**
 * Circumscribes an ADERDG Kernel
 * 
 * For a description how to add new variants consult SolverFactory.
 */
public class ADERDGKernel {
  
   
  /**
   * Configuration parameter: id of the options
   */
  public static final String GENERIC_OPTION_ID           = "generic";
  public static final String OPTIMISED_OPTION_ID         = "optimised";

  public static final String LINEAR_OPTION_ID            = "linear";
  public static final String NONLINEAR_OPTION_ID         = "nonlinear";
  public static final String USER_DEFINED_OPTION_ID      = "user";
  public static final String LEGENDRE_OPTION_ID          = "gausslegendre";
  public static final String LOBATTO_OPTION_ID           = "gausslobatto";
    
  public static final String FLUX_OPTION_ID              = "flux";
  public static final String SOURCE_OPTION_ID            = "source";
  public static final String NCP_OPTION_ID               = "ncp";
  public static final String POINTSOURCE_OPTION_ID       = "pointsources";
  public static final String MATERIALPARAMETER_OPTION_ID = "materialparameters";

  public static final String NO_TIME_AVG_OPTION_ID       = "notimeavg";
  public static final String PATCHWISE_ADJUST_OPTION_ID  = "patchwiseadjust";
  
  private Set<String> type;
  private Map<String, Integer> terms;
  private Set<String> optimisation;
  
  public ADERDGKernel(PSolver solver) throws IllegalArgumentException {
    if(solver instanceof AAderdgSolver) {
      type = parseIds(((AAderdgSolver) solver).getKernelType());
      terms = parseIdsToMap(((AAderdgSolver) solver).getKernelTerms());
      optimisation = parseIds(((AAderdgSolver) solver).getKernelOpt());
    } else if(solver instanceof ALimitingAderdgSolver) {
      type = parseIds(((ALimitingAderdgSolver) solver).getKernelType());
      terms = parseIdsToMap(((ALimitingAderdgSolver) solver).getKernelTerms());
      optimisation = parseIds(((ALimitingAderdgSolver) solver).getKernelOpt());
    } else {
      throw new IllegalArgumentException("No kernel definition found");
    }
    
    validate();
  }
  
  //return null on error, use only after the program should already have failed with invalid kernel
  public static ADERDGKernel noExceptionContructor(PSolver solver) {
    try {
      return new ADERDGKernel(solver);
    } catch(IllegalArgumentException e) {
      return null;
    }
  }
  
  private static Set<String> parseIds(PIds idsRaw) {
    return ((AIds)idsRaw).getId().stream().map(e -> (e instanceof AIdentifierId) ? ((AIdentifierId)e).getValue().getText() : ((AIdentifierWithMultId)e).getValue().getText()).collect(Collectors.toSet());
  }
  
  private static Map<String, Integer> parseIdsToMap(PIds idsRaw) {
    return ((AIds)idsRaw).getId().stream().collect(Collectors.toMap(
      (e -> (e instanceof AIdentifierId) ? ((AIdentifierId)e).getValue().getText() : ((AIdentifierWithMultId)e).getValue().getText()),
      (e -> (e instanceof AIdentifierId) ? Integer.valueOf(-1) : Integer.valueOf(((AIdentifierWithMultId)e).getMultiplicity().getText()))
    ));
  }
  
  private void validate() throws IllegalArgumentException {
    if(!type.contains(LINEAR_OPTION_ID) ^ type.contains(NONLINEAR_OPTION_ID)) {//should be only one
      throw new IllegalArgumentException("nonlinear or linear not specified or both specified in the kernel type");
    }
    if(usePointSource() && getNumberOfPointSources() < 0) {
      throw new IllegalArgumentException("point sources used but number not specified! In the specification file, use "+POINTSOURCE_OPTION_ID+":X, with X the number of point sources.");
    }
  }

  public enum KernelType {
    GenericADERDG,
    OptimisedADERDG,
    UserDefined,
    Unknown
  }

  public boolean isLinear() throws IllegalArgumentException {
    return type.contains(LINEAR_OPTION_ID);
  }

  public KernelType getKernelType() {
    if (optimisation.contains(OPTIMISED_OPTION_ID)) {
     return KernelType.OptimisedADERDG;
   }
      
   // default kernel - must be last   
   if(optimisation.contains(GENERIC_OPTION_ID)) {
      return KernelType.GenericADERDG;
  }
  
    return  KernelType.Unknown;
  }

  public boolean usesOptimisedKernels() {
    // assert: !optimisation.contains(GENERIC_OPTION_ID)
    return optimisation.contains(OPTIMISED_OPTION_ID);
  }

  public boolean useFlux() {
    return terms.containsKey(FLUX_OPTION_ID);
  }
  
  public boolean useSource() {
    return terms.containsKey(SOURCE_OPTION_ID);
  }
  
  public boolean useNCP() {
    return terms.containsKey(NCP_OPTION_ID);
  }
  
  public boolean usePointSource() {
    return terms.containsKey(POINTSOURCE_OPTION_ID);
  }
  
  public boolean useMaterialParameterMatrix() {
    return terms.containsKey(MATERIALPARAMETER_OPTION_ID);
  }
  
  public boolean noTimeAveraging() {
    return optimisation.contains(NO_TIME_AVG_OPTION_ID);
  }
  
  public boolean patchwiseAdjust() {
    return optimisation.contains(PATCHWISE_ADJUST_OPTION_ID);
  }
  
  public int getNumberOfPointSources() {
    if(usePointSource()) {
      return terms.get(POINTSOURCE_OPTION_ID);
    }
    
    return -1;
  }
  
  //(type: [...], terms: [...], opt: [...])
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("(type: [");
    for(String s : type) {
      sb.append(s);
      sb.append(", ");
    }
    sb.deleteCharAt(sb.length()-2);
    sb.append("], terms: [");
    for(String s : terms.keySet()) {
      sb.append(s);
      sb.append(", ");
    }
    sb.deleteCharAt(sb.length()-2);
    sb.append("], opt: [");
    for(String s : optimisation) {
      sb.append(s);
      sb.append(", ");
    }
    sb.deleteCharAt(sb.length()-2);
    sb.append("])");
    
    return sb.toString();
  }
  
}