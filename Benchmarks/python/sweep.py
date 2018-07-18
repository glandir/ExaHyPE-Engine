#!/usr/bin/env python3
"""
.. module:: sweep
  :platform: Unix, Windows, Mac
  :synopsis: Generate benchmark suites for ExaHyPE.
   
.. moduleauthor:: Dominic Etienne Charrier <dominic.e.charrier@durham.ac.uk>, 

:synopsis: Generate benchmark suites for ExaHyPE.
"""
def parseArgument(argv,i):
    if i<len(argv):
        return argv[i]
    else:
        return None

def haveToPrintHelpMessage(argv):
    """
    Check if we have to print a help message.
    """
    result = parseArgument(argv,2) not in subprograms or \
             parseArgument(argv,1)==None
    for arg in argv:
        result = result or ( arg=="-help" or arg=="-h" )
    return result

def dictProduct(dicts):
    """
    Computes the Cartesian product of a dictionary of lists as 
    a list of dictionaries.
    
    Gladly copied this code from:
    https://stackoverflow.com/questions/5228158/cartesian-product-of-a-dictionary-of-lists
    
    Example input:
    options = {"number": [1,2,3], "color": ["orange","blue"] }
    
    Example output:
    [ {"number": 1, "color": "orange"},
      {"number": 1, "color": "blue"},
      {"number": 2, "color": "orange"},
      {"number": 2, "color": "blue"},
      {"number": 3, "color": "orange"},
      {"number": 3, "color": "blue"}
    ]
    """
    return (dict(zip(dicts, x)) for x in itertools.product(*dicts.values()))

def hashDictionary(dictionary):
    """
    Hash a dictionary. Sort the dictionary according
    to the keys beforehand.
    """
    chain = ""
    for key,value in sorted(dictionary.items()):
        chain += key+","+value+";"
    
    return hashlib.md5(chain.encode()).hexdigest()[0:8]

def clean(subFolder=""):
    """
    Clean the complete output folder or just a subfolder
    if specified.
    """
    folder = exahypeRoot+"/"+outputPath+"/"+subFolder
    print("rm -r "+folder)
    subprocess.call("rm -r "+folder, shell=True)


def renderSpecFile(templateBody,parameterDict,ranksPerNode,coresPerRank,consumerTasks):
    renderedFile = templateBody
    
    context = dict(parameterDict)
    
    consistent = True
    # verify mandatory options file parameters can be found in template file
    global createdFirstSpecFile
    if not createdFirstSpecFile:
        keysInTemplate = [m.group(2) for m in re.finditer("(\{\{((\w|-)+)\}\})",templateBody)]
        for key in context:
            if key not in keysInTemplate:
                consistent = False
                print("ERROR: parameter '{{"+key+"}}' not found in spec file template!",file=sys.stderr) 
    # optional parameters
    context["ranksPerNode"]  = ranksPerNode
    context["coresPerRank"]  = coresPerRank
    context["consumerTasks"] = consumerTasks

    if not createdFirstSpecFile:
        for key in context:
            if key not in keysInTemplate:
                print("WARNING: parameter '{{"+key+"}}' not found in spec file template!",file=sys.stderr)
        # verify template parameters are defined in options file
        for key in keysInTemplate:
            if key not in context:
                consistent = False
                print("ERROR: specification file template parameter '{{"+key+"}}' not defined in sweep options file!",file=sys.stderr)
        if not consistent:
            print("ERROR: subprogram aborted as specification file template and sweep options file are inconsistent.",file=sys.stderr)
            sys.exit()
        createdFirstSpecFile = True
    
    for key,value in context.items():
        renderedFile = renderedFile.replace("{{"+key+"}}", value)
    
    return renderedFile

def verifyLogFilterExists(justWarn=False):
    foundLogFilter = False
    for file in os.listdir(exahypeRoot + "/" + projectPath):
        foundLogFilter = foundLogFilter or file.endswith(".log-filter")

    messageType = "ERROR"
    if justWarn:
        messageTypeV = "WARNING"
    if not foundLogFilter:
        print(messageType+": no 'exahype.log-filter' file could be found in the project folder",file=sys.stderr)
        if not justWarn:
            sys.exit()

def verifyEnvironmentIsCorrect(justWarn=False):
    environmentIsCorrect = True
    
    messageType = "ERROR"
    if justWarn:
        messageType = "WARNING"
    
    cpus = jobs["num_cpus"]
    for environmentDict in dictProduct(environmentSpace):
        for key,value in environmentDict.items():
            os.environ[key]=value
        for config in ranksNodesCoreCounts:
            ranks = config.ranks
            nodes = config.nodes
            if (os.environ["DISTRIBUTEDMEM"].strip() not in ["MPI"]) and int(ranks)>1:
                print(messageType+": DISTRIBUTEDMEM environment variable set to "+environmentDict["DISTRIBUTEDMEM"]+" and ranks is set to "+ranks+" > 1",file=sys.stderr)
                environmentIsCorrect = False
            for coreCount in config.coreCounts:
                cores = coreCount.cores
                if (os.environ["SHAREDMEM"].strip() not in ["TBB","CPP14","OMP","TBBInvade"]) and int(cores)>1:
                    print(messageType+": SHAREDMEM environment variable set to "+environmentDict["SHAREDMEM"]+" and cores set to "+cores+" > 1",file=sys.stderr)
                    environmentIsCorrect = False
                        
    if not justWarn and not environmentIsCorrect:
        print("ERROR: subprogram failed as environment variables are not chosen setup correctly. Please adopt your options file according to the error messages.\n" + \
              "       Then rerun the subprogram.",file=sys.stderr)
        sys.exit()

def verifyAllRequiredParametersAreGiven(specFileTemplate):
    if "order" not in parameterSpace and "patchSize" not in parameterSpace:
        print("ERROR: Neither 'order' nor 'patchSize' found in section 'parameters' or 'parameters_grouped'.",file=sys.stderr)
        sys.exit()
    elif "dimension" not in parameterSpace:
        print("ERROR: 'dimension' not found in section 'parameters' or section 'parameters_grouped'.",file=sys.stderr)
        sys.exit()
    elif "optimisation" not in parameterSpace:
        print("ERROR: 'optimisation' not found in section 'parameters' or section 'parameters_grouped'.",file=sys.stderr)
        sys.exit()
    elif "architecture" not in parameterSpace:
        print("ERROR: 'architecture' not found in section 'parameters' or section 'parameters_grouped'.",file=sys.stderr)
        sys.exit()

    foundLimitingADERDG = "limiter-type" in specFileTemplate
    if foundLimitingADERDG:
        if "limiterType" not in parameterSpace:
            print("ERROR: 'limiterType' not found in section 'parameters' or section 'parameters_grouped'.",file=sys.stderr)
            sys.exit()
        elif "limiterOptimisation" not in parameterSpace:
            print("ERROR: 'limiterOptimisation' not found in section 'parameters' or section 'parameters_grouped'.",file=sys.stderr)
            sys.exit()

    return foundLimitingADERDG

def build(buildOnlyMissing=False, skipMakeClean=False):
    """
    Build the executables.
    
    Args:
    buildOnlyMissing(bool):
       Build only missing executables.
    
    skipMakeClean(bool):
        Do not run make clean, only clean application folder
    """
    print("currently loaded modules:")
    subprocess.call("module list",shell=True)
    print("")
    print("ExaHyPE build environment (unmodified):")
    exahypeEnv = ["COMPILER", "MODE", "SHAREDMEM", "DISTRIBUTEDMEM", "EXAHYPE_CC", "PROJECT_C_FLAGS", "PROJECT_L_FLAGS", "COMPILER_CFLAGS", "COMPILER_LFLAGS", "FCOMPILER_CFLAGS", "FCOMPILER_LFLAGS"]
    for variable in exahypeEnv:
        if variable in os.environ:
            print(variable+"="+os.environ[variable])
    print("")
    
    if not os.path.exists(buildFolderPath):
        print("create directory "+buildFolderPath)
        os.makedirs(buildFolderPath)
    
    verifyLogFilterExists(justWarn=True)        
    verifyEnvironmentIsCorrect(justWarn=True)
    
    architectures = parameterSpace["architecture"]
    optimisations = parameterSpace["optimisation"]
    dimensions    = parameterSpace["dimension"]
    
    orderKey = "order"
    if "order" in parameterSpace:
        orderKey = "order"
    elif "patchSize" in parameterSpace:
        orderKey = "patchSize"
    else:
       print("ERROR: Neither 'order' nor 'patchSize' found in section 'parameters' and 'parameters_ungrouped'.",file=sys.stderr)
       sys.exit()
    ordersOrPatchSizes = parameterSpace[orderKey]

    limiterTypes         = [None] 
    limiterOptimisations = [None]
    if foundLimitingADERDG:
        limiterTypes         = parameterSpace["limiterType"]
        limiterOptimisations = parameterSpace["limiterOptimisation"]

    buildParameterDict = list(dictProduct(parameterSpace))[0]
        
    firstIteration = True
    executables=0
    for environmentDict in dictProduct(environmentSpace):
        for key,value in environmentDict.items():
            os.environ[key]=value
        environmentDictHash = hashDictionary(environmentDict)
        
        for architecture in architectures:
            for dimension in dimensions:
                if not firstIteration and not skipMakeClean:
                    command = "make clean"
                    print(command)
                    process = subprocess.Popen([command], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                    (output, err) = process.communicate()
                    process.wait()
                for optimisation in optimisations:
                    for order in ordersOrPatchSizes:
                        for limiterType in limiterTypes:
                            for limiterOptimisation in limiterOptimisations:
                                oldExecutable = exahypeRoot + "/" + projectPath+"/ExaHyPE-"+projectName
                                suffix = architecture+"-d" + dimension + "-" + optimisation+ "-p" + order
                                if foundLimitingADERDG:
                                    suffix += "-"+limiterType+"-"+limiterOptimisation
                                
                                executable = buildFolderPath + "/ExaHyPE-"+projectName+"-"+environmentDictHash+"-"+suffix
                                
                                if not os.path.exists(executable) or not buildOnlyMissing:
                                    buildParameterDict["optimisation"] = optimisation
                                    buildParameterDict["architecture"] = architecture
                                    buildParameterDict["dimension"]    = dimension
                                    buildParameterDict[orderKey]       = order
                                        
                                    buildSpecFileBody = renderSpecFile(specFileTemplate,buildParameterDict,"1","1","1")
                                        
                                    buildSpecFilePath = outputPath+"/"+buildFolder+"/"+projectName+"-"+suffix+".exahype"
                                        
                                    with open(exahypeRoot + "/" + buildSpecFilePath, "w") as buildspecFile:
                                        buildspecFile.write(buildSpecFileBody)
                                        
                                    print("building executable with: \n" + \
                                          "- environment="+str(environmentDict) + "\n"\
                                          "- architecture='"+architecture + "'\n"\
                                          "- dimension="+dimension + "\n"\
                                          "- optimisation='"+optimisation + "'\n"\
                                          "- "+orderKey+"="+order + "\n"\
                                          "- limiterType='"+str(limiterType) + "'\n"\
                                          "- limiterOptimisation='"+str(limiterOptimisation)+"'")
                                    
                                    # clean application folder only
                                    command = "rm -r *.o cipofiles.mk cfiles.mk ffiles.mk kernels"
                                    print(command)
                                    process = subprocess.Popen([command], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                                    process.communicate()
                                    process.wait()
                                    
                                    # run toolkit
                                    toolkitCommand = "(cd "+exahypeRoot+" && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive "+buildSpecFilePath+")"
                                    print(toolkitCommand,end="",flush=True)
                                    process = subprocess.Popen([toolkitCommand], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                                    (output, toolkitErr) = process.communicate()
                                    process.wait()
                                    if "setup build environment ... ok" in str(output):
                                        print(" [OK]")
                                    else:
                                        print(" [FAILED]")
                                        print("toolkit output=\n"+output.decode('UTF-8'),file=sys.stderr)
                                        print("toolkit errors/warnings=\n"+toolkitErr.decode('UTF-8'),file=sys.stderr)
                                        sys.exit()
                                    
                                    if firstIteration and not skipMakeClean:
                                        command = "make clean"
                                        print(command)
                                        process = subprocess.Popen([command], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                                        (output, err) = process.communicate()
                                        process.wait()
                                    
                                    # call make
                                    make_threads=general["make_threads"]
                                    makeCommand="make -j"+make_threads
                                    print(makeCommand,end="",flush=True)
                                    process = subprocess.Popen([makeCommand], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                                    (output, makeErr) = process.communicate()
                                    process.wait()
                                    if "build of ExaHyPE successful" in str(output):
                                        print(" [OK]")
                                    else:
                                        print(" [FAILED]",file=sys.stderr)
                                        print("make errors/warnings=\n"+makeErr.decode('UTF-8'),file=sys.stderr)
                                        sys.exit()

                                    if not os.path.exists(oldExecutable):
                                        print("ERROR: could not find built executable '"+oldExecutable+"'. The parameter 'project' in your configuration file is probably wrong." ,file=sys.stderr)
                                        sys.exit()
                                    os.rename(oldExecutable,executable)
                                    print("created executable:"+executable)
                                     
                                    print("SUCCESS!")
                                    print("--------------------------------------------------------------------------------")
                                    print("toolkit errors/warnings=\n"+toolkitErr.decode('UTF-8'),file=sys.stderr)
                                    print("make errors/warnings=\n"+makeErr.decode('UTF-8'),file=sys.stderr)
                                    print("--------------------------------------------------------------------------------")
                                    executables+=1
                                    
                                    firstIteration = False
                                else:
                                    print("skipped building of '"+executable+"' as it already exists.")

    print("built executables: "+str(executables))


def renderJobScript(jobScriptTemplate,jobScriptBody,jobs,
                    jobName,jobScriptFilePath,outputFileName,errorFileName,
                    ranks,nodes,ranksPerNode,cores): # cores still necessary?
    """
    Render a job script.
    """
    renderedFile = jobScriptTemplate
    
    context = {}
    # mandatory
    context["ranks"]       = ranks
    context["nodes"]       = nodes
    context["output_file"] = outputFileName
    context["error_file"]  = errorFileName
    context["job_name"]    = jobName 
    
    context["body"] = jobScriptBody   
    
    consistent = True
    # verify all mandatory(!) sweep options are defined in template
    global createdFirstJobScript
    if not createdFirstJobScript:
        keysInTemplate = [m.group(2) for m in re.finditer("(\{\{((\w|-)+)\}\})",jobScriptTemplate)]
        for key in context:
            if key not in keysInTemplate:
                consistent = False
                print("ERROR: parameter '{{"+key+"}}' not found in job script template!",file=sys.stderr)
    
    # put optional sweep options in context
    context["mail"]         = jobs["mail"]
    context["ranksPerNode"]        = ranksPerNode
    context["time"]         = jobs["time"]
    context["class"]        = jobClass
    context["islands"]      = islands
    context["coresPerRank"] = str( int ( int(jobs["num_cpus"]) / int(ranksPerNode) ) )
    
    # now verify template parameters are defined in options file
    if not createdFirstJobScript:
        for key in keysInTemplate:
            if key not in context:
                consistent = False
                print("ERROR: job script template parameter '{{"+key+"}}' not defined in sweep options file!",file=sys.stderr)
        if not consistent:
            print("ERROR: subprogram aborted as job script template and sweep options file are inconsistent.",file=sys.stderr)
            sys.exit()
    
    context["body"] = jobScriptBody 
 
    for key,value in context.items():
        renderedFile = renderedFile.replace("{{"+key+"}}", value)
    
    return renderedFile

def verifyAllExecutablesExist(justWarn=False):
    """
    Verify that all executables exist.
    """    
    architectures = parameterSpace["architecture"]
    optimisations = parameterSpace["optimisation"]
    dimensions    = parameterSpace["dimension"]
    orders        = parameterSpace["order"]
    
    limiterTypes         = [None] 
    limiterOptimisations = [None]
    if foundLimitingADERDG:
        limiterTypes         = parameterSpace["limiterType"]
        limiterOptimisations = parameterSpace["limiterOptimisation"]
    
    messageType = "ERROR"
    if justWarn:
      messageType = "WARNING"
    
    if not justWarn and not os.path.exists(buildFolderPath):
        print("ERROR: build folder '"+buildFolderPath+"' doesn't exist! Please run subprogram 'build' beforehand.",file=sys.stderr)
        sys.exit()
    
    allExecutablesExist = True
    for environmentDict in dictProduct(environmentSpace):
        environmentDictHash = hashDictionary(environmentDict)
        for architecture in architectures:
            for dimension in dimensions:
                for optimisation in optimisations:
                    for order in orders:
                        for limiterType in limiterTypes:
                            for limiterOptimisation in limiterOptimisations:
                                suffix = architecture+"-d" + dimension + "-" + optimisation+ "-p" + order
                                if foundLimitingADERDG:
                                    suffix += "-"+limiterType+"-"+limiterOptimisation
                                
                                executable = buildFolderPath + "/ExaHyPE-"+projectName+"-"+environmentDictHash+"-"+suffix
                
                                if not os.path.exists(executable):
                                    allExecutablesExist = False
                                    print(messageType+ ": application for " + \
                                          "environment="+str(environmentDict) + \
                                          ", dimension="+dimension + \
                                          ", order="+order + \
                                          " does not exist! ('"+executable+"')",file=sys.stderr)
            
    if not justWarn and not allExecutablesExist:
        print("ERROR: subprogram failed as not all executables exist. Please adopt your options file according to the error messages.\n" + \
              "       Then rerun the 'build' subprogram.",file=sys.stderr)
        sys.exit()

def verifySweepAgreesWithHistoricalExperiments():
    """
    If there are any previous experiments ensure that the sweep 
    parameter spaces contain the same axes.
    """
    if os.path.exists(historyFolderPath):
        previousSweeps = [f for f in os.listdir(historyFolderPath) if f.endswith(".ini")]
        print(previousSweeps)
    
        for f in previousSweeps:
            otherOptions = sweep_options.parseOptionsFile(historyFolderPath+"/"+f)
        
            otherEnvironmentSpace = otherOptions.environmentSpace
            otherParameterSpace   = otherOptions.parameterSpace
        
            environmentSpaceIntersection = set(environmentSpace.keys()).intersection(otherEnvironmentSpace.keys())
            parameterSpaceIntersection   = set(parameterSpace.keys()).intersection(otherParameterSpace.keys())
            if len(set(environmentSpace.keys()))!=len(environmentSpaceIntersection):
                print("ERROR: subprogram failed as environment variables differ from previous experiments found in the output folder.",file=sys.stderr)
                print("environment variables found for CURRENT experiment: " + ", ".join(sorted(environmentSpace.keys())))
                print("environment variables used in PREVIOUS experiment:  " + ", ".join(sorted(otherEnvironmentSpace.keys())))
                sys.exit()
            if len(set(parameterSpace.keys()))!=len(parameterSpaceIntersection):
                print("ERROR: subprogram failed as parameters differ from previous experiments found in the output folder.",file=sys.stderr)
                print("parameters found for CURRENT experiment: "+ ", ".join(sorted(parameterSpace.keys())))
                print("parameters used in PREVIOUS experiment:  "+ ", ".join(sorted(otherParameterSpace.keys())))
                sys.exit()

def getSpecFilePath(parameterDictHash,ranksPerNode,cores,consumers):
   return scriptsFolderPath + "/" + projectName + "-" + \
          parameterDictHash + "-t"+ranksPerNode+"-c"+cores+"-b"+consumers+".exahype"

def getJobNameAndFilePaths(environmentDictHash,ungroupedParameterDictHash,configId,ranks,nodes,run):
   JobInfo = collections.namedtuple("ranksNodesCoreCounts", \
                                    "jobName jobScriptFilePath jobOutputFilePath jobErrorFilePath")
   jobName = projectName + "-" + environmentDictHash + "-" + ungroupedParameterDictHash + \
             "-C"+str(configId)+ "-n" + ranks + "-N" + nodes + "-r"+run
   jobScriptFilePath = scriptsFolderPath + "/" + jobName + ".job"
   jobOutputFilePath = resultsFolderPath + "/" + jobName + ".job_out"
   jobErrorFilePath  = resultsFolderPath + "/" + jobName + ".job_err"
   return JobInfo(jobName,jobScriptFilePath,jobOutputFilePath,jobErrorFilePath) 

def generateScripts():
    """
    Generate spec files and job scripts.
    """
    cpus       = jobs["num_cpus"]
        
    if not os.path.exists(scriptsFolderPath):
        print("create directory "+scriptsFolderPath)
        os.makedirs(scriptsFolderPath)
    
    # spec files
    specFiles=0
    for parameterDict in dictProduct(parameterSpace):
        parameterDictHash = hashDictionary(parameterDict)
        
        for config in ranksNodesCoreCounts:
            ranks = config.ranks
            nodes = config.nodes
            ranksPerNode = str( math.ceil(float(ranks)/float(nodes)) )
            for coreCounts in config.coreCounts:
                cores     = coreCounts.cores
                consumers = coreCounts.consumers
                
                specFilePath = getSpecFilePath(parameterDictHash,ranksPerNode,cores,consumers)
                
                specFileBody = renderSpecFile(specFileTemplate,parameterDict,ranksPerNode,cores,consumers)
                with open(specFilePath, "w") as specFile:
                    specFile.write(specFileBody)
                specFiles+=1
    
    print("generated specification files: "+str(specFiles))
    
    # check if required executables exist
    verifyEnvironmentIsCorrect(justWarn=True)
    verifyLogFilterExists(justWarn=True)        
    verifyAllExecutablesExist(justWarn=True)
    
    # generate job scrips
    jobScripts = 0
    for run in runNumbers:
        for configId,config in enumerate(ranksNodesCoreCounts):
            ranks = config.ranks
            nodes = config.nodes
            ranksPerNode = str( math.ceil(float(ranks)/float(nodes)) )
            for environmentDict in dictProduct(environmentSpace):
                environmentDictHash = hashDictionary(environmentDict)
                for ungroupedParameterDict in dictProduct(ungroupedParameterSpace):
                    ungroupedParameterDictHash = hashDictionary(ungroupedParameterDict)
                    
                    jobInfo = getJobNameAndFilePaths(environmentDictHash,ungroupedParameterDictHash,\
                                                     configId,ranks,nodes,run)

                    # aggregate the job script body
                    jobScriptBody = ""
                    for coreCount in config.coreCounts:
                        cores     = coreCount.cores
                        consumers = coreCount.consumers
                        for runGrouped in runNumbersGrouped:
                            myRun = run
                            if run=="+":
                                myRun = runGrouped
                            for groupedParameterDict in dictProduct(groupedParameterSpace):
                                parameterDict     = {}
                                parameterDict.update(ungroupedParameterDict)
                                parameterDict.update(groupedParameterDict)
                                parameterDict.pop(None) # ensure we do not hash a dummy None key
                                parameterDictHash = hashDictionary(parameterDict)
                                
                                architecture = parameterDict["architecture"]
                                optimisation = parameterDict["optimisation"]
                                dimension    = parameterDict["dimension"]
                                order        = parameterDict["order"]
                                
                                suffix = architecture+"-d" + dimension + "-" + optimisation+ "-p" + order
                                if foundLimitingADERDG:
                                    limiterType         = parameterDict["limiterType"]
                                    limiterOptimisation = parameterDict["limiterOptimisation"]
                                    suffix += "-"+limiterType+"-"+limiterOptimisation
                                
                                executable     = buildFolderPath + "/ExaHyPE-"+projectName+"-"+environmentDictHash+"-"+suffix
                                
                                specFilePath   = getSpecFilePath(parameterDictHash,ranksPerNode,cores,consumers)
                                                 
                                outputFileName = projectName + "-" + environmentDictHash + "-" + parameterDictHash + \
                                                 "-n" + ranks + "-N" + nodes + "-t"+ranksPerNode+"-c"+cores+"-b"+consumers+"-r"+myRun+".out"
                                outputFilePath = resultsFolderPath + "/" + outputFileName 

                                if "preamble" in jobs:
                                    renderedPreamble = jobs["preamble"].strip("\"")
                                    context = dict(parameterDict)
                                    context["ranks"]=ranks
                                    context["nodes"]=nodes
                                    context["ranksPerNode"]=ranksPerNode
                                    context["coresPerRank"]=cores
                                    context["backgroundTasks"]=consumers
                                    context["run"]=myRun
                                    for key,value in context.items():
                                        renderedPreamble = renderedPreamble.replace("{{"+key+"}}", value)
                                    jobScriptBody += "\n\n" + renderedPreamble + "\n\n"
                                
                                # pipe some information into output file
                                jobScriptBody += "echo \"Timestamp (YYYY/MM/dd:hh:mm:ss): `date +%Y/%m/%d:%H:%M:%S`\" > "+outputFilePath+"\n"
                                jobScriptBody += "echo \"\" >> "+outputFilePath+"\n" 
                                jobScriptBody += "module list >> "+outputFilePath+"\n"
                                jobScriptBody += "echo \"\" >> "+outputFilePath+"\n" 
                                jobScriptBody += "printenv >> "+outputFilePath+"\n"
                                jobScriptBody += "echo \"\" >> "+outputFilePath+"\n" 
                                jobScriptBody += "echo \""+jobInfo.jobScriptFilePath+":\" >> "+outputFilePath+"\n" 
                                jobScriptBody += "cat \""+jobInfo.jobScriptFilePath+"\" >> "+outputFilePath+"\n"  
                                jobScriptBody += "echo \"\" >> "+outputFilePath+"\n" 
                                jobScriptBody += "echo \""+specFilePath+":\" >> "+outputFilePath+"\n" 
                                jobScriptBody += "cat \""+specFilePath+"\" >> "+outputFilePath+"\n"
                                jobScriptBody += "echo \"\" >> "+outputFilePath+"\n" 
                                # pipe environment and parameter dicts into output file
                                jobScriptBody += "echo \"sweep/environment="+json.dumps(environmentDict).replace("\"","\\\"")+"\" >> "+outputFilePath+"\n"
                                jobScriptBody += "echo \"sweep/parameters="+json.dumps(parameterDict).replace("\"","\\\"")   +"\" >> "+outputFilePath+"\n"
                                # pipe the commands into the output file
                                runCommand = general["run_command"].replace("\"","")
                                runCommand = runCommand.replace("{{ranks}}",ranks);
                                runCommand = runCommand.replace("{{nodes}}",nodes);
                                runCommand = runCommand.replace("{{ranksPerNode}}",ranksPerNode);
                                runCommand = runCommand.replace("{{cores}}",cores);
                                if "./"==runCommand.strip():
                                    runCommand = runCommand.strip()
                                else:
                                    runCommand += " "
                                jobScriptBody += runCommand+executable+" "+specFilePath+" >> "+outputFilePath+"\n" # no whitespace after runCommand
                                
                                if "likwid" in general:
                                    groups = sweep_options.parseList(general["likwid"])
                                    likwid_pin_arg = int(cores) - 1
                                    if likwid_pin_arg == 0:
                                        likwid_pin_arg = str(likwid_pin_arg)
                                    else:
                                        likwid_pin_arg = "0-"+str(likwid_pin_arg) # e.g. 0-7 pins likwid to cores 0-7 
                                    for group in groups:
                                        likwidCommand = "{runCommand} likwid-perfctr -f -C {likwid_pin_arg} -g {group} {executable} {specFilePath} >> {outputFilePath}.likwid\n"
                                        jobScriptBody += likwidCommand.format(
                                            runCommand=runCommand,
                                            likwid_pin_arg=likwid_pin_arg,
                                            group=group,
                                            executable=executable,
                                            specFilePath=specFilePath,
                                            outputFilePath=outputFilePath
                                        )
                                jobScriptBody += "\n" 
                    
                    # write job file
                    renderedJobScript = renderJobScript(\
                                            jobScriptTemplate,jobScriptBody,jobs,
                                            jobInfo.jobName,jobInfo.jobScriptFilePath,jobInfo.jobOutputFilePath,jobInfo.jobErrorFilePath,
                                            ranks,nodes,ranksPerNode,cores)
                    with open(jobInfo.jobScriptFilePath, "w") as jobScriptFile:
                        jobScriptFile.write(renderedJobScript)
                    
                    jobScripts+=1

    print("generated job scripts: "+str(jobScripts))
                             
def verifyAllJobScriptsExist():
    """
    Verify that all job scripts exist.
    """
    cpus       = jobs["num_cpus"]
    
    if not os.path.exists(scriptsFolderPath):
        print("ERROR: job script folder '"+scriptsFolderPath+"' doesn't exist! Please run subprogram 'scripts' beforehand.",file=sys.stderr)
        sys.exit()
    
    allJobScriptsExist = True
    for run in runNumbers:
        for configId,config in enumerate(ranksNodesCoreCounts):
            ranks = config.ranks
            nodes = config.nodes
            ranksPerNode = str( math.ceil(float(ranks)/float(nodes)) )
            for environmentDict in dictProduct(environmentSpace):
                environmentDictHash = hashDictionary(environmentDict)
                for ungroupedParameterDict in dictProduct(ungroupedParameterSpace):
                    ungroupedParameterDictHash = hashDictionary(ungroupedParameterDict)
                    
                    jobInfo = getJobNameAndFilePaths(environmentDictHash,ungroupedParameterDictHash,\
                                                     configId,ranks,nodes,run)

                    if not os.path.exists(jobInfo.jobScriptFilePath):
                        allJobScriptsExist = False
                        print("ERROR: job script for " + \
                              "environment="+str(environmentDict)+ \
                              ", (ungrouped)parameters="+str(ungroupedParameterDict) + \
                              ", nodes="+nodes + \
                              ", ranksPerNode="+ranksPerNode + \
                              ", cores="+cores + \
                              ", run="+run + \
                              " does not exist! ('"+jobInfo.jobScriptFilePath+"')",file=sys.stderr)
    if not allJobScriptsExist:
        print("ERROR: subprogram failed! Please adopt your sweep options file according to the error messages.\n" + \
              "       Then rerun the 'scripts' subprogram.")
        sys.exit()

def verifyAllSpecFilesExist():
    """
    Verify that all ExaHyPE specification files exist.
    """
    cpus       = jobs["num_cpus"]
    
    if not os.path.exists(scriptsFolderPath):
        print("ERROR: job script folder '"+scriptsFolderPath+"' doesn't exist! Please run subprogram 'scripts' beforehand.",file=sys.stderr)
        sys.exit()
 
    allSpecFilesExist = True
    for parameterDict in dictProduct(parameterSpace):
        parameterDictHash = hashDictionary(parameterDict)
        
        for config in ranksNodesCoreCounts:
            ranks = config.ranks
            nodes = config.nodes
            ranksPerNode = str( math.ceil(float(ranks)/float(nodes)) )
            for coreCounts in config.coreCounts:
                cores     = coreCounts.cores
                consumers = coreCounts.consumers
                
                specFilePath = getSpecFilePath(parameterDictHash,ranksPerNode,cores,consumers)
              
                if not os.path.exists(specFilePath):
                     allSpecFilesExist = False
                     print("ERROR: specification file for \n" + \
                           "parameters="+str(parameterDict) + \
                           ", ranksPerNode="+ranksPerNode + \
                           ", cores="+cores + \
                           " does not exist! ('"+specFilePath+"')",file=sys.stderr)
    
    if not allSpecFilesExist:
        print("ERROR: subprogram failed! Please adopt your sweep options file according to the error messages.\n" + \
              "       Then rerun the 'scripts' subprogram.")
        sys.exit()

def hashSweep():
    chain = ""
    for config in ranksNodesCoreCounts:
        chain += config.ranks+";"+config.nodes+";"
        for coreCount in config.coreCounts:
            chain += coreCount.cores+";"+coreCount.consumers
    for value in runNumbers:
        chain += str(value)+";"
    for value in runNumbersGrouped:
        chain += str(value)+";"  
        
    for environmentDict in dictProduct(environmentSpace):
        chain += hashDictionary(environmentDict)
    for parameterDict in dictProduct(parameterSpace):
        chain += hashDictionary(parameterDict)
        
    return hashlib.md5(chain.encode()).hexdigest()[0:8]

def extractJobId(processOutput):
    jobId = "unknown"
    lines = processOutput.split("\n")
    for line in lines:
        # SLURM
        # hamilton: 'Submitted batch job 67586'
        # coolmuc:  'Submitted batch job 67586 on cluster mpp3'
        # LOAD-LEVELER:
        # supermuc: 'llsubmit: The job "srv24ib.840220" has been submitted.'
        if "Submitted batch job " in line:
            jobId = line.strip().split(" ")[3]
        if "llsubmit: The job " in line:
            jobId = line.strip().split("\"")[1]
    return jobId

def submitJobs():
    """
    Submit all jobs spanned by the options.
    """
    jobSubmissionTool    = general["job_submission"].replace("\"","")
    
    cpus = jobs["num_cpus"]
    
    # verify everything is fine
    verifyEnvironmentIsCorrect()
    verifyLogFilterExists(justWarn=True)        
    verifyAllExecutablesExist()
    verifyAllJobScriptsExist()
    verifyAllSpecFilesExist()
    
    if not os.path.exists(resultsFolderPath):
        print("create directory "+resultsFolderPath)
        os.makedirs(resultsFolderPath)
    
    # loop over job scrips
    jobIds = []
    for run in runNumbers:
        for configId,config in enumerate(ranksNodesCoreCounts):
            ranks = config.ranks
            nodes = config.nodes
            ranksPerNode = str( math.ceil(float(ranks)/float(nodes)) )
            for environmentDict in dictProduct(environmentSpace):
                environmentDictHash = hashDictionary(environmentDict)
                for ungroupedParameterDict in dictProduct(ungroupedParameterSpace):
                    ungroupedParameterDictHash = hashDictionary(ungroupedParameterDict)
                    
                    jobInfo = getJobNameAndFilePaths(environmentDictHash,ungroupedParameterDictHash,\
                                                     configId,ranks,nodes,run)
                            
                    command=jobSubmissionTool + " " + jobInfo.jobScriptFilePath
                    print(command)
                    process = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
                    (output, err) = process.communicate()
                    process.wait()
                    jobIds.append(extractJobId(output.decode("UTF_8")))
    
    if not os.path.exists(historyFolderPath):
        print("create directory "+historyFolderPath)
        os.makedirs(historyFolderPath)
    
    submittedJobsPath = historyFolderPath + "/" + hashSweep() + ".submitted"
    with open(submittedJobsPath, "w") as submittedJobsFile:
        submittedJobsFile.write(json.dumps(jobIds))
    
    print("submitted "+str(len(jobIds))+" jobs")
    print("job ids are memorised in: "+submittedJobsPath)
    command="cp "+optionsFile+" "+submittedJobsPath.replace(".submitted",".ini")
    print(command)
    process = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
    (output, err) = process.communicate()
    process.wait()

def cancelJobs():
    """
    Cancel submitted jobs.
    """
    jobCancellationTool = general["job_cancellation"].replace("\"","")
    
    submittedJobsPath   = historyFolderPath + "/" + hashSweep() + ".submitted"
    
    jobIds = None
    try:
        with open(submittedJobsPath, "r") as submittedJobsFile:
            jobIds = json.loads(submittedJobsFile.read())
    except IOError:
        print("ERROR: couldn't find any submitted jobs for current sweep ('"+submittedJobsPath+"').")
        sys.exit()

    for jobId in jobIds:
        command = jobCancellationTool + " " + jobId
        print(command)
        subprocess.call(command,shell=True)
    print("cancelled "+str(len(jobIds))+" jobs")    

    command = "rm "+submittedJobsPath
    print(command)
    subprocess.call(command,shell=True)

if __name__ == "__main__":
    import sys,os
    import subprocess
    import itertools
    import hashlib
    import json
    import re
    import math
    import collections   
 
    import sweep_analysis
    import sweep_options
    
    subprograms = ["build","buildMissing","buildLocally","scripts","submit","cancel","parseAdapters","parseTotalTimes","parseTimeStepTimes","parseMetrics","cleanBuild", "cleanScripts","cleanResults","cleanHistory","cleanAll"]
    
    if haveToPrintHelpMessage(sys.argv):
        info = \
"""sweep.py:

1) run:

./sweep.py myoptions.ini <subprogram>

available subprograms:

* build              - build all executables
* buildMissing       - build only missing executables
* buildLocally       - rebuild only the local application folder (no make clean)
* scripts            - submit the generated jobs
* cancel             - cancel the submitted jobs
* parseAdapters      - read the job output and parse adapter times
* parseTotalTimes    - read the adapter times table and:
                       Per configuration, calculate the total simulation time 
                       and minimise over all runs.
* parseTimeStepTimes - read the adapter times table and:
                       Per configuration, calculate the time spent per time step
                       and minimise over all runs.
* parseMetrics       - read the job output and parse likwid metrics
* cleanAll           - remove the whole sweep benchmark suite
* cleanBuild         - remove the build subfolder
* cleanScripts       - remove the scripts subfolder
* cleanResults       - remove the results subfolder
* cleanHistory       - clean the submission history

2) typical workflow:

./sweep.py myoptions.ini build
./sweep.py myoptions.ini scripts
./sweep.py myoptions.ini submit

(after jobs have finished)

./sweep.py myoptions.ini parseAdapters [--compress]
./sweep.py myoptions.ini parseTotalTimes
./sweep.py myoptions.ini parseTimeStepTimes

note: with --compress you remove table columns containing the
same value in every row.

(if likwid measurements have been performed)

./sweep.py myoptions.ini parseMetrics [--compress]


3) general options file (*.ini) structure

the options file must contain the following sections:

* [general]     - general settings such as the ExaHyPE-Engine path, the number of make threads, how to submit and cancel jobs, and so on
* [jobs]        - settings for the jobs to be run, e.g. the number of nodes, ranksPerNode, cores, and so on
* [environment] - environments for the executables to be built, e.g. the MODE, SHAREDMEM, DISTRIBUTEDMEM environment variables can be specified here

It must further contain at least one of the following sections:

* [parameters]         - The specification file parameters. Most of them are optional but some are required such as 'dimension', 'order', and so on. 
* [parameters_grouped] - More specificiation file parameters. But the corresponding runs are placed all
                         on the same job script.

"""
        print(info) # correctly indented
        sys.exit()
    
    optionsFile = parseArgument(sys.argv,1)
    subprogram  = parseArgument(sys.argv,2)
    
    compressTable = parseArgument(sys.argv,3)=="--compress"
    
    options = sweep_options.parseOptionsFile(optionsFile)
    
    general                 = options.general
    jobs                    = options.jobs
    environmentSpace        = options.environmentSpace
    parameterSpace          = options.parameterSpace
    ungroupedParameterSpace = options.ungroupedParameterSpace
    groupedParameterSpace   = options.groupedParameterSpace    
 
    exahypeRoot      = options.exahypeRoot
    outputPath       = options.outputPath
    projectPath      = options.projectPath
    projectName      = options.projectName
   
    buildFolder       = options.buildFolder
    buildFolderPath   = options.buildFolderPath
    scriptsFolderPath = options.scriptsFolderPath
    resultsFolderPath = options.resultsFolderPath
    historyFolderPath = options.historyFolderPath
    
    jobClass             = options.jobClass
    islands              = options.islands
    ranksNodesCoreCounts = options.ranksNodesCoreCounts
    runNumbers           = options.runNumbers
    runNumbersGrouped    = options.runNumbersGrouped   

    createdFirstSpecFile  = False
    createdFirstJobScript = False
 
    verifySweepAgreesWithHistoricalExperiments()
    
    specFileTemplatePath = exahypeRoot+"/"+general["spec_template"]
    specFileTemplate     = None
    try:
        with open(specFileTemplatePath, "r") as specFileTemplateFile:
            specFileTemplate=specFileTemplateFile.read()
    except IOError:
        print("ERROR: couldn\'t open specification file template file: "+specFileTemplatePath,file=sys.stderr)
        sys.exit()
        
    jobScriptTemplatePath = exahypeRoot+"/"+general["job_template"]    
    jobScriptTemplate = None
    try:
        with open(jobScriptTemplatePath, "r") as jobScriptTemplateFile:
            jobScriptTemplate=jobScriptTemplateFile.read()
    except IOError:
        print("ERROR: couldn\'t open job script template file: "+jobScriptTemplatePath,file=sys.stderr)
        sys.exit()
        
    foundLimitingADERDG = verifyAllRequiredParametersAreGiven(specFileTemplate)
 
    # select subprogram
    if subprogram == "cleanAll":
        clean()
    elif subprogram == "cleanBuild":
        clean("build")
    elif subprogram == "cleanScripts":
        clean("scripts")
    elif subprogram == "cleanResults":
        clean("results")
    elif subprogram == "cleanHistory":
        clean("history")
    elif subprogram == "build":
        build(buildOnlyMissing=False, skipMakeClean=False)
    elif subprogram == "buildMissing":
        build(buildOnlyMissing=True, skipMakeClean=False)
    elif subprogram == "buildLocally":
        build(buildOnlyMissing=False, skipMakeClean=True)
    elif subprogram == "scripts":
        generateScripts()
    elif subprogram == "submit":
        submitJobs()
    elif subprogram == "cancel":
        cancelJobs()
    elif subprogram == "parseAdapters":
        sweep_analysis.parseAdapterTimes(resultsFolderPath,projectName,compressTable)
    elif subprogram == "parseTotalTimes":
        sweep_analysis.parseSummedTimes(resultsFolderPath,projectName)
    elif subprogram == "parseTimeStepTimes":
        sweep_analysis.parseSummedTimes(resultsFolderPath,projectName,timePerTimeStep=True)
    elif subprogram == "parseMetrics":
        sweep_analysis.parseMetrics(resultsFolderPath,projectName,compressTable)
