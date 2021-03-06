import argparse
from argparse import RawTextHelpFormatter

import pylab
import runtimeParser


Colors           = ['#ff0000','#00ff00','#0000ff','#ffbb00','#ff00ff','#00ffff']
Markerfacecolors = ['None','None','None','None','None','None','k','k','k']
Markers          = ['o','^','s','o','^','s']


xDataMax = 0
yDataMin = 65636
yDataMax = 0
symbolCounter = 0
plotCounter   = 0

def addData(table,normalisation,plotLabels,experimentSetCounter,label):
  global xDataMax
  global yDataMin
  global yDataMax
  global symbolCounter
  global plotCounter
  
  symbolCounter = plotCounter / int(args.changecolour)
  plotCounter   = plotCounter + 1
  experimentSetCounter = experimentSetCounter - plotCounter / int(args.changecolour)
  
  xdata         = runtimeParser.readColumnFromTable(table,0)

  for adap in args.adapter:
    sub_adapters = adap.split('+')
    
    totalTime    = []
    count        = 0
    i = 0
    for sub_adap in sub_adapters:
      sub_adap_time  = runtimeParser.readColumnFromTable(table, runtimeParser.getAdapterRuntimeColumnFromTable(table,sub_adap) )
      count          = max(count ,runtimeParser.readColumnFromTable(table, runtimeParser.getAdapterCountColumnFromTable(table,sub_adap) ));
      
      if (len(totalTime)==0):
        totalTime = sub_adap_time
      else:
        for i in range(0,len(sub_adap_time)):
          totalTime[i] += sub_adap_time[i]
    
    ydata = []
    for i in range(0,len(totalTime)):
      if count[i]==0:
        print "ERROR " + str(i) + "th entry in " + table + "'s " + str(runtimeParser.getAdapterCountColumnFromTable(table,adap)) + "th column (adapter " + adap + ") equals 0"
        ydata = []
        break 
      ydata.append( totalTime[i]*normalisation/count[i])
    if len(ydata)==0:
      print "WARNING: file " + table + " seems to be empty for adapter " + adap
    else:
      for i in range(0,len(ydata)):
        EvaluateMinMaxAt = 1
        EvaluateMinMaxAt = int(args.singlecore)
        if ydata[i]>yDataMax and xdata[i]==EvaluateMinMaxAt:
          yDataMax = ydata[i]
        if ydata[i]<yDataMin and xdata[i]==EvaluateMinMaxAt:
          yDataMin = ydata[i]
      if xdata[-1]>xDataMax:
        xDataMax = xdata[-1]
      if (plotLabels):
        pylab.plot(xdata,ydata,markersize=experimentSetCounter+4,label=adap[:20],color=Colors[symbolCounter],marker=Markers[symbolCounter],markerfacecolor=Markerfacecolors[symbolCounter],markevery=1,lw=1.2) 
      else:
        pylab.plot(xdata,ydata,markersize=experimentSetCounter+4,color=Colors[symbolCounter],marker=Markers[symbolCounter],markerfacecolor=Markerfacecolors[symbolCounter],markevery=1,lw=1.2) 
    symbolCounter = symbolCounter + 1

    if adap==args.adapter[-1] and len(xdata)>0:
     pylab.text(xdata[-1],ydata[-1],label )



def plotLinearSpeedupCurves():
  global xDataMax
  global yDataMin
  global yDataMax

  speedupCurveYdata = [yDataMin/(x-int(args.singlecore)+1) for x in range(int(args.singlecore),int(xDataMax))]
  pylab.plot(range(int(args.singlecore),int(xDataMax)),speedupCurveYdata,markersize=4,markevery=1,lw=1.2,linestyle='dashed',color='grey',label="linear") 

  speedupCurveYdata = [yDataMax/(x-int(args.singlecore)+1) for x in range(int(args.singlecore),int(xDataMax))]
  pylab.plot(range(int(args.singlecore),int(xDataMax)),speedupCurveYdata,markersize=4,markevery=1,lw=1.2,linestyle='dashed',color='grey') 

  if int(args.singlecore)>1:
    speedupCurveYdata = [yDataMin*int(args.singlecore)/x for x in range(int(args.singlecore),int(xDataMax))]
    pylab.plot(range(int(args.singlecore),int(xDataMax)),speedupCurveYdata,markersize=4,markevery=1,lw=1.2,linestyle='dotted',color='black',label="linear w.o. multicore") 

    speedupCurveYdata = [yDataMax*int(args.singlecore)/x for x in range(int(args.singlecore),int(xDataMax))]
    pylab.plot(range(int(args.singlecore),int(xDataMax)),speedupCurveYdata,markersize=4,markevery=1,lw=1.2,linestyle='dotted',color='black') 




  
def initGlobalPlotterSettings():
  pylab.xlabel( args.xaxislabel )
  pylab.grid(True)
  pylab.legend(loc='best',fontsize='%d' % int(args.fontsize))  
  #fig, ax = pylab.subplots()
  #ax.yaxis.set_major_formatter(pylab.FormatStrFormatter('%.2e'))

  pylab.ticklabel_format(axis='y',style='sci', scilimits=(-2,2))
  #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))


def switchToLogScales():
  global xDataMax
  
  pylab.loglog( basex=2, basey=2 )
  XTicks  = [ int(args.singlecore) ]
  XLabels = [ "serial" ]
  #XTicks  = [ int(args.singlecore), int(args.singlecore) ]
  #XLabels = [ "no MPI", "serial" ]
  #XTicks  = [ int(args.singlecore) ]
  #XLabels = [ "no MPI" ]
  for i in range(1,int(xDataMax)+2):
    if i>1 and ((i & (i - 1)) == 0) and i>int(args.singlecore):
      XTicks.append( i )
      XLabels.append( str(i) ) 
  pylab.xticks(XTicks,XLabels)
  pylab.xlim(int(args.singlecore)-1,xDataMax)

  #and i>int(args.singlecore)

  #for i in range(1,len(xdata)):
  #  if xdata[i]==args.singlecore:
  #    tSerial = yData[i]


########################################################################
# START OF THE PROGRAM
########################################################################
# For the meaning of the command line arguments, see the documentation
# of the function 'plot_multithreading_adapter_scaling' above.
help = '''
Creates a speedup plot based on Peano output files with specific file naming pattern.
If multiple adapters are specified, then the cumulative user times are used to compute the speedup.
\n\n
Sample usages:\n
python ../plotRuntime.py -experimentdescription '' -table x4-1.results.table -output experiment1 -xaxislabel "MPI Ranks" -dimension 2 \n
python ../plotRuntime.py -experimentdescription '' -table x4-1.results.table -adapter Total ADERDGTimeStep PredictorRerun -output experiment1 -xaxislabel "MPI Ranks" -dimension 2 \n
python ../plotRuntime.py -experimentdescription 'depth 6' 'depth 7' 'depth 8' -table x4-1.results.table x4-2.results.table x4-3.results.table -adapter ADERDGTimeStep PredictorRerun -output experimentx4 -xaxislabel "MPI Ranks" -dimension 2 \n
python ../plotRuntime.py -experimentdescription 'depth 4' 'depth 5' 'depth 6' 'depth 7' 'depth 8' 'depth 9' -table x16-0.results.table x16-1.results.table x16-2.results.table x16-3.results.table x16-4.results.table x16-5.results.table -adapter ADERDGTimeStep PredictorRerun -output experimentx16 -xaxislabel "MPI Ranks" -dimension 2 \n
'''

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument('-table',nargs='+',required=True,help="Tables that are to be read.")
parser.add_argument('-output',required=True,help="Output file (extensions pdf and png are added).")
parser.add_argument('-adapter',default=['Total'],nargs='+',help="Name of the adapters separated by a \'+\'. (Use \'Total\' for the  cumulative time of all specified adapters. Does not make sense with \'per_iteration\' switched on.)")
parser.add_argument('-dimension',required=True,help="Dimension of problem. Either 2 or 3.")
parser.add_argument('-xaxislabel',required=True,help="Label of x axis.")
parser.add_argument('-fontsize',default=10,required=False,help="Font size of the legend and tick labels. Axis labels are computed by ceiling the font size times a factor 1.2.")
parser.add_argument('-experimentdescription',nargs='+',required=True,help="Per table entry, one experiment desription is required")
parser.add_argument('-singlecore',required=True,help="Which core count is the sequential run time")
parser.add_argument('-meshsize',required=False,nargs='+',help="Specify per table how big the mesh size is")
parser.add_argument('-changecolour',required=False,default=65536,help="Tell plotter after how many plots to change colouring scheme (colour furthermore is changed for each adapter)")
args   = parser.parse_args();

dim = int(args.dimension)



#
# Raw runtime
#
outputFile = args.output + "-raw-runtime"
pylab.clf()
pylab.xlabel('cores')

yDataMin = 65636
yDataMax = 0
plotCounter   = 0

experimentSetCounter =  0
for (table,label) in zip(args.table,args.experimentdescription):
  print "read " + table
  maxLevel = runtimeParser.readColumnFromTable(table,1)   
  addData(table,1.0,True,experimentSetCounter,label)
  experimentSetCounter = experimentSetCounter + 1

plotLinearSpeedupCurves()

initGlobalPlotterSettings()

pylab.ylabel( "time per time step [t]=s" )

pylab.savefig( outputFile + ".png" )
pylab.savefig( outputFile + ".pdf" )
switchToLogScales()
pylab.savefig( outputFile + "-log.png" )
pylab.savefig( outputFile + "-log.pdf" )
print "written " + outputFile




#
# Runtime scaled by number of cells. Does not yield reasonable output if grid statistics are switched off
#
outputFile = args.output + "-scaled-by-cells"
pylab.clf()
pylab.xlabel('cores')

yDataMin = 65636
yDataMax = 0
plotCounter   = 0

experimentSetCounter =  0
for (table,label) in zip(args.table,args.experimentdescription):
  print "read " + table
  maxLevel      = runtimeParser.readColumnFromTable(table,1)
  # normalisation = 1.0
  if len(maxLevel)>0 and maxLevel[-1]==0:
    print "WARNING: max level seems not to have been plotted (perhaps code had been translated without -DTrackGridStatistics). Assume increase by one per table"
  normalisation = 1.0/(3**(experimentSetCounter*dim))
  #print "normalise with " + str(normalisation)
  if len(maxLevel)>0  and maxLevel[-1]>0:
    normalisation = 1.0/(3**(maxLevel[-1]*dim))
    print "reset normalisation to " + str(normalisation)

  addData(table,normalisation,table==args.table[-1],experimentSetCounter,label)

  experimentSetCounter = experimentSetCounter + 1

plotLinearSpeedupCurves()

initGlobalPlotterSettings()

pylab.ylabel( "time per cell per time step [t]=s" )

pylab.savefig( outputFile + ".png" )
pylab.savefig( outputFile + ".pdf" )
switchToLogScales()
pylab.savefig( outputFile + "-log.png" )
pylab.savefig( outputFile + "-log.pdf" )
print "written " + outputFile
    

    
  
#
# Runtime scaled by regular grid of max depth
#
outputFile = args.output + "-scaled-by-regular-grid"
pylab.clf()
pylab.xlabel('cores')

yDataMin = 65636
yDataMax = 0
plotCounter   = 0

if args.meshsize is not None:

 experimentSetCounter =  0
 for (table,label,h) in zip(args.table,args.experimentdescription,args.meshsize):
  print "read " + table
  depth = 1
  while 3**(-depth*dim)>float(h):
    depth = depth+1
  normalisation = 3**(-depth*dim)

  addData(table,normalisation,table==args.table[-1],experimentSetCounter,label)
  experimentSetCounter = experimentSetCounter + 1

 plotLinearSpeedupCurves()
 
 initGlobalPlotterSettings()

 pylab.ylabel( "time per cell per time step [t]=s" )

 pylab.savefig( outputFile + ".png" )
 pylab.savefig( outputFile + ".pdf" )
 switchToLogScales()
 pylab.savefig( outputFile + "-log.png" )
 pylab.savefig( outputFile + "-log.pdf" )
 print "written " + outputFile
    
    
