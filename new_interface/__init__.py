## import the required files
# example: import maddm_interface as maddm_interface # local file
#          import madgraph.various.cluster as cluster #MG5 distribution file
# Three types of functionality are allowed in a plugin
#   1. new output mode
#   2. new cluster support
#   3. new interface

import maddump_interface as maddump_interface
import maddump as mdd

# 1. Define new output mode
#    example: new_output = {'myformat': MYCLASS}
#    madgraph will then allow the command "output myformat PATH"
#    MYCLASS should inherated of the class madgraph.iolibs.export_v4.VirtualExporter 
new_output = {'maddump': mdd.MadDump}

# 2. Define new way to handle the cluster.
#    example new_cluster = {'mycluster': MYCLUSTERCLASS}
#    allow "set cluster_type mycluster" in madgraph
#    MYCLUSTERCLASS should inherated from madgraph.various.cluster.Cluster
new_cluster = {}


 # 3. Define a new interface (allow to add/modify MG5 command)
#    This can be activated via ./bin/mg5_aMC --mode=PLUGINNAME
## Put None if no dedicated command are required
new_interface = maddump_interface.MadDump_interface
 
########################## CONTROL VARIABLE ####################################
__author__ = 'Luca Buonocore'
__email__ = 'lbuono@na.infn.it'
__version__ = (1,0,0)
minimal_mg5amcnlo_version = (2,6,0) 
maximal_mg5amcnlo_version = (1000,1000,1000)
latest_validated_version = (2,6,0)
