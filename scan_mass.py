#! /usr/bin/env python

if '__main__' == __name__:
    import sys
    sys.path.append('./')
    sys.path.append('madgraph/various')

import os,shutil
import logging
#import math
import lhe_parser
import madgraph.various.banner as banner_mod
import numpy as np
import subprocess
#from meshfitter2D import *

pjoin = os.path.join
pwd_path=os.getcwd()
#logging.basicConfig()


def ChangeParameter(paramcard, flag, value):
    #-----------------------------------------------------------------------#
    #                                                                        #
    #  Changes the parameter in the param card.                                #
    #  flag is a string name of the parameter. Value is the numerical value #
    #  that the parameter should be set to.                                        #
    #                                                                              #
    #-----------------------------------------------------------------------#

    try:
        foundflag = False
        #Open the param card to read in all the lines
        finput = open(paramcard, 'r')
        lines = finput.readlines()
        finput.close()
        #Open it again to rewrite it. Necessary in order to
        #not leave junk at the end.
        finput = open(paramcard, 'w')
        for line in lines:#fileinput.FileInput(self._paramcard,inplace=1):
            line_list = line.split()
            if flag in line_list:
                #Make sure that the flag is not a part of another variable name
                foundflag = True
                index_of_value = line_list.index(flag)-2
                newline=''
                #Find where the value is in the string
                #Put the line back together
                for j in range(0, index_of_value):
                    if j==0:
                        newline=newline+line_list[j]
                    else:
                        newline=newline+' '+line_list[j]
                newline=newline+' '+str(value)
                for j in range(index_of_value+1, len(line_list)):
                    newline=newline+' '+line_list[j]
                newline=newline+'\n'
                finput.write(newline)
                #print "New Line: "+newline
            else:
            	#print "Old Line: "+line
                finput.write(line)

        finput.close()
        if not foundflag:
            print "ERROR: Parameter "+flag+" not found!"

    except OSError, error:
        logging.debug(error)
        print "ERROR: Could not open the Param Card file!"
        sys.exit(1)

prod_dir = 'prod_pp400GeV'
detec_dir = 'SHiP_detection_DIS'
paramcard_path = '/Cards/param_card.dat'

#input parameters production/detection
npot = 10**20         #number of proton on target
A_Mo = 96             #Molibdeno mass number 
Na = 6.022*10**23     #Avogadro's number 
d = 1804.75           #distance target-detector in cm
md = 10**7            #mass of the detector in g
th= 0.05              #aperture angle of particle hitting the detector in rad
xsec_pp = 40*10**9    #cross section pp in pb    
S = np.pi*(d*np.sin(th))**2*10**36 #detector surface area hit by particles in pb^-1

fac = npot*A_Mo/A_Mo**(.7)/xsec_pp*Na*md/S
#range info in GeV
min_mass=3.#1      
max_mass=4 #10     
bin_size=0.5
nstep = (max_mass-min_mass)/bin_size

mass = [min_mass+bin_size*i for i in range(int(nstep)+1)] 

out = open('mass_scan.dat','w')

for m in mass:
    ChangeParameter(prod_dir+paramcard_path,'MZB',m)
    ChangeParameter(detec_dir+paramcard_path,'MZB',m)

    run_dir = 'ZBmass_'+str(m)+'GeV'
    os.chdir(prod_dir)
    os.system('./bin/generate_events < input')
    rundirname=subprocess.check_output("ls Events -trl | tail -n 1 | awk '{print $9}'",shell=True)
    
    try:
        shutil.rmtree('Events/'+run_dir)
    except OSError:
        pass

    os.system('mv Events/'+rundirname.rstrip()+' Events/'+run_dir)


    os.chdir('../'+detec_dir+'/Cards')

    try:
        os.remove('unweighted_events.lhe.gz')
    except OSError:
        pass

    os.system('ln -s ../../'+prod_dir+'/Events/'+run_dir+'/unweighted_events.lhe.gz')
    os.system('./lhe-meshfitter.py')
    os.chdir('../')
    os.system('./bin/generate_events < input')    
    rundirname=subprocess.check_output("ls Events -trl | tail -n 1 | awk '{print $9}'",shell=True)

    try:
        shutil.rmtree('Events/'+run_dir)
    except OSError:
        pass

    os.system('mv Events/'+rundirname.rstrip()+' Events/'+run_dir)
    events_lhefile = 'Events/'+run_dir+'/unweighted_events.lhe.gz'
    lhe = lhe_parser.EventFile(events_lhefile)
    xsec = lhe.cross    
    nevt=fac*xsec

    os.chdir('../')
    out.write(str(m) + '\t' + str(xsec) + '\t' + str(nevt) + '\n' )
    
out.close()
