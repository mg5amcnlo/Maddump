import logging
import os

import maddump_run_interface as maddump_run_interface
#import madgraph.core.diagram_generation as diagram_generation

import madgraph.interface.master_interface as master_interface
import madgraph.interface.madgraph_interface as madgraph_interface
import madgraph.various.misc as misc

import madgraph.interface.common_run_interface as common_run

import models.check_param_card as check_param_card
import models.model_reader as model_reader
import madgraph.core.helas_objects as helas_objects

from madgraph.iolibs.files import cp, ln, mv

import re
pjoin = os.path.join

logger = logging.getLogger('madgraph.plugin.maddump')

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    GRAY = '\033[90m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

class DMError(Exception): pass

# Root path
MDumpDIR = os.path.dirname(os.path.realpath( __file__ ))

    
class MadDump_interface(master_interface.MasterCmd):

    # intro_banner= 'MadDump plugin'
    intro_banner=\
  "            ====================================================\n"+\
  "            |                  "+bcolors.OKBLUE+"  Maddump v1.0 (beta)            "+bcolors.ENDC+"|\n"\
  "            ====================================================\n"+\
  "            ====================================================\n"+\
  "            |           "+bcolors.OKBLUE+" A MadGraph5_aMC@NLO plugin.            "+bcolors.ENDC+"|\n"\
  "            ====================================================\n"+\
  "%s" 


    _define_options = ['darkmatter']
    
    # process number to distinguish the different type of matrix element
    process_tag = {'prod': 100,
                   'decay': 200,
                   'DIS': 1100,
                   'electron': 1200,
                   'ENS': 1300}
    
    # eff_operators_SI = {1:'SIEFFS', 2:'SIEFFF', 3:'SIEFFV'}
    # eff_operators_SD = {1:False, 2:'SDEFFF', 3:'SDEFFV'} 
    
    def __init__(self, *args, **opts):
        
        super(MadDump_interface, self).__init__(*args, **opts)
        self._dm_candidate = []
        self._param_card = None
        self._out_dir = ''
        
################################################################################        
# DEFINE COMMAND
################################################################################        
    def do_define(self, line, **opts):
        """pass"""
        args = self.split_arg(line)
        if len(args) and args[0] in MadDump_interface._define_options:
            if args[0] == 'darkmatter':
                if len(args)==2:
                    self._dm_candidate = [self._curr_model.get_particle(args[1])]
                    if not self._dm_candidate[0]:
                        raise DMError, '%s is not a valid particle for the model.' % args[1] 
                    # self.update_model_with_EFT()
                else:
                    raise DMError, 'The user must supply a dark matter candidate!'
                
        else:
            return super(MadDump_interface, self).do_define(line,**opts)
        
    def complete_define(self, text, line, begidx, endidx, formatting=True):
        """Complete particle information + handle maddump options"""

        out = {}
        args = self.split_arg(line[0:begidx])
        if len(args) == 1:
            out['maddump options'] = self.list_completion(text, self._define_options, line)
        out['particles name']  = self.model_completion(text, line[6:], line)
         
        return self.deal_multiple_categories(out, formatting)

    def help_define(self):
        """ """
        logger.info("**************** MADDUMP OPTION ***************************")
        logger.info("syntax: define darkmatter [OPTIONS]", '$MG:color:BLUE')
        logger.info(" -- define the current (set) of darkmatter.")
        logger.info(" OPTIONS:", '$MG:color:BLACK')
        logger.info("    - You can specify the darkmatter by specifying their name/pdg code")
        logger.info("     o example: define darkmatter n1 n2",'$MG:color:GREEN')
        logger.info("    - You can remove some particle from the search by prefixing the particle name by '/'.")
        logger.info("     o example: define darkmatter /n1", '$MG:color:GREEN')                    
        logger.info("")        
        logger.info("**************** MG5AMC OPTION ***************************")
        super(MadDump_interface, self).help_define()

    # To modify!
    
    # def do_import(self, line,*args, **opts):
    #     """normal import but perform a cleaning for MadDM  variables"""
        
    #     lineargs = self.split_arg(line)
    #     self.check_import(lineargs)
    #     if lineargs and lineargs[0].startswith('model'):
    #         #reset DM information
    #         self._param_card = None
    #         self._dm_candidate = []
    #         self._coannihilation = []
            
    #     return super(MadDumpM_interface, self).do_import(line, *args, **opts)
      

    def do_add(self, line):
        """ """
        
        args = self.split_arg(line)
        args.pop(0)
        ''' command: generate production process 
            (process is mandatory)
        ''' 

        if len(args) and args[0] == 'process':
            args.pop(0)

        if len(args) and args[0] in['production']:
            args.pop(0)
            line = 'process'
            for x in args:
                line += ' ' + x
            line += ' ' + '@prod' 
            self.do_add(line)

        #  command: generate interaction @channel 
        #  (channel is mandatory)
        #  it is assumed that the user has defined the DM candidate
        #  handle the interaction of the DM candidate in the
        #  a) DIS 
        #  b) elastic electron scattering 
        #  c) elastic nucleon scattering (WIP)
        #  channels. The matrix elements are the same as in madevent
        #  while the output format must be supplied by the maddump
        #  plugin.  
        elif len(args) and args[0] in ['interaction']:
            args.pop(0)            
            if not '@' in line:
                raise DMError, 'The user must supply a valid interaction channel!'
            tag = re.search('(?<=@)\w+', line)
            print(tag.group(0))
            dm_candidate = self._dm_candidate[-1]['name']
            int_proc = {'DIS' : 'process ' +  dm_candidate + ' p > ' + dm_candidate + ' p @DIS',
                        'electron' : 'process ' + dm_candidate + ' e- > ' + dm_candidate + ' e- @electron',
                        'ENS' : 'work in progress'}
            self.do_add(int_proc[tag.group(0)])
            
        else:
            if '@' in line:
                line = re.sub(r'''(?<=@)(%s\b)''' % '\\b|'.join(self.process_tag), 
                              lambda x: `self.process_tag[x.group(0)]`, line)
            return super(MadDump_interface, self).do_add(line)

        
    def complete_generate(self, text, line, begidx, endidx, formatting=True):
        """Complete particle information + handle maddump options"""

        args = self.split_arg(line[:begidx])
        out = super(MadDump_interface, self).complete_generate(text, line, begidx, endidx,formatting=False)
        if not isinstance(out, dict):
            out = {"standard options": out}
        
        if len(args) == 1:
            options = ['production', 'interaction']
            out['maddump options'] = self.list_completion(text, options , line)
        
        return self.deal_multiple_categories(out, formatting)

    
    def complete_add(self, text, line, begidx, endidx, formatting=True):
        """Complete particle information"""

        args = self.split_arg(line[:begidx])
        
        return super(MadDump_interface, self).complete_add(text, line, begidx, endidx,formatting=False)

    
    def do_output(self, line):
        """ """
        
        if not self._curr_amps:
            raise DMError, 'No valid process to output!'
        
        args = self.split_arg(line)

        if len(args)>1:
            raise DMError, 'Output command error: too many options!'
        
        for arg in args:
            if arg in self._export_formats:
                raise DMError, 'Standard madgraph options are not allowed inside the maddump interface!'

        if not args:
            not_created = True
            i = 1
            while not_created:
                self._out_dir = 'maddump_'+str(i)
                try:
                    os.makedirs(self._out_dir)
                except:
                    i+=1
                    continue
                os.chdir(self._out_dir)
                not_created = False                

        else:
            self._out_dir = args[0]
            try:
                os.makedirs(self._out_dir)
            except:
                raise DMError, 'Output command error: the directory %s already exists!' % self._out_dir 
            os.chdir(self._out_dir)

        #create a directory for common cards: param_card and fit2D_card
        #stored in the list cards
        os.makedirs('Cards')
        cards_dir = 'Cards/'
        cards = ['param_card.dat','param_card_default.dat',
                 'fit2D_card.dat','fit2D_card_default.dat']
        
        for proc in self._curr_proc_defs:
            
            if proc['id'] < 1000:
                self._curr_matrix_elements = helas_objects.HelasMultiProcess()
                amps = [ amp for amp in self._curr_amps if amp.get('process')['id'] == proc['id'] ]
                line = 'production'
                with misc.TMP_variable(self, ['_curr_proc_defs','_curr_amps','_done_export'], [proc,amps,None]):                    
                    super(MadDump_interface, self).do_output(line)
                # check whether the param_card is already in the common Cards dir
                # otherwise it is copied from the current process 
                if cards[0] not in os.listdir(cards_dir):
                    for i in range(2):
                        cp('production/Cards/'+cards[i],
                           cards_dir+cards[i])
                        os.remove('production/Cards/'+cards[i])
                        os.symlink(pjoin('..','..',cards_dir,cards[i]), pjoin('production','Cards',cards[i]))

                # put the run_card in the common Cards dir and create a symbolic link
                # in the production/Cards dir
                cp('production/Cards/'+'run_card.dat', cards_dir+'run_card.dat')
                cp('production/Cards/'+'run_card_default.dat', cards_dir+'run_card_default.dat')
                try:
                    for card in cards:
                        os.remove('production/Cards/'+'run_card.dat')
                        os.remove('production/Cards/'+'run_card_default.dat')
                except OSError:
                    pass
                
                os.symlink(pjoin('..','..',cards_dir,'run_card.dat'), pjoin('production','Cards','run_card.dat'))
                os.symlink(pjoin('..','..',cards_dir,'run_card_default.dat'), pjoin('production','Cards/','run_card_default.dat'))
            else:
                self._curr_matrix_elements = helas_objects.HelasMultiProcess()
                processes_list = self.process_tag.keys()
                for process in processes_list:
                    if self.process_tag[process] == proc['id']:
                        channel = process
                        break
                amps = [ amp for amp in self._curr_amps if amp.get('process')['id'] == proc['id'] ]
                
                line = 'maddump interaction_' + channel 
                with misc.TMP_variable(self, ['_curr_proc_defs','_curr_amps','_done_export'], [proc,amps,False]):                    
                    super(MadDump_interface, self).do_output(line)

                current_dir = 'interaction_' + channel + '/'
                if cards[0] not in os.listdir(cards_dir):
                    for i in range(2):
                        cp(current_dir+'Cards/'+cards[i],
                           cards_dir+cards[i])
                if cards[2] not in os.listdir(cards_dir):
                    for i in range(2,4):
                        cp(current_dir+'Cards/'+cards[i],
                           cards_dir+cards[i])
                try:
                    for card in cards:
                        os.remove(current_dir+'Cards/'+card)
                except OSError:
                    pass
                # os.symlink(os.path.join('production',paramcard_path), os.path.join(current_dir, paramcard_path))
                for card in cards:
                    os.symlink(pjoin('..','..',cards_dir,card), pjoin(current_dir,'Cards',card))
                    #ln(cards_dir + card, current_dir + card, log=False)
                # if 'proc_characteristics' not in os.listdir(cards_dir):
                #     ln(current_dir+'SubProcesses/proc_characteristics', cards_dir + 'proc_characteristics', log=False)
                    
        os.chdir('../')
                
    # def find_output_type(self, path):
    #     if os.path.exists(pjoin(path,'matrix_elements','proc_characteristics')):
    #         return 'maddm'
    #     else:
    #         return super(MadDM_interface, self).find_output_type(path)

    
    def do_launch(self, line):
        ''' '''
        args = self.split_arg(line)
        (options, args) = madgraph_interface._launch_parser.parse_args(args)
        #with misc.TMP_variable(self, '_export_formats', self._export_formats + ['maddm']):
        #self.check_launch(args, options)
        options = options.__dict__

        if not args:
            self._MDUMP = maddump_run_interface.MADDUMPRunCmd(dir_path = self._out_dir, options = self.options)
        else:
            self._MDUMP = maddump_run_interface.MADDUMPRunCmd(dir_path = line, options = self.options)
            
        self._MDUMP.exec_cmd('launch')
        # if args[0] != 'maddm':
        #     return super(MadDM_interface, self).do_launch(line)
        # else:
        #     self._MDM = maddm_run_interface.MADDMRunCmd(dir_path=args[1], options=self.options)
        #     self._done_export = (args[1], 'plugin')
            
        #     if options['interactive']:
        #         return self.define_child_cmd_interface(self._MDM)
        #     else:
        #         self.define_child_cmd_interface(self._MDM,  interface=False)
        #         try:
        #             self._MDM.exec_cmd('launch ' + line.replace(args[1], ''))
        #         except:
        #             self._MDM.exec_cmd('quit')
        #             raise
        #         else:
        #             self._MDM.exec_cmd('quit')
        #         self._done_export = (args[1], 'plugin')
        #         return

            
    def define_multiparticles(self, label, list_of_particles):
        """define a new multiparticle from a particle list (add both particle and anti-particle)"""
        
        pdg_list = []
        for p in list_of_particles:
            if p.get('name') == p.get('antiname'):
                pdg_list.append(p.get('pdg_code'))
            else:
                pdg_list.append(p.get('pdg_code'))
                pdg_list.append(-1*p.get('pdg_code'))

        self.optimize_order(pdg_list)
        self._multiparticles[label] = pdg_list
