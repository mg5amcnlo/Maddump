import logging
import os

import maddump_run_interface as maddump_run_interface
#import madgraph.core.diagram_generation as diagram_generation

import madgraph.interface.master_interface as master_interface
import madgraph.interface.madgraph_interface as madgraph_interface
import madgraph.various.misc as misc

import madgraph.interface.common_run_interface as common_run

import models.write_param_card as param_writer
import models.check_param_card as check_param_card
import models.model_reader as model_reader
import madgraph.core.helas_objects as helas_objects

from madgraph.iolibs.files import cp, ln, mv
from .. import fit2D_card as fit2D

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


    _define_options = ['darkmatter','decay_channel']
    _importevts_options = ['decay','interaction']
    
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
        self._decay_channel = []
        self._param_card = None
        self._out_dir = ''
        self._evts_inputfile_todecay = [] 
        self._evts_inputfile_tointeract = [] 
        self._decay_list = []
        self._particle_todisplaced=[]
        self._multiparticles = {}
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
                
            # elif args[0] == 'decay_channel':
            #     if len(args)==2:
            #         self._decay_channel = [self._curr_model.get_particle(args[1])]
            #         if not self._decay_channel[0]:
            #             raise DMError, '%s is not a valid particle for the model.' % args[1] 
            #         # self.update_model_with_EFT()
                # else:
                #     raise DMError, 'The user must supply a decay channel (daughter particle)!'
                
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

    
    def do_import_events(self, line):
        args = self.split_arg(line)
        # if '.lhe' or '.hepmc' not in args:
        #     raise DMError, 'Events file name must contain .lhe or .hepmc extension!'
        if 'decay' in args:
            args.pop(0)
            if not os.path.isfile(args[0]):
                raise DMError, 'Invalid path: the file does not exist!'                
            if not any( [ext in args[0] for ext in ['hepmc','lhe']]):
                raise DMError, 'Invalid events file name: it must contain the hepmc or lhe extension!'
            if not self._evts_inputfile_todecay:
                self._evts_inputfile_todecay.append(args[0])
            else:
                while(1):
                    answ = raw_input ('Events file to decay already loaded. Do you want to overwrite it? [y/n] ')
                    if answ == 'y':
                        self._evts_inputfile_todecay[0] = args[0]
                        break 
                    elif answ == 'n':
                        break
                    else:
                        print("Please, answer with 'y' or 'n'.")                        
        if 'interaction' in args:
            args.pop(0)
            if not os.path.isfile(args[0]):
                raise DMError, 'Invalid path: the file does not exist!'                
            if not any( [ext in args[0] for ext in ['hepmc','lhe']]):
                raise DMError, 'Invalid events file name: it must contain the hepmc or lhe extension!'
            if not self._evts_inputfile_tointeract:
                self._evts_inputfile_tointeract.append(args[0])
            else:
                while(1):
                    answ = raw_input ('Events file to decay already loaded. Do you want to overwrite it? [y/n] ')
                    if answ == 'y':
                        self._evts_inputfile_tointeract[0] = args[0]
                        break 
                    elif answ == 'n':
                        break
                    else:
                        print("Please, answer with 'y' or 'n'.")

                        
    def complete_import_events(self, text, line, begidx, endidx, formatting=True):
        "Complete the import_events command"
        
        out = {}
        args = self.split_arg(line[0:begidx])

        if len(args) == 1:
            out['maddump options'] = self.list_completion(text, self._importevts_options, line)
            return self.deal_multiple_categories(out, formatting)        
        if len(args) == 2:
            base_dir = '.'
        else:
            base_dir = args[2]
        
        return self.path_completion(text, base_dir)
        
        # Directory continuation
        if os.path.sep in args[-1] + text:
            return self.path_completion(text,
                                    pjoin(*[a for a in args if \
                                                      a.endswith(os.path.sep)]))

        
    def do_add(self, line):
        """ """
        
        args = self.split_arg(line)
        args.pop(0)
        ''' command: generate production process 
            (process is mandatory)
        ''' 

        if len(args) and args[0] == 'process':
            args.pop(0)

        if len(args) and args[0] == 'production':
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
        #  c) elastic nucleon scattering (not yet implemented)
        #  channels. The matrix elements are the same as in madevent
        #  while the output format must be supplied by the maddump
        #  plugin.  
        elif len(args) and args[0] == 'interaction':
            args.pop(0)            
            if not '@' in line:
                raise DMError, 'The user must supply a valid interaction channel!'
            tag = re.search('(?<=@)\w+', line)
            #print(tag.group(0))
            excl = re.search('(?<=/)\w+', line)
            if excl:
                excluded = '/'+excl.group(0)
            else:
                excluded = ''
            try:
                dm_candidate = self._dm_candidate[-1]['name']
            except:
                raise DMError, 'Please define a valid dark matter candidate!'
            int_proc = {'DIS' : 'process ' +  dm_candidate + ' p > ' + dm_candidate + ' p ' + excluded + ' @DIS',
                        'electron' : 'process ' + dm_candidate + ' e- > ' + dm_candidate + ' e- ' + excluded + ' @electron',
                        'ENS' : 'work in progress'}
            self.do_add(int_proc[tag.group(0)])

        elif len(args) and args[0] == 'displaced_decay':
            if len(args)==2:
                self._particle_todisplaced = [self._curr_model.get_particle(args[1])]
                if not self._particle_todisplaced[0]:
                    raise DMError, '%s is not a valid particle for the model.' % args[1] 
            else:
                raise DMError, 'The user must supply a valid particle to displaced!'
            return
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
            options = ['production', 'interaction', 'displaced_decay']
            out['maddump options'] = self.list_completion(text, options , line)
        
        return self.deal_multiple_categories(out, formatting)

    
    def complete_add(self, text, line, begidx, endidx, formatting=True):
        """Complete particle information"""

        args = self.split_arg(line[:begidx])
        
        return super(MadDump_interface, self).complete_add(text, line, begidx, endidx,formatting=False)


    def do_decay(self, line):
        self._decay_list.append(line)

    # def do_displaced_decay(self, line):
    #     args = self.split_arg(line)
    #     self._particle_todisplaced = [self._curr_model.get_particle(args[0])]
    #     if not self._particle_todisplaced[0]:
    #         raise DMError, '%s is not a valid particle for the model.' % args[1] 
        
        
    def complete_decay(self, text, line, begidx, endidx, formatting=True):
        """Complete particle information"""

        args = self.split_arg(line[:begidx])
        out = super(MadDump_interface, self).complete_generate(text, line, begidx, endidx,formatting=False)
        if not isinstance(out, dict):
            out = {"standard options": out}
        return self.deal_multiple_categories(out, formatting)
        
    
    def do_output(self, line):
        """ """

        if not self._curr_amps:
            if not self._particle_todisplaced:
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
                not_created = False                

        else:
            self._out_dir = args[0]
            try:
                os.makedirs(self._out_dir)
            except:
                raise DMError, 'Output command error: the directory %s already exists!' % self._out_dir 


        #create a directory for common cards: param_card and fit2D_card
        #stored in the list cards
        cards_dir = pjoin(self._out_dir, 'Cards')
        os.makedirs(cards_dir)
        cards = ['param_card.dat','param_card_default.dat',
                 'fit2D_card.dat','fit2D_card_default.dat']
        #Decay
        if self._evts_inputfile_todecay:
            evts_todecay_dir = pjoin(self._out_dir, 'Events_to_decay')
            os.makedirs(evts_todecay_dir)
            ln(self._evts_inputfile_todecay[0], evts_todecay_dir,)
            self.write_madspin_card(cards_dir,evts_todecay_dir)
        elif self._decay_list:
            decay_dir = pjoin(self._out_dir, 'Events_to_decay')
            os.makedirs(decay_dir)
            self._evts_inputfile_todecay.append('unweighted_events.lhe.gz')
            self.write_madspin_card(cards_dir,decay_dir)
        
        #interaction_only
        if self._evts_inputfile_tointeract:
            evts_tointeract_dir = pjoin(self._out_dir, 'Events_to_interact')
            os.makedirs(evts_tointeract_dir)
            ln(self._evts_inputfile_tointeract[0], evts_tointeract_dir,)
        #Displaced_decay
        if self._particle_todisplaced:
            evts_todisplaced_dir = pjoin(self._out_dir, 'Displaced_Decay')
            os.makedirs(evts_todisplaced_dir)
            proc_characteristics = open(pjoin(evts_todisplaced_dir,'proc_characteristics'),'w')
            proc_characteristics.write('grouped_matrix = True' + '\n')
            proc_characteristics.write('pdg_mother = ' + str(self._particle_todisplaced[0]['pdg_code'])+'\n')
#            proc_characteristics.write('pdg_daughter = '+ str(self._decay_channel[0]['pdg_code'])+'\n')
            proc_characteristics.write('BSM_model = ' + str(self._curr_model.get('modelpath'))+'\n')
            proc_characteristics.close()
            self.create_fit2D_card(cards_dir)
            if not self._curr_amps:
                self.create_param_card(cards_dir)
                return
            
        for proc in self._curr_proc_defs:
            
            if proc['id'] < 1000:
                self._curr_matrix_elements = helas_objects.HelasMultiProcess()
                amps = [ amp for amp in self._curr_amps if amp.get('process')['id'] == proc['id'] ]
                line = pjoin(self._out_dir,'production')
                with misc.TMP_variable(self, ['_curr_proc_defs','_curr_amps','_done_export'], [proc,amps,None]):                    
                    super(MadDump_interface, self).do_output(line)
                # check whether the param_card is already in the common Cards dir
                # otherwise it is copied from the current process 
                if cards[0] not in os.listdir(cards_dir):
                    for i in range(2):
                        pcard = pjoin(self._out_dir,'production', 'Cards')
                        cp(pjoin(pcard, cards[i]), cards_dir)
                        os.remove(pjoin(pcard, cards[i]))
                        ln(pjoin(cards_dir, cards[i]), pcard, ) 


                # put the run_card in the common Cards dir and create a symbolic link
                # in the production/Cards dir
                pcards_dir = pjoin(self._out_dir, 'production', 'Cards') 
                cp(pjoin(pcards_dir, 'run_card.dat'), pjoin(cards_dir))
                cp(pjoin(pcards_dir, 'run_card_default.dat'), pjoin(cards_dir))
                
                try:
                    #for card in cards: #remove this pointless loop
                    os.remove(pjoin(pcards_dir, 'run_card.dat'))
                    os.remove(pjoin(pcards_dir, 'run_card_default.dat'))
                except OSError:
                    pass
                for name in ['run_card.dat', 'run_card_default']:
                    ln(pjoin(cards_dir, name), pcards_dir)
                    
            else:
                self._curr_matrix_elements = helas_objects.HelasMultiProcess()
                processes_list = self.process_tag.keys()
                for process in processes_list:
                    if self.process_tag[process] == proc['id']:
                        channel = process
                        break
                amps = [ amp for amp in self._curr_amps if amp.get('process')['id'] == proc['id'] ]
                
                line = 'maddump %s' % pjoin(self._out_dir, 'interaction_' + channel) 
                with misc.TMP_variable(self, ['_curr_proc_defs','_curr_amps','_done_export'], [proc,amps,False]):                    
                    super(MadDump_interface, self).do_output(line)

                current_dir = pjoin(self._out_dir, 'interaction_' + channel)
                if cards[0] not in os.listdir(cards_dir):
                    for i in range(2):
                        cp(pjoin(current_dir, 'Cards', cards[i]), pjoin(cards_dir, cards[i]))
                if cards[2] not in os.listdir(cards_dir):
                    for i in range(2,4):
                        cp(pjoin(current_dir, 'Cards', cards[i]), pjoin(cards_dir, cards[i]))

                try:
                    for card in cards:
                        os.remove(pjoin(current_dir,'Cards', card))
                except OSError:
                    pass
                # os.symlink(os.path.join('production',paramcard_path), os.path.join(current_dir, paramcard_path))
                for card in cards:
                    ln(pjoin(cards_dir, card),  pjoin(current_dir,'Cards'))
                
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
        self.define_child_cmd_interface(self._MDUMP,  interface=False)
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


    def write_madspin_card(self,cpath,epath):
        """ """
        
        evts_file = pjoin(epath,os.path.basename(self._evts_inputfile_todecay[0]))
        
        out = open(pjoin(cpath,'madspin_card.dat'),'w')
        s = 'set spinmode none\nset cross_section {0:1.0}\nset new_wgt BR\n'
        out.write(s)
        if 'lhe' in evts_file:
            out.write('set input_format lhe_no_banner\n')
        elif 'hepmc' in evts_file:
            out.write('set input_format hepmc\n')

        s = 'import ' + evts_file + '\n' 
        s += 'import model ' + self._curr_model.get("modelpath") + ' ' + pjoin(cpath,'param_card.dat') + '\n'   
        out.write(s)

        if self._multiparticles:
            for tag in self._multiparticles:
                s='define '+str(tag)+' = '
                for pid in self._multiparticles[tag]:
                    s+=str(pid)+' '
                s+='\n'
                out.write(s)
                
        for decay in self._decay_list:
            out.write('decay '+ decay + '\n')
        out.write('launch\n')
                        
        
    #===========================================================================
    #  create the fit2D_card 
    #===========================================================================
    def create_fit2D_card(self,path):
        """  """
        fit2D_card = fit2D.Fit2DCard()
        
        fit2D_card.write(pjoin(path, 'fit2D_card_default.dat'))
        fit2D_card.write(pjoin(path, 'fit2D_card.dat'))

#TODO for decay_displaced

    @staticmethod
    def create_param_card_static(model, output_path, rule_card_path=False,
                                 mssm_convert=True):
        """ create the param_card.dat for a givent model --static method-- """
        #1. Check if a default param_card is present:
        done = False
        if hasattr(model, 'restrict_card') and isinstance(model.restrict_card, str):
            restrict_name = os.path.basename(model.restrict_card)[9:-4]
            model_path = model.get('modelpath')
            if os.path.exists(pjoin(model_path,'paramcard_%s.dat' % restrict_name)):
                done = True
                files.cp(pjoin(model_path,'paramcard_%s.dat' % restrict_name),
                         output_path)
        if not done:
            param_writer.ParamCardWriter(model, output_path)
         
        if rule_card_path:   
            if hasattr(model, 'rule_card'):
                model.rule_card.write_file(rule_card_path)
        
        if mssm_convert:
            model_name = model.get('name')
            # IF MSSM convert the card to SLAH1
            if model_name == 'mssm' or model_name.startswith('mssm-'):
                import models.check_param_card as translator    
                # Check the format of the param_card for Pythia and make it correct
                if rule_card_path:
                    translator.make_valid_param_card(output_path, rule_card_path)
                translator.convert_to_slha1(output_path)        
    
    def create_param_card(self,path):
        """ create the param_card.dat """

        rule_card = pjoin(path, 'param_card_rule.dat')
        model = self._curr_model
        if not hasattr(model, 'rule_card'):
            rule_card=False
        self.create_param_card_static(model=model, 
                                      output_path=pjoin(path, 'param_card.dat'), 
                                      rule_card_path=rule_card, 
                                      mssm_convert=True)

    #===========================================================================
    #  create the param_card 
    #===========================================================================
    # def create_param_card(self,path):
    #     """  """
    #     fit2D_card = fit2D.Fit2DCard()
        
    #     fit2D_card.write(pjoin(path, 'fit2D_card_default.dat'))
    #     fit2D_card.write(pjoin(path, 'fit2D_card.dat'))


    
