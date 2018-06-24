import os

import madgraph.various.misc as misc
import madgraph.interface.common_run_interface as common_run
import madgraph.interface.extended_cmd as cmd
import madgraph.interface.launch_ext_program as launch_ext
import madgraph.iolibs.files as files
import models.check_param_card as param_card_mod
import subprocess
import logging
#import MadSpin.decay as decay
import madgraph.iolibs.save_load_object as save_load_object
import MadSpin.interface_madspin as interface_madspin
import lhe_to_pythia_hadron_std as lheToPythia
import displaced_decay as displ_decay

pjoin = os.path.join
logger = logging.getLogger('madgraph.plugin.maddump')

from .. import fit2D_card as fit2D
from .. import  meshfitter2D as meshfitter

# LAUNCH PROGRAM
_launch_usage = "launch [DIRPATH] [options]\n" + \
         "-- execute the madevent/standalone/standalone_cpp/pythia8/NLO output present in DIRPATH\n" + \
         "   By default DIRPATH is the latest created directory \n" + \
         "   (for pythia8, it should be the Pythia 8 main directory) \n" + \
         "   Example: launch PROC_sm_1 --name=run2 \n" + \
         "   Example: launch ../pythia8 \n"
_launch_parser = misc.OptionParser(usage=_launch_usage)
_launch_parser.add_option("-f", "--force", default=False, action='store_true',
                                help="Use the card present in the directory in order to launch the different program")
_launch_parser.add_option("-n", "--name", default='', type='str',
                                help="Provide a name to the run (for madevent run)")
_launch_parser.add_option("-c", "--cluster", default=False, action='store_true',
                                help="submit the job on the cluster")
_launch_parser.add_option("-m", "--multicore", default=False, action='store_true',
                                help="submit the job on multicore core")

_launch_parser.add_option("-i", "--interactive", default=False, action='store_true',
                                help="Use Interactive Console [if available]")
_launch_parser.add_option("-s", "--laststep", default='',
                                help="last program run in MadEvent run. [auto|parton|pythia|pgs|delphes]")
_launch_parser.add_option("-R", "--reweight", default=False, action='store_true',
                            help="Run the reweight module (reweighting by different model parameter")
_launch_parser.add_option("-M", "--madspin", default=False, action='store_true',
                            help="Run the madspin package")


class paramCardIterator(param_card_mod.ParamCardIterator):

    def write_summary(self, path, order=None, lastline=False, nbcol=20):
        """ """
        predir,name = os.path.split(path)
        
        for t in ['DIS', 'electron', 'ENS']:
            if os.path.exists(pjoin(predir, 'interaction_%s' % t)):
                path = pjoin(predir, 'interaction_%s' % t, name)
                super(paramCardIterator, self).write_summary(path, order)
        


#===============================================================================
# CommonRunCmd
#===============================================================================
class MADDUMPRunCmd(cmd.CmdShell):
    
    intro_banner= 'MadDump plugin'
    
    def __init__(self, dir_path, options, *args, **opts):
        """common"""
  
        self.in_scan_mode = False # please modify this value only with "with" statement
        cmd.Cmd.__init__(self, *args, **opts)
        # Define current MadEvent directory
        if dir_path is None and not MG5MODE:
            dir_path = root_path

        self.dir_path = dir_path
        self.me_dir = dir_path
        self.param_card_iterator = [] #a placeholder containing a generator of paramcard for scanning
        #self.store_results = {}
        self.results = {}
        #A shared param card object to the interface.
        #All parts of the code should read this card.
        #Set at the beginning of launch()
        #self.param_card = None
        
        self.auto_width = set() # keep track of width set on Auto

        self.options = options
        
    # def check_param_card(self, path, run=True):
    #     """
    #     1) Check that no scan parameters are present
    #     2) Check that all the width are defined in the param_card.
    #     - If a scan parameter is defined. create the iterator and recall this function 
    #       on the first element.
    #     - If some widths are set on 'Auto', call the computation tools.
    #     - keep track of the width of auto if some present"""
        
    #     pattern_scan = re.compile(r'''^(decay)?[\s\d]*scan''', re.I+re.M)  
    #     pattern_width = re.compile(r'''decay\s+(\+?\-?\d+)\s+auto(@NLO|)''',re.I)
    #     text = open(path).read()
               
    #     if pattern_scan.search(text):
    #         if not isinstance(self, cmd.CmdShell):
    #             # we are in web mode => forbid scan due to security risk
    #             raise Exception, "Scans are not allowed in web mode"
    #         # at least one scan parameter found. create an iterator to go trough the cards
    #         main_card = param_card_mod.ParamCardIterator(text)
    #         self.param_card_iterator = main_card
    #         first_card = main_card.next(autostart=True)
    #         first_card.write(path)
    #         return self.check_param_card(path, run)
        
    #     pdg_info = pattern_width.findall(text)
    #     if pdg_info:
    #         pdg = [pdg for pdg,nlo in pdg_info]
    #         self.auto_width = set(pdg)
    #         if run:
    #             if not self.in_scan_mode:
    #                 logger.info('Computing the width set on auto in the param_card.dat')
    #             has_nlo = any(nlo.lower()=="@nlo" for _,nlo in pdg_info)
    #             if not has_nlo:
    #                 self.do_compute_widths('%s --path=%s' % (' '.join(pdg), path))
    #                 #self.run_mg5(['compute_widths %s --path=%s' % (' '.join(pdg), path)])
    #             else:
    #                 self.run_mg5(['compute_widths %s --path=%s --nlo' % (' '.join(pdg), path)])
    #         else:
    #             logger.info('''Some widths are on Auto in the card.
    # Those will be computed as soon as you have finish editing the cards.
    # If you want to force the computation right now and re-edit
    # the cards afterwards, you can type \"compute_wdiths\".''')

    def get_proc_characteristics(self,path):
        proc_file = open(path,'r')
        proc_characteristics = {} 
        for line in proc_file:
            if '=' not in line:
                continue
            else:
                args = line.split()
                proc_characteristics[args[0]] = args[2]
        return proc_characteristics
    
    def get_model(self):
        
        model_name =  self.proc_characteristics['BSM_model']
        print(model_name)
        self.mother.exec_cmd('import model %s' % model_name, child=False)
        return self.mother._curr_model

        
    def do_compute_widths(self, line):
        """normal fct but ensure that self.maddm_card is up-to-date"""
        
        # try:
        #     self.mother_interface.maddm_card = self.maddm
        # except Exception,error:
        #     logger.error("Invalid command: %s " % error)
        #     return
        return super(MadDumpSelector, self).do_compute_widths(line)

    def help_compute_widths(self, line):
        
        return self.run_mg5([' help compute_widths ' + line])    

    ############################################################################
    def do_launch(self, line):
        """run the code"""

        args = line.split()
        if '-f' in args or '--force' in args:
            force = True
        else:
            force = False
        # check argument validity and normalise argument
        (options, args) = _launch_parser.parse_args(args)
        #self.check_launch(args, options)
        options = options.__dict__
                    
        interaction_dir_all = ['interaction_DIS','interaction_electron','interaction_ENS']
        listdir = subprocess.check_output("ls %s"%self.dir_path,shell=True).split()

        self.displaced = False
        if 'production' in listdir:
            if 'Events_to_decay' in listdir:
                self.DMmode='production_decay'
                self.displaced = True
            else:
                self.DMmode='production_interaction'
        elif 'Events_to_decay' in listdir:
            if 'Displaced_Decay' in listdir:
                self.DMmode='decay_displaced'
                self.displaced = True
            else:
                self.DMmode = 'decay_interaction'
        elif 'Events_to_interact' in listdir:
            self.DMmode = 'interaction_only'
        else:
            print('Corrupted directory: no production/Events_to_decay/Events_to_interact directory ')
            exit(-1)
            
        self.interaction_dir = []
        for dir in listdir:
            if dir in interaction_dir_all:
                self.interaction_dir.append(dir)

        if options['name']:
            self.run_name = options['name']
        else:
            i = 1
            while os.path.exists(pjoin(self.dir_path, 'production', 'Events', 'run_%02d' %i)) or\
                  os.path.exists(pjoin(self.dir_path, 'Events_to_decay', 'Events', 'run_%02d' %i)) or \
                  os.path.exists(pjoin(self.dir_path, 'Events_to_interact', 'Events', 'run_%02d' %i)) or \
                  any([os.path.exists(pjoin(self.dir_path, dir, 'Events', 'run_%02d' %i)) for dir in self.interaction_dir]) or \
                  os.path.exists(pjoin(self.dir_path, 'Displaced_Decay', 'Events', 'run_%02d' %i)):
                i += 1
            self.run_name = 'run_%02d' % i

        if self.interaction_dir:
            self.proc_characteristics = self.get_proc_characteristics(pjoin(self.dir_path,self.interaction_dir[0],'SubProcesses','proc_characteristics'))

        if self.displaced:
            self.proc_characteristics = self.get_proc_characteristics(pjoin(self.dir_path,'Displaced_Decay','proc_characteristics'))
            
        self.ask_run_configuration(mode=[], force=force)
        self.run_launch()
        
    @common_run.scanparamcardhandling(result_path=lambda obj:  pjoin(obj.me_dir, 'scan_%s.txt'))
    def run_launch(self):
        
        if self.DMmode in ['production_interaction','production_decay']:
            misc.call(['./bin/generate_events', self.run_name, '-f'],
                      cwd=pjoin(self.dir_path, 'production'))
            evts_path = pjoin(self.dir_path, 'production', 'Events', self.run_name, 'unweighted_events.lhe.gz') 
            results = self.load_results_db(pjoin(self.dir_path,'production'),self.run_name)
            label = 'xsec_prod (pb)'
            self.results[label] = results['cross']
            if self.DMmode == 'production_decay':
                decay_dir = pjoin(self.dir_path,'Events_to_decay')
                files.ln(evts_path, decay_dir)
           
        if self.DMmode in ['decay_interaction','production_decay','decay_displaced']:
            decay_dir = pjoin(self.dir_path,'Events_to_decay')
            listdir = subprocess.check_output("ls %s"%decay_dir,shell=True).split()
            for file in listdir:
                if any( [ext in file for ext in ['hepmc','lhe']]):
                    evts_file = file.replace('.hepmc','')
                    evts_file = evts_file.replace('.lhe','')
                    evts_file = evts_file.replace('.gz','')

            try:
                os.makedirs(pjoin(decay_dir,'Events'))
            except:
                pass
            
            madspin_cmd = interface_madspin.MadSpinInterface()
            # pass current options to the interface
            madspin_cmd.mg5cmd.options.update(self.options)
            madspin_cmd.cluster = None #self.cluster
        
            madspin_cmd.update_status = lambda *x,**opt: self.update_status(*x, level='madspin',**opt)
            path = pjoin(self.me_dir, 'Cards', 'madspin_card.dat')

            madspin_cmd.import_command_file(path)

            # create a new run_name directory for this output
            i = 1
            while os.path.exists(pjoin(decay_dir,'Events', '%s_decayed_%i' % (self.run_name,i))):
                i+=1
            new_run = '%s_decayed_%i' % (self.run_name,i)
            evt_dir = pjoin(decay_dir, 'Events',new_run)

            os.makedirs(evt_dir)
            listdir = subprocess.check_output("ls %s"%decay_dir,shell=True).split()

            decayed_file = pjoin(decay_dir,os.path.basename(evts_file)+'_decayed.lhe.gz')
            try:
                files.mv(decayed_file,pjoin(evt_dir,'unweighted_events.lhe.gz'))
            except:
                logger.error('MadSpin fails to create any decayed file.')
                         
            # current_file = evts_file.replace('.lhe', '_decayed.lhe')
            # new_file = pjoin(evt_dir, new_run, )
            # if not os.path.exists(current_file):
            #     if os.path.exists(current_file+'.gz'):
            #         current_file += '.gz'
            #         new_file += '.gz'
            #     elif current_file.endswith('.gz') and os.path.exists(current_file[:-3]):
            #         current_file = current_file[:-3]
            #         new_file = new_file[:-3]
            #     else:

            #         return
            evts_path = pjoin(evt_dir,'unweighted_events.lhe.gz')
            
        if self.DMmode == 'interaction_only':
            interactionevts_dir = pjoin(self.dir_path,'Events_to_interact')
            listdir = subprocess.check_output("ls %s"%interactionevts_dir,shell=True).split()
            for file in listdir:
                if any( [ext in file for ext in ['hepmc','lhe']]):
                    evts_file = file            
            evts_path = pjoin(interactionevts_dir,evts_file)


        if self.displaced:
            displ_dir = pjoin(self.dir_path,'Displaced_Decay')
            try:
                os.makedirs(pjoin(displ_dir,'Events'))
            except:
                pass
            
            run_dir_displ = pjoin(displ_dir,'Events',self.run_name)
            os.makedirs(run_dir_displ)
            files.ln(evts_path, run_dir_displ)
            displacement = displ_decay.displaced_decay(pjoin(run_dir_displ,'unweighted_events.lhe.gz'), pjoin(self.dir_path,'Cards','param_card.dat'),int(self.proc_characteristics['pdg_displ'])) 
            displacement.finalize_output(pjoin(run_dir_displ,'evt_displaced.lhe'))
            label = '#displaced_events'
            self.results[label] = displacement.total_events
            
            
        for dir in self.interaction_dir:
            cpath = pjoin(self.dir_path, dir, 'Cards')

            fit2D_card = fit2D.Fit2DCard(pjoin(cpath, 'fit2D_card.dat'))
            fit2D_card.write_include_file(pjoin(self.dir_path, dir, 'Source'))           
            
            try:
                os.remove('unweighted_events.lhe.gz')
            except OSError:
                pass
            
            if 'DIS' in dir:
                interaction_channel = 'DIS'
            elif 'electron' in dir:
                interaction_channel = 'electron'
            else:
                interaction_channel = None
                
            files.ln(evts_path, cpath)
            
            with misc.chdir(cpath):
                hist2D_energy_angle = meshfitter.fit2D_energy_theta(self.proc_characteristics, \
                                                'unweighted_events.lhe.gz',interaction_channel)
                if hist2D_energy_angle.npass < 1000:
                    raise Exception, "Error: numbers of events entering the detector too small!"
                hist2D_energy_angle.do_fit()

            misc.call(['./bin/generate_events', self.run_name, '-f'],
                      cwd=pjoin(self.dir_path, dir))
            
            run_dir = pjoin(self.dir_path, dir,'Events',self.run_name)

            # store info in the corresponding run_dir 
            try:
                output_tp = ['ehist.dat','cell_fortran.dat','in_DM.dat','mesh2D.png','fit1D.dat']
                for file in output_tp:
                    files.mv(pjoin(self.dir_path,dir,'Cards',file), pjoin(run_dir,file))
            except:
                pass

            # store results
            results = self.load_results_db(pjoin(self.dir_path,dir),self.run_name)
            label = 'nevts_' + interaction_channel
            self.results[label] = results['cross']

            # apply channel tag to events file name
            try:
                files.mv(pjoin(run_dir,'unweighted_events.lhe.gz'), pjoin(run_dir,'unweighted_events'+'_'+interaction_channel+'.lhe.gz'))
            except:
                raise Exception, 'Error: events file not generated!'
            
            # for DIS, generate the LHE events to be showered by Pythia
            if interaction_channel == 'DIS': 
                pythialhe = lheToPythia.LHEtoPYTHIAHadronSTD(pjoin(run_dir,'unweighted_events_DIS.lhe.gz'))
                pythialhe.write_PYTHIA_input(pjoin(run_dir,'pythiainput_DIS.lhe'))


    def load_results_db(self,path,run_name):
        """load the current results status"""
        print(path)
        # load the current status of the directory
        if os.path.exists(pjoin(path,'HTML','results.pkl')):
            try:
                results = save_load_object.load_from_file(pjoin(path,'HTML','results.pkl'))
                return results[run_name][0]
            except:
                raise Exception, "Error: no results!"
        
    def store_scan_result(self):
        return self.results#{'cross (pb)': 1e-5} #self.store_results

    def set_run_name(self, name):
        self.run_name = name
    
    def ask_run_configuration(self, mode=None, force=False):
        """ask the question about card edition """
        
        if self.DMmode == 'production_interaction': 
            cards = ['param_card.dat','fit2D_card.dat','run_card.dat']
        elif self.DMmode == 'production_decay':
            cards = ['param_card.dat','run_card.dat','madspin_card.dat']
        elif self.DMmode in ['decay_interaction','interaction_only']:
            cards = ['param_card.dat','fit2D_card.dat','madspin_card.dat']
        elif self.DMmode == 'decay_displaced':
            cards = ['param_card.dat','madspin_card.dat']
            
        self.ask_edit_cards(cards,plot=False)

        self.auto_width = set() #ensure to reset auto_width! at the 
        # self.check_param_card(pjoin(self.dir_path, 'Cards', 'param_card.dat'))
        # self.param_card = check_param_card.ParamCard(pjoin(self.dir_path, 'Cards', 'param_card.dat'))

        
    def ask_edit_cards(self, cards, mode='fixed', plot=True, first_cmd=None):
        """ """
        self.ask_edit_card_static(cards, mode, False, 60,
                                  self.ask, first_cmd=first_cmd)
        

    @staticmethod
    def ask_edit_card_static(cards, mode='fixed', plot=True,
                             timeout=0, ask=None, **opt):
        if not ask:
            ask = CommonRunCmd.ask

        def path2name(path):
            if '_card' in path:
                return path.split('_card')[0]
            elif path == 'delphes_trigger.dat':
                return 'trigger'
            elif path == 'input.lhco':
                return 'lhco'
            elif path == 'MadLoopParams.dat':
                return 'MadLoopParams'
            else:
                raise Exception, 'Unknow cards name %s' % path

        # Ask the user if he wants to edit any of the files
        #First create the asking text
        question = """Do you want to edit a card (press enter to bypass editing)?\n"""
        possible_answer = ['0', 'done']
        card = {0:'done'}
        
        indent = max(len(path2name(card_name)) for card_name in cards)
        question += '/'+'-'*60+'\\\n'
        for i, card_name in enumerate(cards):
            imode = path2name(card_name)
            possible_answer.append(i+1)
            possible_answer.append(imode)
            #question += '| %-77s|\n'%((' \x1b[31m%%s\x1b[0m. %%-%ds : \x1b[32m%%s\x1b[0m'%indent)%(i+1, imode, card_name))
            if card_name != 'run_card.dat':
                question += '| %-77s|\n'%((' \x1b[31m%%s\x1b[0m. %%-%ds : \x1b[32m%%s\x1b[0m'%indent)%(i+1, imode, card_name))
            else:
                question += '| %-77s|\n'%((' \x1b[31m%%s\x1b[0m. %%-%ds : \x1b[32m%%s\x1b[0m (production)'%indent)%(i+1, imode, card_name))
            card[i+1] = imode
            
        if plot and not 'plot_card.dat' in cards:
            question += '| %-77s|\n'%((' \x1b[31m9\x1b[0m. %%-%ds : \x1b[32mplot_card.dat\x1b[0m'%indent) % 'plot')
            possible_answer.append(9)
            possible_answer.append('plot')
            card[9] = 'plot'

        question += '\\'+'-'*60+'/\n'

        if 'param_card.dat' in cards:
            # Add the path options
            question += ' you can also\n'
            question += '   - enter the path to a valid card or banner.\n'
            question += '   - use the \'set\' command to modify a parameter directly.\n'
            question += '     The set option works only for param_card and run_card.\n'
            question += '     Type \'help set\' for more information on this command.\n'
            question += '   - call an external program (ASperGE/MadWidth/...).\n'
            question += '     Type \'help\' for the list of available command\n'
        else:
            question += ' you can also\n'
            question += '   - enter the path to a valid card.\n'
        if 'transfer_card.dat' in cards:
            question += '   - use the \'change_tf\' command to set a transfer functions.\n'

        out = 'to_run'
        while out not in ['0', 'done']:
            #if 'fit2D_card.dat' in cards:
            out = ask(question, '0', possible_answer, timeout=60,
                      path_msg='enter path', ask_class = MadDumpSelector,
                      cards=cards, mode=mode, **opt)
            # else:
            #     out = ask(question, '0', possible_answer, timeout=60,
            #               path_msg='enter path', ask_class = common_run.AskforEditCard,
            #               cards=cards, mode=mode, **opt)

    ############################################################################
    def do_open(self, line):
        """Open a text file/ eps file / html file"""

        args = self.split_arg(line)
        # Check Argument validity and modify argument to be the real path
        #self.check_open(args)
        file_path = args[0]

        misc.open_file(file_path)

    # def check_open(self, args):
    #     """ check the validity of the line """

    #     if len(args) != 1:
    #         self.help_open()
    #         raise self.InvalidCmd('OPEN command requires exactly one argument')

    #     if args[0].startswith('./'):
    #         if not os.path.isfile(args[0]):
    #             raise self.InvalidCmd('%s: not such file' % args[0])
    #         return True

    #     # if special : create the path.
    #     if not self.dir_path:
    #         if not os.path.isfile(args[0]):
    #             self.help_open()
    #             raise self.InvalidCmd('No MadEvent path defined. Unable to associate this name to a file')
    #         else:
    #             return True

    #     path = self.dir_path
    #     if os.path.isfile(os.path.join(path,args[0])):
    #         args[0] = os.path.join(path,args[0])
    #     elif os.path.isfile(os.path.join(path,'Cards',args[0])):
    #         args[0] = os.path.join(path,'Cards',args[0])
    #     # special for card with _default define: copy the default and open it
    #     elif '_card.dat' in args[0]:
    #         name = args[0].replace('_card.dat','_card_default.dat')
    #         if os.path.isfile(os.path.join(path,'Cards', name)):
    #             files.cp(os.path.join(path,'Cards', name), os.path.join(path,'Cards', args[0]))
    #             args[0] = os.path.join(path,'Cards', args[0])
    #         else:
    #             raise self.InvalidCmd('No default path for this file')
    #     elif os.path.isfile(os.path.join(path, '%s_card.dat' % args[0])):
    #         args[0] = os.path.join(path, '%s_card.dat' % args[0])
    #     elif os.path.isfile(os.path.join(path, 'Cards' ,'%s_card.dat' % args[0])):
    #         args[0] = os.path.join(path, 'Cards', '%s_card.dat' % args[0])
    #     elif not os.path.isfile(args[0]):
    #         raise self.InvalidCmd('No default path for this file')

    
class MadDumpSelector(common_run.AskforEditCard):
    """ """

    to_init_card = ['param', 'run', 'fit2D', 'madspin']

    def __init__(self, question, *args, **opts):
        self.me_dir = opts['mother_interface'].dir_path        
        common_run.AskforEditCard.__init__(self, question,
                                            *args, **opts)
        
    def write_switch(self, path=None):
        """store value of the switch for the default at next run"""
        
        if not path:
            path = pjoin(self.me_dir, 'Cards', '.lastrun')
        fsock = open(path,'w')
        for key, value in self.answer.items():
            fsock.write('%s %s\n' % (key,value))
        fsock.close()
           
           
    def detect_card_type(self, path):
        """detect card type """
        
        output = super(MadDumpSelector, self).detect_card_type(path)
        return output

    
    def init_fit2D(self, path):
        """ initialize cards for the reading/writing of maddump"""

        self.maddump_def = fit2D.Fit2DCard(self.paths['fit2D_default'], consistency=False)
        try:
            self.maddump = fit2D.Fit2DCard(self.paths['fit2D'], consistency=False)
        except Exception as e:
            logger.error('Current fit2D_card is not valid. We are going to use the default one.')
            logger.error('problem detected: %s' % e) 
            files.cp(self.paths['maddump_default'], 'Cards')
            self.maddump = fit2D.Fit2DCard('Cards')
            
        self.maddump_set = list(set(self.maddump_def.keys() + self.maddump_def.hidden_param))
        return self.maddump.keys() 

    # def get_cardcmd(self):
    #     """ return the list of command that need to be run to have a consistent 
    #         set of cards with the switch value choosen """
        
    #     cmd = super(MadDMSelector,self).get_cardcmd()
    #     for c in cmd:
    #         self.exec_cmd(c)        

    #     return cmd 

    
    def define_paths(self, **opt):
        
        super(MadDumpSelector, self).define_paths(**opt)
        self.paths['fit2D'] = pjoin(self.me_dir,'Cards','fit2D_card.dat')
        self.paths['fit2D_default'] = pjoin(self.me_dir,'Cards','fit2D_card_default.dat')
        self.paths['madspin'] = pjoin(self.me_dir,'Cards','madspin_card.dat')
        
      
    def do_update(self, line, timer=0):
        """syntax: update dependent: Change the mass/width of particles which are not free parameter for the model.
        update missing:   add to the current param_card missing blocks/parameters.
        Bypass dependent mode but if request by the user
        """
         
        if timer == 0:
            return super(MadDumpSelector, self).do_update(line)
        else: 
            args = self.split_arg(line)
            if args[0] == 'dependent':
                return
            else:
                return super(MadDumpSelector, self).do_update(line)
         
    
    def check_card_consistency(self):
        super(MadDumpSelector, self).check_card_consistency()
        
        
    def reload_card(self, path):
        """ensure that maddump object are kept in sync"""
        
        if path == self.paths['fit2D']:
            try:
                self.maddump = fit2D.Fit2DCard(path) 
            except Exception as e:
                logger.error('Current fit2D_card is not valid. We are going to use the default one.')
                logger.error('problem detected: %s' % e)
                logger.error('Please re-open the file and fix the problem.')
                logger.warning('using the \'set\' command without opening the file will discard all your manual change')
        else:
            return super(MadDumpSelector,self).reload_card(path)
        
        
    def complete_set(self, text, line, begidx, endidx, formatting=True):
        """ Complete the set command"""
        possibilities = super(MadDumpSelector,self).complete_set(text, line, begidx, endidx, formatting=False)
        args = self.split_arg(line[0:begidx])
        if len(args)>1 and args[1] == 'fit2d':
            start = 2
        else:
            start = 1 
            if len(args) ==1:
                possibilities['category of parameter (optional)'] += \
                                self.list_completion(text, ['fit2d'], line)
                
        if len(args)==start:
            correct = self.maddump_set + ['default']
            possibilities['fit2D Card'] = self.list_completion(text, correct, line)   
        elif len(args)==start+1:
            allowed_for_run = []
            if args[-1].lower() in self.maddump.allowed_value:
                allowed_for_run = self.maddump.allowed_value[args[-1].lower()]
                if '*' in allowed_for_run: 
                    allowed_for_run.remove('*')
            elif isinstance(self.maddump[args[-1]], bool):
                allowed_for_run = ['True', 'False']
            opts = [str(i) for i in  allowed_for_run]
            possibilities['fit2D Card'] = self.list_completion(text, opts)
        
        return self.deal_multiple_categories(possibilities, formatting)

    
    def do_set(self, line):
        """ edit the value of one parameter in the card"""
        
        args = self.split_arg(line)
        # fix some formatting problem
        if '=' in args[-1]:
            arg1, arg2 = args.pop(-1).split('=')
            args += [arg1, arg2]
        if '=' in args:
            args.remove('=')
        args = [ a.lower() for a in args]

        if args[0] in ['fit2d']:
            mode = 'fit2d'
            start = 1
            if args[1] == 'default':
                logging.info('replace %s by the default card' % args[0])
                self.maddump = fit2D.Fit2DCard(self.paths['fit2D_default'])
        else:
            start = 0
            mode = 'unknow'
             
        if args[start] in self.maddump_set:
            if args[start+1] == 'default':
                default = self.maddump_def[args[start]]
                self.setfit2D(args[start], default)
            else:
                if args[start] in self.maddump.list_parameter or \
                       args[start] in self.maddump.dict_parameter:
                    val = ' '.join(args[start+1:])
                    val = val.split('#')[0]
                    self.setfit2D(args[start], val)
                else:
                    self.setfit2D(args[start], args[start+1:][0])
        elif mode == 'unknow':
            return super(MadDumpSelector, self).do_set(line)
        else:
            logger.warning('%s not recognised as fit2D parameter. Please retry!', args[start+1])
            return 
        
        #write the new file
        self.maddump.write(self.paths['fit2D'])
        
    def setfit2D(self, name, value, loglevel=20):
        logger.log(loglevel,'modify parameter %s of the fit2D_card.dat to %s' % (name, value), '$MG:BOLD')
        self.maddump.set(name, value, user=True)
