import os

import madgraph.various.misc as misc
import madgraph.interface.common_run_interface as common_run
import madgraph.interface.extended_cmd as cmd
import madgraph.interface.launch_ext_program as launch_ext
import models.check_param_card as param_card_mod
import subprocess
import logging

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
        
        path1= pjoin(predir, 'DIS', name)
        path2= pjoin(predir, 'electron', name)
        path3= pjoin(predir, 'ENS', name)
        
        super(paramCardIterator, self).write_summary(path1, order)
        super(paramCardIterator, self).write_summary(path2, order)
        super(paramCardIterator, self).write_summary(path3, order)


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
        return self.proc_characteristics['BSM_model']
        # self.mg5.exec_cmd('import model name OPTIONS')
        # return self.mg5._curr_model

        
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

        if options['name']:
            self.run_name = options['name']
        else:
            i = 1
            while os.path.exists(pjoin(self.dir_path, 'production', 'Events', 'run_%02d' %i)) or\
                  os.path.exists(pjoin(self.dir_path, 'production', 'Events', 'run_%02d_01' %i)):
                i += 1
            self.run_name = 'run_%02d' % i

        #self.run_name = 'run_%02d' % 3
        
            
        interaction_dir_all = ['interaction_DIS','interaction_electron','interaction_ENS']
        listdir = subprocess.check_output("ls %s"%self.dir_path,shell=True).split()

        self.interaction_dir = []
        for dir in listdir:
            if dir in interaction_dir_all:
                self.interaction_dir.append(dir)

        self.proc_characteristics = self.get_proc_characteristics(pjoin(self.dir_path,self.interaction_dir[0],'SubProcesses','proc_characteristics'))
        
        # param_card_iterator.write(card_path) #-> this is done by the with statement
        #         name = misc.get_scan_name(orig_name, next_name)
        #         path = result_path(obj) % name 
        #         logger.info("write all cross-section results in %s" % path ,'$MG:BOLD')
        #         order = summaryorder(obj)()
        #         param_card_iterator.write_summary(path, order=order)
            
        # # determine run_name for name in the output directory:
        # if '-n' in args:
        #     self.run_name = args[args.index('-n')+1]
        #     if os.path.exists(pjoin(self.dir_path, 'output', self.run_name)):
        #         shutil.rmtree(pjoin(self.dir_path, 'output', self.run_name))
        #         try:
        #             shutil.rmtree(pjoin(self.dir_path, 'Indirect','Events', self.run_name))
        #         except Exception:
        #             pass
        # else:
        #     i = 1
        #     while os.path.exists(pjoin(self.dir_path, 'output', 'run_%02d' %i)) or\
        #           os.path.exists(pjoin(self.dir_path, 'output', 'run_%02d_01' %i)):
        #         i += 1
        #     self.run_name = 'run_%02d' % i
        
        self.ask_run_configuration(mode=[], force=force)
        self.run_launch()
        
    @common_run.scanparamcardhandling(iteratorclass=paramCardIterator)
    def run_launch(self):
        prod_dir = 'production'
        # generate events: production
        os.chdir(pjoin(self.dir_path,prod_dir))
        os.system('./bin/generate_events %s -f'%self.run_name) 
        
        misc.sprint(os.getcwd())
        os.chdir('..')
        misc.sprint(os.getcwd())
        
        for dir in self.interaction_dir:
            os.chdir(pjoin(dir,'Cards'))
            fit2D_card = fit2D.Fit2DCard(pjoin('fit2D_card.dat'))
            fit2D_card.write_include_file('../Source')            
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
                
            os.system('ln -s ../../'+prod_dir+'/Events/'+self.run_name+'/unweighted_events.lhe.gz')
            hist2D_energy_angle = meshfitter.fit2D_energy_theta(self.proc_characteristics, \
                                                'unweighted_events.lhe.gz',interaction_channel)
            hist2D_energy_angle.do_fit()
            os.chdir('../')
            os.system('./bin/generate_events %s -f'%self.run_name)
            os.chdir('../')
        os.chdir('../')

        
    def store_for_scan(self):
        return {}

    def set_run_name(self, name):
        self.run_name = name
                        
        # if self.param_card_iterator:
        #     self.run_name += '_01'

        # # create output directory.
        # os.mkdir(pjoin(self.dir_path, 'output', self.run_name))
        
                 
#         output = pjoin(self.dir_path, 'output', self.run_name, 'maddm.out') 
#         misc.call(['./maddm.x', pjoin('output', self.run_name, 'maddm.out')], cwd =self.dir_path)
#         #Here we read out the results which the FORTRAN module dumped into a file
#         #called 'maddm.out'. The format is such that the first line is always relic density
#         # , second line is the nucleon scattering cross section (SI) for proton, third is SI
#         # nucleon cross section for the neutron etc. If a quantity was not calculated, we output -1


#         # Define a dictionary holding the results
#         result = {}
#         result['GeV2pb*pb2cm2']   = GeV2pb*pb2cm2 # conversion factor                                                                                                           

#         mdm = self.param_card.get_value('mass', self.proc_characteristics['dm_candidate'][0])

#         result['tot_SM_xsec'] = -1
#         sigv_indirect = 0.

#         for line in open(pjoin(self.dir_path, output)):
      
#                 splitline = line.split()
#                 #If capture rate is calculated.
#                 if 'ccap' in line:
#                     oname = splitline[0].strip(':')+'_'+splitline[1]
#                     val = splitline[2]
#                     result[oname.split(':')[0] ] = val

#                 else:
#                     if self._two2twoLO:
#                         if 'sigma*v' in line: 
#                             sigv_temp = float(splitline[1])
#                             oname = splitline[0].split(':',1)[1] .split('_')
#                             oname = oname[0]+'_'+oname[1] #To eliminate the annoying suffix coming from the '/' notation
#                             result['taacsID#'     + oname] = sigv_temp 
#                             result['err_taacsID#' + oname] = 0 
#                             if oname.split('_')[1] in ['uux','ddx','ssx']:
#                                 result['lim_taacsID#'+oname] = self.limits.ID_max(mdm, 'qqx')
#                             elif oname.split('_')[1] in self.limits._allowed_final_states:
#                                 result['lim_taacsID#'+oname] = self.limits.ID_max(mdm, oname.split('_')[1]) 
#                             elif oname.split('_')[1] not in self.limits._allowed_final_states:
#                                 result['lim_taacsID#'+oname] = -1
#                             sigv_indirect += sigv_temp
#                     result[splitline[0].split(':')[0]] = float(splitline[1])

                    
#         result['taacsID'] = sigv_indirect
                            
#         np_names = ['g','nue','numu','nutau']

#         if str(self.mode['indirect']).startswith('flux'):
#             for chan in np_names:
#                 result['flux_%s' % chan] = -1.0

#         # Calculating the xsi factor. Set == 1 if relic is not evaluated
#         if self.mode['relic']:
#            if result['Omegah^2'] < 0:
#                result['xsi'] = 1.0   
#            elif result['Omegah^2'] < self.limits._oh2_planck and result['Omegah^2'] > 0:
#                result['xsi'] = result['Omegah^2'] / self.limits._oh2_planck
#            else: result['xsi'] = 1.0

#         else: result['xsi'] = 1.0

#         if self.mode['direct']:
#             result['lim_sigmaN_SI_n'] = self.limits.SI_max(mdm)
#             result['lim_sigmaN_SI_p'] = self.limits.SI_max(mdm)
#             result['lim_sigmaN_SD_p'] = self.limits.SD_max(mdm, 'p')
#             result['lim_sigmaN_SD_n'] = self.limits.SD_max(mdm, 'n')

#         self.last_results = result
#         self.last_results['run'] = self.run_name

                          
# #        if self.mode['indirect'] and not self._two2twoLO:
# #            with misc.MuteLogger(names=['madevent','madgraph'],levels=[50,50]):
# #                self.launch_indirect(force)

#         if self.mode['indirect']:
#             self.launch_indirect(force)
        
#         if not self.in_scan_mode and not self.multinest_running:
#             self.print_results()

#         # Saving the output for single point
#         if not self.param_card_iterator:
#             self.save_remove_output(scan = False)

#         # --------------------------------------------------------------------#
#         #   THIS PART IS FOR MULTINEST SCANS
#         # --------------------------------------------------------------------#

#         #multinest can be launched only after one launch has been executed
#         if self.mode['nestscan'] and not self.multinest_running:
#             self.multinest_running = True
#             self.launch_multinest()
#             self.multinest_running = False


#         # --------------------------------------------------------------------#
#         #   THIS PART IS FOR SEQUENTIAL SCANS
#         # --------------------------------------------------------------------#
        
#         if self.param_card_iterator:
            
#             param_card_iterator = self.param_card_iterator

#             parameters, values =  param_card_iterator.param_order , param_card_iterator.itertag
#             self.param_card_iterator = []

#             self.save_remove_output(scan = True) ## this is to remove or save spectra, not the scan summary file!

#             # *** Initialize a list containing the desired variables in the summary output of the scan
#             order = ['run']

#             # *** saving the name of the iterated parameters, and the values in the results dictionary
#             for par,val in zip(parameters, values):
#                 order.append(par)
#                 self.last_results[par] = val

#             # *** Relic density
#             if self.mode['relic']:
#                 order += ['Omegah^2','x_f', 'sigmav(xf)']
#             order.append('xsi')

#             # *** Direct Detection
#             if self.mode['direct'] :
#                 order += ['sigmaN_SI_p', 'lim_sigmaN_SI_p', 
#                           'sigmaN_SI_n', 'lim_sigmaN_SI_n',
#                           'sigmaN_SD_p', 'lim_sigmaN_SD_p',
#                           'sigmaN_SD_n', 'lim_sigmaN_SD_n']

           
#             if self.mode['direct'] == 'directional':
#                 order += ['Nevents', 'smearing']

#             if self.mode['capture']:
#                 detailled_keys = [k for k in self.last_results if k.startswith('ccap_') and '#' not in k]
#                 for key in detailled_keys:
#                     order += [key]
 
#             # *** Indirect detection
#             if self.mode['indirect']:
#                 halo_vel = self.maddm_card['vave_indirect']
#                 if halo_vel > (3*10**(-6)) and halo_vel < ( 1.4*10**(-4) ):  # cannot calculate Fermi limits
                 
#                 #if not self._two2twoLO:
#                     detailled_keys = [k for k in self.last_results if k.startswith('taacsID#') ]
#                     if len(detailled_keys)>1:
#                         for key in detailled_keys:
#                             clean_key_list = key.split("_")
#                             clean_key = clean_key_list[0]+"_"+clean_key_list[1]

#                             order +=[clean_key]
#                             order +=['lim_'+clean_key]

#                     order.append('taacsID')
#                     order.append('tot_SM_xsec')
#                     order.append('Fermi_sigmav')

#                     if self.last_results['xsi'] >0 and self.last_results['xsi'] <1: # thermal and non thermal case 
#                        order = order + ['pvalue_th','like_th','pvalue_nonth','like_nonth']
#                     else:
#                        order = order + ['pvalue_nonth','like_nonth']

#                     if self.mode['indirect'].startswith('flux'):
#                        #for channel in self.Spectra.spectra.keys(): #['neutrinos_e', 'neutrinos_mu' , 'neutrinos_tau']
#                        for channel in ['gammas','neutrinos_e', 'neutrinos_mu' , 'neutrinos_tau']:
#                            if 'antip' in channel or 'pos' in channel: continue
#                            order.append('flux_%s' % channel)

#             for elem in order:
#                 if elem not in self.last_results.keys():
#                     order.remove(elem)

#             run_name = str(self.run_name).rsplit('_',1)[0]
#             summary_file = pjoin(self.dir_path, 'output','scan_%s.txt' % run_name)
#             self.write_scan_output(out_path = summary_file , keys = order, header = True )
#             self.write_scan_output(out_path = summary_file , keys = order )
#             with misc.TMP_variable(self, 'in_scan_mode', True):
#                 with misc.MuteLogger(names=['cmdprint','madevent','madgraph','madgraph.plugin'],levels=[50,50,50,20]):
                                        
#                     for i,card in enumerate(param_card_iterator):
#                         card.write(pjoin(self.dir_path,'Cards','param_card.dat'))
#                         self.exec_cmd("launch -f -n %s_%02d" % (run_name, i+2),
#                                        precmd=True, postcmd=True, errorhandling=False)
#                         for par,val in zip(param_card_iterator.param_order, param_card_iterator.itertag):
#                             self.last_results[par] = val
#                         self.write_scan_output(out_path = summary_file , keys = order, header = False)

#                         ### the following three lines are added by chiara to check the widht = auto function 
#                         # self._param_card = param_card_mod.ParamCard('/Users/arina/Documents/physics/software/maddm_dev2/test_width/Cards/param_card.dat')
#                         # width = self.param_card.get_value('width', 5000000)
#                         # logger.warning('--> try again WY0: %.2e' % width)
#                         #<=-------------- Mihailo commented out max_col = 10
#                         #logger.info('Results for the point \n' + param_card_iterator.write_summary(None, order, lastline=True,nbcol=10)[:-1])#, max_col=10)[:-1])
#                         self.save_remove_output(scan = True)



#             param_card_iterator.write(pjoin(self.dir_path,'Cards','param_card.dat'))

    def ask_run_configuration(self, mode=None, force=False):
        """ask the question about card edition """
        
        # if not hasattr(self, 'mode') or not force:
            # self.mode, cmd_quest = self.ask('', '0', mode=mode, 
            #             ask_class=MadDumpSelector, timeout=60, path_msg=' ',
            #             return_instance=True, force=force)
        # out = self.ask('', '0', mode=mode, 
        #                ask_class=MadDumpSelector, timeout=60, path_msg=' ',
        #                return_instance=True, force=force)

        cards = ['param_card.dat','fit2D_card.dat','run_card.dat']
        self.ask_edit_cards(cards,plot=False)

        self.auto_width = set() #ensure to reset auto_width! at the 
        # self.check_param_card(pjoin(self.dir_path, 'Cards', 'param_card.dat'))
        # self.param_card = check_param_card.ParamCard(pjoin(self.dir_path, 'Cards', 'param_card.dat'))

        
    def ask_edit_cards(self, cards, mode='fixed', plot=True, first_cmd=None):
        """ """
        # if not self.options['madanalysis_path']:
        #     plot = False

        self.ask_edit_card_static(cards, mode, False, 60,
                                  self.ask, first_cmd=first_cmd)
        
        # for c in cards:
        #     if not os.path.isabs(c):
        #         c = pjoin(self.me_dir, c) 
        #     if not os.path.exists(c):
        #         default = c.replace('dat', '_default.dat')
        #         if os.path.exists(default):
        #             files.cp(default, c)

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
            out = ask(question, '0', possible_answer, timeout=60,
                              path_msg='enter path', ask_class = MadDumpSelector,
                              cards=cards, mode=mode, **opt)

        
            # # automatically switch to keep_wgt option
            # #edit the maddm_card to be consistent with self.mode
            # cmd_quest.get_cardcmd()
            # # write Cards/.lastmode to recover 
            # cmd_quest.write_switch()
    
            # self.maddump_card = cmd_quest.maddump
            # for key, value in self.mode.items():
            #     if value == 'ON' or value is True:
            #         self.mode[key] = True
    
            #     elif value == 'OFF':
            #         self.mode[key] = False
            # self.mode['capture'] = False

            # # create the inc file for maddm
            # logger.debug('2to2 in ask_run_configuration: %s' % self._two2twoLO)

        
        # #set fortran switch and write include file
        # if not self.in_scan_mode:
        #     # create the inc file for maddm
        #     self.maddm_card.set('do_relic_density', self.mode['relic'], user=False)
        #     self.maddm_card.set('do_direct_detection', True if self.mode['direct'] else False, user=False)
        #     self.maddm_card.set('do_directional_detection', self.mode['direct'] == 'directional', user=False)
        #     self.maddm_card.set('do_capture', self.mode['capture'], user=False)
        #     self.maddm_card.set('do_indirect_detection', True if self.mode['indirect'] else False, user=False)
        #     self.maddm_card.set('do_flux', True if (self.mode['indirect'] and self.mode['indirect'] != 'sigmav') else False, user=False)
        #     self.maddm_card.set('only2to2lo', self._two2twoLO, user=False)
        #     #self.maddm_card.set('run_multinest', self.mode['run_multinest'], user=False)

        # if not self.in_scan_mode and not self.mode['nestscan']:
        #     logger.info("Start computing %s" % ','.join([name for name, value in self.mode.items() if value]))
        # return self.mode


    # def compile(self):
    #     """compile the code"""

    #     #logger.info(self.mode)

    #     if self.in_scan_mode:
    #         return

    #     #self.maddump_card.write_include_file(pjoin(self.dir_path,'include'))
    #     misc.compile(['all'],cwd=self.dir_path)

    #     # if self.mode['relic'] and self.mode['direct']:
    #     #     misc.compile(['all'],cwd=self.dir_path)
    #     # elif self.mode['relic'] and not self.mode['direct']:
    #     #     misc.compile(['relic_density'],cwd=self.dir_path)
    #     # elif self.mode['direct'] and not self.mode['relic']:
    #     #     misc.compile(['direct_detection'],cwd=self.dir_path)
    #     # elif self.mode['indirect'] and not self.mode['relic']:
    #     #     misc.compile(['relic_density'],cwd=self.dir_path)
    #     # else:
    #     #     raise Exception, "No computation requested. End the computation"
    
    #     # if os.path.exists(pjoin(self.dir_path, 'src', 'maddm.x')) or os.path.exists(pjoin(self.dir_path, 'maddm.x')):
    #     #     logger.info("compilation done")
    #     # else:
    #     #     raise Exception, 'Compilation of maddm failed'

        
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

    to_init_card = ['param', 'run', 'fit2D']

    def __init__(self, question, *args, **opts):

        self.me_dir = opts['mother_interface'].dir_path
        #self.availmode = opts.pop('data', collections.defaultdict(bool))
        
        # param_card_path = pjoin(opts['mother_interface'].dir_path, 'Cards', 'param_card.dat')
        # fit2D_card_path = pjoin(opts['mother_interface'].dir_path, 'Cards', 'fit2D_card.dat')
        # pythia8_card_path = pjoin(opts['mother_interface'].dir_path, 'Cards', 'pythia8_card.dat')
        
        # cards = [param_card_path, fit2D_card_path]

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
        
#     # TODO HERE!
#     def default(self, line):
#         """Default action if line is not recognized"""
        
#         try:
#             return cmd.ControlSwitch.default(self, line, raise_error=True)
#         except cmd.NotValidInput:
#             return common_run.AskforEditCard.default(self, line)     
        
            
#     def do_help(self, line, conflict_raise=False, banner=True):
#         """proxy for do_help"""
             
#         if line:
#             if banner:                      
#                 logger.info('*** HELP MESSAGE ***', '$MG:BOLD')
#             card = common_run.AskforEditCard.do_help(self, line, conflict_raise=conflict_raise, banner=False)
#         else:
#             if banner:
#                 logger.info('*** HELP MESSAGE FOR CARD EDITION ***', '$MG:BOLD')
#             card = common_run.AskforEditCard.do_help(self, line, conflict_raise=conflict_raise, banner=False)
#             logger.info('*** HELP MESSAGE FOR CHANGING SWITCH ***', '$MG:BOLD')
#             card = cmd.ControlSwitch.do_help(self, line, list_command=False)
#             if banner:                      
#                 logger.info('*** END HELP ***', '$MG:BOLD')
            
#             logger_tuto.info("""
# This question allows you BOTH to define what you are going to run (via the value at the top).
# But also to edit the content of the various file defining the  run/benchmark (bottom).

# To bypass the computation of relic density you can do
# > relic=OFF            
# to make the compuatation of the  directional detection
# > direct=directional
                
# You can also edit the card referenced.
# Note that you can 
#    1) edit any parameter like this:
#        > set mxd 10
#        [use auto completion if you need to search a name]
#    2) run a scan over parameter space
#         > set mxd scan:[10, 20, 40]
#         > set mxd scan:[10**i for i in range(5)]
#    3) ask to compute the width automatically for a particle
#         > set my0 Auto  
        
# When you are done with such edition, just press enter (or write 'done' or '0')          
# """) 

#             return
            
#         args = self.split_arg(line)
        
#         if not args:
#             args =['']
        
#         start = 0
        
#         if args[start] in ['maddm', 'maddm_card']:
#             card = 'maddm'
#             start += 1
            
#         #### MADDM CARD 
#         if args[start] in [l.lower() for l in self.maddm.keys()] and card in ['', 'maddm']:
#             if args[start] not in self.maddm_set:
#                 args[start] = [l for l in self.maddm_set if l.lower() == args[start]][0]

#             if args[start] in self.conflict and not conflict_raise:
#                 conflict_raise = True
#                 logger.info('**   AMBIGUOUS NAME: %s **', args[start], '$MG:BOLD')
#                 if card == '':
#                     logger.info('**   If not explicitely speficy this parameter  will modif the maddm_card file', '$MG:BOLD')

#             self.maddm.do_help(args[start])
            
#         ### CHECK if a help_xxx exist (only for those wihtout a do_xxx)
#         if len(args) == 1:
#             if not hasattr(self, 'do_%s' % args[0]) and hasattr(self, 'help_%s' % args[0]):
#                 getattr(self, 'help_%s' %args[0])()
            
#         if banner:                      
#             logger.info('*** END HELP ***', '$MG:BOLD')  
#         return card    
    
      
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
         
    # def update_to_full(self, line):
    #     """ trigger via update to_full LINE"""
        
    #     logger.info("update the maddm_card by including all the hidden parameter")
    #     self.maddm.full_template = True
    #     self.maddm.write(self.paths['maddm'], write_hidden=True)

    
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

        misc.sprint(args[0])
        if args[0] in ['fit2d']:
            start = 1
            if args[1] == 'default':
                logging.info('replace %s by the default card' % args[0])
                self.maddump = fit2D.Fit2DCard(self.paths['fit2D_default'])
                return
        else:
            return super(MadDumpSelector, self).do_set(line)
        
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
        #write the new file

        self.maddump.write(self.paths['fit2D'])
        
    def setfit2D(self, name, value, loglevel=20):
        logger.log(loglevel,'modify parameter %s of the fit2D_card.dat to %s' % (name, value), '$MG:BOLD')
        self.maddump.set(name, value, user=True)


        
