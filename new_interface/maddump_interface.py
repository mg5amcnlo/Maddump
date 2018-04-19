import logging
import os

#import maddm_run_interface as maddm_run_interface

import madgraph.core.diagram_generation as diagram_generation
import madgraph.interface.master_interface as master_interface
import madgraph.interface.madgraph_interface as madgraph_interface
import madgraph.various.misc as misc

import madgraph.interface.common_run_interface as common_run

import models.check_param_card as check_param_card
import models.model_reader as model_reader
import madgraph.core.helas_objects as helas_objects

import re
pjoin = os.path.join

logger = logging.getLogger('madgraph.plugin.maddump')

# class bcolors:
#     HEADER = '\033[95m'
#     OKBLUE = '\033[94m'
#     GRAY = '\033[90m'
#     OKGREEN = '\033[92m'
#     WARNING = '\033[93m'
#     FAIL = '\033[91m'
#     ENDC = '\033[0m'
#     BOLD = '\033[1m'
#     UNDERLINE = '\033[4m'

class DMError(Exception): pass

# Root path
MDumpDIR = os.path.dirname(os.path.realpath( __file__ ))

    
class MadDump_interface(master_interface.MasterCmd):

    intro_banner= 'MadDump plugin' 

    _define_options = ['darkmatter']
    
    # process number to distinguish the different type of matrix element
    process_tag = {'prod': 900, 
                   'DIS': 1100,
                   'electron': 1200,
                   'ENS': 1300}
    
    # eff_operators_SI = {1:'SIEFFS', 2:'SIEFFF', 3:'SIEFFV'}
    # eff_operators_SD = {1:False, 2:'SDEFFF', 3:'SDEFFV'} 
    
    def __init__(self, *args, **opts):
        
        super(MadDump_interface, self).__init__(*args, **opts)
        self._dm_candidate = []
        self._param_card = None
                
        
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
                dirname = 'maddump_'+str(i)
                try:
                    os.makedirs(dirname)
                except:
                    i+=1
                    continue
                os.chdir(dirname)
                not_created = False                

        else:
            dirname = args[0]
            try:
                os.makedirs(dirname)
            except:
                raise DMError, 'Output command error: the directory %s already exists!' % dirname 
            os.chdir(dirname)


        # print(len(self._curr_proc_defs))
        # print(len(self._curr_matrix_elements))
        # print(self._curr_matrix_elements)
        # print(len(self._curr_amps))
        # # print(self._curr_amps)
        # print(self._done_export)
        
        # exit(-1)
              
        for proc in self._curr_proc_defs:
            
            if proc['id'] < 1000:
                self._curr_matrix_elements = helas_objects.HelasMultiProcess()
                amps = [ amp for amp in self._curr_amps if amp.get('process')['id'] == proc['id'] ]
                line = 'production'
                with misc.TMP_variable(self, ['_curr_proc_defs','_curr_amps','_done_export'], [proc,amps,False]):                    
                    super(MadDump_interface, self).do_output(line)
            else:
                self._curr_matrix_elements = helas_objects.HelasMultiProcess()
                processes_list = self.process_tag.keys()
                for process in processes_list:
                    if self.process_tag[process] == proc['id']:
                        channel = process
                        break
                print(proc['id'])
                print(len(self._curr_amps))
                amps = [ amp for amp in self._curr_amps if amp.get('process')['id'] == proc['id'] ]
                print(len(amps))
                
                line = 'maddump interaction_' + channel 
                with misc.TMP_variable(self, ['_curr_proc_defs','_curr_amps','_done_export'], [proc,amps,False]):                    
                    super(MadDump_interface, self).do_output(line)
                current_dir = 'interaction_' + channel
                paramcard_path = '/Cards/param_card.dat'
                try:
                    os.remove(current_dir+paramcard_path)
                except OSError:
                    pass
                # os.symlink(os.path.join('production',paramcard_path), os.path.join(current_dir, paramcard_path))
                os.symlink('../../production'+paramcard_path, current_dir+paramcard_path)
                    
    # def find_output_type(self, path):
    #     if os.path.exists(pjoin(path,'matrix_elements','proc_characteristics')):
    #         return 'maddm'
    #     else:
    #         return super(MadDM_interface, self).find_output_type(path)

    def do_launch(self, line):
        
        
        args = self.split_arg(line)
        (options, args) = madgraph_interface._launch_parser.parse_args(args)
        with misc.TMP_variable(self, '_export_formats', self._export_formats + ['maddm']):
            self.check_launch(args, options)
        options = options.__dict__        
        
        
        if args[0] != 'maddm':
            return super(MadDM_interface, self).do_launch(line)
        else:
            MDM = maddm_run_interface.MADDMRunCmd(dir_path=args[1], options=self.options) 
            if options['interactive']:              
                return self.define_child_cmd_interface(MDM)
            else:
                self.define_child_cmd_interface(MDM,  interface=False)
                MDM.exec_cmd('launch ' + line.replace(args[1], ''))
                return


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






        
        
    # def generate_relic(self, excluded_particles=[]):
    #     """--------------------------------------------------------------------#
    #     #                                                                      #
    #     #  This routine generates all the 2 -> 2 annihilation matrix elements. #
    #     #  The SM particles, BSM particles in self._bsm_final_states, and the  #
    #     #  DM particles are used as final states.  The initial states are      #
    #     #  looped over all the different combinations of DM particles.         #
    #     #--------------------------------------------------------------------"""

    #     #Here we define bsm multiparticle as a special case for exlusion only in the
    #     #final state. This is different than the standard / notation which will exclude
    #     #the particles from the propagators as well.
    #     #Since we use the exluded_particles everywhere, we just use the presence of
    #     #the bsm in the array to initialize a flag and then remove it from the excluded
    #     #particles array so as not to conflict with the use in other places.
    #     self.define_multiparticles('bsm',[p for p in self._curr_model.get('particles')\
    #                      if abs(p.get('pdg_code')) > 25])
    #     if 'bsm' in excluded_particles:
    #         exclude_bsm = True
    #         while 'bsm' in excluded_particles:
    #             excluded_particles.remove('bsm')
    #     else:
    #         exclude_bsm = False


    #     if not self._dm_candidate:
    #         self.search_dm_candidate(excluded_particles)
    #         if not self._dm_candidate:
    #             return
    #         self.search_coannihilator(excluded=excluded_particles)


    #     # Tabulates all the BSM particles so they can be included in the
    #     # final state particles along with the SM particles.
        
    #     #dm_mass = abs(self._curr_model.get_mass(self._dm_candidate[0]))
    #     #is_lower_mass = lambda p : True #abs(self._curr_model.get_mass(p)) < (1-self.coannihilation_diff)*dm_mass
        
    #     ids_veto = [p.get('pdg_code') for p in self._dm_candidate + self._coannihilation]     
    #     bsm_final_states = [p for p in self._curr_model.get('particles') 
    #                      if abs(p.get('pdg_code')) > 25 and \
    #                      (p.get('name') not in excluded_particles or p.get('antiname') not in excluded_particles) and\
    #                      p.get('width') != 'ZERO' and
    #                      abs(p.get('pdg_code')) not in ids_veto]
    #     if not exclude_bsm:
    #         if bsm_final_states:
    #             logger.info("DM is allowed to annihilate into the following BSM particles: %s",
    #                      ' '.join([p.get('name') for p in bsm_final_states]))
    #             logger.info("use generate relic_density / X to forbid the decay the annihilation to X")
    #             logger.info("if you want to forbid DM to annihilate into BSM particles do relic_density / bsm")
    #     else:
    #         bsm_final_states=[]

    #     # Set up the initial state multiparticles that contain the particle 
    #     #and antiparticle
    #     for i,dm in enumerate(self._dm_candidate + self._coannihilation):
    #         self.define_multiparticles('dm_particle_%s' % i,  [dm])
    #         #self.do_display('multiparticles')

    #     self.define_multiparticles('dm_particles', self._dm_candidate+self._coannihilation)

    #     self.do_display('multiparticles')
    #     sm_pdgs = range(1, 7) + range(11, 17) + range(21, 26) #quarks/leptons/bosons
    #     sm_part = [self._curr_model.get_particle(pdg) for pdg in sm_pdgs]

    #     self.define_multiparticles('fs_particles', sm_part +bsm_final_states)

    #     # generate annihilation diagram
    #     coupling = "QED<=4 SIEFFS=0 SIEFFF=0 SIEFFV=0 SDEFFF=0 SDEFFV=0"
    #     if excluded_particles:
    #         proc = "dm_particles dm_particles > fs_particles fs_particles / %s %s  @DM2SM" \
    #                 % (' '.join(excluded_particles), coupling)
    #     else:
    #         proc = "dm_particles dm_particles > fs_particles fs_particles %s @DM2SM" \
    #                %(coupling)  #changed here
                   
    #     self.do_generate(proc)

    #     # Generate all the DM -> DM processes
    #     nb_dm = len(self._dm_candidate + self._coannihilation)
    #     for i in xrange(nb_dm):
    #         for j in xrange(i,nb_dm):
    #             if  excluded_particles:
    #                 proc = "DM_particle_%s DM_particle_%s > dm_particles dm_particles / %s %s @DM2DM"\
    #                    % (i,j,' '.join(excluded_particles), coupling)
    #             else:
    #                 proc = "DM_particle_%s DM_particle_%s > dm_particles dm_particles %s @DM2DM"\
    #                    % (i,j, coupling)
    #             try:
    #                 self.do_add('process %s' % proc)
    #             except (self.InvalidCmd,diagram_generation.NoDiagramException) :
    #                 continue
    #     # Get the matrix elements and make sure that we don't have any pure 
    #     #scattering processes
    #     for amp in self._curr_amps[:]:
    #         if amp.get('process').get('id') != self.process_tag['DM2DM']:
    #             continue
    #         if set(amp.get('process').get_initial_ids()) == (set(amp.get('process').get_final_ids())):
    #             self._curr_amps.remove(amp)


        # Generate all the DM particles scattering off of the thermal background and
        # change to a different DM species (again, no pure scatterings)

        # for i in xrange(nb_dm):
        #     for j in xrange(nb_dm):
        #         if i == j:
        #             continue
        #         if excluded_particles:
        #             proc = "DM_particle_%s sm_particles > DM_particle_%s sm_particles / %s %s @DMSM"\
        #                % (i,j,' '.join(excluded_particles), coupling)
        #         else:
        #             proc = "DM_particle_%s sm_particles > DM_particle_%s sm_particles %s @DMSM"\
        #                % (i,j, coupling)
        #         try:
        #             self.do_add('process %s' % proc)
        #         except (self.InvalidCmd,diagram_generation.NoDiagramException) :
        #             continue

    # def generate_interactions(self, excluded_particles=[]):
    #     """User level function which performs direct detection functions        
    #     """

    #     if not self._dm_candidate:
    #         self.search_dm_candidate(excluded_particles)
    #         if not self._dm_candidate:
    #             return

    #     if len(self._dm_candidate) > 1:
    #         logger.warning("More than one DM candidate. Can not run Direct Detection.")
    #         return 
  
   
        
    #     #Now figure out the label of the effective vertex to use. The convention is:
    #     #<SI or SD>EFF<F, S, or V>, i.e. SIEFFV for Spin Independent Effective vertex for Vector Current.
    #     dm_spin = int(self._dm_candidate[0]['spin'])
    #     eff_operators_SI = self.eff_operators_SI[dm_spin]
    #     eff_operators_SD = self.eff_operators_SD[dm_spin]
        
    #     logger.info("Generating X Nucleon > X Nucleon diagrams from the full lagrangian...")
    #     has_direct = self.DiagramsDD(eff_operators_SI, eff_operators_SD, 'QED')

    #     if not has_direct:
    #         logger.warning("No Direct Detection Feynman Diagram")
    #         return
        
    #     logger.info("Generating X Nucleon > X Nucleon diagrams from the effective lagrangian...")
    #     #ONLY EFFECTIVE LAGRANGIAN
    #     self.DiagramsDD(eff_operators_SI, eff_operators_SD, 'SI')

    #     logger.info("INFO: Generating X Nucleon > X Nucleon diagrams from the effective+full lagrangian...")
    #     #EFFECTIVE + FULL
    #     self.DiagramsDD(eff_operators_SI, eff_operators_SD, 'SI+QED')
        
    #     if (eff_operators_SD != False):
    #         logger.info("Doing the spin dependent part...")
    #         logger.info("Generating X Nucleon > X Nucleon diagrams from the effective lagrangian...")

    #         self.DiagramsDD(eff_operators_SI, eff_operators_SD, 'SD')
    #         #EFFECTIVE + FULL
    #         logger.info("Generating X Nucleon > X Nucleon diagrams from the effective + full lagrangian...")
    #         self.DiagramsDD(eff_operators_SI, eff_operators_SD, 'SD+QED')
        

        
        
    # #-----------------------------------------------------------------------#
    # def DiagramsDD(self, SI_name, SD_name, type, excluded=[]):
    #     """Generates direct detection diagrams. i_dm is the index of DM part. 
    #              Whether spin dependent or spin independent diagrams should be                 
    #              calculated. XX_order parameters determine maximum order of  a
    #              coupling. For instance, if you only want effective vertices for
    #              spin independent coupling you would set SI_order = 2 and all others
    #              to zero. If you want the spin independent full lagrangian + eff.
    #              then you need to set SI_order=2 and QED_order=2...
    #              WARNING: This function is a proxy to be used inside
    #              GenerateDiagramsDirDetect() function and no place else!
    #     """                                                             
        
    #     quarks = range(1,7) # ['d', 'u', 's', 'c', 'b','t']
    #     antiquarks = [-1*pdg for pdg in quarks] # ['d~', 'u~', 's~', 'c~', 'b~','t~']
        
    #     if type == 'SI':
    #         orders = '%s==2' % SI_name
    #     elif type == 'SD':
    #         orders = '%s==2' % SD_name
    #     elif type == 'QED':
    #         orders = '%s=0 %s=0' % (SD_name, SI_name)
    #     elif type == 'SI+QED':
    #         orders = '%s=0 %s<=99' % (SD_name, SI_name)
    #     elif type == 'SD+QED':
    #         orders = '%s<=99 %s=0' % (SD_name, SI_name)
            
    #     #loop over quarks
    #     has_diagram = False
    #     for i in quarks + antiquarks:
    #         proc = ' %(DM)s %(P)s > %(DM)s %(P)s %(excluded)s %(orders)s @DD' %\
    #                 {'DM': self._dm_candidate[0].get('name'),
    #                  'P': i,
    #                  'excluded': ('/ %s' % ' '.join(excluded) if excluded else ''),
    #                  'orders': orders
    #                  }
            
    #         try:
    #             self.do_add('process %s' % proc)
    #         except (self.InvalidCmd, diagram_generation.NoDiagramException), error:
    #             logger.debug(error)
    #             continue # no diagram generated
    #         has_diagram = True
    #     return has_diagram

    # def generate_indirect(self, argument):
    #     """User level function which performs direct detection functions        
    #        Generates the DM - q,g scattering matrix elements for spin dependent 
    #        and spin independent direct detection cross section                   
    #        Currently works only with canonical mode and one dm candidate.        
    #        The function also merges the dark matter model with the effective        
    #        vertex model.          
    #     related to syntax: generate indirect a g / n3  
    #     """

    #     if not self._dm_candidate:
    #         self.search_dm_candidate(excluded_particles)
    #         if not self._dm_candidate:
    #             return
        
    #     # separate final state particle from excluded particles
    #     if '/' in argument:
    #         ind = argument.find('/')
    #         particles, excluded = argument[:ind], argument[ind+1:]
    #     elif any(a.startswith('/') for a in argument):
    #         line = ' '.join(argument)
    #         particles, excluded = line.split('/',1)
    #         particles = particles.split()
    #         excluded = excluded.replace('/','').split()
    #     else:
    #         particles = argument
    #         excluded = []
        
        
    #     #handling the final state to ensure model independant support
    #     if 'v' in particles:
    #         particles.remove('v')
    #         particles += ['12', '14', '16']
        
    #     antiparticles = []
    #     for i,p in enumerate(particles):
    #         if p.isdigit():
    #             if p not in ['22','21','12','14','16']:
    #                 raise self.InvalidCmd, '%s is not a valid final state for indirect detection' % p
    #             antiparticles.append(str(-1*int(p)))   
    #         else:
    #             if p in ['a', 'g', 've','vm','vt']:
    #                 p ={'a':22, 'g':21, 've':12,'vm':14,'vt':16}[p]
    #             misc.sprint(p)
    #             part = self._curr_model.get_particle(p)
    #             misc.sprint(part)
    #             if part.get('pdg_code') not in [22,21,12,14,16]:
    #                 raise self.InvalidCmd, '%s is not a valid final state for indirect detection' % p
    #             particles[i] = part.get('name') 
    #             antiparticles.append(part.get('antiname'))
        
    #     # First try LO matrix-element
    #     done= []
    #     for dm in self._dm_candidate:
    #         name = dm.get('name')
    #         antiname = dm.get('name')
    #         if name in done:
    #             continue
    #         done += [name, antiname]
    #         for p, antip in zip(particles,antiparticles):
    #             proc = '%s %s > %s %s @ID' % (name, antiname, p,antip)
    #             try:
    #                 self.do_add('process %s' % proc)
    #             except (self.InvalidCmd, diagram_generation.NoDiagramException), error:
    #                 proc = '%s %s > %s %s [virt=ALL] @ID' % (name, antiname, p,antip)
    #                 self.do_add('process %s' % proc)        

      
#     def update_model_with_EFT(self):
#         """ """
#         eff_operators_SD = {1:False, 2:'SDEFFF', 3:'SDEFFV'}
#         eff_model_dm_names = {1:'~sdm', 2:'~fdm', 3:'~vdm'}

#         DM = self._dm_candidate[0]
        
#         if self._dm_candidate[0]['self_antipart']: 
#             EFT = 'REAL'
#         else:
#             EFT = 'COMPLEX'

#         eff_dm_name = eff_model_dm_names[int(DM['spin'])]
#         mg5_command = 'add model %s %s=%s --recreate --keep_decay' % (pjoin(MDMDIR, 'EffOperators', EFT),
#                                                      eff_dm_name,DM.get('name'))
        
#         # We want to preserve the following variable while updating the model
#         backup_amp = self._curr_amps
#         backup_param_card = self._param_card
#         backup_dm_candidate = self._dm_candidate
#         backup_coannihilation = self._coannihilation
        
#         self.exec_cmd(mg5_command)
#         self._curr_amps = backup_amp
#         self._param_card = backup_param_card 
#         self._dm_candidate = backup_dm_candidate
#         self._coannihilation = backup_coannihilation
        
#         # update the param_card value
#         txt = self._curr_model.write_param_card()
#         param_card = check_param_card.ParamCard(self._curr_model.write_param_card())

#         if self._param_card:
#             for block in self._param_card:
#                 if block not in param_card:
#                     logger.debug('%s not valid block entry in the card' ,block)
#                     continue   
#                 for param in self._param_card[block]:
#                     try:
#                         slhaparam = param_card[block].get(param.lhacode)
#                     except KeyError:
#                         logger.debug('%s %s not valid entry in the card', block,param.lhacode)
# #                        misc.sprint(block,param.lhacode, 'fail')
#                     else:    
#                         slhaparam.value =  param.value
 
#         self._param_card = param_card        
#         if not isinstance(self._curr_model, model_reader.ModelReader):
#             self._curr_model = model_reader.ModelReader(self._curr_model) 
#         self._curr_model.set_parameters_and_couplings(self._param_card) 
