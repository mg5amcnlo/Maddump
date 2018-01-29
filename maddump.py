from distutils import dir_util
import os
from madgraph.iolibs.files import cp, ln, mv

import madgraph.iolibs.export_v4 as export_v4
import madgraph.iolibs.file_writers as writers
import madgraph.various.misc as misc
import madgraph.various.banner as banner_mod
import madgraph.iolibs.group_subprocs as group_subprocs

pjoin = os.path.join

class MadDump(export_v4.ProcessExporterFortranMEGroup):
    
    check = True
    # check status of the directory. Remove it if already exists
    # Language type: 'v4' for f77/ 'cpp' for C++ output
    exporter = 'v4'
    # Output type:
    #[Template/dir/None] copy the Template, just create dir  or do nothing 
    output = 'Template' 
    # Decide which type of merging to used [madevent/madweight]
    grouped_mode = 'madevent'    
    # if no grouping on can decide to merge uu~ and u~u anyway:
    sa_symmetry = False 

    def __init__(self, *args, **opts):
        misc.sprint("Initialise MadDump")
        return super(MadDump, self).__init__(*args, **opts)

    def copy_template(self, *args, **opts):
        
        maddump_dir = pjoin(self.mgme_dir, 'PLUGIN/maddump/')
        misc.sprint("copy the associate template")
        super(MadDump, self).copy_template(*args, **opts)

        dir_util.copy_tree(pjoin(maddump_dir, 'FITPACK'), 
                           pjoin(self.dir_path, 'Source/FITPACK'))
        
        #cp(maddump_dir + 'dummy_fct.f',
        #       self.dir_path + '/SubProcesses/dummy_fct.f')
        #cp(maddump_dir + 'unwgt.f',
        #       self.dir_path + '/SubProcesses/unwgt.f')

        cp(maddump_dir + 'meshfitter2D.py',
               self.dir_path + '/Cards/meshfitter2D.py')
        cp(maddump_dir + 'lhe-meshfitter.py',
               self.dir_path + '/Cards/lhe-meshfitter.py')
        cp(maddump_dir + 'README_MADDUMP',
               self.dir_path + '/Cards/README_MADDUMP')
        
        with open(self.dir_path+'/SubProcesses/makefile') as f:
            text = []
            for line in f:
                if 'LINKLIBS' in line:
                    text.append(line)
                    text.append('LINKLIBS += -lfitpack\n')
                else:
                    text.append(line)
                    
        with open(self.dir_path+'/SubProcesses/makefile', "w") as f:
            for line in text:
                f.write(line)
        
        files = ["dummy_fct.f","unwgt.f"]
        remove_list = [['get_dummy_x1','get_dummy_x1_x2'],["write_leshouche"]]
        for name, to_rm in zip(files, remove_list):
            template = open(pjoin(self.dir_path, "SubProcesses", name),"r").read()
            plugin = open(pjoin(self.mgme_dir, "PLUGIN", "maddump", name),"r").read()
            misc.sprint(pjoin(self.mgme_dir, "PLUGIN", "maddump", name))
            ff = writers.FortranWriter(pjoin(self.dir_path, "SubProcesses", name))
            ff.remove_routine(template, to_rm, formatting=False)
            ff.writelines(plugin, formatting=False)
            ff.close()


    # def get_source_libraries_list(self):
    #     set_of_lib = super(MadDump, self).get_source_libraries_list()
    #     set_of_lib.append('$(LIBDIR)libfitpack.$(libext)')
    #     return set_of_lib


    def write_source_makefile(self, writer):
        """Write the nexternal.inc file for MG4"""
        
        replace_dict = super(MadDump, self).write_source_makefile(None)

        libfitpack_line='''$(LIBDIR)libfitpack.$(libext):\n\tcd FITPACK; make\n\nlifitpack: $(LIBDIR)libfitpack.$(libext)\n'''


        replace_dict['libraries'] += ' $(LIBDIR)libfitpack.$(libext)'
        replace_dict['additional_dependencies'] += libfitpack_line
        
        if writer:
            path = pjoin(self.mgme_dir, 'madgraph', 'iolibs','template_files','madevent_makefile_source')
            text = open(path).read() % replace_dict
            writer.write(text)
            writer.write('\tcd FITPACK; make clean; cd ..\n')
            
        return replace_dict
    
    def pass_information_from_cmd(self, cmd):
        """pass information from the command interface to the exporter.
           Please do not modify any object of the interface from the exporter.
        """
        self.cmd = cmd
        return super(MadDump, self).pass_information_from_cmd(cmd)
        
    def finalize(self,*args, **opts):
        
        from madgraph import MG5DIR    
        filename = os.path.join(self.cmd._export_dir, 'Cards', 'me5_configuration.txt')
        self.cmd.do_save('options %s' % filename.replace(' ', '\ '), check=False,
                         to_keep={'mg5_path':MG5DIR})
        
        return super(MadDump, self).finalize(*args, **opts)

    #===========================================================================
    #  create the run_card 
    #===========================================================================
    def create_run_card(self, matrix_elements, history):
        """ Create the run card under Cards directory """
        run_card = banner_mod.RunCard()
        run_card.remove_all_cut()
        run_card['lpp1'] = 9

        # to do: distinguish the case between a collision with a proton or an electron in the detector(it needs to be fixed, lpp2 is not update automatically in the two different cases)
        # if run_card['lpp2']  == 1:
        #     run_card['ebeam2'] = 0.938
        # else:
        #     run_card['ebeam2'] = 0.000511

        run_card.write(pjoin(self.dir_path, 'Cards', 'run_card_default.dat'))
        run_card.write(pjoin(self.dir_path, 'Cards', 'run_card.dat'))
        

    def link_files_in_SubProcess(self, Ppath):
        """ Create the necessary links in the P* directory path Ppath"""

        #import ehist.dat and cell_frotran.dat into Subprocesses/P
        ln(self.dir_path + '/Cards/ehist.dat', cwd=Ppath)
        ln(self.dir_path + '/Cards/cell_fortran.dat', cwd=Ppath)

        super(MadDump, self).link_files_in_SubProcess(Ppath)
    
