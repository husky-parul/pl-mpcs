#!/usr/bin/env python                                            
#                                                            _
# mpcs ds app
#
# (c) 2016-2019 Fetal-Neonatal Neuroimaging & Developmental Science Center
#                   Boston Children's Hospital
#
#              http://childrenshospital.org/FNNDSC/
#                        dev@babyMRI.org
#

import  sys
import  os
import  random
import  copy
import  shutil
import  csv
import  numpy               as      np

# import the Chris app superclass
from    chrisapp.base       import ChrisApp

Gstr_title = """
 _ __ ___  _ __   ___ ___ 
| '_ ` _ \| '_ \ / __/ __|
| | | | | | |_) | (__\__ \\
|_| |_| |_| .__/ \___|___/
          | |             
          |_|             
"""

Gstr_synopsis = """

    NAME

       mpcs.py 

    SYNOPSIS

        python mpcs.py                                                  \\
            [--random] [--seed <seed>]                                  \\
            [-p <f_posRange>] [--posRange <f_posRange>]                 \\
            [-n <f_negRange>] [--negRange <f_negRange>]                 \\
            [-z <zFile>] [--zFile <zFile>]                              \\
            [-v <level>] [--verbosity <level>]                          \\
            [--version]                                                 \\
            [--man]                                                     \\
            [--meta]                                                    \\
            <inputDir>                                                  \\
            <outputDir> 

    BRIEF EXAMPLE

        * Bare bones execution

            mkdir in out && chmod 777 out
            python mpcs.py  in    out

    DESCRIPTION

        `mpcs.py` simulates a call to a remote Multi-Party Compute (MPC)
        in the context of a FreeSurfer workflow.

        This particular application simply returns a z-score file to be
        consumed by a downstream plugin, typciall ``pl-z2labelmap``.

        NOTE: The <inputDir> is largely ignored by this plugin.

    ARGS

        [--random] [--seed <seed>]
        If specified, generate a z-score file based on <posRange> and 
        <negRange>. In addition, if a further optional <seed> is passed,
        then initialize the random generator with that seed, otherwise
        system time is used.

        [-p <f_posRange>] [--posRange <f_posRange>]
        Positive range for random max deviation generation.

        [-n <f_negRange>] [--negRange <f_negRange>]
        Negative range for random max deviation generation.

        [-z <zFile>] [--zFile <zFile>]
        z-score file to save in output directory. Defaults to 'zfile.csv'.

        [-v <level>] [--verbosity <level>]
        Verbosity level for app. Not used currently.

        [--version]
        If specified, print version number. 
        
        [--man]
        If specified, print (this) man page.

        [--meta]
        If specified, print plugin meta data.

    EXAMPLES

    Create a z-file with values between -3.0 and +3.0

        mpcs.py     --random --seed 1                   \\
                    --posRange 3.0 --negRange -3.0      \\
                    in out

"""

class Mpcs(ChrisApp):
    """
    This app simulates an MPC compute call..
    """
    AUTHORS                 = 'Rudolph Pienaar (rudolph.pienaar@gmail.com)'
    SELFPATH                = os.path.dirname(os.path.abspath(__file__))
    SELFEXEC                = os.path.basename(__file__)
    EXECSHELL               = 'python3'
    TITLE                   = 'A simple exemplar multi-party-compute (MPC) simulator plugin.'
    CATEGORY                = 'MPC'
    TYPE                    = 'ds'
    DESCRIPTION             = 'This app simulates an MPC compute call.'
    DOCUMENTATION           = 'https://github.com/FNNDSC/pl-mpcs'
    VERSION                 = '1.0.3'
    ICON                    = '' # url of an icon image
    LICENSE                 = 'Opensource (MIT)'
    MAX_NUMBER_OF_WORKERS   = 1  # Override with integer value
    MIN_NUMBER_OF_WORKERS   = 1  # Override with integer value
    MAX_CPU_LIMIT           = '' # Override with millicore value as string, e.g. '2000m'
    MIN_CPU_LIMIT           = '' # Override with millicore value as string, e.g. '2000m'
    MAX_MEMORY_LIMIT        = '' # Override with string, e.g. '1Gi', '2000Mi'
    MIN_MEMORY_LIMIT        = '' # Override with string, e.g. '1Gi', '2000Mi'
    MIN_GPU_LIMIT           = 0  # Override with the minimum number of GPUs, as an integer, for your plugin
    MAX_GPU_LIMIT           = 0  # Override with the maximum number of GPUs, as an integer, for your plugin

    # Use this dictionary structure to provide key-value output descriptive information
    # that may be useful for the next downstream plugin. For example:
    #
    # {
    #   "finalOutputFile":  "final/file.out",
    #   "viewer":           "genericTextViewer",
    # }
    #
    # The above dictinary is saved when plugin is called with a ``--saveoutputmeta`` 
    # flag. Note also that all file paths are relative to the system specified
    # output directory.
    OUTPUT_META_DICT = {}
 
    def manPage_show(self):
        """
        Print some quick help.
        """
        print(Gstr_synopsis)

    def metaData_show(self):
        """
        Print the plugin meta data
        """
        l_metaData  = dir(self)
        l_classVar  = [x for x in l_metaData if x.isupper() ]
        for str_var in l_classVar:
            str_val = getattr(self, str_var)
            print("%20s: %s" % (str_var, str_val))

    def a2009sStructList_define(self):
        """
        The list of structures in the a2009s cortical parcellation
        """

        self.l_a2009s = [
            'G_and_S_frontomargin',
            'G_and_S_occipital_inf',
            'G_and_S_paracentral',
            'G_and_S_subcentral',
            'G_and_S_transv_frontopol',
            'G_and_S_cingul-Ant',
            'G_and_S_cingul-Mid-Ant',
            'G_and_S_cingul-Mid-Post',
            'G_cingul-Post-dorsal',
            'G_cingul-Post-ventral',
            'G_cuneus',
            'G_front_inf-Opercular',
            'G_front_inf-Orbital',
            'G_front_inf-Triangul',
            'G_front_middle',
            'G_front_sup',
            'G_Ins_lg_and_S_cent_ins',
            'G_insular_short',
            'G_occipital_middle',
            'G_occipital_sup',
            'G_oc-temp_lat-fusifor',
            'G_oc-temp_med-Lingual',
            'G_oc-temp_med-Parahip',
            'G_orbital',
            'G_pariet_inf-Angular',
            'G_pariet_inf-Supramar',
            'G_parietal_sup',
            'G_postcentral',
            'G_precentral',
            'G_precuneus',
            'G_rectus',
            'G_subcallosal',
            'G_temp_sup-G_T_transv',
            'G_temp_sup-Lateral',
            'G_temp_sup-Plan_polar',
            'G_temp_sup-Plan_tempo',
            'G_temporal_inf',
            'G_temporal_middle',
            'Lat_Fis-ant-Horizont',
            'Lat_Fis-ant-Vertical',
            'Lat_Fis-post',
            'Pole_occipital',
            'Pole_temporal',
            'S_calcarine',
            'S_central',
            'S_cingul-Marginalis',
            'S_circular_insula_ant',
            'S_circular_insula_inf',
            'S_circular_insula_sup',
            'S_collat_transv_ant',
            'S_collat_transv_post',
            'S_front_inf',
            'S_front_middle',
            'S_front_sup',
            'S_interm_prim-Jensen',
            'S_intrapariet_and_P_trans',
            'S_oc_middle_and_Lunatus',
            'S_oc_sup_and_transversal',
            'S_occipital_ant',
            'S_oc-temp_lat',
            'S_oc-temp_med_and_Lingual',
            'S_orbital_lateral',
            'S_orbital_med-olfact',
            'S_orbital-H_Shaped',
            'S_parieto_occipital',
            'S_pericallosal',
            'S_postcentral',
            'S_precentral-inf-part',
            'S_precentral-sup-part',
            'S_suborbital',
            'S_subparietal',
            'S_temporal_inf',
            'S_temporal_sup',
            'S_temporal_transverse'
        ]
        maxlen          = len(max(self.l_a2009s, key = len))
        # Right pad the structure names with spaces to set all names to same length.
        self.l_a2009s   = [x + (' ' * (maxlen - len(x))) for x in self.l_a2009s]
        return self.l_a2009s

    def randomZscoreFile_generate(self, astr_parcellation):
        """
        Generate a "random" z-score file, based on the range given in the
        --random  --posRange and --negRange  command line flags. 
        
        This file has three columns,

            <structName> <leftHemisphere-z-score> <rightHemisphere-z-score>

        Save file to both input and output directories.

        """

        def file_write(f):
            writer  = csv.writer(f)
            for row in self.rows_zscore:
                writer.writerow(row) 

        l_parc  = self.d_parcellation[astr_parcellation]['structureNames']

        self.d_parcellation[astr_parcellation]['lh']['zScore'] =   np.random.uniform(  
                                            low     = self.options.f_negRange, 
                                            high    = self.options.f_posRange, 
                                            size    = (len(l_parc,))
                                ).tolist()
        self.d_parcellation[astr_parcellation]['rh']['zScore'] =   np.random.uniform(  
                                            low     = self.options.f_negRange, 
                                            high    = self.options.f_posRange, 
                                            size    = (len(l_parc,))
                                ).tolist()

        self.rows_zscore = zip( 
                        self.d_parcellation[astr_parcellation]['structureNames'], 
                        self.d_parcellation[astr_parcellation]['lh']['zScore'],
                        self.d_parcellation[astr_parcellation]['rh']['zScore']
                        )
        with open('%s/%s' % (self.options.outputdir, self.options.zFile), 
                    "w", newline = '') as f:
            file_write(f)

        return {
            'status':   True
        }

    def define_parameters(self):
        """
        Define the CLI arguments accepted by this plugin app.
        """

        self.add_argument("-p", "--posRange",
                            help        = "positive range for random max deviation generation",
                            type        = float,
                            dest        = 'f_posRange',
                            optional    = True,
                            default     = 1.0)
        self.add_argument("-n", "--negRange",
                            help        = "negative range for random max deviation generation",
                            type        = float,
                            dest        = 'f_negRange',
                            optional    = True,
                            default     = -1.0)
        self.add_argument("-z", "--zFile",
                            help        = "z-score file to read (relative to input directory)",
                            type        = str,
                            dest        = 'zFile',
                            optional    = True,
                            default     = 'zfile.csv')
        self.add_argument('--random',
                            help        = 'if specified, generate a z-score file',
                            type        = bool,
                            dest        = 'b_random',
                            action      = 'store_true',
                            optional    = True,
                            default     = False)
        self.add_argument("-d", "--seed",
                            help        = "random number seed",
                            type        = str,
                            dest        = 'seed',
                            optional    = True,
                            default     = '')
        self.add_argument("-v", "--verbosity",
                            help        = "verbosity level for app",
                            type        = str,
                            dest        = 'verbosity',
                            optional    = True,
                            default     = "0")
        self.add_argument('--man',
                            help        = 'if specified, print man page',
                            type        = bool,
                            dest        = 'b_man',
                            action      = 'store_true',
                            optional    = True,
                            default     = False)
        self.add_argument('--meta',
                            help        = 'if specified, print plugin meta data',
                            type        = bool,
                            dest        = 'b_meta',
                            action      = 'store_true',
                            optional    = True,
                            default     = False)
        self.add_argument('--version',
                            help        = 'if specified, print version number',
                            type        = bool,
                            dest        = 'b_version',
                            action      = 'store_true',
                            optional    = True,
                            default     = False)

    def run(self, options):
        """
        Define the code to be run by this plugin app.
        """

        b_zFileProcessed    = False

        self.options        = options
        self.d_hemiStats    = {
            'zScore':       [],
            'stats':        {'min': 0.0, 'max': 0.0},
            'posNorm':      [],
            'negNorm':      []
        }
        self.d_core         = {
            'structureNames':   [],
            'lh':               copy.deepcopy(self.d_hemiStats),
            'rh':               copy.deepcopy(self.d_hemiStats),
            'zScoreFile':       "",
            'labelMapFile':     "",
            'daM_color':        None
        }
        self.d_parcellation = {
            'a2009s':   copy.deepcopy(self.d_core),
            'DKatlas':  copy.deepcopy(self.d_core),
            'default':  copy.deepcopy(self.d_core)
        }

        if options.b_man:
            self.manPage_show()
            sys.exit(0)

        if options.b_meta:
            self.metaData_show()
            sys.exit(0)

        if options.b_version:
            print('Plugin Version: %s' % Mpcs.VERSION)
            sys.exit(0)

        print(Gstr_title)
        print('Version: %s' % Mpcs.VERSION)

        self.d_parcellation['a2009s']['structureNames']     = \
                 self.a2009sStructList_define()
        self.d_parcellation['a2009s']['zScoreFile']         = '%s/%s' % \
                (self.options.inputdir, self.options.zFile)
        self.d_parcellation['a2009s']['labelMapFile']       = '%s/%s' % \
                (self.options.outputdir, 'aparc.annot.a2009s.ctab')

        if len(options.seed):
            # Use user specified seed
            random.seed(options.seed)
        else:
            # else use system time
            random.seed()
        print("Generating %s/%s..." % (options.outputdir, options.zFile))
        self.randomZscoreFile_generate('a2009s')
        b_zFileProcessed    = True


# ENTRYPOINT
if __name__ == "__main__":
    app = Mpcs()
    app.launch()
