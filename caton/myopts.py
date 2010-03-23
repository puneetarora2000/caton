from optparse import OptionParser,Option

def set_fast(*args):
        import caton.core
        caton.core.extract_wave = caton.core.extract_wave_simple

def set_fast2(*args):
        import caton.subset_sorting
        caton.subset_sorting.MINCLUSTERS = 5
        caton.subset_sorting.MAXCLUSTERS = 5

probe = Option("-p","--probe",action="store",dest="probe",
                  help="Specify the location of the probe file. By default, program looks for a .probe file in .dat directory")
max_spikes = Option("-n",type=int,dest="max_spikes",
                  help="Extract and cluster only the first n spikes.")
output = Option("-o","--output",action="store",dest="output",
                help="""Directory where the output directory 'basename/' ends up. By default, output goes next to .dat file.""")                  
fast = Option("--fast",action="callback",callback=set_fast,
                  help="Makes pre-clustering go much faster by skipping slowest step: interpolation around the inferred peak.")
fast2 = Option("--fast2",action="callback",callback=set_fast2,
                  help="Makes clustering go fast by doing EM algorithm only once with 5 clusters, instead of all numbers between 3 and 14")


