#!/usr/bin/python -W ignore::Warning

import sys,os

from optparse import OptionParser
from caton.core import classify_from_raw_data
import caton.myopts

usage = """
This script is used if you have run cluster_from_raw_data on a subset of your data, \
and you want to extend this sort to the rest of your data.

%prog your_dat_file.dat directory_with_clu_file [options]
%prog -h displays help"""

parser = OptionParser(usage)
parser.add_options([caton.myopts.probe,caton.myopts.max_spikes,caton.myopts.output,caton.myopts.fast,caton.myopts.fast2])


if __name__ == '__main__':

    (opts,args) = parser.parse_args()
    if len(args) == 0:
        raise Exception("Must specify a dat file")
    DatFileName = args[0]
    clu_dir = args[1]
    if not os.path.exists(DatFileName):
        raise Exception("Dat file not found: %s"%DatFileName)
    if not os.path.exists(clu_dir):
        raise Exception("Directory not found: %s"%clu_dir)

    classify_from_raw_data("generalize",DatFileName,opts.probe,max_spikes=opts.max_spikes,output_dir=opts.output,clu_dir=clu_dir)
    
