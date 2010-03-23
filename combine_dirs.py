#!/usr/bin/python -W ignore::Warning

from caton.core import combine_h5s
from optparse import OptionParser
import caton.myopts


usage = """%prog [dir1 [dir2 ...]]"""
parser = OptionParser(usage)
parser.add_option(caton.myopts.output)

(opts,args) = parser.parse_args()
combine_h5s(*args)
