from __future__ import with_statement
import os,sys,numpy

# Give scripts the correct first line.
scripts = ["check_crosstalk.py","cluster_from_raw_data.py","extract_intra_spikes.py","generalize_from_raw_data.py",
           "make_features_from_spk.py","plot_features.py","plot_probe.py","combine_dirs.py"]

from distutils.core import setup,Extension
from Cython.Distutils import build_ext

    
setup(name="caton",
      scripts=scripts,
      version="0.2.1-beta",
      author="John Schulman",
      author_email="joschu@caltech.edu",
      description="Spike sorting for multi-site probes",
      license="GPL3",
      url="http://caton.googlecode.com",
      packages=["caton"],
      package_data={"caton":["data/*.txt"]},
      cmdclass = {'build_ext': build_ext},      
      include_dirs = [numpy.get_include(),'.'],
      ext_modules = [Extension("caton.CEM_extensions", sources=["caton/CEM_extensions.c"])])
