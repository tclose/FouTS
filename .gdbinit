python
import sys

#Set pretty printers system paths.
sys.path.insert(0, '/home/tclose/Code/python/gdb_printers/python')
sys.path.insert(0, '/home/tclose/git/BaFTrS/gdb_printers/python')
sys.path.insert(0, '/home/tclose/git/MRtrix/gdb_printers/python')

#Import and register 'std::*' pretty printers.
from libstdcxx.v6.printers import register_libstdcxx_printers
register_libstdcxx_printers (None)

#Import and register 'BTS::Fibre::*' pretty printers.
from bts.fibre.printers import register_bts_fibre_printers
register_bts_fibre_printers(None)

#Import and register 'BTS::Image::*' pretty printers.
from bts.image.printers import register_bts_image_printers
register_bts_image_printers(None)

#Import and register 'BTS::MCMC::*' pretty printers.
from bts.mcmc.printers import register_bts_mcmc_printers
register_bts_mcmc_printers(None)

#Import MR::Math* pretty printers.
from mr.math.printers import register_mr_math_printers
register_mr_math_printers(None)

