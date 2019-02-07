#! /usr/bin/env python

import numpy  as np
import scipy.integrate   as     integrate
import math


# Include the python wrapper for the cpp functions
import CRXS.XS_wrapper as XS

#
#   Some examples how to call functions from the wrapper:
#

print('Evaluate: "XS.dE_AA_pbar_LAB_incNbarAndHyperon(  100        , 3           )"  automaticaly assumes pp scattering and KORSMEIER_II parametrization')
print(            XS.dE_AA_pbar_LAB_incNbarAndHyperon(  100        , 3           ) )


print('Evaluate: "XS.dE_AA_p_LAB                     (  100        , 3           )"  automaticaly assumes pp scattering and ANDERSON parametrization')
print(            XS.dE_AA_p_LAB                     (  100        , 3           ) )


# You may also vectorize the function:
_dE_AA_pbar_LAB_incNbarAndHyperon = np.vectorize(XS.dE_AA_pbar_LAB_incNbarAndHyperon)

print('Evaluate: "  _dE_AA_pbar_LAB_incNbarAndHyperon( [100., 10.]   , [30., 3.] )" ')
print(              _dE_AA_pbar_LAB_incNbarAndHyperon( [100., 10.]   , [30., 3.] ) )
