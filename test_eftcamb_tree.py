###############################################################################
# import modules:
###############################################################################

import sys, platform, os
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import integrate
import copy

###############################################################################
# import CAMB:
###############################################################################

#here = os.path.dirname(os.path.abspath(__file__))
here = './'
camb_path = os.path.realpath(os.path.join(os.getcwd(),here))
sys.path.insert(0,camb_path)
import camb
camb.set_feedback_level(10)
from camb import model, initialpower
from camb.baseconfig import CAMBError

#Function to check parameters (it is a prototype still ... but it gets all the model with a double loop)

def check_params(eftcamb_params):

    temp = {}
    flag = ['stability','time','level','priors']

    check = False

    modelflag = ['PureEFTmodel','AltParEFTmodel','DesignerEFTmodel','FullMappingEFTmodel']

    while True :

        # I couldn't find a better way than these flags to come out from the loops (may be improved)

        if temp and temp.get('DesignerEFTmodel',0) == 0 :
            for KEY in temp:
                if ( 'model' in KEY and temp[KEY] >= 4):
                        check = True

        if (temp and (temp.get('DesignerEFTmodel',0) != 0 or temp.get('FullMappingEFTmodel',0) != 0)):
            check = True

        if not temp:
            temp = eftcamb_params.copy()

        if check :
            print(temp)
            return temp
            break


        pars = camb.set_params(H0=67.3,
                           EFTCAMB_params=temp,eft_header=False)
        read_par = pars.EFTCAMB.read_parameters()
        #print(temp)

        # This for loop controls the filling of the temp dictionary and assure that only parameter flag are updated
        # in the loop.

        for key in read_par:
            if (flag[0] not in key and flag[1] not in key and flag[2] not in key and flag[3] not in key) :
                    if ('EFTwDE' != key and type(read_par[key]) is not bool and 'EFTflag' != key and 'RPHwDE' != key\
                    and (modelflag[0] != key and modelflag[1] != key and modelflag[2] != key and modelflag[3] != key)):
                        #print(key,read_par[key],read_par[key] == 0 and 'EFTwDE' != key )
                        temp.update({key : read_par[key]+1})
                    else:
                        temp.update({key : read_par[key]})



#Now we just use python error handling to obtain the whole tree

modelflag = ['PureEFTmodel','AltParEFTmodel','DesignerEFTmodel','FullMappingEFTmodel']

for j,i in enumerate(modelflag):
    for xx in range(1,10):
        try :
            params = check_params(eftcamb_params = {'EFTflag':j+1,i:xx})
        except ValueError:
            break

# we need a function to get the used arguments i.e. what EFTCAMB looked for in the dictionary
# define helper that tries to initialize and returns the parameters that are read:
def helper_initialize(pars):
    try:
        temp = camb.set_params(As=2.12605e-9, ns=0.96, H0=67.,
                               ombh2=0.022445, omch2=0.120557, mnu=0.06, tau=0.079,
                               num_massive_neutrinos=1, nnu=3.046,
                               EFTCAMB_params=pars)
        read_params = temp.EFTCAMB.read_parameters()
    except:
        read_params = None
    #
    return read_params

# firt call to have a lay of the land:

print(helper_initialize({'EFTflag':0}))
print(helper_initialize({'EFTflag':1}))


pass
