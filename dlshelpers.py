# -*- coding: utf-8 -*-
# utils.py

import numpy as np
import pandas as pd

def angleToQ(deg, refrac, wavelen):
    return 4.*np.pi*refrac/wavelen*np.sin(deg*np.pi/360.)

# constant for the experiment at each angle
def calcDLSconst(q, temp, visc):
    kB = 1.38066E-23 # Boltzman constant, [kB] = Newton m / Kelvin
    return kB*temp/(6*np.pi*visc)*q**2

def getDLSgammaSi(angle, refrac, wavelen, temp, visc):
    return calcDLSconst(angleToQ(angle, refrac, wavelen), temp, visc)

def parseValue(val):
        val = val.strip('" ')
        try:
            val = float(val)
        except ValueError:
            pass
        return val

def getDLSFileMeta(filenameOrBuffer):
        # read the measurement settings (temperature, viscosity, refractive index and wavelength)
        # only for 'ALV-7004 CGS-8F Data' data files
        meta = pd.read_csv(filenameOrBuffer, r'\s*:\s+', skiprows=1, nrows=36, encoding='cp1250',
                           names=['name', 'value'], index_col='name', engine='python')
        meta = {key: parseValue(value)
                for (key, value) in meta.to_dict()['value'].items()}
        return meta
