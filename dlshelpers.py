# -*- coding: utf-8 -*-
# utils.py

import os
from pathlib import Path
from dateutil.parser import parse as parseDateTime
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

def getDLSFileMeta(filenameOrBuffer, encoding='cp1250'):
        # read the measurement settings (temperature, viscosity, refractive index and wavelength)
        # only for 'ALV-7004 CGS-8F Data' data files
        meta = pd.read_csv(filenameOrBuffer, sep=r'\s*:\s+', skiprows=1, nrows=36, encoding=encoding,
                           names=['name', 'value'], index_col='name', engine='python')
        meta = {key: parseValue(value)
                for (key, value) in meta.to_dict()['value'].items()}
        return meta

def getDLSFileData(filename, showProgress=False, encoding='cp1250'):
    if showProgress:
        print('.', end="") # some progress output
    data = dict(filename=Path(filename).resolve())
    header = getDLSFileMeta(str(filename))
    data.update(sampleName=header["Samplename"])
    data.update(timestamp=parseDateTime(header['Date']+' '+header['Time']))
    memostr = "".join([value for key,value in header.items() if key.startswith("SampMemo")])
    # try to get the concentration from the memo field
    #print("memostr", memostr)
    memofields = [field.strip(' ,;') for field in memostr.split()]
    #print("memofields", memofields)
    try:
        concentration = [float(field.split(':')[1]) for field in memofields if ':' in field][0]
        concentration = 1/float(concentration)
    except (IndexError, ValueError):
        concentration = 1
    data.update(concentration=concentration)
    angles = [value for key,value in header.items() if key.startswith("Angle")]
    data.update(angles=angles)
    # ATTENTION! MeanCR is stored in reversed order in these files for
    # the acquisition software used here: "ALV MultiAngle" Version 3.1.4.3
    meancr = np.flip(np.array([value for key,value in header.items() if key.startswith("MeanCR")]))
    data.update(meancr=pd.DataFrame(meancr.reshape(1,-1), columns=angles))
    for name in "Temperature", "Viscosity", "Refractive Index", "Wavelength":
        for key,value in header.items():
            if key.startswith(name):
                data[key] = value
    with open(filename, encoding='cp1250') as fd:
        lines = fd.readlines()
        #[print(ln.strip()) for ln in lines[17:25]]
        value = [line for line in lines if line.startswith('Monitor Diode')]
        if len(value): # empty for '_averaged' files
            value = float(value[0].split()[-1])
            data.update(monitorDiode=value)
        #print("data", data)
        def findMatchingLine(lineLst, pattern):
            idx = [idx for idx, ln in enumerate(lineLst) if pattern in ln]
            return int(idx[0]) if len(idx) else len(lineLst)-1
        corrstart = findMatchingLine(lines, 'Correlation')
        crstart   = findMatchingLine(lines, 'Count Rate')
        crend     = findMatchingLine(lines, 'Monitor Diode')
    # read tabular data after file was closed by leaving scope of 'with'
    if corrstart > 0:
        corr = pd.read_csv(filename, skiprows=corrstart+1, nrows=crstart-1-corrstart,
                           sep=r'\s+', names=["tau"]+angles, index_col=0, encoding=encoding)
        data.update(correlation=corr)
    if (crend-2-crstart) > 0:
        cr   = pd.read_csv(filename, skiprows=crstart+1, nrows=crend-2-crstart,
                           sep=r'\s+', names=["time"]+angles, index_col=0, encoding=encoding)
        data.update(countrate=cr)
    return data
