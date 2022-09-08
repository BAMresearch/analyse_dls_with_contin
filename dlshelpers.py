# -*- coding: utf-8 -*-
# utils.py

import os, re
from warnings import warn
from pathlib import Path
from dateutil.parser import parse as parseDateTime
import numpy as np
import pandas as pd

# transmission levels of ALV-CGS4F/8F given in measurement software, angles/detectors 1-4
transmissionLevels = np.array(((1., .1, .03, .01),)+((1., .3, .1, .03),)*3)

def angleToQ(deg, refrac, wavelen):
    return 4.*np.pi*refrac/wavelen*np.sin(deg*np.pi/360.)

# constant for the experiment at each angle
def calcDLSconst(q, temp, visc):
    kB = 1.38066E-23 # Boltzman constant, [kB] = Newton m / Kelvin
    return kB*temp/(6*np.pi*visc)*q**2

def getDLSgammaSi(angle, refrac, wavelen, temp, visc):
    return calcDLSconst(angleToQ(angle, refrac, wavelen), temp, visc)

def parseValue(val):
    val = val.strip('"')
    try:
        val = float(val)
    except ValueError:
        pass
    return val

def getDLSFileMeta(filenameOrBuffer, encoding='cp1250'):
    # read the measurement settings (temperature, viscosity, refractive index and wavelength)
    # only for 'ALV-7004 CGS-8F Data' data files
    meta = pd.read_csv(filenameOrBuffer, sep=':', skiprows=1, nrows=39, encoding=encoding,
                       index_col=0,
                       names=range(5), # also read fields containing two colons
                       na_filter=False, # don't replace empty fields with NaN
                       header=None,
                       #quotechar='"',
                       escapechar="\t", # filters any TAB whitespace
                       skipinitialspace=True # strips whitespace from fields outside of quotes
                      )
    meta = {name.strip(): parseValue(':'.join(field for field in meta.T[name]
                                              if len(field)).strip())
            for name in meta.index}
    return meta

def getDLSFileData(filename, showProgress=False, encoding='cp1250',
                   attenKeys={'key': 'abgeschwächt', 'detectorKey': 'detektor', 'levelKey': 'stufe'}):
    if showProgress:
        print('.', end="") # some progress output
    data = dict(filename=Path(filename).resolve())
    header = getDLSFileMeta(str(filename))
    data.update(sampleName=header["Samplename"])
    data.update(timestamp=parseDateTime(header['Date']+' '+header['Time']))
    # parsing the scattering angles
    angles = [value for key,value in header.items() if key.startswith("Angle")]
    data.update(angles=angles)
    # try to get the concentration from the memo field
    memostr = "".join([value for key,value in header.items() if key.startswith("SampMemo")])
    #print(f"memostr: '{memostr}'")
    memofields = [field.strip(' ,;') for field in memostr.split()]
    try:
        concentration = [float(field.split(':')[1]) for field in memofields if ':' in field][0]
        concentration = 1/float(concentration)
    except (IndexError, ValueError):
        concentration = 1
    data.update(concentration=concentration)
    # parsing attenuation values
    atten = [value for key,value in header.items() if key.startswith("Atten")]
    atten += [1.0] * (len(angles) - len(atten)) # pad attenuation list to have same amount as angles (only 3 stored in file)
    if attenKeys['key'] in memostr:
        attenManual = getAttenuationFromMemo(memostr, len(angles), **attenKeys)
        if atten != attenManual:
            warn("\n"f"  Attenuation values specified in sample description: {attenManual}\n"
                 f"  ('{memostr}'),\n"
                 f"  differ from the values stored by software:          {atten}!\n"
                 f"  -> Using values from sample description:            {attenManual}.", # 'cause they may contain 4 instead of 3 values
                 None, stacklevel=2)
            atten = attenManual
    data.update(attenuation=pd.DataFrame(np.array(atten).reshape(1,-1), columns=angles))
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

def getAttenuationFromMemo(memo, count=4, detectorKey='detektor', levelKey='stufe', key=None):
    """Extract manually specified attenuation of individual scattering angles from a given
    measurement data files (ASC) *SampMemo* field.
    Different notations and numbers of detectors are accepted (test cases):

    >>> getAttenuationFromMemo("Toluolmessung für Statik,, Detektor 1 und 2, 1 Stufe abgeschwächt")
    [0.1, 0.3, 1.0, 1.0]

    >>> getAttenuationFromMemo("Toluolmessung für Statik,, Detektor 1 und 2, 1 Stufe abgeschwächt, Detektor 3 , 2 Stufen abgeschwächt")
    [0.1, 0.3, 0.1, 1.0]

    >>> getAttenuationFromMemo("Toluolmessung für Statik,, Detektor 1 und 2, 1 Stufe abgeschwächt, Detektor 3 und 4, 2 Stufen abgeschwächt")
    [0.1, 0.3, 0.1, 0.1]

    >>> getAttenuationFromMemo("Toluolmessung für Statik,, Detektor 1, 1 Stufe abgeschwächt, Detektor 3 , 2 Stufen abgeschwächt")
    [0.1, 1.0, 0.1, 1.0]

    >>> getAttenuationFromMemo("1:100 verdünnt mit Wasser, ungefiltert,, Detektor 1(10%) und 2(30%), 1 Stufe abgeschwächt")
    [0.1, 0.3, 1.0, 1.0]

    Table of transmission levels 1-4 (columns) as defined in the program for detectors 1-4 (rows):

    >>> print(transmissionLevels)
    [[1.   0.1  0.03 0.01]
     [1.   0.3  0.1  0.03]
     [1.   0.3  0.1  0.03]
     [1.   0.3  0.1  0.03]]
    """
    def getnumbers(text, idxoffset=0):
        return ([(int(idx)+idxoffset, int(inner)/100 if len(inner) else 1.)
                for idx, _, inner in re.findall(r"\b(\d)\b(\((\d+)%\))?", text)])
    atten, detector = list(np.ones(count)), ()
    for field in memo.split(','):
        field = field.strip().lower()
        #print(f"-> field '{field}':")
        # extract detectors indices first, may contain transmission % already
        if detectorKey in field:
            detector = getnumbers(field, idxoffset=-1)
            #print(detector)
        # extract transmission level next, level indices map to hard coded transmission factor table
        if levelKey in field:
            level = getnumbers(field)
            detector = [(d, transmissionLevels[d,level[d%len(level)][0]]) for d, _ in detector]
            #print(detector)
        for d, val in detector:
            atten[d] = val
    return atten

if __name__ == "__main__":
    import doctest
    doctest.testmod()
