# -*- coding: utf-8 -*-
# utils.py

import os, re, io, codecs
from warnings import warn
from pathlib import Path
from dateutil.parser import parse as parseDateTime
import zipfile
import bson
import pprint
from uuid import UUID

import numpy as np
import pandas as pd
from jupyter_analysis_tools.utils import isList

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

def bufferFromLines(ln):
    assert isList(ln)
    return io.StringIO("".join(ln))

def readDLSMetaASC(lines, encoding='cp1250'):
    """Reads the measurement settings (temperature, viscosity, refractive index and wavelength).
    The *MeanCR* field is ignored because they were found to be either inconsistent with the
    recorded count rate in the same file (type `ALV-7004/USB`) or they have the wrong ordering
    compared to the recorded count rate in the same file (type `ALV-7004 CGS-8F Data`). The
    means can be calculated easily from the count rate read by *readDLSDataASC()*.
    **For examples**, please see the tests for this function in `tests/test_dlshelpers.py`.
    """
    assert isList(lines)
    meta = {'filetype': lines[0].strip()}
    # section names to look for and their keys in the metadata dict
    sections = [('corrStartLn', '"Correlation"'), ('crStartLn', '"Count Rate"'),
                ('monitorDiodeLn', 'Monitor Diode'), ('cum1StartLn', '"Cumulant 1.Order"'),
                ('cum2StartLn', '"Cumulant 2.Order"'), ('cum3StartLn', '"Cumulant 3.Order"'),
                ('stddevStartLn', '"StandardDeviation"')]
    meta['sections'] = {key: -1 for key, _ in sections}
    key, pattern = sections.pop(0)
    for idx, line in enumerate(lines):
        if pattern in line:
            meta['sections'][key] = idx
            if len(sections):
                key, pattern = sections.pop(0)
            else:
                break
    assert meta['sections']['corrStartLn'] > 0, \
            'The mandatory "Correlation" section was not found!'
    # read the recorded environment data, device settings and sample notes
    df = pd.read_csv(bufferFromLines(lines), sep=':', encoding=encoding,
                       skiprows=1, nrows=meta['sections']['corrStartLn']-2,
                       index_col=0,
                       names=range(5), # also read lines containing more colons = multiple fields
                       na_filter=False, # don't replace empty fields with NaN
                       header=None,
                       #quotechar='"',
                       escapechar="\t", # filters any TAB whitespace
                       skipinitialspace=True # strips whitespace from fields outside of quotes
                      )
    # parse floats and put datetimes back together
    meta.update({name.strip():
                 parseValue(':'.join(field for field in df.T[name] if len(field)).strip())
                 for name in df.index
                 if "MeanCR" not in name
                })
    if meta['sections']['monitorDiodeLn'] > 0:
        meta['monitorDiode'] = float(lines[meta['sections']['monitorDiodeLn']].split()[-1])
    # return the sorted dict
    return meta

def readDLSDataASC(filename, showProgress=False, encoding='cp1250',
                   attenKeys={'key': 'abgeschwächt', 'detectorKey': 'detektor', 'levelKey': 'stufe'}):
    """
    *Mode*: The mode describes which channels (count rate columns) were used for auto- or
    cross-correlation, as given in
    [the specs sheet by ALV](https://www.alvgmbh.de/dwnload/alv7000specs.pdf).
    *corr*: Data with `C-CH0/1+1/0` mode contain the cross-correlation of two channels
            in both directions in the first and second column. Both columns are averaged
            on load and the averaged column is prepended and indexed with the measurement
            angle, to be used for further processing (improves data quality, see issue #1).
    **For examples**, please see the tests for this function in `tests/test_dlshelpers.py`.
    """
    if showProgress:
        print('.', end="") # some progress output
    filename = Path(filename)
    filelines = []
    with codecs.open(filename, encoding=encoding) as fh:
        filelines = [ln for ln in fh.readlines() if len(ln.strip())] # filter empty lines
    data = readDLSMetaASC(filelines)
    data['filename'] = filename
    data['timestamp'] = parseDateTime(data['Date']+' '+data['Time'])

    # parsing the scattering angles
    angles = [value for key,value in data.items() if key.startswith("Angle")]
    data['angles'] = angles

    # try to get the concentration from the memo field, TODO: write a test for this
    data['memo'] = "".join([value for key,value in data.items() if key.startswith("SampMemo")])
    #print(f"memostr: '{data['memo']}'")
    data['memofields'] = [field.strip(' ,;') for field in data['memo'].split()]
    try:
        concentration = [float(field.split(':')[1])
                         for field in data['memofields'] if ':' in field][0]
        concentration = 1/float(concentration)
    except (IndexError, ValueError):
        concentration = 1
    data['concentration'] = concentration

    # parse attenuation values
    atten = [value for key,value in data.items() if key.startswith("Atten")]
    atten += [1.0] * (len(angles) - len(atten)) # pad attenuation list to have same amount as angles (only 3 stored in file)
    if attenKeys['key'] in data['memo']:
        attenManual = getAttenuationFromMemo(data['memo'], len(angles), **attenKeys)
        if atten != attenManual:
            warn("\n"f"  Attenuation values specified in sample description: {attenManual}\n"
                 f"  ('{memostr}'),\n"
                 f"  differ from the values stored by software:          {atten}!\n"
                 f"  -> Using values from sample description:            {attenManual}.", # 'cause they may contain 4 instead of 3 values
                 None, stacklevel=2)
            atten = attenManual
    data['attenuation'] = pd.DataFrame(np.array(atten).reshape(1,-1), columns=angles)

    data = dict(sorted(data.items(), key=lambda x: x[0].lower()))

    # read a section of data columns from header until begin of next section, or file end
    lineIndices = np.array([idx for key, idx in data['sections'].items()], dtype=int)
    def sectionEndLn(startLn): # return the index of any section directly following the given index
        followingSections = lineIndices[lineIndices > startLn]
        return followingSections.min() if followingSections.size else len(filelines)

    # read correlation data columns (mandatory)
    startln, endln = data['sections']['corrStartLn']+1, sectionEndLn(data['sections']['corrStartLn'])
    corr = pd.read_csv(bufferFromLines(filelines), encoding=encoding, sep=r'\s+',
                       index_col=0, header=None,
                       skiprows=startln, nrows=endln-startln)
    if data['Mode'].strip() == "C-CH0/1+1/0": # cross-correlation of two channels in both directions
        corr.insert(0, 0, corr.iloc[:,:2].mean(axis=1))

    # adjust column names if more columns than angles were found
    columns = angles
    if len(corr.columns) > len(columns):
        # repeat the last angle with suffix if there are more columns
        columns = columns + ([columns[-1],] * (len(corr.columns) - len(columns)))
        columns = [f"{int(angle)}_{idx}" if idx else angle for idx, angle in enumerate(columns)]
    corr.rename(columns=dict(zip(corr.columns, columns)), inplace = True)
    corr.index.names = ["tau"]
    data['correlation'] = corr
    
    # read the count rate data columns
    if data['sections']['crStartLn'] > 0: # does not exist in an 'averaged' data file
        startln, endln = data['sections']['crStartLn']+1, sectionEndLn(data['sections']['crStartLn'])
        cr = pd.read_csv(bufferFromLines(filelines), encoding=encoding, sep=r'\s+',
                         #names=["time"]+columns,
                         header=None,
                         index_col=0,
                         skiprows=startln, nrows=endln-startln)
        if data['Mode'].strip() == "C-CH0/1+1/0": # cross-correlation of two channels in both directions
            cr.insert(0, 0, cr.iloc[:,:2].mean(axis=1))
        cr.rename(columns=dict(zip(cr.columns, columns)), inplace = True)
        cr.index.names = ["time"]
        data['countrate'] = cr

    # read cumulants from data
    cumulants = []
    for num in range(3):
        cumulants.append(None)
        secname = f"cum{num+1}StartLn"
        if data['sections'][secname] <= 0:
            continue
        startln, endln = data['sections'][secname]+1, sectionEndLn(data['sections'][secname])
        secname = filelines[startln-1].strip().strip('"')
        cumulants[num] = pd.read_csv(bufferFromLines(filelines), encoding=encoding, sep=r'\t',
                                     names=("value",), index_col=0, engine='python',
                                     skiprows=startln, nrows=endln-startln)
    if any([cumu is not None for cumu in cumulants]):
        data['cumulants'] = cumulants

    # read the standard deviation
    stddevStartLn = data['sections']['stddevStartLn']
    if stddevStartLn > 0:
        startln, endln = stddevStartLn+1, sectionEndLn(stddevStartLn)
        data['stddev'] = pd.read_csv(bufferFromLines(filelines), encoding=encoding, sep=r'\s+',
                                     names=["tau", columns[0]], index_col=0,
                                     skiprows=startln, nrows=endln-startln)
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

def convertAPKWentries(data):
    """Gets a data dictionary loaded from an APKW file and adds or converts field
    values for compatibility with existing code."""
    data['timestamp'] = data['MetaInfo']['StartTime']
    data['Samplename'] = data['InputParameter']['Name']
    data['memo'] = data['InputParameter']['Comment']
    data['memofields'] = [field.strip(' ,;') for field in data['memo'].split()] # same a in getDLSFileData()
    data['concentration'] = 1
    # use the ref. index value of the solvent
    data['Refractive Index'] = data['InputParameter']['Solvent']['RefractiveIndex']
    if data['Refractive Index'] == 0.:
        # find the ref. index in look-up table if not set
        data['Refractive Index'] = [entry['Value']
             for entry in data['InputParameter']['Solvent']['RefractiveIndices']
             if     entry['Key']['Temperature'] == data['InputParameter']['Solvent']['Temperature']
                and entry['Key']['Wavelength'] == data['InputParameter']['Solvent']['Wavelength']][0]
    data['Wavelength [nm]'] = data['InputParameter']['Solvent']['Wavelength']*1e9
    data['Temperature [K]'] = data['InputParameter']['Solvent']['Temperature']
    data['Viscosity [cp]'] = data['InputParameter']['Solvent']['Viscosity']*1e3
    if data['MetaInfo']['InstrumentType'] == 'Litesizer 100':
        data['Angle [°]'] = 175.
        data['angles'] = [data['Angle [°]']]
    data['correlation'] = pd.DataFrame(data['RawData'][0]['CorrelationDataScaled'],
                                       columns=data['angles'],
                                       index=data['RawData'][0]['CorrelationDelayTimes'],)
    data['correlation'] -= 1.
    data['correlation'].index *= 1e3
    data['correlation'].index.names = ["tau"]
    data['countrate'] = pd.DataFrame(data['RawData'][0]['IntensityTrace'], columns=data['angles'])
    del data['RawData'] # not needed anymore, duplicate data
    data['Duration [s]'] = data['InputParameter']['MeasurementTime']
    data['countrate'].index = np.linspace(0, data['Duration [s]'], data['countrate'][175].size)
    data['attenuation'] = pd.DataFrame({data['Angle [°]']: (data['Attenuation'],)})
    return data

def readDLSDataAPKW(filename):
    """Yields none or more unique dicts, one for each measurement,
    stored in the given APKW file name."""
    #print(filename)
    filename = Path(filename)
    if zipfile.is_zipfile(filename):
        with zipfile.ZipFile(filename, 'r') as zf:
            for zi in zf.infolist():
                if zi.filename.startswith('measurement'):
                    #print(zi)
                    data = bson.loads(zf.read(zi))
                    #pprint.pprint(data)
                    if data['StorageStatus'] == 3 and data['UsedAngle'] > 0:
                        data['filename'] = filename / str(data['Id']) # remember it here
                        yield convertAPKWentries(data)

def readDLSData(files, *args, **kwargs):
    """Read DLS measurements from provided files, *.ASC or *.apkw are supported."""
    if not isinstance(files, list) and not isinstance(files, tuple):
        files = (files,)
    # gather all measurement data from APKW files first
    dirData = list({datadict['Id']: datadict
                    for fn in files
                    for datadict in readDLSDataAPKW(fn)}.values())
    # try to read ASC files
    if not len(dirData):
        dirData = [readDLSDataASC(fn) for fn in files]
    return dirData

if __name__ == "__main__":
    import doctest
    doctest.testmod()
