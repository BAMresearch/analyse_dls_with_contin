# -*- coding: utf-8 -*-
# utils.py

from dlshelpers import getDLSgammaSi, getDLSFileMeta

import os, shutil, subprocess
import time
import urllib
import requests
from pathlib import Path
from bs4 import BeautifulSoup
import pandas as pd
import numpy as np
from jupyter_analysis_tools.utils import isWindows, isMac, pushd

def getContinForWindows(targetPath):
    binaryName = "contin-windows.exe"
    baseurl="http://www.s-provencher.com"
    html_page = urllib.request.urlopen(baseurl+"/contin.shtml").read()
    soup = BeautifulSoup(html_page)
    binurl = [link.get('href') for link in soup.findAll('a') if link.text.strip() == binaryName]
    binurl = '/'.join((baseurl,binurl[0]))
    binary = requests.get(binurl, allow_redirects=True)
    open(targetPath, 'wb').write(binary.content)
    return targetPath

def getContinPath():
    if isMac():
        # get local path to the CONTIN executable
        return Path.home() / "code" / "cntb2" / "bin" / "contin OSX"
    elif isWindows():
        continCmd = Path(os.getenv('APPDATA')) / "contin" / "contin.exe"
        #print(continCmd)
        if continCmd.is_file():
            return continCmd
        continCmd.parent.mkdir(parents=True, exist_ok=True)
        getContinForWindows(continCmd)
        if continCmd.is_file():
            print(f"Installed CONTIN at '{continCmd}'.")
        else:
            print("Failed to retrieve CONTIN executable!")
    raise NotImplementedError("Don't know how to retrieve the CONTIN executable!")

def genContinInput(filename, **continConfig):
    IWT = 5 if continConfig['weighResiduals'] else 1
    # transform data? Trd=0: no transform
    Trd = -1 # Trd=1: initial g(2), input sqrt[g(2)-1]; Trd=-1: initial g(2)-1, input sqrt[g(2)-1]
    # read the measurement settings (temperature, viscosity, refractive index and wavelength)
    meta = getDLSFileMeta(filename)
    # select the measurement angle, make sure it's in the file
    angle = continConfig['angle']
    anglesFromFile = [meta[key] for key in meta if key.startswith("Angle")]
    assert angle in anglesFromFile, \
        f"Given angle ({angle}) not found in file '{filename}': {anglesFromFile}"
    # get environment values for storage in contin file
    temp    = meta['Temperature [K]']
    visc    = meta['Viscosity [cp]']
    refrac  = meta['Refractive Index']
    wavelen = meta['Wavelength [nm]']
    gamma = getDLSgammaSi(angle, refrac, wavelen*1e-9, temp, visc*1e-3)
    fitmin  = min(continConfig['fitRangeM'])/gamma
    fitmax  = max(continConfig['fitRangeM'])/gamma
    name = "{} a{:03.0f}".format(meta['Samplename'], angle)
    Im, dIm = 1, 0
    # get measured correlation data and tau
    dataStartEnd = 0, 0
    with open(filename, encoding='cp1250') as fd:
        content = fd.readlines()
        dataStartEnd = [idx for idx, line in enumerate(content)
                        if line.strip() in ('"Correlation"', '"Count Rate"')]
        if len(dataStartEnd) == 1: # no end found in 'averaged' data files
            dataStartEnd.append(len(content))
    dlsData = pd.read_csv(filename, r'\s+', skiprows=dataStartEnd[0]+1,
                       nrows=dataStartEnd[1]-dataStartEnd[0]-1, encoding='cp1250', 
                       names=['tau']+anglesFromFile)
    dlsData.tau *= 1e-3 # convert to seconds
    # restrict data to given range
    tmin, tmax = min(continConfig['ptRangeSec']), max(continConfig['ptRangeSec'])
    tmask = np.logical_and(tmin <= dlsData.tau, dlsData.tau <= tmax)
    tauCropped = dlsData.tau[tmask]
    corCropped = dlsData[angle][tmask]
    a2s_kwargs = dict(floatmode='fixed', sign=' ', max_line_width=80,
                      formatter={'float_kind': '{0: .5E}'.format})
    tauStr  = np.array2string(tauCropped.values, **a2s_kwargs)[1:-1]
    corStr  = np.array2string(corCropped.values, **a2s_kwargs)[1:-1]
    npts = len(tauCropped)
    # generate CONTIN input file
    content = f"""{name}
 IFORMY    0    .00
 (6E13.7)
 IFORMT    0    .00
 (6E13.7)
 IFORMW    0    .00
 (6E13.7)
 NINTT     0   -1.00
 DOUSNQ    0    1.00
 IQUAD          1.00
 PRWT      0    1.00
 PRY       0    1.00
 IPLRES    1    3.00
 IPLRES    2    3.00
 IPRINT    1    0.00
 IPRINT    2    2.00
 IPLFIT    1    0.00
 IPLFIT    2    0.00
 NONNEG    1    1.00
 IWT       0    {IWT:.2f}
 NQPROG    1    5.00
 NQPROG    2    5.00
 GMNMX     1    {fitmin:.3E}
 GMNMX     2    {fitmax:.3E}
 RUSER    10    {Trd:.2f}
 NG        0    {{gridpts:.2f}}
 NLINF     0    {{freeBaseline:.2f}}
 IUSER    10    4.00
 RUSER    21    1.0
 RUSER    22   -1.0
 RUSER    23    0.0
 RUSER    18         {temp:.5f}
 RUSER    17         {angle:.5f}
 RUSER    19         {visc:.5f}
 RUSER    15         {refrac:.5f}
 RUSER    16         {wavelen:.5f}
 RUSER    25         {Im}
 RUSER    26         {dIm}
 END       0    0.00
 NY{npts: >9d}
 {tauStr}
 {corStr}
""".format(**continConfig
    )
    return content.encode('ascii')

def getContinOutputDirname(angle):
    return f"contin_{angle:03.0f}"

def runContin(filenameAndConfigTuple, recalc=True):
    """Starts a single CONTIN process for the given DLS DataSet
    (which should contain a single angle only)."""
    continCmd = getContinPath()
    assert continCmd.is_file(), "CONTIN executable not found!"
    filename, continConfig = filenameAndConfigTuple
    filename = Path(filename)
    try:
        continInData = genContinInput(filename, **continConfig)
    except AssertionError:
        print(f"Scattering angle {continConfig['angle']} not found!\n"
              f"    Skipping '{filename.name}'.")
        return
    workDir = filename.parent
    #ts = datetime.datetime.now().strftime("%Y%m%d-%H%M%S") # timestamp
    tmpDir = workDir / (getContinOutputDirname(continConfig['angle'])+' '+filename.stem)
    if tmpDir.is_dir(): # deleting old results
        if not recalc:
            return tmpDir
        shutil.rmtree(tmpDir)
    os.mkdir(tmpDir)
    continInDataPath  = tmpDir / continConfig['infn']
    continOutDataPath = tmpDir / continConfig['outfn']
    # Store input data
    with open(continInDataPath, 'wb') as fd:
        fd.write(continInData)
    with pushd(tmpDir):
        proc = subprocess.run([str(continCmd)], input=continInData,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if len(proc.stderr):
            print(str(continInDataPath)+':')
            print(proc.stderr.decode())
    # Store output data
    with open(continOutDataPath, 'wb') as fd:
        fd.write(proc.stdout)
    return tmpDir

def processFiles(fnLst, config, nthreads=None):
    start = time.time()
    if nthreads == 1:
        resultDirs = [runContin((dn, config), recalc=True) for dn in fnLst]
    else: # Using multiple CPU cores if available
        import multiprocessing
        if not nthreads:
            nthreads = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes = nthreads)
        resultDirs = pool.map(runContin, [(dn, config) for dn in fnLst])
        pool.close()
        pool.join()
    print()
    print(f"CONTIN analysis with {nthreads} thread{'s' if nthreads > 1 else ''} took {time.time()-start:.1f}s.")
    return [rd for rd in resultDirs if rd is not None]
