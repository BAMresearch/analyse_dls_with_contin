# -*- coding: utf-8 -*-
# utils.py

import os, shutil, subprocess
import io
import time
import urllib
import requests
import itertools
from pathlib import Path
from bs4 import BeautifulSoup
import pandas as pd
import numpy as np
from .jupyter_analysis_tools.utils import isWindows, isMac, isList, pushd, grouper
from .dlshelpers import getDLSgammaSi, getDLSFileData

InputFn = "contin_in.txt"
OutputFn = "contin_out.txt"

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
    data = getDLSFileData(filename)
    # select the measurement angle, make sure it's in the file
    angle = continConfig['angle']
    assert angle in data['angles'], \
        f"Given angle ({angle}) not found in file '{filename}': {data['angles']}"
    # get environment values for storage in contin file
    temp    = data['Temperature [K]']
    visc    = data['Viscosity [cp]']
    refrac  = data['Refractive Index']
    wavelen = data['Wavelength [nm]']
    gamma = getDLSgammaSi(angle, refrac, wavelen*1e-9, temp, visc*1e-3)
    fitmin  = min(continConfig['fitRangeM'])/gamma
    fitmax  = max(continConfig['fitRangeM'])/gamma
    Im, dIm = 1, 0
    # get measured correlation data and tau
    dlsData = data["correlation"]
    dlsData.reset_index(inplace=True)
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
    content = f"""{filename.name}
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

def runContin(filenameAndConfigTuple):
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
        if not continConfig.get("recalc", True):
            return tmpDir
        shutil.rmtree(tmpDir)
    os.mkdir(tmpDir)
    continInDataPath  = tmpDir / InputFn
    continOutDataPath = tmpDir / OutputFn
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

def runContinOverFiles(fnLst, config, nthreads=None):
    start = time.time()
    if not isList(config):
        config = (config,)
    if nthreads == 1:
        resultDirs = [runContin((dn, cfg)) for dn in fnLst for cfg in config]
    else: # Using multiple CPU cores if available
        import multiprocessing
        if not nthreads:
            nthreads = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes = nthreads)
        resultDirs = pool.map(runContin, [(dn, cfg) for dn in fnLst for cfg in config])
        pool.close()
        pool.join()
    print()
    print(f"CONTIN analysis with {nthreads} thread{'s' if nthreads > 1 else ''} took {time.time()-start:.1f}s.")
    return [rd for rd in resultDirs if rd is not None]

def convertContinResultsToSizes(lines, df):
    # user variables for environmental values as set by CNTb scripts
    varmap = dict(temp="RUSER    18", angle="RUSER    17", visc="RUSER    19",
                  refrac="RUSER    15", wavelen="RUSER    16")
    indices = getLineNumber(lines, list(varmap.values()))
    for name, idx in zip(varmap.keys(), indices):
        varmap[name] = float(lines[idx].split()[-1])
    # convert to SI units
    varmap["visc"] *= 1e-3
    varmap["wavelen"] *= 1e-9
    gamma = getDLSgammaSi(varmap["angle"], varmap["refrac"], varmap["wavelen"], varmap["temp"], varmap["visc"])
    # unclear if gamma needs doubled due to G2(t)-1 = g1(t)^2 = exp(-t*gamma)^2
    df["abscissa"] *= gamma
    df.rename(columns={'abscissa': 'radius', 'ordinate': 'distrib', 'error': 'err'},
              inplace=True)
    return df, varmap

def getLineNumber(lines, phrases, debug=False):
    """Returns the line numbers containing the provided phrases after searching
    for the previous phrases sequentially. Search starts with the first phrase,
    once it is found, search starts with the 2nd phrase from that line,
    until the last phrase is found. Ignores early matches of the final phrase."""
    nums = []
    for i, line in enumerate(lines):
        if phrases[len(nums)] in line:
            if debug:
                print("found '{}' on line {}.".format(phrases[len(nums)], i))
            nums.append(i)
            if len(phrases) == len(nums):
                return nums
    return nums

def getContinInputCurve(inputAbsPath):
    assert inputAbsPath.is_file()
    # read in line by line, some adjustments required for parsing floats
    startLine, count = 0, 0
    with open(inputAbsPath) as fd:
        startLine, count = [(idx, int(line.split()[-1]))
                            for idx, line in enumerate(fd) if "NY" in line][0]
    lines = []
    with open(inputAbsPath) as fd:
        lines = fd.readlines()
    return [float(f) for line in lines[startLine+1:] for f in line.split()][count:]

def getContinResults(sampleDir, angle=None):
    """*sampleDir*: A pathlib Path of the location where the CONTIN results can be found."""
    sampleDir = Path(sampleDir)
    # check first if there was any CONTIN output generated
    resultsFile = sampleDir / OutputFn
    if not resultsFile.is_file():
        # try the subdir first
        assert angle is not None, "An angle in degrees has to be provided"
        resultsFile = sampleDir / getContinOutputDirname(float(angle)) / OutputFn
        if not resultsFile.is_file():
            print("No distribution found in\n '{}'!".format(resultsFile.parent))
            return None, None
    # read in line by line, some adjustments required for parsing floats
    lines = []
    with open(resultsFile) as fd:
        lines = fd.readlines()
    # find the beginning and end of the fitted correlation curve
    startLines = getLineNumber(lines, ["T            Y", "0PRECIS"])
    if not len(startLines):
        print(f"Fitted curve not found in CONTIN output!\n ({resultsFile})")
        return None, None
    dfStart, dfEnd = startLines[-2]+1, startLines[-1]
    dfFit = pd.DataFrame([f for line in lines[dfStart:dfEnd] for f in grouper(line.split(), 2)],
                         columns=('tau', 'corrFit'), dtype=float)
    dfFit.corrFit = dfFit.corrFit**2 # to be compared with measured data
    # get input correlation curve first, to be added to fitted correlation curve
    dfFit['corrIn'] = getContinInputCurve(sampleDir/InputFn)
    # find the beginning and end of the distribution data
    startLines = getLineNumber(lines, ["CHOSEN SOLUTION", "ORDINATE", "LINEAR COEFFICIENTS"])
    if not len(startLines):
        print(f"Distribution data not found in CONTIN output!\n ({resultsFile})")
        return None, None
    dfStart, dfEnd = startLines[-2]+1, startLines[-1]
    lineEnd = lines[dfStart].index("X")
    # convert CONTIN output distrib to parseable data for pandas
    dfDistrib = pd.read_csv(io.StringIO("\n".join([line.replace("D", "E")[:lineEnd]
                                for line in lines[dfStart:dfEnd]])),
                            delim_whitespace=True, names=("ordinate", "error", "abscissa"))
    dfDistrib, varmap = convertContinResultsToSizes(lines, dfDistrib)
    # parse original input data as well, if available
    infn = sampleDir.parent / lines[0][52:].strip()
    varmap['dataFilename'] = infn
    return dfDistrib, dfFit, varmap
