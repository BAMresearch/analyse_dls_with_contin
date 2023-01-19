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
import sys
cwd = Path().resolve()
if cwd not in sys.path:
    sys.path.insert(0, str(cwd))
from jupyter_analysis_tools.utils import isWindows, isMac, isList, pushd, grouper, updatedDict
from jupyter_analysis_tools.analysis import getModZScore
from dlshelpers import getDLSgammaSi, readDLSData

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

def genContinInput(filedata, **continConfig):
    """Expects a dictionary of file data created by readDLSData()."""
    IWT = 5 if continConfig['weighResiduals'] else 1
    # transform data? Trd=0: no transform
    Trd = -1 # Trd=1: initial g(2), input sqrt[g(2)-1]; Trd=-1: initial g(2)-1, input sqrt[g(2)-1]
    # select the measurement angle, make sure it's in the file
    angle = continConfig['angle']
    assert angle in filedata['angles'], \
        f"Given angle ({angle}) not found in file '{filedata['filename']}': {filedata['angles']}"
    # get environment values for storage in contin file
    temp    = filedata['Temperature [K]']
    visc    = filedata['Viscosity [cp]']
    refrac  = filedata['Refractive Index']
    wavelen = filedata['Wavelength [nm]']
    gamma = getDLSgammaSi(angle, refrac, wavelen*1e-9, temp, visc*1e-3)
    fitmin  = min(continConfig['fitRangeM'])/gamma
    fitmax  = max(continConfig['fitRangeM'])/gamma
    Im, dIm = 1, 0
    # get measured correlation data and tau
    dlsData = filedata["correlation"].reset_index()
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
    # store the score for this data if it was determined
    scoreRecord = f"\n RUSER    11    {filedata['score'][angle]:.3E}" if 'score' in filedata else ""
    storedFn = filedata['filename'].name
    if not filedata['filename'].is_file():
        storedFn = filedata['filename'].parent.suffix + '/' + storedFn
    # generate CONTIN input file
    content = f"""{storedFn}
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
 LINEPG    0   50.
 NONNEG    1    1.00
 IWT       0    {IWT:.2f}
 NQPROG    1    5.00
 NQPROG    2    5.00
 GMNMX     1    {fitmin:.3E}
 GMNMX     2    {fitmax:.3E}
 RUSER    10    {Trd:.2f}
 NG        0    {{gridpts:.2f}}
 NLINF     0    {{baselineCoeffs:.2f}}
 IUSER    10    4.00{scoreRecord}
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

def workerInit(_queue):
    """Initializes a queue for log messages in each worker process during multiprocessing.
    Queue object is global within each process only, not in the parent."""
    global queue
    queue = _queue

def runContin(filedata, continConfig, useQueue=True):
    """Starts a single CONTIN process for the given DLS DataSet
    (which should contain a single angle only)."""
    continCmd = getContinPath()
    assert continCmd.is_file(), "CONTIN executable not found!"
    logPrefix = f"{filedata['filename'].name}@{continConfig['angle']}°: "
    workDir = filedata['filename'].parent
    if workDir.is_file():
        logPrefix = workDir.stem + '/' + logPrefix
        workDir = workDir.parent / workDir.stem
    logFunc = queue.put if useQueue else print
    def log(text):
        #print("="+logPrefix+text)
        logFunc(" "+logPrefix+text)
    try:
        continInData = genContinInput(filedata, **continConfig)
    except AssertionError:
        log(f"Scattering angle {continConfig['angle']} not found! "
            f"Skipping…")
        return
    #ts = datetime.datetime.now().strftime("%Y%m%d-%H%M%S") # timestamp
    tmpDir = workDir / (getContinOutputDirname(continConfig['angle'])+' '+filedata['filename'].stem)
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
            log(proc.stderr.decode().strip())
    # Store output data
    with open(continOutDataPath, 'wb') as fd:
        fd.write(proc.stdout)
    return tmpDir

def readData(fnLst, configLst):
    angles = [cfg['angle'] for cfg in configLst]
    dataLst = readDLSData(fnLst)
    # calc modified Z-Score based on median absolute deviation
    # for each count rate at the same angle
    for angle in set(angles):
        try:
            idx, cr = zip(*((i, filedata['countrate'][angle].values)
                             for i, filedata in enumerate(dataLst)
                                 if angle in filedata['countrate']))
        except (ValueError, KeyError): # filedata without 'countrate', most probably
            continue
        #print(angle, idx, getModZScore(np.stack(cr)))
        for i, score in zip(idx, getModZScore(np.stack(cr))):
            if 'score' not in dataLst[i]:
                dataLst[i]['score'] = dict()
            dataLst[i]['score'][angle] = score
    return dataLst

def runContinOverFiles(fnLst, configLst, nthreads=None, outputCallback=None):
    """*fnLst*: List of file paths to .ASC files
       *configLst*: List of parameters, one dict for each file, such as
            {'recalc': True, 'gridpts': 200, 'transformData': True,
             'ptRangeSec': (3e-07, 1.0), 'fitRangeM': (7e-10, 3.9e-07),
              'baselineCoeffs': 0, 'weighResiduals': True}
       *nthreads*: number of parallel CONTIN processes to launch,
           1: Sequential processing, one file after another
           None: number of processes equals the number of computing cores
       *outputCallback*: A function with one argument to called repeatedly (0.5s)
           with new output messages combined from all CONTIN processes.
    """
    start = time.time()
    if not isList(configLst):
        configLst = (configLst,)
    dataLst = readData(fnLst, configLst)
    # get all combinations of CONTIN parameters and data files
    dataNConfig = [(data, cfg) for data in dataLst for cfg in configLst]
    if nthreads == 1:
        resultDirs = [runContin(data, cfg, False) for data, cfg in dataNConfig]
    else: # Using multiple CPU cores if available
        import multiprocessing
        if not nthreads:
            nthreads = multiprocessing.cpu_count()
        from multiprocessing import Queue as MPQueue
        # use a queue to collect stdout messages from subprocesses
        logQueue = MPQueue()
        pool = multiprocessing.Pool(processes=nthreads, initializer=workerInit, initargs=(logQueue,))
        resultDirs = pool.starmap_async(runContin, dataNConfig)
        pool.close()
        def resultReady(asyncResult): # checks if the overall result is ready
            try:
                asyncResult.successful()
            except ValueError:
                return False
            return True
        outputBuffer = [] # buffer to store output messages from queue in,
                          # for sorting, for deterministic testing
        while not resultReady(resultDirs):
            time.sleep(.5) # update interval of output
            while not logQueue.empty():
                newOutput = logQueue.get_nowait()
                if not outputCallback:
                    print(newOutput) # the traditional way
                else:
                    outputBuffer.extend(newOutput.splitlines())
            if callable(outputCallback):
                # use a custom callback to handle the output from subprocesses
                outputCallback("\n".join(sorted(outputBuffer)))
        #print("READY!")
        resultDirs = resultDirs.get()

    summary = f"CONTIN analysis with {nthreads} thread{'s' if nthreads > 1 else ''} took {time.time()-start:.1f}s."
    return [rd for rd in resultDirs if rd is not None], summary

def getValueDictFromLines(lines, **kwargs):
    """Searches the given list of lines for the keys of the given arguments
    and converts the values to float. Returns the completed dict."""
    # search begin of common variables
    lstart = [idx for idx, line in enumerate(lines)
                  if "INPUT DATA FOR CHANGES TO COMMON VARIABLES" in line]
    # search end of common variables section (where the next begins)
    lend   = [idx for idx, line in enumerate(lines)
                  if "FINAL VALUES OF CONTROL VARIABLES" in line]
    result = dict()
    if len(lstart) and len(lend):
        result = {key: float(line.split()[-1])
                    for line in lines[lstart[0]:lend[0]]
                    for key, pattern in kwargs.items() if pattern in line}
    return result

def getContinUserVars(lines):
    """Extract previously set user variables for environmental values
    from CONTIN output data.
    *lines*: List of lines of CONTIN output data."""
    varmap = getValueDictFromLines(lines,
                temp="RUSER    18", angle="RUSER    17", visc="RUSER    19",
                refrac="RUSER    15", wavelen="RUSER    16", score="RUSER    11")
    # convert to SI units
    varmap["visc"] *= 1e-3
    varmap["wavelen"] *= 1e-9
    varmap["gamma"] = getDLSgammaSi(varmap["angle"], varmap["refrac"], varmap["wavelen"],
                                    varmap["temp"], varmap["visc"])
    return varmap

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
            return None, None, None
    # read in line by line, some adjustments required for parsing floats
    lines = []
    with open(resultsFile) as fd:
        lines = fd.readlines()
    # find the beginning and end of the fitted correlation curve
    startLines = getLineNumber(lines, ["T            Y", "0PRECIS"])
    if not len(startLines):
        print(f"Fitted curve not found in CONTIN output!\n ({resultsFile})")
        return None, None, None
    dfStart, dfEnd = startLines[-2]+1, startLines[-1]
    dfFit = pd.DataFrame([f for line in lines[dfStart:dfEnd] for f in grouper(line.split(), 2)],
                         columns=('tau', 'corrFit'), dtype=float)
    dfFit.corrFit = dfFit.corrFit**2 # to be compared with measured data
    # get input correlation curve first, to be added to fitted correlation curve
    dfFit['corrIn'] = getContinInputCurve(sampleDir/InputFn)
    # find the beginning and end of the distribution data
    startLines = getLineNumber(lines, ["CHOSEN SOLUTION", "ORDINATE"])
    if not len(startLines):
        print(f"Distribution data not found in CONTIN output!\n ({resultsFile})")
        return None, None, None
    gridSize = int(getValueDictFromLines(lines, distribSize="NG        0").get('distribSize',0))
    dfStart = startLines[1]+1
    lineEnd = 31 # do not parse floats beyond this column
    # convert CONTIN output distrib to parseable data for pandas
    fixedFloatFmt = io.StringIO("\n".join([line[:lineEnd].replace("D", "E")
                                for line in lines[dfStart:dfStart+gridSize]]))
    dfDistrib = pd.read_csv(fixedFloatFmt, delim_whitespace=True,
                            names=("distrib", "err", "decay"))
    dfDistrib = dfDistrib[["decay", "distrib", "err"]] # reorder to (x,y,u)
    # update x/abscissa with values from another section of the output
    # to avoid duplicates due to low precision in solution output parsed above
    startLines = getLineNumber(lines, ["GRID POINT"])
    if len(startLines):
        dfStart = startLines[0]+1
        decayNew = np.fromiter([line.split()[0] for line in lines[dfStart:dfStart+gridSize]], float)
        dfDistrib.decay = decayNew
    varmap = getContinUserVars(lines)
    # parse original input data filename as well, if available
    # see genContinInput() for creating storedFn
    storedFn = lines[0][52:].strip()
    infn = sampleDir.parent / storedFn
    if storedFn[0] == '.': # starts with a dot
        infn = sampleDir.parent.parent / (sampleDir.parent.name + storedFn)
    varmap['dataFilename'] = infn
    return dfDistrib, dfFit, varmap
