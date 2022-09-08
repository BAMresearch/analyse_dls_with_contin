# Test data for Dynamic Light Scattering analysis

## Measurement conditions

This folder contains raw DLS measurements acquired with the DLS instrument *ALV/CGS-3 Compact Goniometer System* bei ALV GmbH (Langen, Germany) measuring at 8 scattering angles with the control software *ALV-MultiAngle Version 3.1.4.3* on Windows 10.

The main environmental conditions of each measurement are stored in each data file by the control software: measurement mode, number of runs, measurement duration, ambient temperature, viscosity, refractive index, laser wavelength, date and time as well as a short sample description provided by the operator.

## Materials measured

- Water, the buffer medium surrounding the particles  
  https://en.wikipedia.org/wiki/Water
- Toluene, used for calibrating the instrument  
  https://en.wikipedia.org/wiki/Toluene
- MW002-03: Polypropylene particles in water with a radius of approximately 64 nm.
  Their synthesization was adapted from the methods described in https://doi.org/10.1039/C8NA00210J.  
  Cumulant analysis according to ISO 13321 gave:
  ```
  Mean hydrodynamic radius             <Rh>     = 63.982 +/- 2.003 nm 
  Extrapolation to q^2 = 0,            R_h(0)   = 61.558 +/- 0.688 nm 
  Extrapolated R_h at 90 deg,          R_h(90)  = 64.066 +/- 1.006 nm 
  Extrapolated R_h at ZetaSizer angle, R_h(173) = 66.767 +/- 1.673 nm 
  Mean PDI                             <PDI>    = 0.0799 +/- 0.0310 
  Extrapolated PDI at 90 deg,          PDI(90)  = 0.0720 +/- 0.0062 
  Extrapolated PDI at ZetaSizer angle, PDI(173) = 0.0717 +/- 0.0062 
  Mean sigma(Rh)                             <sigma>    = 12.1139 +/- 2.3710 
  Extrapolated sigma(Rh) at 90 deg,          sigma(90)  = 11.8678 +/- 0.4446 
  Extrapolated sigma(Rh) at ZetaSizer angle, sigma(173) = 11.8624 +/- 0.0692 
  
  Apparent diffusion constant extrapolated to q^2 = 0: 
                           D_0 = 3.463e-08 +/- 3.898e-10 cm^2/s 
  Hydrodynamic radius from D_0 :
                           Rh(D_0) = 61.524 +/- 0.692 nm
  ```
- *cgs3*: Measurements with a ALV/CGS-3 Goniometer device in a slightly different file format as provided by Vincent in issue #1. At 13 different angles, three measurements per angle were aquired plus an averaged file computed from these measurements.

## Intended use
The main purposes of this data are functional and regression tests of the analysis Python code, respectively reading the data format correctly and reproducing analysis results during iterative software development steps consistently.

### Disclaimer
The data are not meant to verify the analysis algorithms (and results) for correctness since the materials measured here are **no** certified reference materials.
