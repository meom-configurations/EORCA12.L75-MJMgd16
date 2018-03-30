# EORCA12.L75-MJMgd16

EORCA12.L75-MJMgd16 is the NEMO configuration set up and run during the Grand Challenge 2016 on OCCIGEN super computer at CINES.

This repository hold the code, and configuration files usefull for running this configuration.

## REFERENCE CODE : 
 This configuration is based on the NEMO configuration from rev 7046 of the NEMO trunk. It is used together with the XIOS server at rev. 924 of the XIOS trunk. Both codes can be downloaded from the IPSL forge using the following statements :
### NEMO
 ```svn co -r 7046 http://forge.ipsl.jussieu.fr/nemo/svn/trunk/NEMOGCM ```
### XIOS
 ```svn co -r 924 http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/trunk xios_2.0_rev_924 ```
 
## BRIEF DESCRIPTION:
### Overview
   This global configuration  uses the eORCA12 grid, with the standard DRAKKAR 75 levels. The ORCA12 domain was extendend  southward in order to be able to take into account the Ice Shelves (eORCA12). In this particular realisation, we do not resolve explicitely the under ice-shelf ocean circulation but uses the Mathiot parameterizations with prescribed fresh water fluxes due to ice-shelf melting.  For production reasons, the domain for this configuration is reduced to open ocean only grid points (*i.e* the domain corresponding to ice-shelves has been removed). This is why the configuration is called EORCA12, and not eORCA12.
   
   This simulation was performed in the frame of the 2016 *grands defis* on CINES/OCCIGEN2 super computer. The project was entitled **SINGOUT** for SImulation Nouvelle Generation 
   
###  Parameterizations:
 1. use non linear free surface (VVL)
 2. use UBS advection scheme for dynamics
 3. use FCT4 advection scheme for tracers.
 4. no explicit lateral viscosity  ( UBS scheme does it)
 5. use laplacian isopycnal diffusivity for tracers.
 
 ### Forcing:
  1. Atmospheric forcing is DFS5.2, with ```CORE bulk formulae``` 
  2. SSS restoring toward WOA9 with a piston velocity of 167 mm/day ( 60 days/10 meters).
  3. Run-off from Dai-Trenberth as usual, except around Antarctica ( explicit iceberg calving and ice-shelf melting).
  
  ### Antarctic fresh water fluxes:
  1. use of ICB module to represent iceberg calving and melting explicitely (Merino et al )
  2. use ice-shelf parameterization to represent melting ( Mathiot et al.).
  
  ### XIOS output ( all file in netcdf4 with deflation level 1).
  1. Due to vvl use weighted average (e3 ) when relevant.
  2. 1d output (170 Gb/year)
     * **gridTsurf** files :SST, SSS, SST, MXL
     * **gridUsurf** files :SSU, U10m, U30m, U50m 
     * **gridVsurf** files :SSV, V10m, V30m, V50m
  3. 5d output (1.1 Tb/year)
     * **gridT** files : e3t, votemper, vosaline, sossheig
     * **gridU** files : e3u, vozocrtx, vozotaux
     * **gridV** files : e3v, vomecrty, vometauy
     * **gridW** files : e3w, vovecrtz, voavt, voavmu, voavmv 
     * **flxT** files : 17 fluxes/forcing variables
     * **icemod3** files : 22 LIM3 variables.
     * **ICB** files : 16 iceberg related variables.
     
   ### Run time files:
      Most of the run time files are indicated in the namelist files, except for :
      * bathymetry : ```EORCA12_bathymetry_v2.5.nc```
      * coordinates : ```EORCA12_coordinates.nc```
      * bottom friction : ```EORCA12_bfr2d_UKmod.nc ```
