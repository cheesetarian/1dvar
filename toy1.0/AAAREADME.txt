
...................................................................
Last update 1/13/16 (DD, CSU) to ver1.0 
...................................................................

...................................................................
Toy version of full ver1.0---
Instead of reading in a full orbit of L1 data, the toy version takes
 all necessary inputs from the header of the definition file and prints
 out the retrieved variables.
...................................................................

CRTM-based 1DVAR with support for multiple sensors, channel configurations,
 emissivity models, ancillary data sources, and output formats.


last updates:
**** 12/24 'BCMZ'-- updated (v3) eof binary file that allows 
         interpolation between SST bins. posterior errors of TPW and LWP 
	 now correct. also now able to use MERRA data for first guess.
**** 01/13/16 -- V1.0 -- Same as BCMZ but renamed and cleaned up.


contains:
************************************************************************

run_' '.sh -- shell script that calls the preprocessor and calls oeorbit.
              this is where number of channels is defined, as is the
              calibration table used in the preprocessor with no need to 
              recompile the code. emissivity model, method of output, and
	      specification file all contained here.

oeorbit.f -- call routines to read pp file, calculate azimuth angle, import
	     analysis data, read and assign sfc type and sst to each
             pixel, initialize CRTM, then calls the OE, runs quality check,
	     then calls subroutine to output to a binary file

subr_oeorbit.f -- subroutines called by oeorbit. includes code to read pp
             data, read and assign sst and surface type, 
             read lookup tables and analysis data, output retrieved values, etc.

nr_subrs.f -- fortran module containing the subroutines of the 1DVAR

define_oeorbit.f -- fortran module with explicit and dynamic definitions 
		    of common variables used in multiple subroutines

crtm.f90 -- fortran module containing subroutines to initialize CRTM for many
	    sensor/channel configurations, run the CRTM based on input from
	    nr_subrs.f, and destroy the CRTM.

**********************************************************************************************


old versions log:
*** new version 9/30/14 ***
- using ECMWF Interim reanalysis-derived values to create profile of 
   temperature and moisture, a priori and variances for TPW and winds
   based on SST bins. this requires modifying the RT and a fair amount
   of the nr_subrs module.
**** 11/25/14 version using L2A Tbs, RSS SST (not retrieved) ****
**** 02/15/15 'l2_crtm' uses L2A data as input, but is mainly just
        a first attempt at using the CRTM instead of the old Kummerow
	RT model.
**** 03/04/15 'gmi_crtm' is an OE version that uses the CRTM for the fwd
	model, but is specifically tailored to work for GMI. Note: THE
	LUTs BEING USED CURRENTLY FOR PRIOR INFO AREN'T FROM THE RIGHT YEAR!
	also, GMI is tricky due to having two different EIAs for low/high
	frequency channels, taken into account in crtm.f90.
**** 03/11/15 'gmic_p1' follows on from gmi_crtm but utilizes a wholly new
	profile scheme to provide more physical profiles of atmospheric
	water vapor (in noscatter subroutine)
**** 03/23/15 'gmic_pL' uses the same profile scheme as p1/p2 but now allows
	LRH to vary as a retrieved variable
**** 04/01/15 'gmic_pL1' is the first retrieval to read in means and variances
	for the retrieved RH values, as well as lapse rate, covariances,
	and wind information. look in oe/make_luts/rh_analysis and rh_condense
	for what is in the binary files. all data are taken from 6-hrly ECMWF
	Interim reanalysis, then split into monthly chunks, smoothed, etc.
	lapse rates are from 0-10km. vertical levels are at 1km spacing, and
	defined at 15.5km, 10.5km, 4.5km, 0.5 km [middle of 1km levels].
**** 04/03/15 'gmic_pL2' makes many tweaks to pL1, including what is output,
	what can be retrieved, bug fixes, and some cleanup too.
**** 04/21/15 'gmic_pL3' includes the first attempt at a full Sy matrix 
	that comes from recent analysis
**** 05/06/15 'GCL6' is gmic_pL3 but with no sa covariances and some tweaks
	to make it better at running 5/7ch. RH=100 in cloud layers,
	ability to run 7/11ch too, miss_flt for output whether chisq too big
  	or no convergence reached.
**** 05/07/15 'GCL7' changes the vertical tie points from GCL6 to 0, 2, 6km
	reading in vR3 input files. new output field for posteriori error, 
	tb_diff matches available channels. posterior errors for every pixel
	where convergence reached.
	(adding 5/10)-- added functionality for using RSS emissivity model
	instead of CRTM's FASTEM4/5. also improved first guess.
**** 05/11 'GCL9' keeps same tie points as '7' (since '8' with vR4 wasn't 
	better), improves first guess, tries keeping forced RH in cloud and
	bringing back some Sa covariances.
**** 05/15 'BCL1' improves upon GCL9's success. default is to not fix RH
	at 100 in cloud, now have tb offset possible to match Sy analysis,
	using new SY3 channel errors. 'B' signals move toward ability to 
	run AMSR2 or GMI (or other sensors), though not functional yet.
**** 05/29 'BCE1' switches to using EOFs for water vapor profiles, hence
	new LUTs and output structure...
**** 06/05 'BCE2' gets the EOFs right side up, new output structure, 
	assumes constant lapse rate everywhere of 6.26K/km (for now!),
	and a few other tweaks
**** 06/09 'BCM1' moves to pressure levels in fwd model and EOFs of mass
	mixing ratio instead of relative humidity
**** 06/15 'BCM2' changes output structure, now has max of 3 EOFs
**** 06/19 'BCM3' uses new constraints on EOF coefficients from analysis
	of gridded EC files, code found in make_luts
	update 7/8-- a couple 'bug' fixes in OE code itself that greatly
	 limits # iterations needed. no RH constraints/checks now, and
	 moved to outputting a 'normalized' chisq (chisq/nch).
**** 07/08 'BCM4' is bug fixes, back to better first guess from most recent
	 nearby convergence, no RH constraint. set in stone?
**** 07/20 'BCM5' includes new SST-indexed Sy matrix with biases included for 
	 FASTEM5 (default), FASTEM4, and MW/RSS, along with call for emis
	 model built into script functionality. also, ability to read
	 in EC winds with direction, and get FOV azimuthal angle for GMI
	 to give CRTM wind angle. reading in analysis data :/. also now
	 interpolating assumed mr profile between SST bins.
**** 07/23 'BCM6' will retain enhancements from BCM5 but remove SST-indexed
	 Sy matrix functionality since that didn't show any improvment.
	 more representative mean and EOF profiles for high SST bins now
	 incorporated in BCM6.
**** 08/03 'BCM7' experimental setup to read in wind, direction, as well as
	 mixing ratio and temperature profiles maybe from ERA-I as first 
	 guess and a priori
**** 08/11 'BCM8' will contain a switch for reading in ERA analysis data or 
	 not (prior state is either analysis or climatology), as well as 
	 off-diagonal elements of Sa matrix from SYM3 analysis
**** 08/12 'BCM9' builds on BCM8 but includes ability to scale down eof
	 variances if using analysis data
**** 08/27 'BCMS' S for surface since using SLP output from EC now in the 
	 code as well. a few other tweaks learned from SY improvements.
**** 10/04 'BCMT' just a few minor changes: bigger minimum wind speed, moving
	 closer into the coast, clwc/ciwc output in RT file
**** 10/05 'BCMU' now reads EC data directly from binary instead of calling
	 IDL function, gets Sy from tdata1 folder, has experimental way of
	 cutting down number of times Jacobian is calculated to speed it up.
**** 10/14 'BCMV' updates Sy read, adds L1CR output field to RT files, gets
	 rid of symult and replaces with 'syf' for sy fraction
**** 10/28 'BCMW' builds on BCMV, removes Tb offsets (for now, at least in
	 analysis mode), cleans up code, includes ODAS/ODPS binary files,
	 new Sy matrix calls for updated GMI spectral functions in CRTM
**** 11/09 'BCMX' changes nothing in BCMW apart from the output format.
	 sun glint is now read and used for determining a quality flag.
	 secondary convergence criteria slightly modified.
**** 11/16 'BCMY' removes 2nd conv crit (tightens 1st slightly) and using EC
	 LWP as prior, decreases LWP prior from -2 to -3 and increases
	 LWP sigma to 2.0 from 1.6; this solves iterative issues in the
	 tropics of trading off LW/WV. also decreases chisq threshold to
	 4 from 7, and uses x2 for jacobian slope instead of xa2.
	 'BCMY' is a fairly major version update, affecting # iterations
	 and LWP in Tropics especially.
	 more mods 11/17, still BCMY: new output screening procedure
	 (qual_czech) to screen out anomalously high TPW pixels. the issue
	 of TPW vs. LWP persists, but this is a way to mitigate it. a value
	 has also be introduced to the def file to determine the intermediate
	 quality flag via chisq.
	 Crucially, BCMY moves to intrinsic matrix math functions of f90+,
	 speeding up the retrieval slightly and providing greater accuracy.
**** 12/24 'BCMZ'-- updated (v3) eof binary file that allows 
         interpolation between SST bins. posterior errors of TPW and LWP 
	 now correct. also now able to use MERRA data for first guess.
