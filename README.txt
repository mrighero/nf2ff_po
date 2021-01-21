Matlab scripts and data are in the compressed folder code_and_data.zip

The main driver is call_nf2ff_po.m; it sets the different parameters and call the main routine nf2ff_po.p.
The parameters are explained within call_nf2ff.m.
The main routine nf2ff_po calls some routines stored in the folder ./auxiliary_routines_p

The folder ./field/ contains measured fields, stored according to the format described below (*)
The folder ./mesh/ contains meshes used to run the main routine




(*) Measured data stored in a matlab structure with the following fields
- ScanDistance'='NEAR'
- 'freq' = working frequency
- 'dataCoords'= N, x, y, z, u1x, u1y, u1z, ..., uNx, uNy, uNz, w, where:
	- N = number of polarizations (usually two);
	- x, y, z = Cartesian coordinates of the measured point [m]
	- unx, uny, unz	= Cartesian components of the polarization versor “n”, n = 1, ..., N
	- w		= weight if needed to de-emphasize some samples (e.g. around poles)
- 'dataMeas'= E1Re, E1Im, ..., ENRe, ENIm, where EnRe, EnIm are the Real and Imaginary part of the field along the versor “n”, n = 1, ..., N
- 'Comment'= a string with additional comment 
