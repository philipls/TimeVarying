
Matlab/C routines to fit the standard and time-varying diffusion models to the data from the Dutilh et al. (2019) study. To accompany Smith, P. L. & Lilburn, S. D. Vision for the blind: Visual psychophysics and blinded inference for decision models.  Psychonomic Bulletin & Review, 2020. 

The raw data for the Dutilh et al. (2019) study are downloadable from the OSF website at:

https://osf.io/jy7pi

The individual subject data structures in the file osf4.mat were created from the raw data at this site as described below. 

I. Resources (Matlab)
---------------------
1. osfnd1.m - Fitting routine for the standard diffusion model
2. osfvwn1x.m - Fitting routine for the time-varying diffusion model.
3. osfqpf.m - Quantile-probability plot of fit to data
4. osf4.mat - Data for individual subjects

II. Dependencies
-----------------

osfnd1 --> etawiener --> ewiener (all in same file)
osfvwn1x --> vwien400x --> wiensustn400x.c (separate files)

The C-function wiensustn400x.c was compiled under Linux using the GNU scientific library, gsl, and the gslcblas library. These libraries provide the function gsl_sf_gamma_inc_P (the incomplete gamma function). Other libraries may provide this function in other environments. The function must be compiled before use under Matlab with the command:

>>  mex wiensustn400x.c -lgsl -lgslcblas  -lm

III. Data Structures
--------------------

The data structure for an individual subject is a 12 x 14 array. The rows are stimulus conditions; the columns are data for a condition. The first 7 columns are correct responses; the next 7 columns are error responses. The columns are P(R) (response probability); Ni (number of trials) and Q1...Q5, the .1, .3, .5, .7, and .9 RT quantiles in seconds. The rows are 
  Hi F Hard  (Speed)
  Hi F Easy
  Eq F Hard
  Eq F Easy
  Lo F Hard
  Lo F Easy
  Repeat (Accuracy)

Missing data were treated as described in the paper. If there were too few errors in a condition to compute error distribution quantiles, but enough to compute medians (M), then the median was computed and G^2 was computed on two bins. For those conditions, columns 10 to 14 are 0, 0, M, 0, 0. For conditions in which there were too few errors to compute medians, columns 10 to 15 are 0, 0, 0, 0, 0, and G^2 was computed on one bin.  

The file osf4.mat contains the data structures S1, S2, ... S20, and Sg. These are quantile structures for subjects #1 through #20 and quantile-averaged group data for the five subjects for whom there was complete error data. 

IV. Calling Conventions
-----------------------
Help calls: "help osfnd1", "help osfvwn1," etc., give the calling conventions including a list of parameters. For osfnd1 and osfvwn1x the calls are:

 [G2, AIC, BIC, Pred] = osfnd1(Pvar,Pfix, Sel, Dmatrix, {trace});

 [G2, AIC, BIC, Pred] = osfvwn1x(Pvar,Pfix, Sel, Dmatrix, {trace});

These functions take a vector of variable parameters, Pvar, which are estimated during fitting, a vector of fixed parameters, Pfix, which remain fixed, a selector vector, Sel, used to select and assemble parameter vectors, a data matrix, Dmatrix, and an optional trace switch, which gives information about parameter bounds. 

For example, osfvwn1x takes a total of 14 parameters, 

P = [as, aa, vh, ve, c,  pz,  eta, szs, sza, Ter  st, sigma2, beta, m]
      1   2   3  4   5   6     7    8   9    10   11   12    13   14

The time-varying reference model forced four of the parameters (c, szs, stz, and st) to zero. This can be accomplished by defining (say), 

>> P = [0.089 0.22 0.17 0.27 0 -0.02 0.12 0  0 0.19 0 0.06 17.5 4.90];

>> Sel = [1  1  1  1  0  1  1  0  0  1  0  1  1 1];

and calling

>> [g2,a,b,Pred]=osfvwn1x(P(Sel==1), P(Sel==0), Sel, S1)

V. Fitting Sequence
-------------------
The typical fitting sequence in Matlab is as follows. For osfnd1 you need to begin with an 11-element vector of starting parameters, P, and an 11-element selector vector, Sel. For osfvwn1x, P and Sel must both be 14 elements. For Nelder-Mead Simplex fitting:

>> setopt % Creat an options structure and turn on iteration trace
>> pest = fminsearch(@osfnd1, P(Sel==1), options, P(Sel==0), Sel, S1) % Do the fit for subject #1
>> P(Sel==1) = pest % Update the working parameter vector
>> [g2,b,Pred]=osfnd1(P(Sel==1), P(Sel==0), Sel, S1)  % Generate predictions with the new parameters
>> osfqpf(S1, Pred) % Quantile-probability plot of the fit

VI. Miscellaneous Notes and Cautions
------------------------------------

There are a number of context-sensitive flags and variables that are set internally in osfnd1 and osfvwn1x. These include tmax, the maximum time index (currently set to 3.0 s). There is also a combination of hard and soft (quadratically penalized) constraints on the parameters to keep Simplex out of parts of the space where things may go bad numerically. These are ad hoc, but work well enough in practice, so long as you're aware of what the settings are. Setting the trace flag to 1 in the function will give some information about where the current parameter vector is in relation to the constraint set. Both osfnd1 and osfvwn1x save a copy of the working parameter vector to disk each time they are called. These calls add some disk-access overheads but are convenient if you want to exit a fit prematurely or to diagnose where something went bad. Because of the context-sensitive nature of these flags and constraints, these routines should not be treated as black-boxes, but instead as starting points that you should adapt to your own requirements. 

These routines are supplied as is, for nonprofit research purposes only, without any assumed, intended, or implied liability on the part of the author.  










