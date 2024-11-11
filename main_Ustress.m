%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% main_Ustress.m 
%%
%% Ustress (uniform stress) is an inversion for slip on a single fault plane 
%% assuming a uniform stress drop. The rectangular fault plane is subdivided
%% into many small rectangular slip patches. Patches are assigned a binary
%% parameter that designates the patch as either locked or unlocked. Slip
%% on unlocked patches is determined using a boundary element approach
%% assuming uniform stress drop. The inversion uses a Monte Carlo-Metropolis
%% method to solve for fault geometry, distribution of locked and unlocked
%% patches, strike- and dip- components of shear stress drop, and data
%% weighting parameters.
%%
%% This inversion might work for small earthquakes where the uniform stress
%% drop assumption might be appropriate, but probably would not work for large
%% earthquakes where the uniform stress drop assumption might be fail.
%% 
%% First set up the inversion in Input_file_Ustress.m
%% Then run inversion using this program, main_Ustress.
%%
%% At any time you can terminate the inversion before it is finished by pressing 
%% CNTRL-C. Then run close_files.m to properly close up text files
%% containing outputs from the inversion.
%%
%% To restart the inversion run montecarlo_inversion_restart.m.
%%
%% Kaj M. Johnson, Indiana University, December 2008
%% Contributions by Jianbao Sun, Institute of Geology, China Earthquake Administration,
%% for incorporating InSAR data and errors
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUTS
%
numsamples = 10^4;  %number of Monte Carlo samples

%name of inversion (for appending output file names)
inversion_name = 'synthetic';

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath tools



%specify inversion paramters in Input_file_Ustress.m
Input_file_Ustress

%load data
[data,datasig,xysites,XYsites,numdata]=LoadData(filename,datatype,origin);


%do monte carlo inversion
montecarlo_inversion_Ustress


