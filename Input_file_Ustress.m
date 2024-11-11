
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs for Ustress
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Note on data file formats:
%Data should be in .txt files with the following formats:

%for horizontal GPS data
%columns: long, lat, east displacement, north displacement, east standard dev, north standard dev

%for vertical GPS data/leveling data
%columns: long, lat, vertical displacement, standard deviation

%for InSAR data
%columns: long, lat, LOS displacement, standard deviation (if not known, column of ones)

%NOTE: displacements need to be in meters!

%for each data set, enter the following information. You can add as many
%data sets as you have by adding more lines like those below.
%data types
%1 = horizontal GPS
%2 = vertical GPS or leveling
%3 = InSAR

%data set 1
filename{1}='./data/GPS_synthetic.txt';
datatype{1}=1; 
dataname{1}='GPS';
look{1}=nan;  %look direction for InSAR data -- not applicable for GPS data
covname=[];  %no InSAR covariance matrix

%data set 2
filename{2}='./data/leveling_synthetic.txt';
datatype{2}=2; 
dataname{2}='LEV';
look{2}=nan;  %look direction for InSAR data -- not applicable for leveling data
covname=[];  %no InSAR covariance matrix

%data set 2
filename{3}='./data/insar_synthetic.txt';
datatype{3}=3; 
dataname{3}='INSAR';
look{3}=[120 75];  %look direction for InSAR data 
covname=[];  %no InSAR covariance matrix

%origin of coordinate system (lat,long)
origin=[36 -121];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fixed fault parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fixed fault parameters
%(length (km), width (km), number of elements along strike, number of elements along dip)
faults_fixed=[37.5 30 3 3]';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify parameters as fixed (known) or free (unknown)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for each parameter, specify a value if you want the parameter fixed,
%otherwise specify 'nan' and the inversion will solve for this parameter

%fault geometry parameters
%faults_unknown=[*depth, dip(degrees), strike(degrees), *east offset, *north offset]
% *--refers to center of top edge
%faults_unknown=[nan nan nan  nan  nan]';  %must be column vector
faults_unknown=[nan nan nan nan nan]'; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify upper and lower bounds on unknown parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%specify 'inf' and '-inf' if no bounds
faults_upperbound=[inf 80 45 inf inf]';
faults_lowerbound=[0 10 -70 -inf -inf]';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify initial parameters to serve as starting point in inversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Note: results should not be sensitive to initial conditions, but the inversion
%will be more efficient if initial conditions are similar to final values

%fault_initial=[*depth, dip(degrees), strike(degrees), *east offset, *north offset]
% *--refers to center of top edge
faults_initial=[5 30 0  0  0]'; 
%note: initial values of fixed parameters must be equal to the fixed value

%optional starting values for data weights
%leave empty, [], to use default values
Xsig_initial=[3;2;10];  %IMPORTANT NOTE: must be a column vector

%optional starting value for locking
%leave empty, [], to use default values
Xlocked_initial=[]; %IMPORTANT NOTE: must be a column vector

%optional starting value for shear stress components
%leave empty, [], to use default values
%Xstress =  [strike-slip dip-slip]
%positive is left-lateral and reverse sense of shear
Xstress_initial=[-2 5]'; %IMPORTANT NOTE: must be a column vector

%optional stepsize vector for X
%leave empty, [], to use default values
stepsize_X=2*[1 1 1 1 1]'; %IMPORTANT NOTE: must be a column vector

%optional stepsize vector for Xsig
%leave empty, [], to use default values
stepsize_Xsig=[.5;.5;.5]; %IMPORTANT NOTE: must be a column vector

%optional stepsize vector for Xstress
%leave empty, [], to use default values
stepsize_Xstress=[.2;.5]; %IMPORTANT NOTE: must be a column vector

