%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conduct Monte Carlo Metropolis inversion
% NO INPUTS REQUIRED IN THIS SCRIPT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%number of patches
npatch=faults_fixed(3)*faults_fixed(4);

%make Xsig vector
if ~isempty(Xsig_initial)
    Xsig=Xsig_initial;
else
    Xsig=ones(length(numdata),1);
end

%make Xlocked vector
if ~isempty(Xlocked_initial)
    Xlocked=Xlocked_initial;
else
    Xlocked=zeros(npatch,1);
end

%make Xstress vector
if ~isempty(Xstress_initial)
    Xstress=Xstress_initial;
else
    Xstress=[1;1];
end

%make stepsize_X vector
if isempty(stepsize_X)
    for loop=1:size(faults_initial,1)
        stepsize_X=[stepsize_X; ones(5,1)];
    end
end

%make stepsize_Xsig vector
if isempty(stepsize_Xsig)
    stepsize_Xsig=0.1*ones(length(numdata),1);
end

%make stepsize_Xstress vector
if isempty(stepsize_Xstress)
    stepsize_Xstress=0.1;
end

%make X  (vector of fault geometry parameters)
X=faults_initial;
notfixed=isnan(faults_unknown);
notfixed=logical(notfixed);  %parameters that are not fixed and vary in the inversion

invcov=[];

%%Do first calculation using initial values
tic
G=[];stressG=[]; %initialize empty matrices
geom_switch=1;  %need to compute new pm and G and Lapp
[dhat,slip,logDET,logrho,G,stressG]=get_log_prob_Ustress(X,Xsig,Xstress,Xlocked,numdata,look,invcov,...
       data,datasig,datatype,faults_fixed,xysites,XYsites,G,stressG,geom_switch);
t1=toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%determine computation time

Xlocked_temp=round(rand(npatch,1));  %assume about half of patches are locked
tic
for loop=1:npatch
       geom_switch=0;  %don't compute new pm and G and Lapp
       [dummy1,dummy2,dummy3,dummy4,dummy5,dummy6]=get_log_prob_Ustress(X,Xsig,Xstress,Xlocked_temp,numdata,look,invcov,...
       data,datasig,datatype,faults_fixed,xysites,XYsites,G,stressG,geom_switch);
end
t2=toc;

basetime=t1*sum(notfixed)+t2; %time to compute one sample
disp(['It takes ~' num2str(round(basetime*numsamples/60)) ' minutes (' num2str(round(basetime*numsamples/60/24)) ' hours) to compute ' num2str(numsamples) ' samples.']); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%generate files to store results
fid = fopen(['./MCMC_outputs/M_Ustress_' inversion_name '.txt'],'w'); fclose(fid);  %store model parameters in M.txt
fid = fopen(['./MCMC_outputs/Msig_Ustress_' inversion_name '.txt'],'w'); fclose(fid); 
fid = fopen(['./MCMC_outputs/Mstress_Ustress_' inversion_name '.txt'],'w'); fclose(fid); 
fid = fopen(['./MCMC_outputs/Mlocked_Ustress_' inversion_name '.txt'],'w'); fclose(fid); 
fid = fopen(['./MCMC_outputs/Slip_Ustress_' inversion_name '.txt'],'w'); fclose(fid); 
fid = fopen(['./MCMC_outputs/logrho_Ustress_' inversion_name '.txt'],'w'); fclose(fid);
fid = fopen(['./MCMC_outputs/logDET_Ustress_' inversion_name '.txt'],'w'); fclose(fid);
fid = fopen(['./MCMC_outputs/DHAT_Ustress_' inversion_name '.txt'],'w'); fclose(fid);


fidM = fopen(['./MCMC_outputs/M_Ustress_' inversion_name '.txt'],'a');  %store model parameters in M.txt
fidMsig = fopen(['./MCMC_outputs/Msig_Ustress_' inversion_name '.txt'],'a'); 
fidMlocked = fopen(['./MCMC_outputs/Mlocked_Ustress_' inversion_name '.txt'],'a');
fidMstress = fopen(['./MCMC_outputs/Mstress_Ustress_' inversion_name '.txt'],'a'); 
fidSlip = fopen(['./MCMC_outputs/Slip_Ustress_' inversion_name '.txt'],'a'); 
fidrho = fopen(['./MCMC_outputs/logrho_Ustress_' inversion_name '.txt'],'a'); 
fidDET = fopen(['./MCMC_outputs/logDET_Ustress_' inversion_name '.txt'],'a'); 
fidDHAT = fopen(['./MCMC_outputs/DHAT_Ustress_' inversion_name '.txt'],'a'); 



%intialize some parameters
rand('state', sum(100*clock)); %reset random number generator
slipprev=slip;
Xprev=X;
Xsigprev=Xsig;
Xstressprev=Xstress;
Xlockedprev=Xlocked;
logrhoprev=logrho;
logDETprev=logDET;
dhatprev=dhat;

Xcount=0;
allstepsizeX=[stepsize_X;stepsize_Xsig;stepsize_Xstress];
allnotfixed=[notfixed;ones(size(Xsig));ones(size(Xstress))];

num_fault_params=length(allstepsizeX);

%begin Monte Carlo Metropolis walk
countaccept=0;
for iter=1:numsamples*(npatch+num_fault_params)

    nocompute=0;  %need to compute logprob

    if mod(iter,npatch+num_fault_params)<= num_fault_params & mod(iter,npatch+num_fault_params) ~= 0 

        slip_patch=0;
        geom_switch=1;  %need to compute new pm and G and Lapp
        Xcount=Xcount+1;
        
        allX=[X;Xsig;Xstress];
        %take a random step in model space
        r=(-1).^round(rand(1)).*rand(1);
        r=r.*allstepsizeX(Xcount);
        if allnotfixed(Xcount)==1
            allX(Xcount)=allX(Xcount)+r; %new model parameters to evaluate
        else
            nocompute=1;  %no need to compute logprob
        end
        
        
      
    else
        
        Xcount=0;
        geom_switch=0;  %no need to compute new pm and G and Lapp
        slip_patch=slip_patch+1; %cycle through slip patches -- change one at a time
        Xlocked(slip_patch)=mod(Xlocked(slip_patch)+1,2);

    end
    
X=allX(1:length(X));
Xsig=allX(length(X)+1:length(X)+length(Xsig));
Xstress=allX(length(X)+length(Xsig)+1:length(X)+length(Xsig)+length(Xstress));

        
    
%forward model caclulation


if nocompute==1 | sum(Xsig < 0) >0 | sum(X < faults_lowerbound(:))>0 | sum(X > faults_upperbound(:))>0
    accept=0;
else


[dhat2,slip2,logDET2,logrho2,G,stressG]=get_log_prob_Ustress(X,Xsig,Xstress,Xlocked,numdata,look,invcov,...
    data,datasig,datatype,faults_fixed,xysites,XYsites,G,stressG,geom_switch);

   
%use Metropolis rule to decide whether or not to accept model
accept=metropolis_log(logDET,logrho,logDET2,logrho2);

end


%do not accept negative sigmas or locking depths
%X=[GPSsigscale;smoothing;depth;dip;strike;eastpos;northpos]; %vector of unknowns

    
if accept==1
    
    if geom_switch==1 
        countaccept=countaccept+1;
    end


    %if accept==1, keep the model
   dhat=dhat2; 
   slip=slip2;
   logrho=logrho2;
   logDET=logDET2;
  
   dhatprev=dhat; 
   slipprev=slip;
   Xprev=X;
   Xsigprev=Xsig;
   Xstressprev=Xstress;
   Xlockedprev=Xlocked;
   logrhoprev=logrho;
   logDETprev=logDET;

   

else
    %if accept==0, discard this model and retain previous    
    dhat=dhatprev;
    slip=slipprev;
    X=Xprev;
    Xsig=Xsigprev;
    Xstress=Xstressprev;
    Xlocked=Xlockedprev;
    logrho=logrhoprev;
    logDET=logDETprev;
    
end

%store resuls in text files
if mod(iter,npatch+num_fault_params)==1
        
     sample_number=iter/(npatch+num_fault_params);
   
    disp(['Completed sample number ' num2str(ceil(sample_number))...
        '. Acceptance rate for X, Xsig, Xstress is ' num2str(round(countaccept/(sample_number*length(allX))*100)) ' percent.' ])
          
    fprintf(fidM,'\n',' '); fprintf(fidM,'%6.8f\t',X');
    fprintf(fidMsig,'\n',' '); fprintf(fidMsig,'%6.8f\t',Xsig');
    fprintf(fidMstress,'\n',' '); fprintf(fidMstress,'%6.8f\t',Xstress');
    fprintf(fidMlocked,'\n',' '); fprintf(fidMlocked,'%6.8f\t',Xlocked');
    fprintf(fidSlip,'\n',' '); fprintf(fidSlip,'%6.8f\t',slip');
    fprintf(fidrho,'\n',' '); fprintf(fidrho,'%6.8f\t',logrho);
    fprintf(fidDET,'\n',' '); fprintf(fidDET,'%6.8f\t',logDET);
    fprintf(fidDHAT,'\n',' '); fprintf(fidDHAT,'%6.8f\t',dhat');

end





end %iter
    
%close text files
fprintf(fidM,'\n',' ');
fprintf(fidMsig,'\n',' ');
fprintf(fidMstress,'\n',' ');
fprintf(fidMlocked,'\n',' ');
fprintf(fidSlip,'\n',' ');
fprintf(fidrho,'\n',' ');
fprintf(fidDET,'\n',' ');
fprintf(fidDHAT,'\n',' ');

fclose(fidM);
fclose(fidMsig);
fclose(fidMstress);
fclose(fidMlocked);
fclose(fidSlip);
fclose(fidrho);
fclose(fidDET);
fclose(fidDHAT);

