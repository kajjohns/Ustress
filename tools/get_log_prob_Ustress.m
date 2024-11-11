
function [dhat,returnslip,logDET,logrho,G,stressG]=get_log_prob_Ustress(X,Xsig,Xstress,Xlocked,numdata,...
       InSAR_vectors,invcov,data,datasig,datatype,faults_fixed,xysites,XYsites,G,stressG,geom_switch)

       
if geom_switch==1
    
pm=[];

faults=[faults_fixed(1:2)' X'];

%make fault patch matrix, pm, and Laplacian matrix
%specify components of slip to be calculate ([strike-slip,dip-slip,opening]) -- e.g. [0 1 0] means dip slip only
dis_geom1  = [faults, [1 1 0]];
dis_geom = movefault(dis_geom1);  % move the fault so that the coordinates of the midpoint refer to the
                                             % fault bottom as in Okada
%% Create slip patches
nhe=faults_fixed(3);
nve=faults_fixed(4);
pm1=patchfault(dis_geom(1,1:7),nhe,nve);
pm = [pm; pm1];



%build G matrix
G=[];
for loop=1:length(XYsites)

    if datatype{loop}==1
        dispG=make_dispG_novert(pm,XYsites{loop});
    end
    
    if datatype{loop}==2
        dispG=make_dispG_vert(pm,XYsites{loop});
    end
    
    if datatype{loop}==3
        dispG=make_dispG_insar(pm,XYsites{loop},InSAR_vectors{loop});
    end
    
    G=[G;dispG];

end %loop=1:length(XYsites)


stressG=make_stressG(pm);

end %if geom_switch==1

Xlocked=[logical(Xlocked);logical(Xlocked)];
%compute slip assuming uniform stress drop    
stressGp=stressG;stressGp(Xlocked,:)=[];stressGp(:,Xlocked)=[];
Gc=G;Gc(:,Xlocked)=[];

stressdrop=[Xstress(1)*ones(size(stressGp,1)/2,1);Xstress(2)*ones(size(stressGp,1)/2,1)];
%slip=-inv(stressGp)*(stressdrop)*10^6;  %convert from Pa to MPa
slip=-stressGp\stressdrop*10^6;  %convert from Pa to MPa

dhat=Gc*slip;

returnslip=zeros(size(Xlocked));
index=Xlocked==0;
returnslip(index)=slip;


%multiply data sigmas by weights   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
weight=ones(size(datasig));
weight(1:numdata(1))=Xsig(1); %weight of the first data set 
for k=2:length(numdata)
    weight(1+sum(numdata(1:k-1)):sum(numdata(1:k)))=Xsig(k);  %weight of other data sets   %%should be  Xsig(k)^2  ???
end                                                           %smoothing is Xsig(end)

datasig=weight.*datasig;  %We use the full cov, so the datasig is not used any more. 
                           %Only the weight is kept for weighting each data set.  
                           %Difficult to consider the weight scheme of Yuri
                           %Fialko, 2004 JGR due to the full cov matrix of
                           %SAR data.---Jianbao



logDET=0;
for k=1:length(numdata)
    logDET=logDET + .5*(-numdata(k))*log(Xsig(k)^2);
end
%logrho=-.5*(data./datasig-dhat./datasig)'*(data./datasig-dhat./datasig);
logrho=-.5*(data./datasig-dhat./datasig)'*(data./datasig-dhat./datasig);




