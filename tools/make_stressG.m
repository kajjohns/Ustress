function stressG=make_stressG(pm)

mu=3*10^10;
%convert km to m
pm(:,1:3)=1000*pm(:,1:3);
pm(:,6:7)=1000*pm(:,6:7);

npatches=size(pm,1);  %number of patches

%calculate patch centers
angle=90-pm(:,5)-90;
centerx=pm(:,6)-.5*pm(:,2).*cos(pm(:,4)*pi/180).*cos(angle*pi/180);
centery=pm(:,7)-.5*pm(:,2).*cos(pm(:,4)*pi/180).*sin(angle*pi/180);
centerz=pm(:,3)-.5*pm(:,2).*sin(pm(:,4)*pi/180);



Xcenter=[centerx';centery';-centerz'];



%calculate patch normals
strike=90-pm(:,5);
strikevec=[cos(strike*pi/180) sin(strike*pi/180) zeros(size(strike))];
dipvec=[pm(:,6) pm(:,7) -pm(:,3)]-[centerx centery -centerz];
vecnorm=sqrt(dipvec(:,1).^2+dipvec(:,2).^2+dipvec(:,3).^2);
vecnorm=repmat(vecnorm,1,3);
dipvec=dipvec./vecnorm;
normal=cross(strikevec,dipvec,2);


GssT1=[];GssT2=[];GdsT1=[];GdsT2=[];

%calculate coseismic stress changes
for k=1:npatches
   mss=[pm(k,:) 1 0 0]';
   mds=[pm(k,:) 0 1 0]';
   
   [U,D,Sss,flag]=disloc3d(mss,Xcenter,mu,.25);
   [U,D,Sds,flag]=disloc3d(mds,Xcenter,mu,.25);
   
   T1ss=Sss(1,:)'.*normal(:,1) + Sss(2,:)'.*normal(:,2) + Sss(3,:)'.*normal(:,3);
   T2ss=Sss(2,:)'.*normal(:,1) + Sss(4,:)'.*normal(:,2) + Sss(5,:)'.*normal(:,3);
   T3ss=Sss(3,:)'.*normal(:,1) + Sss(5,:)'.*normal(:,2) + Sss(6,:)'.*normal(:,3);
   
   T1ds=Sds(1,:)'.*normal(:,1) + Sds(2,:)'.*normal(:,2) + Sds(3,:)'.*normal(:,3);
   T2ds=Sds(2,:)'.*normal(:,1) + Sds(4,:)'.*normal(:,2) + Sds(5,:)'.*normal(:,3);
   T3ds=Sds(3,:)'.*normal(:,1) + Sds(5,:)'.*normal(:,2) + Sds(6,:)'.*normal(:,3);
   
   %component along strikevec
   GssTsv=T1ss.*strikevec(:,1)+T2ss.*strikevec(:,2)+T3ss.*strikevec(:,3);
   GdsTsv=T1ds.*strikevec(:,1)+T2ds.*strikevec(:,2)+T3ds.*strikevec(:,3);
   %component along dipvec
   GssTdv=T1ss.*dipvec(:,1)+T2ss.*dipvec(:,2)+T3ss.*dipvec(:,3);
   GdsTdv=T1ds.*dipvec(:,1)+T2ds.*dipvec(:,2)+T3ds.*dipvec(:,3);
   
   GssT1=[GssT1 GssTsv];
   GssT2=[GssT2 GssTdv];
   GdsT1=[GdsT1 GdsTsv];
   GdsT2=[GdsT2 GdsTdv];
   
end

stressG=[GssT1 GdsT1;GssT2 GdsT2];
