function dispG=make_dispG(pm,xystats)

xloc=[xystats';zeros(1,size(xystats',2))];

npatches=size(pm,1);

for k=1:npatches
   m1=[pm(k,:) 1 0 0]';
   m2=[pm(k,:) 0 1 0]';
   
   [U1,D,S,flag]=disloc3d(m1,xloc,1,.25);
   [U2,D,S,flag]=disloc3d(m2,xloc,1,.25);
   
   G1(:,k)=U1(:);
   G2(:,k)=U2(:);
   
 
end

dispG=[G1 G2];