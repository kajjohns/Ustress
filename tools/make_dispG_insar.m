function dispG=make_dispG_insar(pm,xystats,insar_look)

%make look-direction vector for inSAR data
insar_look=insar_look*pi/180; %convert to radians
insar_vector(1)=cos(insar_look(2))*cos(pi/2-insar_look(1));
insar_vector(2)=cos(insar_look(2))*sin(pi/2-insar_look(1));
insar_vector(3)=sin(insar_look(2));

%InSAR_vectors: numdata x 3(Ue Un Uu)   Jianbao

xloc=[xystats';zeros(1,size(xystats',2))];

npatches=size(pm,1);

for k=1:npatches
   m1=[pm(k,:) 1 0 0]';   %%SS
   m2=[pm(k,:) 0 1 0]';   %%DS
   
   [U1,D,S,flag]=disloc3d(m1,xloc,1,.25);   
   [U2,D,S,flag]=disloc3d(m2,xloc,1,.25);
   %U1=[east, east, east,... ;
   %    north, north, north,... ;
   %    vertical,vertical, vertical...] (3xn)
%           for j = 1:size(InSAR_vectors,1)
%           tmp1_range(j) =  InSAR_vectors(j,:)*U1(:,j);
%           tmp2_range(j) =  InSAR_vectors(j,:)*U2(:,j);
%           end  
% 	      G1(:,k)=tmp1_range(:);
%           G2(:,k)=tmp2_range(:);   
   G1(:,k)=U1(1,:)'*insar_vector(1) + U1(2,:)'*insar_vector(2) + U1(3,:)'*insar_vector(3);
   G2(:,k)=U2(1,:)'*insar_vector(1) + U2(2,:)'*insar_vector(2) + U2(3,:)'*insar_vector(3);
   
 
end

dispG=[G1 G2];


