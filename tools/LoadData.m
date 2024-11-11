function [data,datasig,xysites,XYsites,numdata]=LoadData(filename,datatype,origin)


%loop over number of data files, load data, and construct data and datasig
%vectors
data=[]; datasig=[]; xysites=[];
for loop=1:length(filename)
    d=load(filename{loop});
    
    if datatype{loop}==1  %horizontal GPS; format is e,n,e,n,...
        
        dtemp=zeros(2*size(d,1),1);
        dtemp(1:2:end)=d(:,3);  %east component 
        dtemp(2:2:end)=d(:,4);  %north component
        
        sigtemp=zeros(2*size(d,1),1);
        sigtemp(1:2:end)=d(:,5); %east component
        sigtemp(2:2:end)=d(:,6); %north component
        
        
        data=[data; dtemp];
        datasig=[datasig; sigtemp];
      
        %coordinates -- convert to local xy coordinates
        llhstats = d(:,1:2);
        xysites = [xysites; llh2local(llhstats', fliplr(origin))'];
        
        numdata(loop)=length(dtemp);
        XYsites{loop}=llh2local(llhstats', fliplr(origin))';
        
    end
    
    if datatype{loop}==2  %vertical GPS or leveling 
        
        data=[data; d(:,3)];
        datasig=[datasig; d(:,4)];
        
        %coordinates -- convert to local xy coordinates
        llhstats = d(:,1:2);
        xysites = [xysites; llh2local(llhstats', fliplr(origin))'];
       
         numdata(loop)=size(d,1);
        XYsites{loop}=llh2local(llhstats', fliplr(origin))';
        
    end
    
     if datatype{loop}==3  %InSAR 

        data=[data; d(:,3)];  %The data is in *METER* unit.
%        data=[data; -d(:,5)];  %The origin data is in radar convention, Jianbao
        
        %The full covariance matrix is used or diagonal mtrix is used, see
        %Input_file.m 
        if size(d,2)>3
            datasig=[datasig; d(:,4)];
        else
            datasig=[datasig; 0.01*ones(size(d(:,3)))];  %diagonal matrix only, assume 1 cm error !!!!!
        end    
       
        %coordinates -- convert to local xy coordinates
         llhstats = d(:,1:2);
        xysites = [xysites; llh2local(llhstats', fliplr(origin))'];
        
         numdata(loop)=size(d,1);
        XYsites{loop}=llh2local(llhstats', fliplr(origin))';
        
       
     end    
end

   
    
    