setup_nctoolbox
clear all
clc
day_yr=[365 366 365 365 365 365 366 365 365 365 366 365];

 file='input climate_projection_lat_lon_common_grids.txt';
          loc=importdata(file);
          grid_row=loc(:,1);
         
tic      
 lon=zeros(31450,1);
 lat=zeros(31450,1);
counter=0;
         for j=1:170
             for i=1:185
                 counter=counter+1;
                 lon(counter,1)=870+j;
                 lat(counter,1)=780-i;
                
             end
         end
%mod = {'MPI-ESM-MR';'bcc-csm1-1';'CESM1-CAM5';'CESM1-BGC';'CMCC-CMS';'EC-Earth';'CanESM2';'GISSE2R';'GFDL-ESM2M';'HadGEM2-ES';'MPI-ESM-LR';'inmcm4';'HadGEM2-CC';'IPSL-CM5A-LR'};
mod = {'bcc-csm1-1';'CESM1-BGC';'CESM1-CAM5';'CMCC-CMS';'CanESM2';'GISSE2R';'EC-Earth';'GFDL-ESM2M';'HadGEM2-ES';'MPI-ESM-LR';'MPI-ESM-MR';'HadGEM2-CC';'inmcm4';'IPSL-CM5A-LR'};
for mm=1:1
  model=mod{mm}


     count_yr=0;
for year=2000:2000
    count_yr=count_yr+1;
    yyyy=num2str(year);
   

filename1=['RCP4.5/',model,'/',yyyy,'/Extraction_pr.nc'];
filename2=['RCP4.5/',model,'/',yyyy,'/Extraction_tasmax.nc'];
filename3=['RCP4.5/',model,'/',yyyy,'/Extraction_tasmin.nc'];
% nc1=cfdataset(filename1);
% nc2=cfdataset(filename2);
% nc3=cfdataset(filename3);
precip=ncread(filename1,'pr');
latData=ncread(filename2,'lat');
lonData=ncread(filename2,'lon');
temp_max=ncread(filename2,'tasmax');
temp_min=ncread(filename3,'tasmin');

%   for member=11:11
%       for leadtime=3:2:9
          clear time1
          clear time2
          clear time3
clear forcing
% first_pcp1=[day member leadtime 25 86];
% last_pcp1=[day member leadtime 54 130];

%first_pcp2=[day member 3 20 65];
%last_pcp2=[day member 34 40 95];

t1=precip(1:242,1:162,1:60);
t1(isnan(t1))=0;
t1=squeeze(t1);
t2=precip(1:242,1:162,335:365);
t2(isnan(t2))=0;
t2=squeeze(t2);
t-cold=t1+t2;

tmp1=temp_max(1:242,1:162,1:60);
tmp1(isnan(tmp1))=0;
tmp1=squeeze(tmp1);
%t1=t1*9/5+32;
tmp2=temp_min(1:242,1:162,1:60);
tmp2(isnan(tmp2))=0;
tmp2=squeeze(tmp2);
temp-cold=tmp1+tmp2/
%t2=t2*9/5+32;
%tt=nc2.data('Total_precipitation',first_pcp2,last_pcp2);
%LL=nc2.data('intValidTime',[day,3],[day,34]);
% xx=zeros(114955,1);
% yy=zeros(114955,1);
forcing1=zeros(185,170);
forcing2=zeros(185,170);
forcing3=zeros(185,170);
forcing4=zeros(185,170);
% Vq=zeros(415,277);
% Tmax=zeros(415,277);
% Tmin=zeros(415,277);
% V=zeros(114955,1);
% T1=zeros(114955,1);
% T2=zeros(114955,1);
Vq_new=zeros(22736,1);
new=zeros(22736,1);

Input_pcp=zeros(185,170);
[X,Y]=meshgrid(67.03:0.0623:82.09,38.91:0.0623:48.97);
[Xq,Yq]=meshgrid(67.03:0.03636:82.09,38.91:0.03636:48.97);

Vq=interp2(X,Y,t,Xq,Yq);
Tmax=interp2(X,Y,t1,Xq,Yq);
Tmin=interp2(X,Y,t2,Xq,Yq);
xx=reshape(Xq,114955,1);
yy=reshape(Yq,114955,1);
V=reshape(Vq,114955,1);
T1=reshape(Tmax,114955,1);
T2=reshape(Tmin,114955,1);
m=0;
count=0;
for j=1:170
    for i=1:185
        count=count+1;
        m=grid_row(count,1);
        forcing1(i,j)=V(m,1);
        forcing2(i,j)=T1(m,1);
        forcing3(i,j)=T2(m,1);
    end
end
forcing4=(forcing2+forcing3)/2;

% hh1='00';hh2='01';hh3='02',hh4='03';hh5='04';hh6='05';hh7='06';hh8='07';hh9='08';hh10='09';hh11='10';hh12='11';
% hh13='12';hh14='13';hh15='14',hh16='15';hh17='16';hh18='17';hh19='18';hh20='19';hh21='20';hh22='21';hh23='22';hh24='23';

ct=0;
for kk=1:1
    hh=ct
%    clear forcing4
%    if hh<8 
%forcing4=forcing3;
%elseif hh>=8 && hh<13
%forcing4=forcing3+(hh-8)*(forcing2-forcing3)/12;
%elseif hh>=13 && hh<16
%forcing4=forcing3+(hh-12)*(forcing2-forcing3)/6;
%elseif hh>=16 && hh<21
%forcing4=forcing2-(hh-16)*(forcing2-forcing3)/3
%else
%forcing4=forcing3
%end
    if hh<10 
        hh1=['0',num2str(hh)];
    else hh1=num2str(ct);
    end
    
    
%pcp=interp2(hrap_x,hrap_y,forcing,lon,lat)
 output_file1=['RCP4.5/',model,'/',yyyy,'/precip/xmrg',mm,dd,yyyy,hh1,'z.asc'];
 fid=fopen(output_file1,'wt');
 fprintf(fid,'ncols 170  \n');
 fprintf(fid,'nrows 185  \n');
 fprintf(fid,'xllcorner 870.000000  \n');
 fprintf(fid,'yllcorner 595.000000  \n');
 fprintf(fid,'cellsize 1.000000 \n');
 fprintf(fid,'NODATA_value -1.000000  \n');
 dlmwrite(output_file1,forcing1,'delimiter',' ','precision','%e','-append');
 fclose(fid);
 

output_file4=['RCP4.5/',model,'/',yyyy,'/temp_new/tair',mm,dd,yyyy,hh1,'z.asc'];
fid=fopen(output_file4,'wt');
fprintf(fid,'ncols 170  \n');
fprintf(fid,'nrows 185  \n');
fprintf(fid,'xllcorner 870.000000  \n');
fprintf(fid,'yllcorner 595.000000  \n');
fprintf(fid,'cellsize 1.000000 \n');
fprintf(fid,'NODATA_value -1.000000  \n');
dlmwrite(output_file4,forcing4,'delimiter',' ','precision','%e','-append');
fclose(fid);


 system(['./asctoxmrg ',output_file1]);
 system(['rm ',output_file1]);
% system(['./asctoxmrg ',output_file2]);
% system(['./asctoxmrg ',output_file3]);
 system(['./asctoxmrg ',output_file4]);
 system(['rm ',output_file4]);
ct=ct+1;
end
end
end

