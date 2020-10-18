% Compositing CONUS precip on nino 3.4
% The Nino 3.4 index was downloaded from a NOAA site, and reformatted into a .mat file for easy use
close all
clear all

%for correlation pnorm=1, for regression pnorm=0 *looped*

meanplot=0  % if 1 plot the monthly climatology
devplot = 0  % if 1 plot the monthly standard deviation
tstatplot = 1  %  if 1 plot the t-statistic
onlytstat = 1   %  if 1 skip every plot but tstat plot

coast = load('coast.mat');  % load coastline data

alphaup=1.0-0.01/2    ;  % 99% T-Statistic two-tailed
tcrit=tinv(alphaup,64)  % Critical t-statistic are there 64 years assumed independent
nn=zeros(2,12);
   
%  Load Precip data
filenme = 'precip.V1.0.mon.mean.nc'
lat = ncread(filenme,'lat');
lon = ncread(filenme,'lon');
lat=double(lat);
lon=double(lon);
latx=length(lat);
lonx=length(lon);
[latgrid, longrid] = meshgrid(lat, lon);
p = ncread(filenme,'precip');
p=p(:,:,1:792);
p(p<0.0)=NaN;   % Matlab plotting ignores NaN, so I set missing data to NaN
p=reshape(p,lonx*latx,12,66);    % 66 years of precip data
pm=double(squeeze(nanmean(p,3)));    % compute monthly climatology 
pd=double(squeeze(nanstd(p,0,3)));   % compute standard deviation of monthly means

% seasons
spring_p=p(:,3,:)+p(:,4,:)+p(:,5,:);
summer_p=p(:,[6,7,8],:);
summer_p=sum(summer_p,2);
autumn_p=p(:,[9,10,11],:);
autumn_p= sum(autumn_p,2);
winter_p=p(:,[12,1,2],:);
winter_p=sum(winter_p,2);
  
% set a seasonal p
pseason(:,1,:)=spring_p;
pseason(:,2,:)=summer_p;
pseason(:,3,:)=autumn_p;
pseason(:,4,:)=winter_p;
p=pseason;

% read nino3.4 data
% this loads month and nino34, which are 816 long
% the precip data are only 792 long, so we need to cut off some of the nino34 and month
load('nino34s_1948.mat');
nino34=nino34(1:792);   % just include the first 66 years of Nino index
month=month(1:792);

nino34=reshape(nino34,12,66);
% seasons
spring_34=nino34([3,4,5],:);
spring_34=sum(spring_34,1);
summer_34=nino34([6,7,8],:);
summer_34=sum(summer_34,1);

autumn_34=nino34([9,10,11],:);
autumn_34=sum(autumn_34,1);

winter_34=nino34([12,1,2],:);
winter_34=sum(winter_34,1);

% create new nino34 --seasonal
nino34=reshape(nino34,12,66);
ninoseason(1,:)=spring_34;
ninoseason(2,:)=summer_34;
ninoseason(3,:)=autumn_34;
ninoseason(4,:)=winter_34;
nino34=ninoseason;

% normalize
nino34=reshape(nino34,264,1);
ninmean=squeeze(nanmean(nino34));
nindev=squeeze(nanstd(nino34,0));
ninnd=(nino34-ninmean)/nindev; 

% remove mean climatology and divide by standard deviation of each quarter
% to facilitate calculation of correlation.
for j=1:66
	for i=1:4
		p(:,i,j)=p(:,i,j)-pm(:,i);
                pn(:,i,j)=p(:,i,j)./pd(:,i); 
end  % quarter loop
end  % year loop

ninnd=reshape(ninnd,4,66);  %reshape for regression by quarter

tfac=sqrt(64);

% Now, I think we can compute correlation maps
  for i=1:4
	  cor(:,i)=squeeze(ninnd(i,:))*squeeze(pn(:,i,:))'/66.;
	  reg(:,i)=squeeze(ninnd(i,:))*squeeze(p(:,i,:))'/66.; 
  tstat(:,i)=tfac*cor(:,i)./sqrt((1.0-cor(:,i).*cor(:,i)));
end  % end of quarter loop
tstatll = reshape(tstat,lonx,latx,4);

% plot correlation
for i=1:4  % Loop on seasons   * * * * * * * * * * * * * * * * * 
corr=double(reshape(cor(:,i),lonx,latx));
regr=double(reshape(reg(:,i),lonx,latx));
istat=zeros(lonx*latx,1);
istat(abs(tstat(:,i))>tcrit)=1;
maxem=squeeze(max(max(abs(corr))))

for pnorm=1:2    %  loop on correlaton or regression * * * * * * * * 
if pnorm ==1
maxe=0.80;
else
maxe=maxem;
end  % end of contour choice
maxe10=maxe/10.;
conts=[-maxe:maxe10:maxe];

contt=[-tcrit tcrit];
pslat1=25.0; pslat2=50.0; 


 figure;
  ax =axesm('MapProjection','robinson','Grid','on','FLineWidth',1.0,'MapLatLimit',[pslat1 pslat2],'MapLonLimit',[235 295]);
[cmap]=buildcmap('rwb');
colormap(cmap);
if pnorm==1
[C, H] = contourfm(lat,lon,corr',conts,'LineWidth',0.01,'LineStyle','none');
else
[C, H] = contourfm(lat,lon,regr',conts,'LineWidth',0.01,'LineStyle','none');
end
% the following line puts in a black line where tstat = tcrit
[C, H] = contourm(lat,lon,double(squeeze(real(tstatll(:,:,i))))',contt,'LineWidth',2.0,'LineColor','black');
caxis([-maxe maxe]);
hold on
% put stipples in
ifstip = 0
if ifstip == 1
for ii=1:lonx*latx
	 if istat(ii)>0.8;
   geoshow(latgrid(ii), longrid(ii), 'DisplayType','point','Marker','o',...
'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',0.5);
end  % end of stipple if
end %of grid loop
end  % of stipple if big

plotm(coast.lat,coast.long,'Color',[0.5 0.5 0.5],'LineWidth',0.5);
%geoshow(coast.lat,coast.long,'DisplayType','polygon','FaceColor','black');
if pnorm==1;
title(['Precip Nino 3.4 Correlation ',' season=',num2str(i),', max=',num2str(maxem)],'FontSize',8);
 else
title(['Precip Nino 3.4 Regression 1 std ',' season=',num2str(i),', max=',num2str(maxem)],'FontSize',8);
end
%set(gca,'FontSize',12);
colorbar;
set(gca,'Color','none');

hold off
end % end of pnorm loop
end  % end of months loop

allskip=1;
if allskip == 0

%First do warm events
	    pluscrit=1.0;    % here we choose anomaly index of pluscrit as the cutoff
	    negcrit=-1.0;
%  Here we make a mask to pick out the data that pass the criteria
ninneg=nino34;  %set up the mask arrays
ninpos=nino34;
ninpos(ninpos>pluscrit)=1.0;
ninpos(ninpos<pluscrit)=0.0;
ninneg(ninneg>negcrit)=0.0;
ninneg(ninneg<negcrit)=1.0;

% Make El Nino composite
ppos=reshape(p, lonx*latx, 792);
for i=1:792
	ppos(:,i)=ppos(:,i).*ninpos(i);  %any months not chosen are set to zero
end
ppos=reshape(ppos,lonx*latx,12,66);

% Start over for La Nina Composite
pneg=reshape(p, lonx*latx, 792);
for i=1:792
	pneg(:,i)=pneg(:,i).*ninneg(i);
end
pneg=reshape(pneg,lonx*latx,12,66);

%  OK, now loop on months   ************************************************
  for i=1:12  % this loop should go to end  (1)
	  pposa=double(sum(squeeze(ppos(:,i,:)),2));  %add up the chosen months
	  pposd=squeeze(ppos(:,i,:)).*squeeze(ppos(:,i,:)); %square chosen months
    ninpos=reshape(ninpos,12,66);
    ninpoa=sum(ninpos(i,:),2);    % count the number of months chosen
    nn(1,i)=ninpoa;  %save number chosen into a 2x12 array
    pposd=double(squeeze((sum(pposd,2))/(ninpoa)));  %this computes average of squares
    pposa=pposa./ninpoa;  %divide sum by sample size to get mean
    pposd=pposd-pposa.*pposa;  %  subtract square of average from sum of squares to get variance
	pposd=pposd*ninpoa/(ninpoa-1);  %unbiased estimate of variance
pposa=reshape(pposa,lonx,latx);


% do it all again for the negative composite
	  pnega=double(sum(squeeze(pneg(:,i,:)),2));
pnegd=squeeze(pneg(:,i,:)).*squeeze(pneg(:,i,:));

ninneg=reshape(ninneg,12,66);
ninnea=sum(ninneg(i,:),2);
nn(2,i)=ninnea;
pnegd=double(squeeze(sum(pnegd,2))/(ninnea));
pnega=pnega./ninnea;
	     pnegd=pnegd-pnega.*pnega; %sum of squares minus square of sums
	     pnegd=pnegd*ninnea/(ninnea-1);  % unbiased variance
pnega=reshape(pnega,lonx,latx);
if ninpoa>5   % if either of the samples is less than 6, don't bother (2)
if ninnea>5   %   (3)

%  Here we compute the t-statistic using (1.35) in the chapter 1 notes
pdifa=pposa-pnega;  % raw difference
%tsig=sqrt((pposd/ninpoa+pnegd/ninnea));
tsig=sqrt((pposd*ninpoa+pnegd*ninnea));
  factor=sqrt(((1.0/ninpoa+1.0/ninnea)./(ninpoa+ninnea-2)));
tsig=tsig*factor;
tsig=reshape(tsig,lonx,latx);
tstat=pdifa./tsig;  %final t-statisic array

dof=ninpoa+ninnea-2;  %degrees of freedom
tcrit= tinv(alphaup,dof);   
%tinv is a function in Matlab to compute the crtical t-statistic

%   OK the rest is just plotting  this is line 96   ***************************

if tstatplot == 1  
% plot t statistic
maxe=squeeze(max(max(abs(tstat))));
maxe10=maxe/10.;
conts=[-maxe:maxe10:maxe];
contt=[-tcrit tcrit];
		   pslat1=25.0; pslat2=50.0;
coast = load('coast.mat');


 figure
  ax =axesm('MapProjection','robinson','Grid','on','FLineWidth',1.0,'MapLatLimit',[pslat1 pslat2],'MapLonLimit',[235 295]);
[cmap]=buildcmap('rwb');
colormap(cmap);
[C, H] = contourfm(lat,lon,tstat',conts,'LineWidth',0.01,'LineStyle','none');
[C, H] = contourm(lat,lon,tstat',contt,'LineWidth',2.0,'LineColor','black');
%[C, H] = contourm(lat,lon,pmc(:,:,i)','LineWidth',0.5);
caxis([-maxe maxe]);
hold on
plotm(coast.lat,coast.long,'Color',[0.5 0.5 0.5],'LineWidth',0.5);
%geoshow(coast.lat,coast.long,'DisplayType','polygon','FaceColor','black');
title(['Precip, t-statistic (1.35) ',num2str(pluscrit),' month=',num2str(i),	  ' nn=',num2str(nn(1,i)),' - ',num2str(nn(2,i)),' dof= ',num2str(dof),' tcrit=',num2str(tcrit)  ...
,' tmx=',num2str(maxe)])
set(gca,'FontSize',12);
colorbar
set(gca,'Color','none');
saveas(gca,['PrecipTstat_',num2str(pluscrit),'-',num2str(i),'.epsc']);
export_fig ['PrecipTstat_',num2str(pluscrit),'-',num2str(i),'.epsc'] -transparent

end  % tstatplot if loop end  (end 4)

  if onlytstat==0   %  (5)

% First plot El Nino
maxe=squeeze(max(max(max(pposa))));
maxe10=maxe/15.;
conts=[maxe10:maxe10:maxe];

		   pslat1=25.0; pslat2=50.0;
coast = load('coast.mat');


 figure
  ax =axesm('MapProjection','robinson','Grid','on','FLineWidth',1.0,'MapLatLimit',[pslat1 pslat2],'MapLonLimit',[235 295]);
[cmap]=buildcmap('rwb');
colormap(cmap);
[C, H] = contourfm(lat,lon,pposa',conts,'LineWidth',0.5,'LineStyle','none');
%[C, H] = contourm(lat,lon,pmc(:,:,i)','LineWidth',0.5);
caxis([0 maxe]);
hold on
plotm(coast.lat,coast.long,'Color',[0.5 0.5 0.5],'LineWidth',0.5);
%geoshow(coast.lat,coast.long,'DisplayType','polygon','FaceColor','black');
title(['GPCP Precip, El Nino ',num2str(pluscrit),' month=',num2str(i),' nn=',num2str(nn(1,i))])
set(gca,'FontSize',12);
colorbar
set(gca,'Color','none');
saveas(gca,['PrecipNino_',num2str(pluscrit),num2str(i),'.epsc']);
export_fig ['PrecipNino_',num2str(pluscrit),num2str(i),'.epsc'] -transparent

% Next plot La Nina  ***********************************************
maxe=squeeze(max(max(max(pnega))));
maxe10=maxe/15.;
conts=[maxe10:maxe10:maxe];

pslat1=25.0; pslat2=50.0;

figure
 ax =axesm('MapProjection','robinson','Grid','on','FLineWidth',1.0,'MapLatLimit',[pslat1 pslat2],'MapLonLimit',[235 295]);
[cmap]=buildcmap('rwb');
colormap(cmap);
[C, H] = contourfm(lat,lon,pnega',conts,'LineWidth',0.5,'LineStyle','none');
%[C, H] = contourm(lat,lon,pmc(:,:,i)','LineWidth',0.5);
caxis([0 maxe]);
hold on
plotm(coast.lat,coast.long,'Color',[0.5 0.5 0.5],'LineWidth',0.5);
%geoshow(coast.lat,coast.long,'DisplayType','polygon','FaceColor','black');
title(['GPCP Precip, Warm ',num2str(pluscrit),' month=',num2str(i)])
title(['GPCP Precip, La Nina ',num2str(pluscrit),' month=',num2str(i),' nn=',num2str(nn(2,i))])
set(gca,'FontSize',12);
colorbar
set(gca,'Color','none');
saveas(gca,['PrecipNina_',num2str(pluscrit),num2str(i),'.epsc']);
export_fig ['PrecipNina_',num2str(pluscrit),num2str(i),'.epsc'] -transparent

% Finally plot El Nino minus La Nina  ***************************************
maxe=squeeze(max(max(max(pdifa))));
maxe10=maxe/15.;
conts=[-maxe:maxe10:maxe];

pslat1=25.0; pslat2=50.0;
coast = load('coast.mat');


 figure
  ax =axesm('MapProjection','robinson','Grid','on','FLineWidth',1.0,'MapLatLimit',[pslat1 pslat2],'MapLonLimit',[235 295]);
[cmap]=buildcmap('rwb');
colormap(cmap);
[C, H] = contourfm(lat,lon,pdifa',conts,'LineWidth',0.5,'LineStyle','none');
%[C, H] = contourm(lat,lon,pma(:,:,i)','LineWidth',0.5);
caxis([-maxe maxe]);
hold on
plotm(coast.lat,coast.long,'Color',[0.5 0.5 0.5],'LineWidth',0.5);
%geoshow(coast.lat,coast.long,'DisplayType','polygon','FaceColor','black');
title(['GPCP Precip, Difference ',num2str(pluscrit),' month=',num2str(i)])
title(['GPCP Precip, El Nino minus La Nina ',num2str(pluscrit),' month=', ...
	    num2str(i),' nn=',num2str(nn(1,i)),' - ',num2str(nn(2,i))])
set(gca,'FontSize',12);
colorbar
set(gca,'Color','none');
saveas(gca,['PrecipDif_',num2str(pluscrit),num2str(i),'.epsc']);
export_fig ['PrecipDif_',num2str(pluscrit),num2str(i),'.epsc'] -transparent

end % end of ninpoa test
end % end of ninnea test


if devplot == 1
%  Next plot standard deviations
maxe=squeeze(max(max(max(pmd))));
maxe10=maxe/15.;
conts=[maxe10:maxe10:maxe];

 figure
  ax =axesm('MapProjection','robinson','Grid','on','FLineWidth',1.0,'MapLatLimit',[pslat1 pslat2],'MapLonLimit',[235 295]);
[cmap]=buildcmap('rwb');
colormap(cmap);
[C, H] = contourfm(lat,lon,pmd(:,:,i)',conts,'LineWidth',0.5,'LineStyle','none');
%[C, H] = contourm(lat,lon,pmc(:,:,i)','LineWidth',0.5);
caxis([0 maxe]);
hold on
plotm(coast.lat,coast.long,'Color',[0.5 0.5 0.5],'LineWidth',0.5);
%geoshow(coast.lat,coast.long,'DisplayType','polygon','FaceColor','black');
title(['GPCP Precip Deviation, month ',num2str(i)])
set(gca,'FontSize',12);
colorbar
set(gca,'Color','none');
saveas(gca,['PrecipSTD_',num2str(i),'.epsc']);
export_fig ['PrecipSTD_',num2str(i),'.epsc'] -transparent

end   % end of devplot if
end   % end of onlytstat if
end   % end of month loop

end  %  of allskip if


