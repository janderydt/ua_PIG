function Figure1

%addpath(genpath('/home/UNN/wchm8/Documents/Projects/IceOceanCoupling_PIG/JGR/Matlab_Functions'));
addpath(genpath('/Volumes/mainJDeRydt/Documents/Projects/IceOceanCoupling_PTDC/JGR/Matlab_Functions'));
%addpath('/home/UNN/wchm8/Documents/UaMITgcm_v2/Matlab_tools');
addpath(genpath('/Volumes/mainJDeRydt/Documents/Projects/IceOceanCoupling_PTDC/Matlab_tools_UaMITgcm_v2'));
%addpath(genpath('/media/wchm8/data1-JDeRydt/Antarctic_datasets'));
addpath(genpath('/Volumes/bckpJDeRydt/Antarctic_datasets_bckp'));

froot = '/Volumes/mainJDeRydt/Ua/cases/PIG_IceShelf_Integrity/final';

% xmin = -2.5e6; xmax = -0.5e6; 
% ymin = -1.5e6; ymax = 1.7e6;

xmin = -1.8e6; xmax = -1.3e6; 
%ymin = -1.0682e6; ymax = 0.6228e6;
ymin = -5.5e5; ymax = 2.269e5;

Basins = DefineBasins;

%% B-SOSE data
% file_sose = '/media/wchm8/data1-JDeRydt/Antarctic_datasets/B-SOSE/bsose_i105_2008to2012_monthly_Theta.nc';
% lonC_sose = ncread(file_sose,'XC');
% latC_sose = ncread(file_sose,'YC');
% [LONC_sose,LATC_sose] = ndgrid(lonC_sose,latC_sose);
% [XC_sose,YC_sose] = ll2psxy(LATC_sose,LONC_sose,-71,0);
% [Xm,Ym] = ndgrid(linspace(xmin,xmax,500),linspace(ymin,ymax,500));
% 
% depth_sose = ncread(file_sose,'Z');
% I = min(find(depth_sose<=-200));
% 
% Theta_sose = ncread(file_sose,'THETA');
% % find maximum average temperature below -200m
% Theta_sose_av = mean(Theta_sose(:,:,I:end,:),4);
% Theta_sose_av_max = max(Theta_sose_av,[],3);
% FTheta = scatteredInterpolant(XC_sose(:),YC_sose(:),Theta_sose_av_max(:));
% Theta_sose = FTheta(Xm,Ym);

%% offshore bathymetry
% addpath('/media/wchm8/data1-JDeRydt/Antarctic_datasets/BedMachine_Antarctica');
% bathy_BMA = interpBedmachineAntarctica(Xm,Ym,'bed','2019-09-05');
% 
%% Antarctic outline
[x,y,~] = bedmap2_data('surface','xy');
mask = bedmap2_data('icemask');
[C,h] = contour(x,y,mask,[127 127]);
ii=1; kk=1;
while ii<size(C,2)
   step = C(2,ii);
   Cont(kk).step = step;
   Cont(kk).x = C(1,ii+1:ii+step);
   Cont(kk).y = C(2,ii+1:ii+step);
   ii = ii + step + 1;
   kk = kk+1;   
end
[Maxstepsize,I] = max([Cont(:).step]);
ContAnt_x = Cont(I).x(:);
ContAnt_y = Cont(I).y(:);
% 
% I = find(inpoly([Xm(:),Ym(:)],[ContAnt_x(:) ContAnt_y(:)]));
% bathy_BMA(I)=NaN;
% Theta_sose(I)=NaN;

%% Measures 96 velocity
%load('/home/UNN/wchm8/Dropbox (Northumbria University)/PIG_IceShelf_Integrity_v2/InputData/Velocities/ERS_1996.mat');
load([froot,'/InputData/Velocities/ERS_1996.mat']);
X_ers = X/1e3; Y_ers = Y/1e3; speed_ers = sqrt(VX.^2+VY.^2);

%% Ua grid

%load('/home/UNN/wchm8/Dropbox (Northumbria University)/PIG_IceShelf_Integrity_v2/MeshFileAdapt_PIG1996_Ele109300_Nod3');
load([froot,'/MeshFileAdapt_PIG1996_Ele109300_Nod3']);
MUA96 = MUA;
%load('/home/UNN/wchm8/Dropbox (Northumbria University)/PIG_IceShelf_Integrity_v2/MeshFileAdapt_PIG2016_Ele107705_Nod3');
load([froot,'/MeshFileAdapt_PIG2016_Ele107705_Nod3']);
MUA16 = MUA;

MUA.coordinates = MUA96.coordinates/1e3;

CtrlVar.PlotNodes=0;
CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;
CtrlVar.MeshColor = [1 0.5 0.5];
CtrlVar.FEmeshPlotTitle='';
CtrlVar.PlotsXaxisLabel = 'psx [km]';
CtrlVar.PlotsYaxisLabel = 'psy [km]';
CtrlVar.PlotXYscale = 1;
% 
% lonCMIT(mask==0)=NaN;
% latCMIT(mask==0)=NaN;

%% load changes in speed and thickness
gs = '25000';
m = '3';

load([froot,'/Inversion_ERS_1996_REMAthickenedfinal_gs',gs,'_m',m,'_InverseRestartFile.mat']);
MUA_1996 = MUA; GF_1996 = GF;
Fh = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.h);
Fspeed = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),sqrt(F.ub.^2+F.vb.^2));

%% flowline Figure 2
%x0 = -1.584e6; y0 = -1.853e5;
x0 = -1.4552e6; y0 = -3.0792e4;
x = MUA.coordinates(:,1); y = MUA.coordinates(:,2);
dx = 5e2; dy = 5e2;
[Xm,Ym] = meshgrid([min(x):dx:max(x)],[min(y):dy:max(y)]);
Fu = scatteredInterpolant(x,y,F.ub);
Fv = scatteredInterpolant(x,y,F.vb);
um = Fu(Xm,Ym);
vm = Fv(Xm,Ym);
XY = stream2(Xm,Ym,um,vm,x0,y0,[0.1 10000]);

xy=XY{1}; x_fl = xy(:,1); y_fl = xy(:,2);

I = find(~inpoly([x_fl y_fl],[MUA.Boundary.x MUA.Boundary.y]));
x_fl(I)=[]; y_fl(I)=[];
L_fl = cumsum(sqrt((x_fl(2:end)-x_fl(1:end-1)).^2+(y_fl(2:end)-y_fl(1:end-1)).^2));
I = find(L_fl<L_fl(end)-120e3);
x_fl(I)=[]; y_fl(I)=[];

load([froot,'/Inversion_MEaSUREs_2015_2016_Bedmachine_gs',gs,'_m',m,'_InverseRestartFile.mat']);
CtrlVar = CtrlVarInRestartFile; CtrlVar.PlotXYscale = 1e3;
speed_2015_2016 = sqrt(F.ub.^2+F.vb.^2);

dspeed = speed_2015_2016 - Fspeed(MUA.coordinates(:,1),MUA.coordinates(:,2));
Fdspeed = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),dspeed);

Fdh = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.h-Fh(MUA.coordinates(:,1),MUA.coordinates(:,2)));

%% CPOM and Paolo datapoints
addpath(genpath([froot,'/InputData/Observed_dhdt']));
years = {'1997_2001','2002_2006','2007_2011','2012_2016'};
[~,X_CPOM,Y_CPOM,~,X_Paolo,Y_Paolo,~] = Interpolate_dhdt(CtrlVar,MUA,years,'PIG_Paolo_trimmedv2');

%% Plotting
H=fig('units','inches','width',80*12/72.27,'fontsize',16,'font','Helvetica');
%set(H,'visible','off');

subplot('position',[0.05 0.05 0.45 0.92]); hold on;
% 
% ax1 = axes; hold on;
% 

% pcolor(ax1,Xm/1e3,Ym/1e3,Theta_sose); shading flat; alpha(alphaVal);
% caxis(ax1,[min(Theta_sose(:)) max(Theta_sose(:))]);
% CM = othercolor('RdBu6',64); CM = flipdim(CM,1); CM1 = CM(15:end,:);
% 
% [C,h]=contour(ax1,Xm/1e3,Ym/1e3,bathy_BMA,[-4500:500:-500],'-','color',[0.5 0.5 0.5]);
% 
% axis(ax1,'equal'); axis(ax1,'tight');
% xlim(ax1,[xmin xmax]/1e3); ylim(ax1,[ymin ymax]/1e3);

pcolor(X_ers,Y_ers,speed_ers); shading flat;
CM2 = [[linspace(0.75,0,64)';linspace(0,0.75,64)'] linspace(0.75,0,128)' linspace(0.75,1,128)'];

alphaVal = 0.3;
X1_mask = [xmin xmax xmax xmin]/1e3; X2_mask = [MUA96.Boundary.x'/1e3];
Y1_mask = [ymin ymin ymax ymax]/1e3; Y2_mask = [MUA96.Boundary.y'/1e3];
pgon = polyshape({X1_mask,X2_mask},{Y1_mask,Y2_mask});
plot(pgon,'FaceColor','w','FaceAlpha',0.6);

PlotLatLonGrid(1e3, 5, 10, 1, [0.5 0.5 0.5], 0.1, 14);

plot(MUA96.Boundary.x/1e3,MUA96.Boundary.y/1e3,'-k','linewidth',2);

PlotGroundingLines(CtrlVar,MUA_1996,GF_1996,[],[],[],'-k','linewidth',1.5);

BasinShape = shaperead('Basins_Antarctica_v02.shp');

for ii=1:length(BasinShape)
     X = double(BasinShape(ii).X)/1e3;
     Y = double(BasinShape(ii).Y)/1e3;
    plot(X,Y,'color',[0.4 0.4 0.4]); 
end

plot(x_fl/1e3,y_fl/1e3,'-w','linewidth',1);

plot(ContAnt_x/1e3,ContAnt_y/1e3,'-k');

% patch([-2400 -2400 -1650 -1650],[-1400 -1250 -1250 -1400],'w','edgecolor','w');
% PlotLengthScale(-2350,-1370,[0 250 500],30);
xmin_scale = xmin/1e3+270;
ymin_scale = ymin/1e3+20;
dx_scale = 220;
dy_scale = 50;
patch([xmin_scale xmin_scale xmin_scale+dx_scale xmin_scale+dx_scale],...
    [ymin_scale ymin_scale+dy_scale ymin_scale+dy_scale ymin_scale],'w','edgecolor','w');
PlotLengthScale(xmin_scale+10,ymin_scale+5,[0 50 100 150],dy_scale/5);

plot([-1715 -1715 -1475 -1475 -1715],[-370 -110 -110 -370 -370],'--k','linewidth',1.5);

xticks(gca,[]); yticks(gca,[]);
axis('equal'); axis('tight');
xlim([xmin xmax]/1e3); ylim([ymin ymax]/1e3);

%% Give each one its own colormap
% colormap(ax1,CM1); caxis(ax1,[0 2]);
% cb=colorbar(ax1,'Position',[0.27 0.7 0.15 0.015],'Location','northoutside','AxisLocation','in',...
%  'Ticks',[0:0.5:2],'TickLabels',{'0','0.5','1','1.5','2'});
% cb.XColor='k';
% cb.YColor='k';
% cb.TickLength=0.04;
% cb.FontSize=12;
% cb.Label.String = {'Maximum potential temperature B-SOSE'};

colormap(gca,CM2); caxis([0 3000]);
cb=colorbar(gca,'Position',[0.08 0.08 0.15 0.015],'Location','northoutside','AxisLocation','in',...
 'Ticks',[0:1000:3000],'TickLabels',{'0','1000','2000','3000'});
cb.XColor='k';
cb.YColor='k';
cb.TickLength=0.04;
cb.FontSize=16;
cb.Label.String = {'$U_{96}$ [m/yr]'};
cb.Label.Interpreter = 'latex';

%% Change in ice thickness
subplot('position',[0.5 0.09 0.43 0.44]); hold on;

PlotNodalBasedQuantities_JDR(gca,MUA.connectivity,MUA.coordinates,F.h-Fh(MUA.coordinates(:,1),MUA.coordinates(:,2)),CtrlVar);

x = MUA.coordinates(:,1); y = MUA.coordinates(:,2);
[Xm,Ym] = ndgrid([min(x):1e3:max(x)],[min(y):1e3:max(y)]);
dh = Fdh(Xm,Ym);
I = find(~inpoly([Xm(:) Ym(:)],[MUA.Boundary.x MUA.Boundary.y])); dh(I)=NaN;
[Ccont,hcont]=contour(Xm/CtrlVar.PlotXYscale,Ym/CtrlVar.PlotXYscale,dh,[-100 -80 -60 -40 -20 20],'EdgeColor',[0.4 0.4 0.4]);
[Ccont,hcont]=contour(Xm/CtrlVar.PlotXYscale,Ym/CtrlVar.PlotXYscale,dh,[0 0],'EdgeColor',[0 0 0],'linewidth',0.8);
%clabel([],hcont,[0 0],'color',[0.4 0.4 0.4],'fontsize',12);

PlotGroundingLines(CtrlVar,MUA_1996,GF_1996,[],[],[],'-k','linewidth',1.5);
plot(MUA_1996.Boundary.x/CtrlVar.PlotXYscale,MUA_1996.Boundary.y/CtrlVar.PlotXYscale,'--b','linewidth',1);
PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'-b','linewidth',1.5);
plot(MUA.Boundary.x/CtrlVar.PlotXYscale,MUA.Boundary.y/CtrlVar.PlotXYscale,'-k','linewidth',1);

plot(x_fl/1e3,y_fl/1e3,'-w','linewidth',1);

% g(1)=plot(double(X_CPOM(:)/1e3),double(Y_CPOM(:)/1e3),'ok','markersize',2,'color',[0.4 0.4 0.4]);
% g(2)=plot(double(X_Paolo(:)/1e3),double(Y_Paolo(:)/1e3),'.k','markersize',4,'color',[1 1 1]);

caxis(gca,[-100 100]);
CM1 = othercolor('PuRd9',32); CM1 = flipdim(CM1,1);
CM2 = othercolor('Blues9',32); %CM = flipdim(CM,1);
CM = [CM1 ; CM2];
colormap(gca,CM);
cb = colorbar(gca,'Position',[0.62 0.44 0.2 0.015],'Location','northoutside','AxisLocation','out',...
     'Ticks',[-100:50:100],'ticklabels',{'-100','','0','','100'});
cb.XColor='k';
cb.YColor='k';
cb.TickLength=0.04;
cb.FontSize=16;
cb.Label.String = {'$\Delta H = H_{16} - H_{96}$ [m]'};
cb.Label.Interpreter = 'latex';

axis equal;
xlim([-1715 -1475]); ylim([-370 -110]);
ylabel(gca,'psy [km]'); xlabel(gca,'psx [km]');
grid on; box on;

%% change in surface speed
subplot('position',[0.5 0.53 0.43 0.44]); hold on;

PlotNodalBasedQuantities_JDR(gca,MUA.connectivity,MUA.coordinates,dspeed,CtrlVar);

caxis(gca,[0 2e3]);
CM = othercolor('Reds8',32); CM(1,:)=[1 1 1];
colormap(gca,CM);
cb = colorbar(gca,'Position',[0.62 0.885 0.2 0.015],'Location','northoutside','AxisLocation','out',...
     'Ticks',[0:500:2000],'ticklabels',{'0','','1000','','2000'});
cb.XColor='k';
cb.YColor='k';
cb.TickLength=0.04;
cb.FontSize=16;
cb.Label.String = {'$\Delta U = U_{16} - U_{96}$ [m/yr]'};
cb.Label.Interpreter = 'latex';

plot(x_fl/1e3,y_fl/1e3,'-w','linewidth',1);

dspeed = Fdspeed(Xm,Ym);
I = find(~inpoly([Xm(:) Ym(:)],[MUA.Boundary.x MUA.Boundary.y])); dspeed(I)=NaN;
[Ccont,hcont]=contour(Xm/CtrlVar.PlotXYscale,Ym/CtrlVar.PlotXYscale,dspeed,[50 250:250:1500],'EdgeColor','m');
clabel([],hcont,'manual','color','m','fontsize',12);

PlotGroundingLines(CtrlVar,MUA_1996,GF_1996,[],[],[],'-k','linewidth',1.5);
plot(MUA_1996.Boundary.x/CtrlVar.PlotXYscale,MUA_1996.Boundary.y/CtrlVar.PlotXYscale,'--b','linewidth',1);
PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'-b','linewidth',1.5);
plot(MUA.Boundary.x/CtrlVar.PlotXYscale,MUA.Boundary.y/CtrlVar.PlotXYscale,'-k','linewidth',1);

axis equal;
xlim([-1715 -1475]); ylim([-370 -110]);
xlabel(gca,''); ylabel(gca,'psy [km]');
set(gca,'xticklabel',{});
grid on; box on;

pos = get(H,'Position');
set(H,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
fname = './Figure1';
print(H,fname,'-dpng','-r400');

end