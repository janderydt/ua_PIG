function FigureA1

froot = '/Volumes/mainJDeRydt/';

%addpath(genpath('/home/UNN/wchm8/Documents/Projects/IceOceanCoupling_PIG/JGR/Matlab_Functions'));
addpath(genpath([froot,'Documents/Projects/IceOceanCoupling_PTDC/JGR/Matlab_Functions']));
%addpath('/home/UNN/wchm8/Documents/UaMITgcm_v2/Matlab_tools');
addpath(genpath([froot,'Documents/Projects/IceOceanCoupling_PTDC/Matlab_tools_UaMITgcm_v2']));
%addpath(genpath('/media/wchm8/data1-JDeRydt/Antarctic_datasets'));
addpath(genpath([froot,'Antarctic_datasets/Bedmachine_Antarctica']));

froot = [froot,'Ua/cases/PIG_IceShelf_Integrity/final'];

% xmin = -2.5e6; xmax = -0.5e6; 
% ymin = -1.5e6; ymax = 1.7e6;

xmin = -1.7e6; xmax = -1.55e6; 
%ymin = -1.0682e6; ymax = 0.6228e6;
ymin = -4e5; ymax = -2e5;

Basins = DefineBasins;

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

load([froot,'/Inversion_MEaSUREs_2015_2016_Bedmachine_gs',gs,'_m',m,'_InverseRestartFile.mat']);
CtrlVar = CtrlVarInRestartFile; CtrlVar.PlotXYscale = 1e3;
CtrlVar.PlotGLs = 0;

Fdh = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.h-Fh(MUA.coordinates(:,1),MUA.coordinates(:,2)));

%% Generate uncorrected GL location
addpath(genpath([froot,'/InputData']));
years = {'1997_2001','2002_2006','2007_2011','2012_2016'};
dhdt = Interpolate_dhdt(CtrlVar,MUA_1996,years,'PIG_Paolo_trimmedv2');
fprintf('Loading Bedmachine B, s, h \n');
[B,s,h,~] = generateBedmachineInterpolants(MUA_1996);
S = 0*s;
h = h - dhdt; % this is the thickness used in the 1996 inversions
b = s - h;
[b,s,h,GF_1996_uncorrected]=Calc_bs_From_hBS(CtrlVar,MUA_1996,h,S,B,917,1027);
[GF_1996_uncorrected,~,~,~]=IceSheetIceShelves(CtrlVar,MUA_1996,GF_1996_uncorrected);
[xGL_1996_uncorrected,yGL_1996_uncorrected,~]=PlotGroundingLines(CtrlVar,MUA_1996,GF_1996,[],[],[]);

%% CPOM and Paolo datapoints
addpath(genpath([froot,'/InputData/Observed_dhdt']));
years = {'1997_2001','2002_2006','2007_2011','2012_2016'};
[~,X_CPOM,Y_CPOM,~,X_Paolo,Y_Paolo,~] = Interpolate_dhdt(CtrlVar,MUA,years,'PIG_Paolo_trimmedv2');

%% Plotting
H=fig('units','inches','width',80*12/72.27,'fontsize',16,'font','Helvetica');
%set(H,'visible','off');

%% Change in ice thickness
subplot('position',[0.1 0.1 0.85 0.85]); hold on;

PlotNodalBasedQuantities_JDR(gca,MUA.connectivity,MUA.coordinates,F.h-Fh(MUA.coordinates(:,1),MUA.coordinates(:,2)),CtrlVar);

%clabel([],hcont,[0 0],'color',[0.4 0.4 0.4],'fontsize',12);
%% 2016 GL
[xGL,yGL,~] = PlotGroundingLines(CtrlVar,MUA,GF,[],[],[]);
g(1) = plot(xGL/1e3,yGL/1e3,'-k','linewidth',1);

%% Plot Insar GL 1996
S = shaperead('InSAR_GL_Antarctica_v02.shp');
for ii=1:length(S)
   date=S(ii).DATE1;
   year = date(1:4);
   [x,y] = ll2psxy(S(ii).Y,S(ii).X,-71,0);
   if ismember(str2double(year),[1992:1996])
        gInsar=plot(x/1e3,y/1e3,'-m','linewidth',1);
   end
end

g(2) = gInsar;

%% Uncorrected model Grounding line
g(3) = plot(xGL_1996_uncorrected/1e3,yGL_1996_uncorrected/1e3,'-b','linewidth',2);

plot(MUA_1996.Boundary.x/CtrlVar.PlotXYscale,MUA_1996.Boundary.y/CtrlVar.PlotXYscale,'--k','linewidth',1);
%PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'-b','linewidth',1.5);
%plot(MUA.Boundary.x/CtrlVar.PlotXYscale,MUA.Boundary.y/CtrlVar.PlotXYscale,'-b','linewidth',1);

g(4)=plot(double(X_CPOM(:)/1e3),double(Y_CPOM(:)/1e3),'ok','markersize',4,'color','k');
g(5)=plot(double(X_Paolo(:)/1e3),double(Y_Paolo(:)/1e3),'.k','markersize',10,'color','k');

caxis(gca,[-100 100]);
CM1 = gray(64);
CM2 = othercolor('Blues8',64);
CM = [CM1(33:64,:);CM2(1:32,:)];
colormap(gca,CM);
cb = colorbar(gca,'Position',[0.35 0.18 0.2 0.015],'Location','northoutside','AxisLocation','out',...
     'Ticks',[-100:50:100],'ticklabels',{'-100','','0','','100'});
cb.XColor='k';
cb.YColor='k';
cb.TickLength=0.04;
cb.FontSize=16;
cb.Label.String = {'$\Delta H = H_{16} - H_{96}$ [m]'};
cb.Label.Interpreter = 'latex';

legend(g([2,3,1,4,5]),{'1992-1996 DInSAR grounding line [Rignot et al., 2014]',...
    '1996 model grounding line','2016 DInSAR and model grounding line',...
    'dh/dt data points CPOM [Sheperd et al., 2016]','dh/dt data points ice shelf'});

axis equal;
xlim([xmin xmax]/1e3); ylim([ymin ymax]/1e3);
ylabel(gca,'psy [km]'); xlabel(gca,'psx [km]');
grid on; box on;

pos = get(H,'Position');
set(H,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
fname = './FigureA1';
print(H,fname,'-dpng','-r400');

end