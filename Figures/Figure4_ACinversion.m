function Figure4_ACinversion

version = 'final';
%% plots difference between 1996 and 2016 AGlen and C

% a. Change in AGlen (constant C), m=3, gamma_s = 50000
% b. Change in C (constant AGlen), m=3, gamma_s = 50000

%addpath(genpath('/home/UNN/wchm8/Documents/Projects/IceOceanCoupling_PIG/JGR/Matlab_Functions'));
addpath(genpath('/Volumes/mainJDeRydt/Documents/Projects/2019_Brunt_Stresses/Matlab_functions'));
addpath('/Volumes/mainJDeRydt/Matlab/Matlab_Functions/arrow/');
%addpath('/home/UNN/wchm8/Documents/UaMITgcm_v2/Matlab_tools');
addpath(genpath('/Volumes/mainJDeRydt/Documents/Projects/IceOceanCoupling_PTDC/Matlab_tools_UaMITgcm_v2'));

froot = ['/Volumes/mainJDeRydt/Ua/cases/PIG_IceShelf_Integrity/',version];
%froot = ['/Volumes/bckpJDeRydt/Ua_bckp/cases/PIG_IceShelf_Integrity/',version];

gs = [25000];
m = 'm3';
%m = 'Optimalm_NoRestrictionOndhdt_Interpolated';

mtitle = 'optimal';
titlestr = {['{\boldmath${\cal E}_{AC}^{',mtitle,'}$}'],...%\,{\rm \,(calving)}
    ['{\boldmath${\cal E}_{AC}^{',mtitle,'}$}']};

load([froot,'/Inversion_ERS_1996_REMAthickenedfinal_gs',num2str(gs),'_',m,'_InverseRestartFile.mat']);
Fspeed = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),sqrt(F.ub.^2+F.vb.^2));
x = MUA.coordinates(:,1); y = MUA.coordinates(:,2);
dx = 5e2; dy = 5e2;
[Xm,Ym] = meshgrid([min(x):dx:max(x)],[min(y):dy:max(y)]);
load([froot,'/Inversion_MEaSUREs_2015_2016_Bedmachine_gs',num2str(gs),'_',m,'_InverseRestartFile.mat']);
speed_2015_2016 = sqrt(F.ub.^2+F.vb.^2);
dspeed = speed_2015_2016 - Fspeed(MUA.coordinates(:,1),MUA.coordinates(:,2));
Fdspeed = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),dspeed);
Fspeed16 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),speed_2015_2016);
FGF16 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),GF.node);

%% 1996 data
load([froot,'/Inversion_ERS_1996_REMAthickenedfinal_gs',num2str(gs),'_',m,'_InverseRestartFile.mat']);
CtrlVar = CtrlVarInRestartFile; CtrlVar.PlotXYscale = 1e3;
FAGlen_96 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.AGlen);
FC_96 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.C);
xmin = min(MUA.coordinates(:,1)); xmax = max(MUA.coordinates(:,1)); dx=1e3;
ymin = min(MUA.coordinates(:,2)); ymax = max(MUA.coordinates(:,2)); dy=1e3;
[Xm,Ym] = ndgrid([xmin:dx:xmax],[ymin:dy:ymax]);
Temp = AGlen_vs_Temp(F.AGlen,[]);
FTemp_96 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),Temp);
Tempm_96 = FTemp_96(Xm,Ym);
I = find(~inpoly([Xm(:) Ym(:)],[MUA.Boundary.x(:) MUA.Boundary.y(:)]));
Tempm_96(I) = NaN;

[tbx_96,tby_96,tb_96] = CalcBasalTraction(CtrlVar,MUA,F.ub,F.vb,F.C,F.m,GF);
Ftb_96 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),tb_96);

%% 2016 data
load([froot,'/Inversion_MEaSUREs_2015_2016_Bedmachine_gs',num2str(gs),'_',m,'_InverseRestartFile.mat']);
MUA_1516 = MUA;
AGlen_1516 = F.AGlen;
Temp_1516 = AGlen_vs_Temp(F.AGlen,[]);
FTemp = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),Temp_1516);
Tempm_1516 = FTemp(Xm,Ym);
I = find(~inpoly([Xm(:) Ym(:)],[MUA.Boundary.x(:) MUA.Boundary.y(:)]));
Tempm_1516(I) = NaN;

Tempm_diff = Tempm_1516 - Tempm_96;
AGlen_diff = AGlen_1516 - FAGlen_96(MUA.coordinates(:,1),MUA.coordinates(:,2));
Temp_diff = Temp_1516 - FTemp_96(MUA.coordinates(:,1),MUA.coordinates(:,2));

C_1516 = F.C;
C_1516(GF.node<1)=NaN;

C_diff = C_1516 - FC_96(MUA.coordinates(:,1),MUA.coordinates(:,2));
C_diff_rel = 100*C_diff./FC_96(MUA.coordinates(:,1),MUA.coordinates(:,2));
FC_diff = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),C_diff);

[tbx_1516,tby_1516,tb_1516] = CalcBasalTraction(CtrlVar,MUA,F.ub,F.vb,F.C,F.m,GF);
tb_diff = tb_1516 - Ftb_96(MUA.coordinates(:,1),MUA.coordinates(:,2));

figure; PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,tb_diff);

return

%% Plotting

H=fig('units','inches','width',70*12/72.27,'fontsize',16,'font','Helvetica');

%% fixed C
subplot('position',[0.07 0.12 0.45 0.85]); hold on;

PlotNodalBasedQuantities_JDR(gca,MUA.connectivity,MUA.coordinates,AGlen_diff,CtrlVar);

nn = 17;
if mod(nn,2)==0
    nn=nn+1;
end
caxis(gca,[-5e-7 5e-7]);
CM = othercolor('RdBu7',nn); CM = flipdim(CM,1);
CM((nn+1)/2,:)=[0.8 0.8 0.8];

nn=35;
CM2 = othercolor('RdBu7',nn); CM2 = flipdim(CM2,1);
CM2((nn+1)/2,:)=[0.8 0.8 0.8];

GreyAGlenMin = -10e-7/(2*nn); GreyAGlenMax = 10e-7/(2*nn);
I = find(AGlen_diff > GreyAGlenMin & AGlen_diff < GreyAGlenMax);
J = find(AGlen_diff < GreyAGlenMin | AGlen_diff > GreyAGlenMax);

%Temp_diff(J)=0;
%figure(999); PlotNodalBasedQuantities_JDR(gca,MUA.connectivity,MUA.coordinates,Temp_diff,CtrlVar);
figure(H); 
%disp(['Grey colours are between ',num2str(min(Temp_diff(I))),' and ',num2str(max(Temp_diff(I)))]);

colormap(gca,CM);
cb=colorbar(gca,'Position',[0.2 0.85 0.2 0.012],'Location','northoutside','AxisLocation','in',...
 'Ticks',[-5e-7 -2.5e-7 0 2.5e-7 5e-7],'TickLabels',{'-5x10^{-7}','','0','','5x10^{-7}'});
cb.XColor='k';
cb.YColor='k';
cb.TickLength=0.055;
cb.FontSize=16;
%cb.Label.String = 'Change in $A$ [yr${}^{-1}$ kPa${}^{-3}$]';
%cb.Label.Interpreter = 'latex';
set(get(cb,'title'),'string','$A_{16} - A_{96}$ [yr${}^{-1}$ kPa${}^{-3}$]','interpreter','latex');

%[ccont,hcont] = contour(Xm/1e3,Ym/1e3,Tempm_diff,[20 20],'color',[0.5 0.5 0.5],'linestyle','-','linewidth',0.5);
%[ccont,hcont] = contour(Xm/1e3,Ym/1e3,Tempm_diff,[20 30 40],'color',[0.6 0.6 0.6]);

%dspeed = Fdspeed(Xm,Ym);
%I = find(~inpoly([Xm(:) Ym(:)],[MUA.Boundary.x MUA.Boundary.y])); dspeed(I)=NaN;
%[Ccont,hcont]=contour(Xm/CtrlVar.PlotXYscale,Ym/CtrlVar.PlotXYscale,dspeed,[50 250:250:1500],'EdgeColor','m');
%clabel([],hcont,'color','m','fontsize',12);
speed16 = Fspeed16(Xm,Ym);
GF16 = FGF16(Xm,Ym);
I = find(~inpoly([Xm(:) Ym(:)],[MUA.Boundary.x MUA.Boundary.y]) | GF16(:)<0.5); speed16(I)=NaN;
[Ccont,hcont]=contour(Xm/CtrlVar.PlotXYscale,Ym/CtrlVar.PlotXYscale,speed16,[100 500:500:4000],'EdgeColor','m');
clabel([],hcont,'manual','color','m','fontsize',12);
%I = find(~inpoly([Xm(:) Ym(:)],[MUA.Boundary.x MUA.Boundary.y]) | GF16==1); dspeed=NaN;
%[Ccont,hcont]=contour(Xm/CtrlVar.PlotXYscale,Ym/CtrlVar.PlotXYscale,dspeed,[50 250:250:1500],'EdgeColor','m');
%clabel([],hcont,'color','m','fontsize',12);

northarrow_psxy(Xm,Ym,-1500e3,-350e3);

plot(MUA.Boundary.x/1e3,MUA.Boundary.y/1e3,'-k','linewidth',1);

PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'-k','linewidth',2);

grid on; box on;
axis equal;
xlim([-1745 -1455]); ylim([-370 -10]);
xlabel(gca,'psx [km]');
ylabel(gca,'psy [km]');

title(titlestr{1},'interpreter','latex');


subplot('position',[0.53 0.12 0.45 0.85]); hold on;

PlotNodalBasedQuantities_JDR(gca,MUA.connectivity,MUA.coordinates,C_diff_rel,CtrlVar);

%caxis(gca,[-5e-2 5e-2]);
%nn=501;
%CM = othercolor('RdBu7',nn); CM = flipdim(CM,1);
%CM = parula(nn);
%CM((nn+1)/2,:)=[0.8 0.8 0.8];
%caxis(gca,[-0.5 0.5]);
caxis(gca,[-300 300]);
colormap(gca,CM2);
%cb=colorbar(gca,'Position',[0.65 0.85 0.2 0.012],'Location','northoutside','AxisLocation','in',...
% 'Ticks',[-5e-2 -2.5e-2 0 2.5e-2 5e-2],'TickLabels',{'-0.05','','0','','0.05'});
cb=colorbar(gca,'Position',[0.65 0.85 0.2 0.012],'Location','northoutside','AxisLocation','in',...
 'Ticks',[-300 0 300],'TickLabels',{'-300','0','300'});
cb.XColor='k';
cb.YColor='k';
cb.TickLength=0.055;
cb.FontSize=16;
%set(get(cb,'title'),'string','$C_{16} - C_{96}$ [m yr${}^{-1}$ kPa${}^{-3}$]','interpreter','latex');
set(get(cb,'title'),'string','$(C_{16} - C_{96})/C_{96}$ [\%]','interpreter','latex');

[Ccont,hcont]=contour(Xm/CtrlVar.PlotXYscale,Ym/CtrlVar.PlotXYscale,speed16,[100 500:500:4000],'EdgeColor','m');

plot(MUA.Boundary.x/1e3,MUA.Boundary.y/1e3,'-k','linewidth',1);

PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'-k','linewidth',2);

grid on; box on;
axis equal;
xlim([-1745 -1455]); ylim([-370 -10]);
xlabel(gca,'psx [km]'); 
yticklabels(gca,{}); ylabel(gca,'');

title(titlestr{2},'interpreter','latex');

pos = get(H,'Position');
set(H,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
fname = ['./Figure4_ACinversion_',version];


print(H,fname,'-dpng','-r400');

%figure; hold on; PlotNodalBasedQuantities_JDR(gca,MUA.connectivity,MUA.coordinates,log10(F.C),CtrlVar);
%figure; hold on; PlotNodalBasedQuantities_JDR(gca,MUA.connectivity,MUA.coordinates,log10(FC_96(MUA.coordinates(:,1),MUA.coordinates(:,2))),CtrlVar);
figure; hold on; PlotNodalBasedQuantities_JDR(gca,MUA.connectivity,MUA.coordinates,C_diff_rel,CtrlVar); caxis([-200 200]);
contour(Xm/1e3,Ym/1e3,FC_diff(Xm,Ym),[-0.5:0.1:0.5],'color','k');
hold on;



