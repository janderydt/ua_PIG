function FigureA2

version = 'final';

%% plots results from 1996 inversion

% a. L-curve
% b. Misfit, gamma_s = 50000
% c. log10(AGlen), gamma_s = 50000
% d. log10(C), gamma_s = 50000

%froot = '/home/UNN/wchm8/Dropbox (Northumbria University)/PIG_IceShelf_Integrity_v2';
%addpath('/media/wchm8/data1-JDeRydt/Documents/Projects/Brunt_Stresses/Matlab_functions/');
%addpath(genpath('/home/UNN/wchm8/Documents/UaMITgcm_v2/Matlab_tools'));
froot = ['/Volumes/mainJDeRydt/Ua/cases/PIG_IceShelf_Integrity/',version];
addpath(genpath('/Volumes/mainJDeRydt/Documents/Projects/Brunt_Stresses/Matlab_functions/'));
addpath(genpath(['/Volumes/mainJDeRydt/Ua/cases/PIG_IceShelf_Integrity/',version,'/Matlab_Tools']));
%addpath(genpath('/Volumes/Untitled/Documents/UaMITgcm_v2/Matlab_tools'));

%gs = [500 1000 2500 5000 7500 10000 12500 25000 50000 75000 100000];
gs = [2500 10000 25000 50000 75000 100000];
%gs_label = {'500','1e3','2.5e3','5e3','7.5e3','1e4','1.25e4','2.5e4','5e4','7.5e4','1e5'};
gs_label = {'2.5e3','1e4','2.5e4','5e4','7.5e4','1e5'};

% loop through files and calculate misfit and regularization
for ii=1:length(gs)
    fname = [froot,'/Inversion_ERS_1996_REMAthickenedfinal_gs',num2str(gs(ii)),'_m3_InverseRestartFile.mat'];
    if ismember(ii,[1])
         fname = [froot,'/Inversion_ERS_1996_REMAthickenedfinal_gs',num2str(gs(ii)),'_fromgs10000_m3_InverseRestartFile.mat'];  
    elseif ismember(ii,[3])
        fname = [froot,'/Inversion_ERS_1996_REMAthickenedfinal_gs',num2str(gs(ii)),'_fromgs50000_m3_InverseRestartFile.mat'];  
     end
    load(fname);  CtrlVarInRestartFile.Report_if_b_less_than_B=0;
    fprintf('Total number of iterations, gs=%s: %s\n',num2str(gs(ii)),num2str(RunInfo.Inverse.Iterations(end)));
    % calculate fractional changes in misfit
    dJ = (RunInfo.Inverse.J(2:end)-RunInfo.Inverse.J(1:end-1))./RunInfo.Inverse.J(1:end-1);
    fprintf('Fractional change after %s iterations: %d\n',num2str(RunInfo.Inverse.Iterations(end)),dJ(end-1));
     CtrlVarInRestartFile.WriteRunInfoFile=0;
    if ~isfield(CtrlVarInRestartFile.Inverse,'dFuvdClambda')
        CtrlVarInRestartFile.Inverse.dFuvdClambda=0;
        CtrlVarInRestartFile.SlidingLaw = "Weertman";
    else
        %CtrlVarInRestartFile.Inverse.dFuvdClambda
    end
    CtrlVar=Ua2D_DefaultParameters();
    CtrlVarInRestartFile.LevelSetMethod = CtrlVar.LevelSetMethod;
    CtrlVarInRestartFile.Enforce_bAboveB = CtrlVar.Enforce_bAboveB;
    CtrlVarInRestartFile.GroundingFloatingMaskContains = CtrlVar.GroundingFloatingMaskContains;
    CtrlVarInRestartFile.uvMinimisationQuantity = CtrlVar.uvMinimisationQuantity;
    CtrlVarInRestartFile.uvExitBackTrackingStepLength = CtrlVar.uvExitBackTrackingStepLength;
    CtrlVarInRestartFile.uvDesiredWorkAndForceTolerances = CtrlVar.uvDesiredWorkAndForceTolerances;
    CtrlVarInRestartFile.uvDesiredWorkOrForceTolerances = CtrlVar.uvDesiredWorkOrForceTolerances;
    UserVarIeenRestartFile=[];
    [UserVar,RunInfo,F,l,drdu,Ruv,Lubvb]= uv(UserVarInRestartFile,RunInfo,CtrlVarInRestartFile,MUA,BCs,F,l);
    [Rgrad,Rampl]=Regularisation_JDR(UserVarInRestartFile,CtrlVarInRestartFile,MUA,BCs,F,l,Priors,Meas,BCsAdjoint,RunInfo);
    [I,dIdp,ddIddp,MisfitOuts]=Misfit(UserVarIeenRestartFile,CtrlVarInRestartFile,MUA,BCs,F,l,Priors,Meas,BCsAdjoint,RunInfo,drdu);
    JM(ii) = I;
    JRgrad(ii) = Rgrad;
    JRampl(ii) = Rampl;
end

load([froot,'/Inversion_ERS_1996_REMAthickenedfinal_gs',num2str(gs(3)),'_fromgs50000_m3_InverseRestartFile.mat']);
[froot,'/Inversion_ERS_1996_REMAthickenedfinal_gs',num2str(gs(3)),'_fromgs50000_m3_InverseRestartFile.mat']
xmin = min(MUA.coordinates(:,1)); xmax = max(MUA.coordinates(:,1)); dx=1e3;
ymin = min(MUA.coordinates(:,2)); ymax = max(MUA.coordinates(:,2)); dy=1e3;
[Xm,Ym] = ndgrid([xmin:dx:xmax],[ymin:dy:ymax]);
Temp = AGlen_vs_Temp(F.AGlen,[]);
FTemp = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),Temp);
Tempm = FTemp(Xm,Ym);
I = find(~inpoly([Xm(:) Ym(:)],[MUA.Boundary.x(:) MUA.Boundary.y(:)]));
Tempm(I) = NaN;
% 
% figure; PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,log10(F.C));
%return

%% Plotting

H=fig('units','inches','width',60*12/72.27,'height',70*12/72.27,'fontsize',16,'font','Helvetica');

%% L-curve
subplot('position',[0.1 0.62 0.42 0.36]); hold on;

plot(JRgrad,JM,'o-k','markersize',5,'markerfacecolor','k','linewidth',2); hold on;
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');

subplot('position',[0.6 0.54 0.38 0.45]); hold on;

plot(JRampl,JM,'o-k','markersize',5,'markerfacecolor','k','linewidth',2); hold on;
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');

return
plot(JR(3)./gs(3).^2,JM(3),'ok','markersize',12,'markerfacecolor','w');
plot(JR(3)./gs(3).^2,JM(3),'or','markersize',8,'markerfacecolor','r');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
text(0.8*JR(1)./gs(1).^2,JM(1),gs_label{1},'HorizontalAlignment','left','VerticalAlignment','top');
for ii=[2 4:length(gs)]
   text(1.1*JR(ii)./gs(ii).^2,JM(ii),gs_label{ii},'HorizontalAlignment','left','VerticalAlignment','bottom');
end
text(1.2*JR(3)./gs(3).^2,JM(3),gs_label{3},'color','r','HorizontalAlignment','left','VerticalAlignment','bottom');
grid on; box on;
xlim([6e-9 4.5e-7]); ylim([54 140]);
xlabel('J_{regularization}');
ylabel('J_{misfit}');

%% Misfit
subplot('position',[0.6 0.54 0.38 0.45]); hold on;

CtrlVar = CtrlVarInRestartFile; CtrlVar.PlotXYscale = 1e3;
I = Meas.us==0;
Meas.us(I) = NaN;
PlotNodalBasedQuantities_JDR(MUA.connectivity,MUA.coordinates,sqrt(F.ub.^2+F.vb.^2)-sqrt(Meas.us.^2+Meas.vs.^2),CtrlVar);

caxis(gca,[-30 30]);
CM = othercolor('RdBu7');
colormap(gca,CM);
cb=colorbar(gca,'Position',[0.69 0.93 0.2 0.012],'Location','northoutside','AxisLocation','in',...
 'Ticks',[-30:15:30],'TickLabels',{'-30','-15','0','15','30'});
cb.XColor='k';
cb.YColor='k';
cb.TickLength=0.055;
cb.FontSize=16;
set(get(cb,'title'),'string','$u_{{\rm model}} - u_{{\rm obs}}\, {\rm [m\, yr^{-1}]}$','Interpreter','latex');

plot(MUA.Boundary.x/1e3,MUA.Boundary.y/1e3,'-k','linewidth',1);

PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'-k','linewidth',1.5);

grid on; box on;
axis equal;
xlim([-1745 -1455]); ylim([-370 20]);
xticklabels(gca,{}); xlabel(gca,'');
ylabel(gca,'psy [km]');

%% AGlen
subplot('position',[0.1 0.08 0.45 0.45]); hold on;

PlotNodalBasedQuantities_JDR(MUA.connectivity,MUA.coordinates,F.AGlen,CtrlVar);

CM2 = [[linspace(0.85,0,64)';linspace(0,1,128)'] [linspace(0.85,0,192)'] [linspace(0.85,1,128)';linspace(1,0,64)']];
TempBins = AGlen_vs_Temp([],[-35 -30 -25 -20 -15 -10 -5 0 5]);
caxis(gca,[TempBins(1) TempBins(end)]);
CM3 = CM2(1:21:end,:);
colormap(gca,CM3);

cb=colorbar(gca,'Position',[0.23 0.48 0.2 0.012],'Location','northoutside','AxisLocation','in',...
 'Ticks',[1e-9 1e-8 1e-7 1e-6],'TickLabels',{'10^{-9}','10^{-8}','10^{-7}','10^{-6}'});
cb.XColor='w';
cb.YColor='w';
cb.TickLength=0.055;
cb.FontSize=16;
set(get(cb,'title'),'string','$A\, {\rm [yr^{-1}\, kPa^{-3}]}$','color','w','interpreter','latex');
% cb.Label.String = {'AGlen'};
% pos = get(cb,'Position');
% cb.Label.Position = [pos(1) pos(2)-0.01];
% cb.Label.Rotation = 0;
%cb.Ruler.Scale='log';
set(gca,'colorscale','log');

[ccont,hcont] = contour(Xm/1e3,Ym/1e3,Tempm,[0 0],'color','k');%[0 0 0]);
%[ccont,hcont] = contour(Xm/1e3,Ym/1e3,Tempm,[-20 -20],'color',[0 0 0 ],'linestyle','--');
%clabel([],hcont,[0 0],'color','w','fontsize',12);

plot(MUA.Boundary.x/1e3,MUA.Boundary.y/1e3,'-k','linewidth',1);

PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'-w','linewidth',2);

grid on; box on;
axis equal;
xlim([-1745 -1455]); ylim([-370 20]);
xlabel(gca,'psx [km]'); ylabel(gca,'psy [km]');
xticks(gca,[-1700:50:-1500]); xticklabels(gca,{'-1700','','-1600','','-1500'});

%% C
subplot('position',[0.6 0.08 0.38 0.45]); hold on;

C = F.C; C(GF.node<1)=NaN;
PlotNodalBasedQuantities_JDR(MUA.connectivity,MUA.coordinates,C,CtrlVar);

CM2 = [[linspace(0.85,0,64)';linspace(0,0.85,64)'] linspace(0.85,0,128)' linspace(0.85,1,128)'];
colormap(gca,CM2);
caxis(gca,[1e-5 10]);
cb=colorbar(gca,'Position',[0.69 0.48 0.2 0.012],'Location','northoutside','AxisLocation','in',...
 'Ticks',[1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1 10],'TickLabels',{'10^{-7}','','10^{-5}','','10^{-3}','','10^{-1}','','10'});
cb.XColor='w';
cb.YColor='w';
cb.TickLength=0.055;
cb.FontSize=16;
set(get(cb,'title'),'string','$C\, {\rm [m\, yr^{-1}\, kPa^{-3}]}$','color','w','interpreter','latex');
%cb.Ruler.Scale='log';
set(gca,'colorscale','log');

plot(MUA.Boundary.x/1e3,MUA.Boundary.y/1e3,'-k','linewidth',1);

PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'-k','linewidth',2);

grid on; box on;
axis equal;
xlim([-1745 -1455]); ylim([-370 20]);
xlabel(gca,'psx [km]');
xticks(gca,[-1700:50:-1500]); xticklabels(gca,{'-1700','','-1600','','-1500'});
yticklabels(gca,{}); ylabel(gca,'');

pos = get(H,'Position');
set(H,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
fname = './FigureA2';
print(H,fname,'-dpng','-r400');

