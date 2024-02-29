function PlotLCurve_ga

version = 'final';

%froot = '/home/UNN/wchm8/Dropbox (Northumbria University)/PIG_IceShelf_Integrity_v2';
%addpath('/media/wchm8/data1-JDeRydt/Documents/Projects/Brunt_Stresses/Matlab_functions/');
%addpath(genpath('/home/UNN/wchm8/Documents/UaMITgcm_v2/Matlab_tools'));
froot = ['/Volumes/mainJDeRydt/Ua/cases/PIG_IceShelf_Integrity/',version];
addpath(genpath('/Volumes/mainJDeRydt/Documents/Projects/Brunt_Stresses/Matlab_functions/'));
addpath(genpath(['/Volumes/mainJDeRydt/Ua/cases/PIG_IceShelf_Integrity/',version,'/Matlab_Tools']));
%addpath(genpath('/Volumes/Untitled/Documents/UaMITgcm_v2/Matlab_tools'));

%gs = [500 1000 2500 5000 7500 10000 12500 25000 50000 75000 100000];
ga = [0.1 1 10 100 1000];
%gs_label = {'500','1e3','2.5e3','5e3','7.5e3','1e4','1.25e4','2.5e4','5e4','7.5e4','1e5'};
ga_label = {'0.1','1','10','1e2','1e3'};

% loop through files and calculate misfit and regularization
for ii=1:length(ga)
    fname = [froot,'/Inversion_ERS_1996_REMAthickenedfinal_gs25000_ga',num2str(ga(ii)),'_m3_InverseRestartFile.mat'];
    if ismember(ii,[1])
         fname = [froot,'/Inversion_ERS_1996_REMAthickenedfinal_gs25000_ga',num2str(ga(ii)),'_m3_InverseRestartFile.mat'];  
    elseif ismember(ii,[3])
        fname = [froot,'/Inversion_ERS_1996_REMAthickenedfinal_gs25000_ga',num2str(ga(ii)),'_m3_InverseRestartFile.mat'];  
    end
    load(fname);  CtrlVarInRestartFile.Report_if_b_less_than_B=0;
    fprintf('Total number of iterations, ga=%s: %s\n',num2str(ga(ii)),num2str(RunInfo.Inverse.Iterations(end)));
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

% 
% figure; PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,log10(F.C));
%return

%% Plotting

H=fig('units','inches','width',50*12/72.27,'fontsize',16,'font','Helvetica');

%% L-curve
subplot('position',[0.15 0.15 0.82 0.82]); hold on;

plot(JRampl./ga.^2,JM,'o-k','markersize',5,'markerfacecolor','k','linewidth',2); hold on;
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
for ii=[1:length(ga)]
   text(1.1*JRampl(ii)./ga(ii).^2,JM(ii),ga_label{ii},'HorizontalAlignment','left','VerticalAlignment','bottom');
end
xlabel('J_{regularization}');
ylabel('J_{misfit}');

xlim([1e-4 300]); ylim([60 1e3]);
grid on;
box on;

pos = get(H,'Position');
set(H,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
fname = './Lcurve_ga';
print(H,fname,'-dpng','-r400');
