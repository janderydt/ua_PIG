function [UserVar,C,m]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)

persistent FC Fm

m = UserVar.SlidingCoefficient;

if ~exist(CtrlVar.NameOfFileForReadingSlipperinessEstimate)
     
    if m>3
        addpath('./Matlab_Tools/');
        
        InputFile = strrep(CtrlVar.Inverse.NameOfRestartOutputFile,['m',num2str(m)],'m3');
        Input = load(InputFile,'F');
        C = EstimateC_for_different_exponents(Input.F.ub,Input.F.vb,Input.F.C,Input.F.m,m,CtrlVar);
        
        %figure; PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,log10(Input.F.C)); drawnow;
        %figure; PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,sqrt(Input.F.ub.^2+Input.F.vb.^2)); drawnow;
        %figure; PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,log10(C)); drawnow;
        fprintf('\n Using rescaled C from %s \n',InputFile);
        
    elseif m==0
        
        load(['Optimal_m_gs',num2str(UserVar.Inverse.gs)]);
        
        I = find(isnan(Optimal_m) | Optimal_m>100);
        
        if strfind(CtrlVar.Experiment,'m3')
            Optimal_m(I) = 3;
            FOptimal_m = scatteredInterpolant(x,y,Optimal_m);
        elseif strfind(CtrlVar.Experiment,'m1')
            Optimal_m(I) = 1;
            FOptimal_m = scatteredInterpolant(x,y,Optimal_m);     
        else
            x(I)=[]; y(I)=[]; Optimal_m(I)=[];
            FOptimal_m = scatteredInterpolant(x,y,Optimal_m);
        end
        
        m = FOptimal_m(MUA.coordinates(:,1),MUA.coordinates(:,2));
        m(m<1)=1; m(m>50)=50;
        
        figure; PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,m);
        
        addpath('./Matlab_Tools/');
        
        InputFile = strrep(CtrlVar.Inverse.NameOfRestartOutputFile,'_m1','');
        InputFile = strrep(InputFile,'_m3','');
        InputFile = strrep(InputFile,'Optimalm','m3');
        Input = load(InputFile,'F');
        C = EstimateC_for_different_exponents(Input.F.ub,Input.F.vb,Input.F.C,Input.F.m,m,CtrlVar);
        
        figure; PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,log10(C));
        
        fprintf('\n Using optimal values for m and rescaled C from %s \n',InputFile);
        
    elseif m==-1

        %optimalmfile = ['Optimal_m_gs',num2str(UserVar.Inverse.gs),'_NoRestrictionOndhdt_v2_Interpolated.mat'];
        optimalmfile = '/mnt/SSD1/Documents/Projects/2023_PIGCalibration/Optimal_m_gs25000.mat';
        
        mopt=load(optimalmfile);
        x = mopt.MUA.coordinates(:,1); y = mopt.MUA.coordinates(:,2);
        Optimal_m = mopt.m_optimal;
        
        FOptimal_m = scatteredInterpolant(x,y,Optimal_m);
        m = FOptimal_m(MUA.coordinates(:,1),MUA.coordinates(:,2));
        m(m<1)=1; m(m>60)=60;

        %figure; PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,m);
        
        %addpath('./Matlab_Tools/');
        
        %InputFile = strrep(CtrlVar.Inverse.NameOfRestartOutputFile,'_OnlyAGlen_Optimalm','_m3');
        InputFile = strrep(CtrlVar.Inverse.NameOfRestartOutputFile,'_onlyAGlen_Optimalm','_m3');
        %InputFile = strrep(CtrlVar.Inverse.NameOfRestartOutputFile,'_Optimalm_NoRestrictionOndhdt_Interpolated','_m3');
        %InputFile = strrep(InputFile,'_Interpolated','');
        if ~exist(InputFile,"file")
            error(InputFile+"does not exist");
        end
        
        %fprintf('\n Using optimal values for m from %s and rescaled C from %s \n',optimalmfile,InputFile);
        Input = load(InputFile,'F');

        %C = EstimateC_for_different_exponents(Input.F.ub,Input.F.vb,Input.F.C,Input.F.m,m,CtrlVar);
        u=sqrt(Input.F.ub.^2+Input.F.vb.^2); 
        %u = 750;
        tau = 80 ; % units meters, year , kPa
        C=s*0+u./tau.^m;

        fprintf('\n Using initial & prior C values based on surface speed and constant basal shear stress');
        %figure; PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,log10(C));
    
    else
        ub=750; tau=80 ; % units meters, year , kPa
        C=s*0+ub/tau^m;
        fprintf('\n Using constant initial C values %s \n',num2str(C(1)));
    end
     
else
    if isempty(FC)
        load(CtrlVar.NameOfFileForReadingSlipperinessEstimate,'xC','yC','C','m');
        FC = scatteredInterpolant(xC,yC,C,'linear');
        Fm = scatteredInterpolant(xC,yC,m,'linear');
        fprintf('\n Read slipperiness from file %s \n', CtrlVar.NameOfFileForReadingSlipperinessEstimate);
    end

    x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);
    C = FC(x,y);
    m = Fm(x,y);

end

end
