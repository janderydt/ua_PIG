function [UserVar,s,b,S,B,alpha]=DefineGeometry(UserVar,CtrlVar,MUA,time,FieldsToBeDefined)

%persistent Fs Fb    
addpath('./InputData/BedMachine_Antarctica');
addpath(genpath('./InputData/Observed_dhdt'));

if nargin<5
    FieldsToBeDefined='sbSB';
end

if contains(CtrlVar.Experiment,{'MEaSUREs_2015_2016_Bedmachine',...
        'PerturbationFinal_1996-GLRetreat',...
        'PerturbationFinal_1996-RetreatChangesInAGlen',...
        'PerturbationFinal_1996-RetreatChangesInC'})
    
    fprintf('Loading Bedmachine B, s, h \n');
   
    [~,s,h,~] = generateBedmachineInterpolants(MUA);
    S = 0*s;
    alpha = 0;
    b = s - h;
    
    if ~exist('./InputData/BedMachine_Antarctica/AdjustedBed_1996GL.mat')
        error('First generate adjusted bed topography to match 1996 GL.');
    else
        x = MUA.coordinates(:,1); y = MUA.coordinates(:,2);
        load('./InputData/BedMachine_Antarctica/AdjustedBed_1996GL.mat');
        B = FB(x,y);
    end
      
elseif contains(CtrlVar.Experiment,{'ERS_1996_REMAthickenedfinal',...
        'NoPerturbationFinal_1996',...
        'PerturbationFinal_1996-Calving',...
        'PerturbationFinal_1996-IceShelfThinning',...
        'PerturbationFinal_1996-ChangesInAGlen',...
        'PerturbationFinal_1996-ChangesInC'})
    
    fprintf('Loading dhdt \n');
    
    years = {'1997_2001','2002_2006','2007_2011','2012_2016'};
    %dhdt = Interpolate_dhdt(CtrlVar,MUA,years,'v01Nov2019_trimmed_3kmGLbuffer');
    dhdt = Interpolate_dhdt(CtrlVar,MUA,years,'PIG_Paolo_trimmedv2');
    
    fprintf('Loading Bedmachine B, s, h \n');
   
    [B,s,h,~] = generateBedmachineInterpolants(MUA);
    
    S = 0*s;
    alpha = 0;
    h = h - dhdt; % this is the thickness used in the 1996 inversions
    b = s - h;
    
    x = MUA.coordinates(:,1); y = MUA.coordinates(:,2);
    
    if ~exist('./InputData/BedMachine_Antarctica/AdjustedBed_1996GL.mat')
        % check GL location and imposes 1996 GL if 
        S_GL = shaperead('./InputData/GroundingLines/InSAR_GL_92-96_Combined.shp');
        I = find(b<=B & inpoly([x y],[[S_GL(1).X(1:end-1)]' [S_GL(1).Y(1:end-1)]']));
        Bold = B;
        B(I) = b(I) - 1;
        FB = scatteredInterpolant(x,y,B);
        save('./InputData/BedMachine_Antarctica/AdjustedBed_1996GL.mat','FB');   
    else
        load('./InputData/BedMachine_Antarctica/AdjustedBed_1996GL.mat');
        B = FB(x,y);
    end
    
    if contains(CtrlVar.Experiment,'PerturbationFinal_1996-IceShelfThinning')
        % now use this geometry to thin only the ice shelf 
        [b,s,h,GF]=Calc_bs_From_hBS(CtrlVar,MUA,h,S,B,917,1027);
        [GF,~,~,~]=IceSheetIceShelves(CtrlVar,MUA,GF);
        h_old = h ;
        h(GF.NodesDownstreamOfGroundingLines) = h(GF.NodesDownstreamOfGroundingLines) + dhdt(GF.NodesDownstreamOfGroundingLines);
        b = s - h;
    end
      
%     figure; PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,h-h_old); hold on;
%     for ii=1:length(S_GL)
%        plot(S_GL(ii).X,S_GL(ii).Y,'-m');
%     end
    
else

   error('Unknown case to define geometry');
    
end

% all s above zero
s(s<0)=0;

% check minimum ice thickness
h = s-b;
I = find(h<=CtrlVar.ThickMin);
s(I) = b(I)+CtrlVar.ThickMin;

if any(isnan(s)|isnan(b)|isnan(B)|isnan(S))
	error('NaN values in s, S, b or B');
end

