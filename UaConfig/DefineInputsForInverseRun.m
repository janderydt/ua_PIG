function [UserVar,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo]=DefineInputsForInverseRun(UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo)

fprintf('Reading inputs for inverse run \n');

persistent Fvx Fvy Ferr;

% if isempty(Fvx)
%     
% 	load(['./InputData/',UserVar.VelocityInterpolants],'Fvx','Fvy','Ferr');
%     
%     fprintf('Done loading Measures interpolants \n');
% 
% end

if isempty(Fvx)
    
   velocities = dir('./InputData/Velocities/*.mat'); 
   %isub = [velocities(:).isdir]; %# returns logical vector
   %velocities(isub)=[];
   for ii=1:length(velocities)
       fprintf('Option %s: %s\n',num2str(ii),velocities(ii).name);
       NameMatch(ii) = contains(CtrlVar.Experiment,velocities(ii).name(1:end-4));
   end
   Option = find(NameMatch==1);
   if isempty(Option)
       error('Velocities do not exist');
   else
       fprintf('Chosen option: %s - %s\n',num2str(Option),velocities(Option).name);
       load([velocities(Option).folder,'/',velocities(Option).name]);
   end
     
   Fvx = griddedInterpolant(X,Y,VX);
   Fvy = griddedInterpolant(X,Y,VY);
   Ferr = griddedInterpolant(X,Y,ERR);
    
end

x = MUA.coordinates(:,1); y = MUA.coordinates(:,2);

if strfind(UserVar.Inverse.Measurements,'dhdt')
    fprintf('Adding dhdt to cost function \n');
    
    dhdtMeas = 0*x;
    dhdtError = 0*x+UserVar.Inverse.dhdtError;
    
    Meas.dhdt=dhdtMeas;
    Meas.dhdtCov=sparse(1:MUA.Nnodes,1:MUA.Nnodes,dhdtError.^2,MUA.Nnodes,MUA.Nnodes);
    if any(isnan(dhdtError))
        error('NaN in dhdtError'); 
    end
end

uMeas = Fvx(x,y);
vMeas = Fvy(x,y);
errMeas = Ferr(x,y); 
%errMeas = 0*uMeas + 1;

%% set larger errors outside area of interest
%% -> does this solve the issue of convergence, which might be related to the
%% steep volcanos to the west of Thwaites?
%I = find(x>-1.35e6 | y<-4.2e5); errMeas(I)=1e2;

% Because I put NaN in where there is no data, this will give NaN at location not surounded by four data points
% I can then afterwards find the NaN and put in some data with very high errors.
% I need values everywhere, so I now set vel to zero where no measurments are available, 
% and the error to a very large value

Mask=find(isnan(uMeas) | isnan(vMeas) );
MaskErr = find(isnan(errMeas));

uMeas(Mask)=0;
vMeas(Mask)=0;
errMeas(Mask)=1e10;
errMeas(MaskErr)=1e10;
errMeas(errMeas==0)=eps;

if strfind(UserVar.Experiment,'onlyAGlen')
    % we only invert for AGlen on the floating parts, so we set large
    % velocity errors for the grounded ice
    I = find(GF.node==1);
    uMeas(I)=0;
    vMeas(I)=0;
    errMeas(I)=1e100;
end

Meas.us = uMeas;
Meas.vs = vMeas;

%figure; PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,sqrt(Meas.us.^2+Meas.vs.^2)); hold on;

usError = errMeas;
vsError = errMeas;

if any(isnan(Meas.us))
	error('NaN in Meas.us');
elseif any(isnan(Meas.vs))
    error('NaN in Meas.vs');
elseif any(isnan(usError))
    error('NaN in usError');
end

Meas.usCov=sparse(1:MUA.Nnodes,1:MUA.Nnodes,usError.^2,MUA.Nnodes,MUA.Nnodes);
Meas.vsCov=sparse(1:MUA.Nnodes,1:MUA.Nnodes,vsError.^2,MUA.Nnodes,MUA.Nnodes);

[UserVar,InvStartValues.C,InvStartValues.m]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,CtrlVar.time,F.s,F.b,F.s-F.b,F.S,F.B,F.rho,F.rhow,GF);
[UserVar,InvStartValues.AGlen,InvStartValues.n]=DefineAGlenDistribution(UserVar,CtrlVar,MUA,CtrlVar.time,F.s,F.b,F.s-F.b,F.S,F.B,F.rho,F.rhow,GF);

% original=[CtrlVar.Experiment,'_InverseResults_Intermediate.mat'] ;
% new=[CtrlVar.Experiment,'_InverseResults_Intermediate_',date,'.mat'];
% if exist(original)==2
%     copyfile(original,new);
%     Frestart=load(original,'F');
%     InvStartValues.C = Frestart.F.C;
%     InvStartValues.AGlen = Frestart.F.AGlen;
% end

listingCC=dir('CC.mat') ; listingCA=dir('CAGlen.mat') ;
CC = []; CAGlen = []; 

%%  Covariance matrices of priors
% 
if CtrlVar.AGlenisElementBased
    CAGlen=sparse(1:MUA.Nele,1:MUA.Nele,1,MUA.Nele,MUA.Nele);
else
    CAGlen=sparse(1:MUA.Nnodes,1:MUA.Nnodes,1,MUA.Nnodes,MUA.Nnodes);
end

if strcmpi(CtrlVar.Inverse.Regularize.Field,'cov')
    Err=1e-2 ; Sigma=1e3 ; DistanceCutoff=10*Sigma;
    
    if CtrlVar.CisElementBased
        [CC]=SparseCovarianceDistanceMatrix(xC,yC,Err,Sigma,DistanceCutoff);
    else
        [CC]=SparseCovarianceDistanceMatrix(xC,yC,Err,Sigma,DistanceCutoff);
    end
    
else
    if CtrlVar.CisElementBased
        CC=sparse(1:MUA.Nele,1:MUA.Nele,1,MUA.Nele,MUA.Nele);
    else
        CC=sparse(1:MUA.Nnodes,1:MUA.Nnodes,1,MUA.Nnodes,MUA.Nnodes);
    end
end

  

%% Define Priors
Priors.B=F.B;
Priors.n=3;

if UserVar.SlidingCoefficient==-1
    Priors.m=InvStartValues.m; 
else
    Priors.m=UserVar.SlidingCoefficient; 
end

if contains(CtrlVar.Experiment,'1996prior')
    load('Inversion_ERS_1996_gs10000_InverseRestartFile.mat','F');
    fprintf('Using prior information from Inversion_ERS_1996_gs10000 with PreMult=I');
    Priors.AGlen=F.AGlen;
    Priors.C=F.C;
else
    [~,Priors.AGlen,~]=DefineAGlenDistribution(UserVar,CtrlVar,MUA,CtrlVar.time,F.s,F.b,F.s-F.b,F.S,F.B,F.rho,F.rhow,GF);
    [~,Priors.C,~]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,CtrlVar.time,F.s,F.b,F.s-F.b,F.S,F.B,F.rho,F.rhow,GF);
    
    %AGlen0 = AGlenVersusTemp(-15);
    %Priors.AGlen=F.s*0+AGlen0;
    %ub=750; tau=80 ; % units meters, year , kPa
    %C0=ub./tau.^Priors.m;
    %Priors.C=F.s*0+C0;
end

Priors.rho=F.rho;
Priors.rhow=F.rhow;
Priors.CovAGlen=CAGlen;
Priors.CovC=CC;


end
