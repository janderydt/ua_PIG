function [UserVar,CtrlVar,MeshBoundaryCoordinates]=Ua2D_InitialUserInput_Inverse(UserVar,CtrlVar)

CtrlVar.Parallel.uvhAssembly.parfor.isOn=0;     % assembly over integration points done in parallel using parfor
CtrlVar.Parallel.uvhAssembly.spmd.isOn=0;       % assembly in parallel using spmd over sub-domain (domain decomposition)  
CtrlVar.Parallel.uvhAssembly.spmd.nWorkers=[];  % If left empty, all workers available are used
CtrlVar.Parallel.isTest=false;

CtrlVar.TimeDependentRun=0;  % {0|1} if true (i.e. set to 1) then the run is a forward transient one, if not
CtrlVar.doInverseStep=1;
CtrlVar.Restart=UserVar.Restart;
CtrlVar.AdaptMesh=0;
CtrlVar.ReadInitialMesh=1;
CtrlVar.dt=1e-3; CtrlVar.TotalNumberOfForwardRunSteps=10e5; CtrlVar.TotalTime=100; 
CtrlVar.time=0;
CtrlVar.ResetTime=0;
CtrlVar.ResetTimeStep=0;
CtrlVar.NRitmax=50;
%CtrlVar.AsymmSolver='AugmentedLagrangian';

%
CtrlVar.TriNodes=3;
CtrlVar.nip = 6;
CtrlVar.niph = 6;

CtrlVar.ATStimeStepTarget=0.05;
CtrlVar.dtmin = 1e-10;
CtrlVar.UaOutputsDt=1;
CtrlVar.ATStimeStepFactorUp=1.5 ;
CtrlVar.ATStimeStepFactorDown=2 ;
CtrlVar.ATSTargetIterations=4;

CtrlVar.InitialDiagnosticStep=1;
CtrlVar.InitialDiagnosticStepAfterRemeshing=0;
CtrlVar.Implicituvh=1;
CtrlVar.TG3=0 ; CtrlVar.Gamma=1;
CtrlVar.uvhTimeSteppingMethod='supg';

CtrlVar.Experiment=UserVar.Experiment;
%CtrlVar.MITgcmDataDirectory=['/data/dataphy/janryd69/Ua_MITgcm/',Experiment,'/MIT_data'];
CtrlVar.UaDataDirectory=pwd;
CtrlVar.logfilename=[CtrlVar.UaDataDirectory,'/',CtrlVar.Experiment,'.log'];

%%
if strfind(UserVar.Experiment,'1996')
    load BoundaryCoordinates_1996 MeshBoundaryCoordinates;
elseif strfind(UserVar.Experiment,'2015_2016')
    load BoundaryCoordinates_2016 MeshBoundaryCoordinates;
end

CtrlVar.WriteRestartFile=1;
CtrlVar.WriteRestartFileInterval=100;

CtrlVar.doplots=1;
CtrlVar.PlotWaitBar=0;
CtrlVar.PlotOceanLakeNodes=0;
CtrlVar.PlotMesh=1;  CtrlVar.PlotBCs=1;


CtrlVar.MeltNodesDefinition='edge-wise';
CtrlVar.MassBalanceGeometryFeedback = 0;
CtrlVar.MeltRateFactor=1;
CtrlVar.MeltReductionTime=Inf;

CtrlVar.NameOfRestartFiletoWrite=[CtrlVar.Experiment,'-RestartFile.mat'];
CtrlVar.NameOfRestartFiletoRead='';

CtrlVar.MeshSizeMax=10e3;
CtrlVar.MeshSize=CtrlVar.MeshSizeMax;
CtrlVar.MeshSizeMin=CtrlVar.MeshSizeMax/5;
%CtrlVar.MeshSizeFastFlow=CtrlVar.MeshSizeMax/10;
%CtrlVar.MeshSizeIceShelves=CtrlVar.MeshSizeMax/10;
CtrlVar.MeshSizeBoundary=CtrlVar.MeshSize;

%     CtrlVar.GmshGeoFileAdditionalInputLines{1}='Field[7] = Box;'; % these lines are added to the gmsh .geo input file each time such a file is created
%     CtrlVar.GmshGeoFileAdditionalInputLines{2}='Field[7].VIn = 10e3;';
%     CtrlVar.GmshGeoFileAdditionalInputLines{3}='Field[7].VOut = 20e3;';
%     CtrlVar.GmshGeoFileAdditionalInputLines{4}='Field[7].XMin = -1700e3;';
%     CtrlVar.GmshGeoFileAdditionalInputLines{5}='Field[7].XMax = -1200e3;';
%     CtrlVar.GmshGeoFileAdditionalInputLines{6}='Field[7].YMin = -600e3;';
%     CtrlVar.GmshGeoFileAdditionalInputLines{7}='Field[7].YMax = -50e3;';
%     CtrlVar.GmshGeoFileAdditionalInputLines{8}='Field[10] = Min;';
%     CtrlVar.GmshGeoFileAdditionalInputLines{9}='Field[10].FieldsList = {7};';
%     CtrlVar.GmshGeoFileAdditionalInputLines{10}='Background Field = 10;';
%

if strfind(UserVar.Experiment,'1996')
    CtrlVar.ReadInitialMeshFileName='./MeshFileAdapt_PIG1996_Ele109300_Nod3.mat';
elseif strfind(UserVar.Experiment,'2015_2016')
    CtrlVar.ReadInitialMeshFileName='./MeshFileAdapt_PIG2016_Ele107705_Nod3.mat';
end
%CtrlVar.ReadInitialMeshFileName='./PIG_500strains_600mshelf_600mgl_5kmstreams_10km_mesh.mat';
CtrlVar.SaveInitialMeshFileName='MeshFile.mat';

CtrlVar.OnlyMeshDomainAndThenStop=0; % if true then only meshing is done and no further calculations. Usefull for checking if mesh is reasonable
CtrlVar.MaxNumberOfElements=150e4;


%% plotting
CtrlVar.PlotXYscale=1;

%%
CtrlVar.InfoLevelNonLinIt=1;


%% adapt mesh

CtrlVar.MeshRefinementMethod='explicit:local:newest vertex bisection';
%CtrlVar.MeshRefinementMethod='explicit:global';    % can have any of these values:
                                                   % 'explicit:global'
                                                   % 'explicit:local'
                                                   % 'implicit:global'  (broken at the moment, do not use)
                                                   % 'implicit:local'   (broken at the moment, do not use)
I=1;
CtrlVar.ExplicitMeshRefinementCriteria(I).Name='effective strain rates';
CtrlVar.ExplicitMeshRefinementCriteria(I).Scale=1e-3;
CtrlVar.ExplicitMeshRefinementCriteria(I).EleMin=[];
CtrlVar.ExplicitMeshRefinementCriteria(I).EleMax=[];
CtrlVar.ExplicitMeshRefinementCriteria(I).p=[];
CtrlVar.ExplicitMeshRefinementCriteria(I).InfoLevel=1000;
CtrlVar.ExplicitMeshRefinementCriteria(I).Use=true;

I=I+1;
CtrlVar.ExplicitMeshRefinementCriteria(I).Name='thickness gradient';
CtrlVar.ExplicitMeshRefinementCriteria(I).Scale=1e-2;
CtrlVar.ExplicitMeshRefinementCriteria(I).EleMin=[];
CtrlVar.ExplicitMeshRefinementCriteria(I).EleMax=[];
CtrlVar.ExplicitMeshRefinementCriteria(I).p=[];
CtrlVar.ExplicitMeshRefinementCriteria(I).InfoLevel=1000;
CtrlVar.ExplicitMeshRefinementCriteria(I).Use=true;

I=I+1;
CtrlVar.ExplicitMeshRefinementCriteria(I).Name='effective strain rates gradient';
CtrlVar.ExplicitMeshRefinementCriteria(I).Scale=5e-7;
CtrlVar.ExplicitMeshRefinementCriteria(I).EleMin=[];
CtrlVar.ExplicitMeshRefinementCriteria(I).EleMax=[];
CtrlVar.ExplicitMeshRefinementCriteria(I).p=[0.7];
CtrlVar.ExplicitMeshRefinementCriteria(I).InfoLevel=1000;
CtrlVar.ExplicitMeshRefinementCriteria(I).Use=true;
                                                   
CtrlVar.MeshAdapt.GLrange=[];%[20000 CtrlVar.MeshSize ; 5000  CtrlVar.MeshSizeMin];

CtrlVar.RefineMeshOnStart=0;
CtrlVar.InfoLevelAdaptiveMeshing=1;                                            
CtrlVar.AdaptMeshInitial=0  ; % remesh in first iteration (Itime=1)  even if mod(Itime,CtrlVar.AdaptMeshInterval)~=0.
CtrlVar.AdaptMeshAndThenStop=0;    % if true, then mesh will be adapted but no further calculations performed
                                   % useful, for example, when trying out different remeshing options (then use CtrlVar.doAdaptMeshPlots=1 to get plots)
%CtrlVar.MeshRefinementMethod='explicit:local';
CtrlVar.AdaptMeshMaxIterations=10;
CtrlVar.SaveAdaptMeshFileName='MeshFileAdapt';    %  file name for saving adapt mesh. If left empty, no file is written
CtrlVar.AdaptMeshRunStepInterval=0 ; % remesh whenever mod(Itime,CtrlVar.AdaptMeshInterval)==0
CtrlVar.doAdaptMeshPlots=0; 


%% Adjoint
CtrlVar.Inverse.MinimisationMethod='MatlabOptimization'; % {,'UaOptimization'}
%CtrlVar.Inverse.MinimisationMethod='UaOptimization';
CtrlVar.CisElementBased=0;   
CtrlVar.AGlenisElementBased=0;
%CtrlVar.Inverse.CalcGradI=true;
%CtrlVar.AGlenmin=AGlenVersusTemp(-15)/1e4;
%CtrlVar.AGlenmax=AGlenVersusTemp(-15)*1e4;
CtrlVar.Cmin=1.0000e-150;
CtrlVar.Cmax=1.0000e+150;
%CtrlVar.Czero = CtrlVar.Cmin;

%CtrlVar.InitTrustRegionRadius = 1;
%CtrlVar.InitBarrierParam = 1e-7;
if CtrlVar.Restart
    RestartData = load(UserVar.Inverse.NameOfRestartInputFile,'RunInfo');
    CtrlVar.Inverse.Iterations = max([0,UserVar.Inverse.Iterations - RestartData.RunInfo.Inverse.Iterations(end)]);
    fprintf('RESTART and %s iterations still to do\n',num2str(CtrlVar.Inverse.Iterations));%
    %CtrlVar.InitTrustRegionRadius = 0.1;
    %CtrlVar.InitBarrierParam = 1e-10;
else 
     CtrlVar.Inverse.Iterations=UserVar.Inverse.Iterations; % Number of inverse iterations
end

% CtrlVar.Inverse.MatlabOptimisationParameters = optimoptions('fminunc',...
%     'Algorithm','quasi-newton',...
%     'MaxIterations',CtrlVar.Inverse.Iterations,...
%     'MaxFunctionEvaluations',1000,...
%     'Display','iter-detailed',...
%     'OutputFcn',@fminuncOutfun,...
%     'Diagnostics','on',...
%     'OptimalityTolerance',1e-20,...
%     'StepTolerance',1e-20,...
%     'PlotFcn',{@optimplotfval,@optimplotstepsize},...
%     'SpecifyObjectiveGradient',CtrlVar.Inverse.CalcGradI);

CtrlVar.Inverse.MatlabOptimisationParameters = optimoptions('fmincon',...
        'Algorithm','interior-point',...
        'CheckGradients',false,...
        'ConstraintTolerance',1e-10,...
        'HonorBounds',true,...
        'Diagnostics','on',...
        'DiffMaxChange',Inf,...
        'DiffMinChange',0,...
        'Display','iter-detailed',...
        'FunValCheck','off',...
        'MaxFunctionEvaluations',1e6,...
        'MaxIterations',CtrlVar.Inverse.Iterations,...,...
        'OptimalityTolerance',1e-20,...
        'OutputFcn',@fminuncOutfun,...
        'PlotFcn',{@optimplotlogfval,@optimplotstepsize},...
        'StepTolerance',1e-30,...
        'FunctionTolerance',1,...
        'UseParallel',true,...
        'HessianApproximation','lbfgs',...%{'lbfgs',250},...
        'HessianFcn',[],...
        'HessianMultiplyFcn',[],...
        'InitBarrierParam',1e-7,...           % On a restart this might have to be reduced if objective function starts to increase
        'ScaleProblem','none',...
        'InitTrustRegionRadius',1,...         % set to smaller value if the forward problem is not converging
        'SpecifyConstraintGradient',false,...
        'SpecifyObjectiveGradient',true,...
        'SubproblemAlgorithm','cg');  % here the options are 'gc' and 'factorization', unclear which one is the better one, 'factorization' is the matlab default
    

CtrlVar.Inverse.WriteRestartFile=1;  % always a good idea to write a restart file. 
CtrlVar.Inverse.NameOfRestartOutputFile=[CtrlVar.Experiment,'_InverseRestartFile.mat'];
CtrlVar.Inverse.NameOfRestartInputFile=UserVar.Inverse.NameOfRestartInputFile;
%CtrlVar.Inverse.NameOfRestartOutputFile;
CtrlVar.NameOfFileForReadingSlipperinessEstimate=UserVar.Inverse.NameOfFileForReadingSlipperinessEstimate;
CtrlVar.NameOfFileForReadingAGlenEstimate=UserVar.Inverse.NameOfFileForReadingAGlenEstimate;
CtrlVar.NameOfFileForSavingSlipperinessEstimate=[CtrlVar.Experiment,'_C-Estimate.mat'];
CtrlVar.NameOfFileForSavingAGlenEstimate=[CtrlVar.Experiment,'_AGlen-Estimate.mat'];

CtrlVar.Inverse.Measurements=UserVar.Inverse.Measurements;

% It is usually better to invert for log(A) and log(C) rather than A and C.
% The default is to invert for log(A) and log(C) simultaneously.
%CtrlVar.Inverse.InvertFor= 'logAGlenlogC'; % {'C','logC','AGlen','logAGlen','logAGlenlogC'}
CtrlVar.Inverse.InvertFor=UserVar.Inverse.InvertFor;

% The gradient of the objective function is calculated using the adjoint method.
% When inverting for C only, one can also use a gradient based on a `FixPoint'
% iteration, which is often a very good initial approach. 
%CtrlVar.Inverse.DataMisfit.GradientCalculation='FixPoint'; % {'Adjoint','FixPoint'}
CtrlVar.Inverse.DataMisfit.GradientCalculation='Adjoint';

% The gradient of the objective function can be premultiplied with the inverse
% of the mass matrix. This creates a `mesh independent' gradient. This has both
% advantages and disadvantages. 
CtrlVar.Inverse.AdjointGradientPreMultiplier='M'; % {'I','M'}

% 
% CtrlVar.Inverse.TestAdjoint.isTrue=0; % If true then perform a brute force calculation 
%                                       % of the directinal derivative of the objective function.  
% %CtrlVar.Inverse.TestAdjoint.FiniteDifferenceType='central' ; % {'central','forward'}
% CtrlVar.Inverse.TestAdjoint.FiniteDifferenceStepSize=1e-8 ;
% CtrlVar.Inverse.TestAdjoint.iRange=[] ;

% Regularisation can be applied on A and C or log(A) and log(C). Also possible
% to use a covariance matrix for A and C. 
%
% Select Bayesian motivated regularisation by setting 
% CtrlVar.Inverse.Regularize.Field='cov' and Tikhonov regularisation
% by setting CtrlVar.Inverse.Regularize.Field to either 'C','logC','AGlen','logAGlen','logAGlenlogC'
%
% Default is Tikhonov regularisation on log(A) and log(C)
CtrlVar.Inverse.Regularize.Field= CtrlVar.Inverse.InvertFor; % {'cov','C','logC','AGlen','logAGlen','logAGlenlogC'}
%CtrlVar.Inverse.Regularize.Field='logC';

% [ -- Parameters specific to Tikhonov regularisation
% See the above definition of R in the case of Tikhonov regularisation.
% The values of these parameters can be expected to be highly problem dependent.
% By default regularisation is switched on, but can the switched off by setting
% the gs and the ga parameters to zero.
CtrlVar.Inverse.Regularize.C.gs=0; 
CtrlVar.Inverse.Regularize.C.ga=0;
CtrlVar.Inverse.Regularize.logC.gs=UserVar.Inverse.gs;
CtrlVar.Inverse.Regularize.logC.ga=1;

CtrlVar.Inverse.Regularize.AGlen.gs=0;
CtrlVar.Inverse.Regularize.AGlen.ga=0;
CtrlVar.Inverse.Regularize.logAGlen.gs=UserVar.Inverse.gs;
CtrlVar.Inverse.Regularize.logAGlen.ga=1e4;

fprintf('CtrlVar.Inverse.Regularize.logC.gs = %s \n',num2str(CtrlVar.Inverse.Regularize.logC.gs));
fprintf('CtrlVar.Inverse.Regularize.logAGlen.gs = %s \n',num2str(CtrlVar.Inverse.Regularize.logAGlen.gs));

%  -]

% I and R are multiplied by these followign DataMisit and Regularisation
% multipliers. This is a convening shortcut of getting rid of either the misfit
% (I) or the regularization term (R) in the objective function (J).
% CtrlVar.Inverse.DataMisfit.Multiplier=1;
% CtrlVar.Inverse.Regularize.Multiplier=1;

% CtrlVar.MUA.MassMatrix=1;
% CtrlVar.MUA.StiffnessMatrix=1;


% [----------  The following parameters are only relevant if using the UaOptimization
% i.e. only if CtrlVar.Inverse.MinimisationMethod='UaOptimization';
% The Ua optimisation is a simple non-linear conjugate-gradient method with automated
% resets, combined with a (one-sided) line search. The reset is done if the angle between
% subsequent steepest decent directions is to far from 90 degrees, or if the
% update parameter becomes negative (only relevant for Polak-Ribiere and
% Hestens-Stiefel).
% CtrlVar.Inverse.GradientUpgradeMethod='ConjGrad' ; %{'SteepestDecent','ConjGrad'}
% CtrlVar.Inverse.InitialLineSearchStepSize=[];
% CtrlVar.Inverse.MinimumAbsoluteLineSearchStepSize=1e-20; % minimum step size in backtracking
% CtrlVar.Inverse.MinimumRelativelLineSearchStepSize=1e-5; % minimum fractional step size relative to initial step size
% CtrlVar.Inverse.MaximumNumberOfLineSeachSteps=50;
% CtrlVar.ConjugatedGradientsRestartThreshold=20 ; % degrees!
% CtrlVar.ConjugatedGradientsUpdate='PR'; % (FR|PR|HS|DY)
%                                         % FR ;Fletcher-Reeves
%                                         % PR :Polak-Ribi\`ere
%                                         % HR: Hestenes-Stiefel
%                                         % DY :Dai-Yan
% end, UaOptimization parameters
% ------------]

% Some less often used parameters related to inversion 
CtrlVar.Inverse.InfoLevel=1;  % Set to 1 to get some basic information, >=2 for additional info on backtrackgin,
                              % >=100 for further info and plots
% In an inversion it it generally better to set other infolevels to a low value. So
% consider setting:
%CtrlVar.InfoLevelNonLinIt=1000; CtrlVar.InfoLevel=1000;                                           
                                                        

%%
CtrlVar.ThicknessConstraints=1;
CtrlVar.ResetThicknessToMinThickness=0;  % change this later on
CtrlVar.ThickMin=1;


end
