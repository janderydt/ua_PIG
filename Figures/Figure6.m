function Figure6

version = 'final';
%froot = ['/Volumes/mainJDeRydt/Ua/cases/PIG_IceShelf_Integrity/',version];
froot = ['/Users/janderydt/Documents/Work/',version];

xmin = -1715; xmax = -1455;
ymin = -385; ymax = 0;

gs = '25000'; m='3';

exp = {'Calving','IceShelfThinning','GLRetreat','GLRetreatCalving',};
%description = {'\begin{tabular}{c} Calving \end{tabular}',...
%    '\begin{tabular}{c} Thinning \\ (ice shelf only)\end{tabular}',...
%    '\begin{tabular}{c} Thinning\\ (grounded and floating)\end{tabular}',...
%    '\begin{tabular}{c} Thinning (grounded and floating)\\ and calving \end{tabular}'};

titlestr = '{\boldmath${\cal E}_{\rm CalvThin}^{\rm optimal}$}';%{'calving','ice shelf thinning','thinning','calving and thinning'};
barstr = {'${\cal E}_{\rm Calv}^{\rm optimal}$','${\cal E}_{\rm ISThin}^{\rm optimal}$ ',...
    '${\cal E}_{\rm Thin}^{\rm optimal}$','${\cal E}_{\rm CalvThin}^{\rm optimal}$'};
%labelstr = {'calving','ice shelf \newline thinning','retreat','calving and \newline retreat'};

addpath(['/Volumes/mainJDeRydt/Ua/cases/PIG_IceShelf_Integrity/',version,'/Matlab_Tools']);

load(['/Volumes/mainJDeRydt/Ua/cases/PIG_IceShelf_Integrity/',version,'/NoPerturbationFinal_1996_gs',gs,'_m',m,'-RestartFile.mat'],'MUA','GF');
MUA_1996 = MUA; GF_1996 = GF;

load([froot,'/NoPerturbationFinal_1996_gs',gs,'_Optimalm_NoRestrictionOndhdt_Interpolated-RestartFile.mat'],'MUA','F','GF','CtrlVarInRestartFile');

CtrlVar = CtrlVarInRestartFile; CtrlVar.PlotXYscale = 1e3;

% fluxgates
dL = 250;
FG(1).xends = [-1.619e6 -1.551e6]; FG(1).yends = [-2.143e5 -2.193e5];
FG(2).xends = [-1.621e6 -1.557e6]; FG(2).yends = [-2.51e5 -2.682e5];
for ii=1:length(FG)
    dx = FG(ii).xends(2)-FG(ii).xends(1);
    dy = FG(ii).yends(2)-FG(ii).yends(1);
    L = sqrt(dx.^2+dy.^2);
    n = round(L/dL);
    FG(ii).x = FG(ii).xends(1) + [0:n-1]*dx/n;
    FG(ii).y = FG(ii).yends(1) + [0:n-1]*dy/n;
    FG(ii).dL = L/n;
    rico = dy/dx; 
    FG(ii).en = [dy/L -dx/L];
    %plot(FG(ii).x,FG(ii).y,'-ok')
    %quiver(FG(ii).xends(1)+dx/2,FG(ii).yends(1)+dy/2,10000*FG(ii).en(1),10000*FG(ii).en(2),'w');
end

% CM=othercolor('RdYlBu11',220);
% CM=flipdim(CM,1);
% CM2 = [repmat(CM(20,:),4,1); CM(100:200,:); repmat(CM(220,:),4,1)];
% CM2(5,:) = [0.9 0.9 0.9];
% CM(end-9:end,:) = repmat([1 0 0],10,1);
% CM(1:10,:) = repmat([0 0 1],10,1);
% CM(110:111,:)=[1 1 1;1 1 1]*0.9;

% CM = [repmat([186 186 283]/360,8,1); othercolor('Reds8',200); repmat([232 124 124]/360,8,1)];
% CM(8:9,:) = [0.9 0.9 0.9;0.9 0.9 0.9];
% CM(208:209,:) = [0.9 0.9 0.9;0.9 0.9 0.9];
%CM = [repmat([286 186 283]/360,8,1); othercolor('Reds8',200); repmat([0.3 0.3 0.3],8,1);];
CM = [repmat([286 186 283]/360,8,1); othercolor('Reds8',200); othercolor('Greys3',200)];
CM(9,:) = [0.95 0.95 0.95];

%gs='50000';
load([froot,'/Inversion_MEaSUREs_2015_2016_Bedmachine_gs',gs,'_Optimalm_NoRestrictionOndhdt_Interpolated_InverseRestartFile.mat'],'MUA','F','GF','CtrlVarInRestartFile');
CtrlVar = CtrlVarInRestartFile; MUA_2015_2016 = MUA; GF_2015_2016=GF;
CtrlVar.PlotXYscale = 1e3;
speed_2015_2016 = sqrt(F.ub.^2+F.vb.^2);
Fspeed_2015_2016 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),speed_2015_2016);
Fu_2015_2016 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.ub);
Fv_2015_2016 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.vb);
Fh_2015_2016 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.h);
GLgeo_2015_2016=GLgeometry(MUA.connectivity,MUA.coordinates,GF,CtrlVar);
[Flux_2015_2016]=FluxAcrossGroundingLine_centralPIG(CtrlVar,MUA,F.ub,F.vb,0*F.ub,0*F.vb,F.h,F.rho,GLgeo_2015_2016);  
data.GLFlux_2015_2016 = sum(Flux_2015_2016(:))/1e9;
fprintf('Grounding line flux 2015-2106, gs=%s: %s Gt/yr \n',gs,num2str(sum(Flux_2015_2016(:))/1e9));

% flux gates
for ff=1:length(FG)
    uFG = Fu_2015_2016(FG(ff).x,FG(ff).y);
    vFG = Fv_2015_2016(FG(ff).x,FG(ff).y);
    speedFG_normal = FG(ff).en*[uFG; vFG];
    data.FGFlux_2015_2016(ff) = sum(FG(ff).dL.*speedFG_normal.*Fh_2015_2016(FG(ff).x,FG(ff).y)*F.rho(1)/1e9);
end

load([froot,'/NoPerturbationFinal_1996_gs',gs,'_Optimalm_NoRestrictionOndhdt_Interpolated-RestartFile.mat'],'MUA','F','GF','CtrlVarInRestartFile');
CtrlVar = CtrlVarInRestartFile; CtrlVar.PlotXYscale = 1e3;
MUA_1996 = MUA; GF_1996=GF;
speed_1996 = sqrt(F.ub.^2+F.vb.^2);
Fspeed_1996 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),speed_1996);
Fu_1996 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.ub);
Fv_1996 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.vb);
Fh_1996 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.h);
[Flux_1996]=FluxAcrossGroundingLine_centralPIG(CtrlVar,MUA,F.ub,F.vb,0*F.ub,0*F.vb,F.h,F.rho,GLgeo_2015_2016);
data.GLFlux_1996 = sum(Flux_1996(:))/1e9;
fprintf('Grounding line flux 1996, gs=%s: %s Gt/yr \n',gs,num2str(sum(Flux_1996(:))/1e9));

% flux gates
for ff=1:length(FG)
    uFG = Fu_1996(FG(ff).x,FG(ff).y);
    vFG = Fv_1996(FG(ff).x,FG(ff).y);
    speedFG_normal = FG(ff).en*[uFG; vFG];
    data.FGFlux_1996(ff) = sum(FG(ff).dL.*speedFG_normal.*Fh_1996(FG(ff).x,FG(ff).y)*F.rho(1)/1e9);
end

data.dGLFlux = data.GLFlux_2015_2016 - data.GLFlux_1996;

for ff=1:length(FG)
    data.dFGFlux(ff) = data.FGFlux_2015_2016(ff) - data.FGFlux_1996(ff);
end

for jj=1:length(exp)

    data.Perturbation(jj).exp = exp{jj};

    load([froot,'/PerturbationFinal_1996-',exp{jj},'_gs',gs,'_Optimalm_NoRestrictionOndhdt_Interpolated-RestartFile.mat'],'MUA','F','GF');
    data.Perturbation(jj).MUA = MUA;

    Fu = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.ub);
    Fv = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.vb);
    Fh = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.h);
    speed = sqrt(F.ub.^2+F.vb.^2);
    dspeed_pert_1996 = speed - Fspeed_1996(MUA.coordinates(:,1),MUA.coordinates(:,2));
    dspeed_2015_1996 = Fspeed_2015_2016(MUA.coordinates(:,1),MUA.coordinates(:,2))-Fspeed_1996(MUA.coordinates(:,1),MUA.coordinates(:,2));
    percspeedup = 100*dspeed_pert_1996./dspeed_2015_1996;

    MaskInd = find(Fspeed_1996(MUA.coordinates(:,1),MUA.coordinates(:,2))<100 | Fspeed_2015_2016(MUA.coordinates(:,1),MUA.coordinates(:,2))<100);
    percspeedup(MaskInd) = 0;
    Mask = 0*percspeedup; Mask(MaskInd)=1;

    data.Perturbation(jj).percspeedup = percspeedup;
    data.Perturbation(jj).Mask = Mask;

    [Flux,~,~,~,~,xGL_central,yGL_central]=FluxAcrossGroundingLine_centralPIG(CtrlVar,MUA,F.ub,F.vb,0*F.ub,0*F.vb,F.h,F.rho,GLgeo_2015_2016);
    data.Perturbation(jj).GLFlux = sum(Flux)/1e9;
    data.Perturbation(jj).dGLFlux = data.Perturbation(jj).GLFlux-data.GLFlux_1996;
    fprintf('Grounding line flux %s, gs=%s: %s Gt/yr -- Explains %s%% of the change \n',...
       exp{jj},gs,num2str(sum(Flux(:))/1e9),num2str(data.Perturbation(jj).dGLFlux/data.dGLFlux*100));

   % flux gates
    for ff=1:length(FG)
        uFG = Fu(FG(ff).x,FG(ff).y);
        vFG = Fv(FG(ff).x,FG(ff).y);
        speedFG_normal = FG(ff).en*[uFG; vFG];
        data.Perturbation(jj).FGFlux(ff) = sum(FG(ff).dL.*speedFG_normal.*Fh(FG(ff).x,FG(ff).y)*F.rho(1)/1e9);
        data.Perturbation(jj).dFGFlux(ff) = data.Perturbation(jj).FGFlux(ff) - data.FGFlux_1996(ff);
    end

end


%% PLOTTING
H = fig('units','inches','width',63*12/72.27,'height',60*12/72.27,'fontsize',16,'font','Helvetica');

subplot('position',[0.1 0.5 0.42 0.47]); hold on;

load(['/Volumes/mainJDeRydt/Ua/cases/PIG_IceShelf_Integrity/',version,'/Matlab_Tools/Optimal_m_gs',gs,'_NoRestrictionOndhdt_v2.mat']);
x = MUA.coordinates(:,1); y = MUA.coordinates(:,2); 
[Xm,Ym] = ndgrid([min(x):3e3:max(x)],[min(y):3e3:max(y)]);
Optimal_m(Optimal_m==Inf)=9999;

I = find(Rsquare<0.9);
Optimal_m(I)=NaN;
% I = find(Rsquare<0.95);
% Optimal_m = FOptimal(MUA.coordinates(:,1),MUA.coordinates(:,2)); 
FOptimalm = scatteredInterpolant(x,y,Optimal_m,'nearest');
Optimal_m_m = FOptimalm(Xm,Ym);
Optimal_m_1996 = FOptimalm(MUA_1996.coordinates(:,1),MUA_1996.coordinates(:,2));
FMask = scatteredInterpolant(MUA_1996.coordinates(:,1),MUA_1996.coordinates(:,2),GF_1996.node); 
GF_1996_m = FMask(Xm,Ym);

INaN_m = find(isnan(Optimal_m_m(:)) & GF_1996_m(:)>0.5 & inpoly([Xm(:) Ym(:)],[MUA_1996.Boundary.x(:) MUA_1996.Boundary.y(:)]));
IInf_m = find(Optimal_m_m(:)==9999 & GF_1996_m(:)>0.5 & inpoly([Xm(:) Ym(:)],[MUA_1996.Boundary.x(:) MUA_1996.Boundary.y(:)]));
INaN = find(isnan(Optimal_m_1996(:)) & GF_1996.node(:)>0.5);

Mask = 0*MUA_1996.coordinates(:,1); Mask(INaN)=0.8;


I = find(~isnan(Optimal_m));
FOptimalm = scatteredInterpolant(x(I),y(I),Optimal_m(I),'nearest');
Optimal_m_1996 = FOptimalm(MUA_1996.coordinates(:,1),MUA_1996.coordinates(:,2));

x = MUA_1996.coordinates(:,1); y = MUA_1996.coordinates(:,2);
Optimal_m = Optimal_m_1996;
fname = ['Optimal_m_gs25000_NoRestrictionOndhdt_Interpolated_',version,'_v2.mat'];
save(fname,'x','y','Optimal_m');

I = find(GF_1996.node<0.5); Optimal_m_1996(I)=NaN;
FOptimalC = scatteredInterpolant(MUA_1996.coordinates(:,1),MUA_1996.coordinates(:,2),Optimal_m_1996(:));
OptimalCm = FOptimalC(Xm,Ym);
Ibound = find(inpoly([Xm(:) Ym(:)],[MUA_1996.Boundary.x(:) MUA_1996.Boundary.y(:)])==0);
OptimalCm(Ibound) = NaN;

% save geotiff
OptimalCm_temp = OptimalCm;
S = shaperead('Optimalm_poly.shp');
Ibound = find(~inpoly([Xm(:) Ym(:)],[S.X(1:end-1)' S.Y(1:end-1)']));
OptimalCm_temp(Ibound) = NaN;
R = makerefmat(Xm(1,1),Ym(1,1),Xm(2,1)-Xm(1,1),Ym(1,2)-Ym(1,1));
geotiffwrite(['./Optimal_m_1996'],OptimalCm_temp',R,'coordRefSysCode','EPSG:3031');

return

PlotNodalBasedQuantities_JDR(MUA_1996.connectivity,MUA_1996.coordinates/1e3,Optimal_m_1996(:)); hold on;
% opengl hardware
% patch('faces',MUA_1996.connectivity,'vertices',MUA_1996.coordinates/1e3,'facecolor','w','EdgeColor','none','FaceAlpha','flat',...
%     'FaceVertexAlphaData',Mask(:),'AlphaDataMapping','none');

[C,h] = contour(Xm/1e3,Ym/1e3,OptimalCm,[5 10 15 20],'EdgeColor',[0.4 0.4 0.4],'linewidth',1);
clabel(C,h,'color',[0.4 0.4 0.4],'fontsize',12);
% [C,h] = contour(Xm/1e3,Ym/1e3,OptimalCm,[3 3],'EdgeColor',[0 0 0],'linewidth',1);
% clabel(C,h,'manual','color',[0 0 0],'fontsize',16);

%plot(contourdata(Imaxcont).x,contourdata(Imaxcont).y,'-m','linewidth',3);
plot(Xm(INaN_m)/1e3,Ym(INaN_m)/1e3,'ow','markersize',2,'markerfacecolor','w');
plot(Xm(IInf_m)/1e3,Ym(IInf_m)/1e3,'ok','markersize',2,'markerfacecolor','k');

plot(MUA_1996.Boundary.x/1e3,MUA_1996.Boundary.y/1e3,'-k');
PlotGroundingLines(CtrlVar,MUA_1996,GF_1996,[],[],[],'k');
PlotGroundingLines(CtrlVar,MUA_2015_2016,GF_2015_2016,[],[],[],'k');

grid on; box on; axis equal;
xlim([xmin xmax]); ylim([ymin ymax]);

xlabel('psx [km]','Fontsize',16); ylabel('psy [km]','Fontsize',16);   
        
CM3 = [[linspace(0.75,0,64)';linspace(0,0.75,64)'] linspace(0.75,0,128)' linspace(0.75,1,128)'];
colormap(gca,CM3); caxis([1 40]);

cb = colorbar(gca,'Ticks',[1 5:5:40],'TickLabels',{'1','','10','','20','','30','','40'});
cb.Label.String = 'Optimal value for the sliding exponent $m$';
cb.Label.Interpreter = 'latex';
cb.FontSize=16;

subplot('position',[0.53 0.5 0.42 0.45]); hold on;

PlotNodalBasedQuantities_JDR(data.Perturbation(jj).MUA.connectivity,...
        data.Perturbation(jj).MUA.coordinates/1e3,data.Perturbation(jj).percspeedup); hold on;

caxis([-4 200]);
colormap(gca,CM); 

if strfind(data.Perturbation(jj).exp,'Calving')
    plot(MUA_1996.Boundary.x/1e3,MUA_1996.Boundary.y/1e3,'--k');   
end
plot(data.Perturbation(jj).MUA.Boundary.x/1e3,data.Perturbation(jj).MUA.Boundary.y/1e3,'-k');
for ff=1:length(FG)
        plot(FG(ff).x/1e3,FG(ff).y/1e3,'-k','linewidth',2)     
end
%plot(xGL_central/1e3,yGL_central/1e3,'-m','linewidth',2);

if strfind(data.Perturbation(jj).exp,'Retreat')
    PlotGroundingLines(CtrlVar,MUA_2015_2016,GF_2015_2016,[],[],[],'-k','linewidth',1);  
else
    PlotGroundingLines(CtrlVar,MUA_1996,GF_1996,[],[],[],'-k','linewidth',1);
end

ylabel(''); yticklabels({''});
xlabel('psx [km]');

text(FG(1).xends(2)/1e3+20,FG(1).yends(2)/1e3,'Gate 1','VerticalAlignment','middle','HorizontalAlignment','left');
text(FG(2).xends(2)/1e3+20,FG(2).yends(2)/1e3,'Gate 2','VerticalAlignment','middle','HorizontalAlignment','left');

% cb=colorbar(gca,'Position',[0.40 0.94 0.24 0.015],'Location','northoutside','AxisLocation','in',...
%     'Ticks',[-100:25:100],'TickLabels',{'-100','','-50','','0','','50','','100'});
% cb.XColor='k';
% cb.YColor='k';
% cb.TickLength=0.04;
cb = colorbar(gca);
cb.FontSize=16;
cb.Label.String = {'$(U_{\rm CalvThin}^{\rm optimal}-U_{96})/(U_{16}-U_{96})\,[\%]$'};
cb.Label.Interpreter = 'latex';

title(titlestr,'interpreter','latex');

grid on; box on; axis equal;
xlim([xmin xmax]); ylim([ymin ymax]);

subplot('position',[0.15 0.07 0.8 0.35]); hold on;

for ff=1:length(FG)
    for jj=1:4
       perc = data.Perturbation(jj).dFGFlux(ff)/data.dFGFlux(ff)*100;
       patch([0 0 perc perc],[(5-jj)-(-1)^ff*0.35 (5-jj) (5-jj) (5-jj)-(-1)^ff*0.35],[0.5 0.5 1]);
       if perc<10
            text(perc/2,5-jj-(-1)^ff*0.175,[num2str(perc,1),'%'],'color','w','HorizontalAlignment','center');       
       elseif perc>=10 & perc<100
            text(perc/2,5-jj-(-1)^ff*0.175,[num2str(perc,2),'%'],'color','w','HorizontalAlignment','center');
       elseif perc>=100
           text(perc/2,5-jj-(-1)^ff*0.175,[num2str(perc,3),'%'],'color','w','HorizontalAlignment','center');
       end
       if jj==1
           text(perc+5,5-jj-(-1)^ff*0.175,['Gate ',num2str(ff)],'color','k','HorizontalAlignment','left');
       end
%        if jj==4
%            patch([perc perc 100 100],[(5-jj)-(-1)^ff*0.35 (5-jj) (5-jj) (5-jj)-(-1)^ff*0.35],[1 0.5 0.5]);
%        end
%        if jj==4 & ff==1
%            text(perc/2+50,5-jj-(-1)^ff*0.175,'unaccounted for','color','k','HorizontalAlignment','center');
%        end
    end
end

plot([0 0],[0 length(exp)+1],'--k');
plot([100 100],[0 length(exp)+1],'--k');
ylim([0.5 length(exp)+0.5]);
yticks(gca,[1:length(exp)]); set(gca,'yticklabel',flipdim(barstr,2),'ticklabelinterpreter','latex');
xlim([-10 110])
xticks(gca,[0 100]);
set(gca,'xticklabels',{'1996','2016'});
xlabel('Flux through gate');
box on;
      
pos = get(H,'Position');
set(H,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
fname = ['Figure6_',version];
print(H,fname,'-dpng','-r400');