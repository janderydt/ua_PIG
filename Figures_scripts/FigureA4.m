function FigureA3

version = 'final';
froot = ['/Volumes/mainJDeRydt/Ua/cases/PIG_IceShelf_Integrity/',version];

xmin = -1715; xmax = -1455;
ymin = -385; ymax = 0;

gs = '25000'; m='3';

barstr = {'${\cal E}_{\rm Calv}^{\rm optimal}$','${\cal E}_{\rm ISThin}^{\rm optimal}$ ',...
    '${\cal E}_{\rm Thin}^{\rm optimal}$','${\cal E}_{\rm CalvThin}^{\rm optimal}$'};

addpath(['/Volumes/mainJDeRydt/Ua/cases/PIG_IceShelf_Integrity/',version,'/Matlab_Tools']);

load([froot,'/Inversion_MEaSUREs_2015_2016_Bedmachine_gs',gs,'_m',m,'_InverseRestartFile.mat']);
MUA_2016 = MUA; GF_2016 = GF;

load(['/Volumes/mainJDeRydt/Ua/cases/PIG_IceShelf_Integrity/',version,'/NoPerturbationFinal_1996_gs',gs,'_m',m,'-RestartFile.mat'],'MUA','GF','CtrlVarInRestartFile');
load(['/Volumes/mainJDeRydt/Ua/cases/PIG_IceShelf_Integrity/',version,'/Matlab_Tools/Optimal_m_gs',gs,'_NoRestrictionOndhdt.mat']);

CtrlVar = CtrlVarInRestartFile; CtrlVar.PlotXYscale = 1e3;
x = MUA.coordinates(:,1); y = MUA.coordinates(:,2); 
[Xm,Ym] = ndgrid([min(x):3e3:max(x)],[min(y):3e3:max(y)]);

%CM = othercolor('Reds8',200);
CM = [[linspace(0.75,0,64)';linspace(0,0.75,64)'] linspace(0.75,0,128)' linspace(0.75,1,128)'];

xpoint = [-1582.0239 -1581.4509 -1602.5465];
ypoint = [-165.4126 -227.3969 -245.3878];

%% PLOTTING
H = fig('units','inches','width',63*12/72.27,'height',50*12/72.27,'fontsize',16,'font','Helvetica');

subplot('position',[0.1 0.1 0.42 0.85]); hold on;

load(['/Volumes/mainJDeRydt/Ua/cases/PIG_IceShelf_Integrity/',version,'/Matlab_Tools/Optimal_m_gs',gs,'_NoRestrictionOndhdt_v2.mat']);
Optimal_m(Optimal_m<1)=1;
Optimal_m(Optimal_m>50)=50;
I = find(Rsquare<0.9);
Optimal_m(I)=NaN;
Optimal_m_orig = Optimal_m;
% I = find(Rsquare<0.95);
% Optimal_m = FOptimal(MUA.coordinates(:,1),MUA.coordinates(:,2)); 
FOptimalm = scatteredInterpolant(x,y,Optimal_m,'nearest');
Optimal_m_m = FOptimalm(Xm,Ym);
Optimal_m = FOptimalm(MUA.coordinates(:,1),MUA.coordinates(:,2));
FMask = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),GF.node); 
GF_m = FMask(Xm,Ym);
INaN_m = find(isnan(Optimal_m_m(:)) & GF_m(:)>0.5 & inpoly([Xm(:) Ym(:)],[MUA.Boundary.x(:) MUA.Boundary.y(:)]));
INaN = find(isnan(Optimal_m(:)) & GF.node(:)>0.5);

Mask = 0*MUA.coordinates(:,1); Mask(INaN)=0.8;

Optimal_m = Optimal_m_orig;
I = find(~isnan(Optimal_m));
FOptimalm = scatteredInterpolant(x(I),y(I),Optimal_m(I),'nearest');
Optimal_m = FOptimalm(MUA_2016.coordinates(:,1),MUA_2016.coordinates(:,2));

I = find(GF_2016.node<0.5); Optimal_m(I)=NaN;
FOptimalC = scatteredInterpolant(MUA_2016.coordinates(:,1),MUA_2016.coordinates(:,2),Optimal_m(:));
OptimalCm = FOptimalC(Xm,Ym);
Ibound = find(inpoly([Xm(:) Ym(:)],[MUA.Boundary.x(:) MUA.Boundary.y(:)])==0);
OptimalCm(Ibound) = NaN;

I1 = find(Rsquare>=0.9 & p1>=100);
I2 = find(Rsquare>=0.9 & p1<100);
I3 = find(Rsquare>=0.9 & p1<0);
J = find(Rsquare<0.9);

%PlotNodalBasedQuantities_JDR(MUA_2016.connectivity,MUA_2016.coordinates/1e3,Optimal_m(:)); hold on;
plot(x(I1)/CtrlVar.PlotXYscale,y(I1)/CtrlVar.PlotXYscale,'.r');
plot(x(I2)/CtrlVar.PlotXYscale,y(I2)/CtrlVar.PlotXYscale,'.k');
plot(x(J)/CtrlVar.PlotXYscale,y(J)/CtrlVar.PlotXYscale,'.','color',[0.7 0.7 0.7]);
g(1) = plot(0,0,'.r','markersize',10);
g(2) = plot(0,0,'.k','markersize',10);
g(3) = plot(0,0,'.','color',[0.7 0.7 0.7],'markersize',10);
%plot(x(I3)/CtrlVar.PlotXYscale,y(I3)/CtrlVar.PlotXYscale,'.k');

%[C,h] = contour(Xm/1e3,Ym/1e3,OptimalCm,[5 10 15 20],'EdgeColor',[0.4 0.4 0.4],'linewidth',1);
%clabel(C,h,'color',[0.4 0.4 0.4],'fontsize',12);

%plot(Xm(INaN_m)/1e3,Ym(INaN_m)/1e3,'ow','markersize',2,'markerfacecolor','w');

plot(MUA.Boundary.x/1e3,MUA.Boundary.y/1e3,'-k');
PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'k');
PlotGroundingLines(CtrlVar,MUA_2016,GF_2016,[],[],[],'k');

for ii=1:length(xpoint)
   plot(xpoint(ii),ypoint(ii),'ow','markersize',7,'markerfacecolor','w');
   text(xpoint(ii)+5,ypoint(ii),num2str(ii),'horizontalalignment','left','color','w');
end

grid on; box on; axis equal;
xlim([xmin xmax]); ylim([ymin ymax]);

xlabel('psx [km]','Fontsize',16); 

colormap(gca,CM); caxis([1 40]);

ylabel('psy [km]','Fontsize',16);

legend(g(:),{'$R^2 \geq 0.9,\, f_1\geq 100$','$R^2 \geq 0.9,\, f_1 < 100$','$R^2 < 0.9$'},'interpreter','latex','location','northwest');

marray = [cellfun(@str2num,m)]';
mfine = [0.1:0.1:40];

for ii=1:length(xpoint)
    
    subplot('position',[0.6 0.07+(3-ii)*0.31 0.38 0.27]); hold on;
    
    I = find(x/1e3>xpoint(ii)-0.1 & x/1e3<xpoint(ii)+0.1 & y/1e3>ypoint(ii)-0.1 & y/1e3<ypoint(ii)+0.1);

    g(1)=plot(marray,percspeedup(I,:),'ob');
    g(2)=plot(mfine,p1(I)*mfine./(mfine+q1(I)),'-k');
    
    xticks([0:5:40]);
    xticklabels({'0','','10','','20','','30','','40'});
    if ii==1
        xticklabels({''});
        mopt = 100*q1(I)/(p1(I)-100);
        plot(mopt,100,'dk','markerfacecolor','k','markersize',10);
        plot([mopt mopt],[0 100],'--k');
        plot([1 mopt],[100 100],'--k');
        text(100*q1(I)/(p1(I)-100)+1.5,100-8,['$m_{\rm optimal}$=',num2str(mopt,2)],'interpreter','latex','horizontalalignment','left');
        legend(g(:),{'model output','best fit to $\frac{f_1 m}{m+f_2}$'},'location','southeast','interpreter','latex');
    elseif ii==2
        ylim([0 110]);
        plot([1 40],[p1(I) p1(I)],'--k');
        text(15,p1(I)+7,'$f_1$: horizontal asymptote','interpreter','latex');
        xticklabels({''});
    elseif ii==3
        xlabel('$m$','interpreter','latex');
        ylim([-45 9]);
    end
    xlim([0 40]);
    title(['\bf\boldmath Location ',num2str(ii),' ($R^2$=',num2str(Rsquare(I),2),', $f_1$=',num2str(p1(I),3),')'],'fontsize',14,'interpreter','latex');

    ylabel('$\Delta U_{\rm CalvThin}/\Delta U$ [\%]','interpreter','latex');

    grid on;
    box on;

end

pos = get(H,'Position');
set(H,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
fname = ['FigureA4_',version];
print(H,fname,'-dpng','-r400');