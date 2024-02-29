function FigureA2

version = 'final';
froot = ['/Volumes/mainJDeRydt/Ua/cases/PIG_IceShelf_Integrity/',version];

xmin = -1715; xmax = -1455;
ymin = -385; ymax = 0;

gs = {'10000','25000','50000'};

addpath(['/Volumes/mainJDeRydt/Ua/cases/PIG_IceShelf_Integrity/',version,'/Matlab_Tools']);

load(['/Volumes/mainJDeRydt/Ua/cases/PIG_IceShelf_Integrity/',version,'/NoPerturbationFinal_1996_gs25000_m3-RestartFile.mat'],'MUA','GF','CtrlVarInRestartFile');
CtrlVar = CtrlVarInRestartFile; CtrlVar.PlotXYscale = 1e3;
x = MUA.coordinates(:,1); y = MUA.coordinates(:,2); 
[Xm,Ym] = ndgrid([min(x):3e3:max(x)],[min(y):3e3:max(y)]);

CM = [[linspace(0.75,0,64)';linspace(0,0.75,64)'] linspace(0.75,0,128)' linspace(0.75,1,128)'];

%% PLOTTING
H = fig('units','inches','width',80*12/72.27,'height',40*12/72.27,'fontsize',16,'font','Helvetica');

for ii=1:length(gs)
    
    subplot('position',[0.08+0.28*(ii-1) 0.1 0.26 0.85]); hold on;
    
    load(['/Volumes/mainJDeRydt/Ua/cases/PIG_IceShelf_Integrity/',version,'/Matlab_Tools/Optimal_m_gs',gs{ii},'_NoRestrictionOndhdt_v2.mat']);
    Optimal_m(Optimal_m<1)=1;
    Optimal_m(Optimal_m==Inf)=9999; 
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
    IInf_m = find(Optimal_m_m(:)==9999 & GF_m(:)>0.5 & inpoly([Xm(:) Ym(:)],[MUA.Boundary.x(:) MUA.Boundary.y(:)]));
    INaN = find(isnan(Optimal_m(:)) & GF.node(:)>0.5);

    Mask = 0*MUA.coordinates(:,1); Mask(INaN)=0.8;

    Optimal_m = Optimal_m_orig;
    I = find(~isnan(Optimal_m));
    FOptimalm = scatteredInterpolant(x(I),y(I),Optimal_m(I),'nearest');
    Optimal_m = FOptimalm(MUA.coordinates(:,1),MUA.coordinates(:,2));

    I = find(GF.node<0.5); Optimal_m(I)=NaN;
    FOptimalC = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),Optimal_m(:));
    OptimalCm = FOptimalC(Xm,Ym);
    Ibound = find(inpoly([Xm(:) Ym(:)],[MUA.Boundary.x(:) MUA.Boundary.y(:)])==0);
    OptimalCm(Ibound) = NaN;

    PlotNodalBasedQuantities_JDR(MUA.connectivity,MUA.coordinates/1e3,Optimal_m(:)); hold on;
    
    [C,h] = contour(Xm/1e3,Ym/1e3,OptimalCm,[5 10 15 20],'EdgeColor',[0.4 0.4 0.4],'linewidth',1);
    clabel(C,h,'color',[0.4 0.4 0.4],'fontsize',12);

    plot(Xm(INaN_m)/1e3,Ym(INaN_m)/1e3,'ow','markersize',2,'markerfacecolor','w');
    plot(Xm(IInf_m)/1e3,Ym(IInf_m)/1e3,'ok','markersize',2,'markerfacecolor','k');
    
plot(MUA.Boundary.x/1e3,MUA.Boundary.y/1e3,'-k');PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'k');

    grid on; box on; axis equal;
    xlim([xmin xmax]); ylim([ymin ymax]);

    xlabel('psx [km]','Fontsize',16); 

    colormap(gca,CM); caxis([1 40]);
    
    if ii==1
        ylabel('psy [km]','Fontsize',16);
    else
        ylabel(''); set(gca,'yticklabel',{''});
    end
    
    if ii==3
        cb=colorbar(gca,'Position',[0.92 0.15 0.015 0.75],'Location','east','AxisLocation','out',...
         'Ticks',[1 10:10:40],'TickLabels',{'1','10','20','30','40'});
        cb.XColor='k';
        cb.YColor='k';
        cb.TickLength=0.04;
        cb.FontSize=16;
        cb.Label.String = 'Optimal value for $m$';%
        cb.Label.Interpreter = 'latex';
    end
    
    title(['$\gamma_s$ = ',gs{ii},' m'],'interpreter','latex');

end
      
pos = get(H,'Position');
set(H,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
fname = ['FigureA3_',version];
print(H,fname,'-dpng','-r400');