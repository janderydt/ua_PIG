function Figure5

version = 'final'; %options: v1, v2, final
froot = ['/Volumes/mainJDeRydt/Ua/cases/PIG_IceShelf_Integrity/',version];
xmin = -1715; xmax = -1440;
ymin = -385; ymax = 0;

gs = {'25000','25000','25000'};
m = {'m1','m7','m13'};

exp = {'GLRetreatCalving',};
%description = {'\begin{tabular}{c} Calving \end{tabular}',...
%    '\begin{tabular}{c} Thinning \\ (ice shelf only)\end{tabular}',...
%    '\begin{tabular}{c} Thinning\\ (grounded and floating)\end{tabular}',...
%    '\begin{tabular}{c} Thinning (grounded and floating)\\ and calving \end{tabular}'};
mstr = {'{\boldmath${\cal E}_{\rm CalvThin}^1$}','{\boldmath${\cal E}_{\rm CalvThin}^7$}','{\boldmath${\cal E}_{\rm CalvThin}^{13}$}'};
%labelstr = {'calving','ice shelf \newline thinning','retreat','calving and \newline retreat'};

addpath(['/Volumes/mainJDeRydt/Ua/cases/PIG_IceShelf_Integrity/',version,'/Matlab_Tools']);

load([froot,'/Inversion_MEaSUREs_2015_2016_Bedmachine_gs',gs{1},'_',m{1},'_InverseRestartFile.mat'],'MUA','F','GF','CtrlVarInRestartFile');

%figure; PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,sqrt(F.ub.^2+F.vb.^2));
%hold on; PlotGroundingLines(CtrlVarInRestartFile,MUA,GF);

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
% CM(end-9:end,:) = repmat([1 0 0],10,1);
% CM(1:10,:) = repmat([0 0 1],10,1);
% CM(110:111,:)=[1 1 1;1 1 1]*0.9;
%CM=othercolor('RdYlBu11',400);
%CM=flipdim(CM,1);
%CM2 = [repmat([186 186 283]/360,4,1); othercolor('YlOrRd6',100); repmat([232 124 124]/360,4,1)];%[0.4039 0 0.051]

% CM = [repmat([186 186 283]/360,8,1); othercolor('Reds8',200); repmat([232 124 124]/360,8,1)];
% CM(8:9,:) = [0.9 0.9 0.9;0.9 0.9 0.9];
% CM(208:209,:) = [0.9 0.9 0.9;0.9 0.9 0.9];
CM = [repmat([286 186 283]/360,8,1); othercolor('Reds8',200); othercolor('Greys3',200)];
CM(8:9,:) = [0.9 0.9 0.9;0.9 0.9 0.9];


if exist([froot,'/Matlab_Tools/FluxData.mat'])
    
    load([froot,'/Matlab_Tools/FluxData.mat']);
    
else
    
    for ii=1:length(m)

        data(ii).gs = gs{ii};
        data(ii).m = m{ii};

        load([froot,'/Inversion_MEaSUREs_2015_2016_Bedmachine_gs',gs{ii},'_',m{ii},'_InverseRestartFile.mat'],'MUA','F','GF','CtrlVarInRestartFile');
        CtrlVar = CtrlVarInRestartFile; MUA_2015_2016 = MUA; GF_2015_2016=GF;
        CtrlVar.PlotXYscale = 1e3;
        speed_2015_2016 = sqrt(F.ub.^2+F.vb.^2);
        Fspeed_2015_2016 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),speed_2015_2016);
        Fu_2015_2016 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.ub);
        Fv_2015_2016 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.vb);
        Fh_2015_2016 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.h);
        GLgeo_2015_2016=GLgeometry(MUA.connectivity,MUA.coordinates,GF,CtrlVar);
        [Flux_2015_2016]=FluxAcrossGroundingLine_centralPIG(CtrlVar,MUA,F.ub,F.vb,0*F.ub,0*F.vb,F.h,F.rho,GLgeo_2015_2016);  
        data(ii).GLFlux_2015_2016 = sum(Flux_2015_2016(:))/1e9;
        fprintf('Grounding line flux 2015-2106, gs=%s: %s Gt/yr \n',gs{ii},num2str(sum(Flux_2015_2016(:))/1e9));
        
        % flux gates
        for ff=1:length(FG)
            uFG = Fu_2015_2016(FG(ff).x,FG(ff).y);
            vFG = Fv_2015_2016(FG(ff).x,FG(ff).y);
            speedFG_normal = FG(ff).en*[uFG; vFG];
            data(ii).FGFlux_2015_2016(ff) = sum(FG(ff).dL.*speedFG_normal.*Fh_2015_2016(FG(ff).x,FG(ff).y)*F.rho(1)/1e9);
        end

        load([froot,'/NoPerturbationFinal_1996_gs',gs{ii},'_',m{ii},'-RestartFile.mat'],'MUA','F','GF','CtrlVarInRestartFile');
        CtrlVar = CtrlVarInRestartFile; CtrlVar.PlotXYscale = 1e3;
        MUA_1996 = MUA; GF_1996=GF;
        speed_1996 = sqrt(F.ub.^2+F.vb.^2);
        Fspeed_1996 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),speed_1996);
        Fu_1996 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.ub);
        Fv_1996 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.vb);
        Fh_1996 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),F.h);
        [Flux_1996]=FluxAcrossGroundingLine_centralPIG(CtrlVar,MUA,F.ub,F.vb,0*F.ub,0*F.vb,F.h,F.rho,GLgeo_2015_2016);
        data(ii).GLFlux_1996 = sum(Flux_1996(:))/1e9;
        fprintf('Grounding line flux 1996, gs=%s: %s Gt/yr \n',gs{ii},num2str(sum(Flux_1996(:))/1e9));

        % flux gates
        for ff=1:length(FG)
            uFG = Fu_1996(FG(ff).x,FG(ff).y);
            vFG = Fv_1996(FG(ff).x,FG(ff).y);
            speedFG_normal = FG(ff).en*[uFG; vFG];
            data(ii).FGFlux_1996(ff) = sum(FG(ff).dL.*speedFG_normal.*Fh_1996(FG(ff).x,FG(ff).y)*F.rho(1)/1e9);
        end
        
        data(ii).dGLFlux = data(ii).GLFlux_2015_2016 - data(ii).GLFlux_1996;
        for ff=1:length(FG)
            data(ii).dFGFlux(ff) = data(ii).FGFlux_2015_2016(ff) - data(ii).FGFlux_1996(ff);
        end
        
        for jj=1:length(exp)

            data(ii).Perturbation(jj).exp = exp{jj};

            load([froot,'/PerturbationFinal_1996-',exp{jj},'_gs',gs{ii},'_',m{ii},'-RestartFile.mat'],'MUA','F','GF');
            data(ii).Perturbation(jj).MUA = MUA;
            
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

            data(ii).Perturbation(jj).percspeedup = percspeedup;
            data(ii).Perturbation(jj).Mask = Mask;

            [Flux,~,~,~,~,xGL_central,yGL_central]=FluxAcrossGroundingLine_centralPIG(CtrlVar,MUA,F.ub,F.vb,0*F.ub,0*F.vb,F.h,F.rho,GLgeo_2015_2016);
            data(ii).Perturbation(jj).GLFlux = sum(Flux)/1e9;
            data(ii).Perturbation(jj).dGLFlux = data(ii).Perturbation(jj).GLFlux-data(ii).GLFlux_1996;
            fprintf('Grounding line flux %s, gs=%s: %s Gt/yr -- Explains %s%% of the change \n',...
               exp{jj},gs{ii},num2str(sum(Flux(:))/1e9),num2str(data(ii).Perturbation(jj).dGLFlux/data(ii).dGLFlux*100));
           
           % flux gates
            for ff=1:length(FG)
                uFG = Fu(FG(ff).x,FG(ff).y);
                vFG = Fv(FG(ff).x,FG(ff).y);
                speedFG_normal = FG(ff).en*[uFG; vFG];
                data(ii).Perturbation(jj).FGFlux(ff) = sum(FG(ff).dL.*speedFG_normal.*Fh(FG(ff).x,FG(ff).y)*F.rho(1)/1e9);
                data(ii).Perturbation(jj).dFGFlux(ff) = data(ii).Perturbation(jj).FGFlux(ff) - data(ii).FGFlux_1996(ff);
            end
           
        end
    end
    
    fname = [froot,'/Matlab_Tools/FluxData.mat'];
    %save(fname,'data');
    
end

%%%%%%%%%%%%%%
%% Plotting %%
%%%%%%%%%%%%%%

%%%
% Perturbation experiments (m=3) and bar chart
%%%

H = fig('units','inches','width',60*12/72.27,'height',56*12/72.27,'fontsize',16,'font','Helvetica');

%% calving, ice shelf thinning, thinning
jj=1;
for ii=1:3
    
    subplot('position',[0.1+0.29*(ii-1) 0.38 0.27 0.54]); hold on;

    PlotNodalBasedQuantities_JDR(data(ii).Perturbation(jj).MUA.connectivity,...
        data(ii).Perturbation(jj).MUA.coordinates/1e3,data(ii).Perturbation(jj).percspeedup); hold on;
%     [C,h]=tricontour(data(ii).Perturbation(jj).MUA.connectivity,data(ii).Perturbation(jj).MUA.coordinates(:,1)/1e3,...
%         data(ii).Perturbation(jj).MUA.coordinates(:,2)/1e3,data(ii).Perturbation(jj).Mask,[0.99 0.89]);
%     for ll=1:length(h)
%         h(ll).EdgeColor = [0.5 0.5 0.5];
%     end

    caxis([-4 200]);
    colormap(CM); 
    
    x = data(ii).Perturbation(jj).MUA.coordinates(:,1); y = data(ii).Perturbation(jj).MUA.coordinates(:,2); 
    [Xm,Ym] = ndgrid([min(x):500:max(x)],[min(y):500:max(y)]);
    Fpercspeedup = scatteredInterpolant(x,y,data(ii).Perturbation(jj).percspeedup);
    percspeedup_m = Fpercspeedup(Xm,Ym); I = find(~inpoly([Xm(:) Ym(:)],[data(ii).Perturbation(jj).MUA.Boundary.x(:) data(ii).Perturbation(jj).MUA.Boundary.y(:)]));
    percspeedup_m(I)=NaN;
    contour(Xm/1e3,Ym/1e3,percspeedup_m,[50 50],'edgecolor','k','linestyle','--');
    contour(Xm/1e3,Ym/1e3,percspeedup_m,[0 100],'edgecolor','w');
        
    if strfind(data(ii).Perturbation(jj).exp,'Calving')
        plot(MUA_1996.Boundary.x/1e3,MUA_1996.Boundary.y/1e3,'--k');   
    end
    plot(data(ii).Perturbation(jj).MUA.Boundary.x/1e3,data(ii).Perturbation(jj).MUA.Boundary.y/1e3,'-k');
    for ff=1:length(FG)
            plot(FG(ff).x/1e3,FG(ff).y/1e3,'-k','linewidth',2)     
    end
    %plot(xGL_central/1e3,yGL_central/1e3,'-m','linewidth',2);

    if strfind(data(ii).Perturbation(jj).exp,'Retreat')
        PlotGroundingLines(CtrlVar,MUA_2015_2016,GF_2015_2016,[],[],[],'-k','linewidth',1);  
    else
        PlotGroundingLines(CtrlVar,MUA_1996,GF_1996,[],[],[],'-k','linewidth',1);
    end

    if ii==1
        ylabel('psy [km]');
        text(FG(1).xends(2)/1e3+20,FG(1).yends(2)/1e3,'Gate 1','VerticalAlignment','middle','HorizontalAlignment','left');
        text(FG(2).xends(2)/1e3+20,FG(2).yends(2)/1e3,'Gate 2','VerticalAlignment','middle','HorizontalAlignment','left');
    else
        yticklabels({});
    end
    xlabel('psx [km]');
    
    if jj==1
%         cb=colorbar(gca,'Position',[0.9 0.62 0.015 0.33],'Location','east','AxisLocation','in',...
%          'Ticks',[-100:25:100],'TickLabels',{'-100','','-50','','0','','50','','100'});
%         cb.XColor='k';
%         cb.YColor='k';
%         cb.TickLength=0.04;
%         cb.FontSize=16;
%         cb.Label.String = {'100(U_{2016}-U_{1996})/U_{2016} [%]'};
        cb=colorbar(gca,'Position',[0.35 0.90 0.34 0.015],'Location','northoutside','AxisLocation','out',...
         'Ticks',[0:50:200],'TickLabels',{'0','50','100','150','200'});
        cb.XColor='k';
        cb.YColor='k';
        cb.TickLength=0.04;
        cb.FontSize=16;
        cb.Label.String = '$(U_{\rm CalvThin}^m-U_{96})/(U_{16}-U_{96}){\rm\, [\%]}$';%
        cb.Label.Interpreter = 'latex';
        %cb.Label.String = {''};

    end
    
    grid on; box on; axis equal;
    xlim([xmin xmax]); ylim([ymin ymax]);

    title(mstr{ii},'interpreter','latex');

end

%cbarrow;

%% bar chart GL flux
subplot('position',[0.12 0.08 0.84 0.29]); hold on;

jj = 1; 
for ff=1:length(FG)
    for ii=1:length(m)
       perc = data(ii).Perturbation(jj).dFGFlux(ff)/data(ii).dFGFlux(ff)*100;
       patch([0 0 perc perc],[(4-ii)-(-1)^ff*0.35 (4-ii) (4-ii) (4-ii)-(-1)^ff*0.35],[0.5 0.5 1]);
       if perc>10
           if perc<100
                text(perc/2,4-ii-(-1)^ff*0.175,[num2str(perc,2),'%'],'color','w','HorizontalAlignment','center');
           else
               text(perc/2,4-ii-(-1)^ff*0.175,[num2str(perc,3),'%'],'color','w','HorizontalAlignment','center');
           end
       end
       if ii==1
           text(-11,4-ii-(-1)^ff*0.175,['Gate ',num2str(ff)],'color','k','HorizontalAlignment','left');
       end      
    end
end

% plot([data(ii).FGFlux_1996(ff) data(ii).FGFlux_1996(ff)],[0 length(exp)+1],'--k');
% plot([data(ii).FGFlux_2015_2016(ff) data(ii).FGFlux_2015_2016(ff)],[0 length(exp)+1],'--k');
plot([0 0],[0 length(m)+1],'--k');
plot([100 100],[0 length(m)+1],'--k');
ylim([0 length(m)+1]);
yticks(gca,[1:length(m)]); set(gca,'yticklabel',flipdim(mstr,2),'ticklabelinterpreter','latex');
%xticks(gca,[data(ii).FGFlux_1996(ff) data(ii).FGFlux_2015_2016(ff)]); 
xlim([-15 105])
xticks(gca,[0 100]);
set(gca,'xticklabels',{'1996','2016'});
xlabel('Flux through gate');
box on;

        
pos = get(H,'Position');
set(H,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
fname = ['Figure5_',version];
print(H,fname,'-dpng','-r400');

return

%%%
% Different m, bar chart
%%%


        PlotNodalBasedQuantities_JDR(MUA.connectivity,MUA.coordinates/1e3,dspeed_pert_1996); hold on;
        caxis([-800 800]);
        xlim([xmin xmax]); ylim([ymin ymax]);
        grid on; box on; axis equal;

        xlabel('psx [km]','Fontsize',16); ylabel('psy [km]','Fontsize',16);
        plot(MUA0.Boundary.x/1e3,MUA0.Boundary.y/1e3,'--k');
        plot(MUA.Boundary.x/1e3,MUA.Boundary.y/1e3,'-k');
        PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'k');

        if strfind(exp{ii},'GLRetreat')
            PlotGroundingLines(CtrlVar,MUA0,GF0,[],[],[],'--k');
        end

        subplot('position',[0.53 0.1 0.42 0.8]); hold on;
        PlotNodalBasedQuantities_JDR(MUA.connectivity,MUA.coordinates/1e3,dspeed_2015_1996); hold on;
        caxis([-800 800]);       
        yticklabels(gca,{}); ylabel(''); xlabel('psx [km]','Fontsize',16); 
        cblabel = 'Change in surface speed [m/yr]';


hold on; plot(MUA_1996.Boundary.x/1e3,MUA_1996.Boundary.y/1e3,'--k');
plot(MUA.Boundary.x/1e3,MUA.Boundary.y/1e3,'-k');

PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'k');

if strfind(exp{jj},'Retreat')
    PlotGroundingLines(CtrlVar,MUA_1996,GF_1996,[],[],[],'--k');
end

htitle = sgtitle(['{\textbf{Response to $\left|',description{jj},'\right|$ (gs = ',gs{ii},', m = ',m{ii},')}}'],'Interpreter','latex');
htitle.FontSize= 20;

grid on; box on; axis equal;
xlim([xmin xmax]); ylim([ymin ymax]);

cb = colorbar;
cb.Label.String = cblabel;
cb.FontSize=16;
colormap(CM);

pos = get(H,'Position');
set(H,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
fname = ['../Figures/ChangeInSpeed_',exp{jj},'_gs',gs{ii},'_m',m{ii},'_fixedC'];
print(H,fname,'-dpng','-r400');


%% bar charts with changes
H = fig('units','inches','width',80*12/72.27,'fontsize',16,'font','Helvetica');
subplot('position',[0.2 0.1 0.75 0.85]); hold on;

for jj=1:length(exp)
patch([sum(Flux_1996) sum(Flux_1996) GLFlux(jj) GLFlux(jj)],[jj-0.25 jj+0.25 jj+0.25 jj-0.25],'b');
text(sum(Flux_1996)+dGLFlux(jj)+1e12,jj,[num2str(dGLFlux(jj)/dFlux*100,2),'%']);
end

plot([sum(Flux_1996) sum(Flux_1996)],[0 length(exp)+1],'--k');
plot([sum(Flux_2015_2016) sum(Flux_2015_2016)],[0 length(exp)+1],'--k');
ylim([0 length(exp)+1]);
yticks(gca,[1:length(exp)]); set(gca,'yticklabel',description,'TickLabelInterpreter','latex');
xticks(gca,[sum(Flux_1996) sum(Flux_2015_2016)]); set(gca,'xticklabels',{'1996','2015-2016'},'TickLabelInterpreter','latex');
xlabel('Grounding line flux','Interpreter','latex');
box on;

title(['GL flux response (gs = ',gs{ii},', m = ',m{ii},', fixed C)'],'FontSize',16);

pos = get(H,'Position');
set(H,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
fname = ['../Figures/Change_Partitions_gs',gs{ii},'_m',m{ii},'_fixedC'];
print(H,fname,'-dpng','-r400');

%% combined bar chart
H = fig('units','inches','width',80*12/72.27,'fontsize',16,'font','Helvetica');
        subplot('position',[0.2 0.1 0.75 0.85]); hold on;
    
for jj=1:length(exp)
   patch([sum(Flux_1996) sum(Flux_1996) GLFlux(jj) GLFlux(jj)],[jj-0.25 jj+0.25 jj+0.25 jj-0.25],'b');
   text(sum(Flux_1996)+dGLFlux(jj)+1e12,jj,[num2str(dGLFlux(jj)/dFlux*100,2),'%']);
end

plot([sum(Flux_1996) sum(Flux_1996)],[0 length(exp)+1],'--k');
plot([sum(Flux_2015_2016) sum(Flux_2015_2016)],[0 length(exp)+1],'--k');
ylim([0 length(exp)+1]);
yticks(gca,[1:length(exp)]); set(gca,'yticklabel',description,'TickLabelInterpreter','latex');
xticks(gca,[sum(Flux_1996) sum(Flux_2015_2016)]); set(gca,'xticklabels',{'1996','2015-2016'},'TickLabelInterpreter','latex');
xlabel('Grounding line flux','Interpreter','latex');
box on;

title(['GL flux response (gs = ',gs{ii},', m = ',m{ii},', fixed C)'],'FontSize',16);

pos = get(H,'Position');
set(H,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
fname = ['../Figures/Change_Partitions_gs',gs{ii},'_m',m{ii},'_fixedC'];
    print(H,fname,'-dpng','-r400');     



%% CHANGES IN C

gs = {'10000','10000','10000'};
m = {'3','5','7'};
exp = {'RetreatChangesInC','ChangesInC','GLRetreat','GLRetreatCalving','IceShelfThinning','Calving'};
description = {'\begin{tabular}{c} Thinning and\\ changes in C\end{tabular}',...
    '\begin{tabular}{c} Changes in C\end{tabular}',...
    '\begin{tabular}{c} Thinning\\ (grounded and floating)\end{tabular}',...
    '\begin{tabular}{c} Thinning (grounded and floating) \\ and calving \end{tabular}',...
    '\begin{tabular}{c} Thinning \\ (ice shelf only)\end{tabular}',...
    '\begin{tabular}{c} Calving \end{tabular}'};
% exp = {'GLRetreat','IceShelfThinning','Calving'};
% description = {'\begin{tabular}{r} Thinning\\ (grounded and floating)\end{tabular}',...
%     '\begin{tabular}{r} Thinning \\ (ice shelf only)\end{tabular}',...
%     '\begin{tabular}{r} Calving\end{tabular}'};

CM=othercolor('RdYlBu11',202);
CM=flipdim(CM,1);
CM(end,:) = [1 0 0];
CM(1,:) = [0 0 1];
CM(101:102,:)=[1 1 1;1 1 1]*0.9;

for ii=1:length(gs)
    
    load(['../Inversion_MEaSUREs_2015_2016_Bedmachine_gs',gs{ii},'_m',m{ii},'_fixedAGlen_InverseRestartFile.mat'],'MUA','F','GF','CtrlVarInRestartFile');
%    load(['../Inversion_MEaSUREs_2015_2016_Bedmachine_gs',gs{ii},'_m3_InverseRestartFile.mat'],'MUA','F','GF','CtrlVarInRestartFile');
    CtrlVar = CtrlVarInRestartFile; 
    speed_2015_2016 = sqrt(F.ub.^2+F.vb.^2);
    Fspeed_2015_2016 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),speed_2015_2016);
    [Flux_2015_2016,~,~,~,~,~,~,~,GLgeo,~,~]=FluxAcrossGroundingLine(CtrlVar,MUA,GF,F.ub,F.vb,0*F.ub,0*F.vb,F.h,F.rho);
    fprintf('Grounding line flux 2015-2106, gs=%s: %s Gt/yr \n',gs{ii},num2str(sum(Flux_2015_2016(:))/1e9));
 
    load(['../NoPerturbationFinal_1996_gs',gs{ii},'_m',m{ii},'-RestartFile.mat'],'MUA','F','GF','CtrlVarInRestartFile');
    CtrlVar = CtrlVarInRestartFile; CtrlVar.PlotXYscale = 1e3;
    MUA_1996 = MUA; GF_1996=GF;
    speed_1996 = sqrt(F.ub.^2+F.vb.^2);
    Fspeed_1996 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),speed_1996);
    [Flux_1996,~,~,~,~,~,~,~,GLgeo,~,~]=FluxAcrossGroundingLine(CtrlVar,MUA,GF,F.ub,F.vb,0*F.ub,0*F.vb,F.h,F.rho);
    fprintf('Grounding line flux 1996, gs=%s: %s Gt/yr \n',gs{ii},num2str(sum(Flux_1996(:))/1e9));

    dFlux = sum(Flux_2015_2016) - sum(Flux_1996);
    
    for jj=1:length(exp)
        
        load(['../PerturbationFinal_1996-',exp{jj},'_gs',gs{ii},'_m',m{ii},'-RestartFile.mat'],'MUA','F','GF');
        speed = sqrt(F.ub.^2+F.vb.^2);
        dspeed_pert_1996 = speed - Fspeed_1996(MUA.coordinates(:,1),MUA.coordinates(:,2));
        dspeed_2015_1996 = Fspeed_2015_2016(MUA.coordinates(:,1),MUA.coordinates(:,2))-Fspeed_1996(MUA.coordinates(:,1),MUA.coordinates(:,2));
        percspeedup = 100*dspeed_pert_1996./dspeed_2015_1996;
        
        MaskInd = find(Fspeed_1996(MUA.coordinates(:,1),MUA.coordinates(:,2))<100 | Fspeed_2015_2016(MUA.coordinates(:,1),MUA.coordinates(:,2))<100);
        percspeedup(MaskInd) = 0;
        Mask = 0*percspeedup; Mask(MaskInd)=1;
        
        [Flux,~,~,~,~,~,~,~,~,~,~]=FluxAcrossGroundingLine(CtrlVar,MUA,GF,F.ub,F.vb,0*F.ub,0*F.vb,F.h,F.rho);
        GLFlux(jj) = sum(Flux);
        dGLFlux(jj) = sum(Flux)-sum(Flux_1996);
        fprintf('Grounding line flux %s, gs=%s: %s Gt/yr -- Explains %s%% of the change \n',...
           exp{jj},gs{ii},num2str(sum(Flux(:))/1e9),num2str(dGLFlux(jj)/dFlux*100));
    
        %% plotting 
        H = fig('units','inches','width',90*12/72.27,'fontsize',16,'font','Helvetica');
        
        
        switch plotoption
            case 'absolutechange'
                
                subplot('position',[0.1 0.1 0.42 0.8]); hold on;
                PlotNodalBasedQuantities_JDR(MUA.connectivity,MUA.coordinates/1e3,dspeed_pert_1996); hold on;
                caxis([-800 800]);
                xlim([xmin xmax]); ylim([ymin ymax]);
                grid on; box on; axis equal;

                xlabel('psx [km]','Fontsize',16); ylabel('psy [km]','Fontsize',16);
                plot(MUA0.Boundary.x/1e3,MUA0.Boundary.y/1e3,'--k');
                plot(MUA.Boundary.x/1e3,MUA.Boundary.y/1e3,'-k');
                PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'k');
        
                if strfind(exp{ii},'GLRetreat')
                    PlotGroundingLines(CtrlVar,MUA0,GF0,[],[],[],'--k');
                end

                subplot('position',[0.53 0.1 0.42 0.8]); hold on;
                PlotNodalBasedQuantities_JDR(MUA.connectivity,MUA.coordinates/1e3,dspeed_2015_1996); hold on;
                caxis([-800 800]);       
                yticklabels(gca,{}); ylabel(''); xlabel('psx [km]','Fontsize',16); 
                cblabel = 'Change in surface speed [m/yr]';
                
            case 'percentagespeedup'
                
                subplot('position',[0.1 0.1 0.85 0.8]); hold on;
                PlotNodalBasedQuantities_JDR(MUA.connectivity,MUA.coordinates/1e3,percspeedup); hold on;
                [C,h]=tricontour(MUA.connectivity,MUA.coordinates(:,1)/1e3,...
                    MUA.coordinates(:,2)/1e3,Mask,[0.99 0.99]);
                for ll=1:length(h)
                    h(ll).EdgeColor = [0.5 0.5 0.5];
                end
                  
                caxis([-101 101]); 
                ylabel('psy [km]','Fontsize',16); xlabel('psx [km]','Fontsize',16); 
                cblabel = '(U_{*}-U_{1996})/(U_{2016}-U_{1996}) [%]';
                
        end
        
        hold on; plot(MUA_1996.Boundary.x/1e3,MUA_1996.Boundary.y/1e3,'--k');
        plot(MUA.Boundary.x/1e3,MUA.Boundary.y/1e3,'-k');
        
        PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'k');
        
        if strfind(exp{jj},'Retreat')
            PlotGroundingLines(CtrlVar,MUA_1996,GF_1996,[],[],[],'--k');
        end
        
        htitle = sgtitle(['{\textbf{Response to $\left|',description{jj},'\right|$ (gs = ',gs{ii},', m = ',m{ii},')}}'],'Interpreter','latex');
        htitle.FontSize= 20;
        
        grid on; box on; axis equal;
        xlim([xmin xmax]); ylim([ymin ymax]);

        cb = colorbar;
        cb.Label.String = cblabel;
        cb.FontSize=16;
        colormap(CM);
      
        pos = get(H,'Position');
        set(H,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
        fname = ['../Figures/ChangeInSpeed_',exp{jj},'_gs',gs{ii},'_m',m{ii},'_fixedAGlen'];
        print(H,fname,'-dpng','-r400');

    end
    
    %% bar charts with changes
    H = fig('units','inches','width',80*12/72.27,'fontsize',16,'font','Helvetica');
        subplot('position',[0.2 0.1 0.75 0.85]); hold on;
    
    for jj=1:length(exp)
       patch([sum(Flux_1996) sum(Flux_1996) GLFlux(jj) GLFlux(jj)],[jj-0.25 jj+0.25 jj+0.25 jj-0.25],'b');
       text(sum(Flux_1996)+dGLFlux(jj)+1e12,jj,[num2str(dGLFlux(jj)/dFlux*100,2),'%']);
    end

    plot([sum(Flux_1996) sum(Flux_1996)],[0 length(exp)+1],'--k');
    plot([sum(Flux_2015_2016) sum(Flux_2015_2016)],[0 length(exp)+1],'--k');
    ylim([0 length(exp)+1]);
    yticks(gca,[1:length(exp)]); set(gca,'yticklabel',description,'TickLabelInterpreter','latex');
    xticks(gca,[sum(Flux_1996) sum(Flux_2015_2016)]); set(gca,'xticklabels',{'1996','2015-2016'},'TickLabelInterpreter','latex');
    xlabel('Grounding line flux','Interpreter','latex');
    box on;
    
    title(['GL flux response (gs = ',gs{ii},', m = ',m{ii},', fixed C)'],'FontSize',16);
    
    pos = get(H,'Position');
    set(H,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
    fname = ['../Figures/Change_Partitions_gs',gs{ii},'_m',m{ii},'_fixedAGlen'];
    print(H,fname,'-dpng','-r400');

end


