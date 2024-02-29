function Figure2

froot = '/Volumes/Untitled/Ua/cases/PIG_IceShelf_Integrity_v2';
xmin = -1715; xmax = -1440;
ymin = -385; ymax = 0;

gs = {'50000','50000','50000'};
m = {'m3','m5','m7'};

exp = {'Calving','IceShelfThinning','GLRetreat','GLRetreatCalving',};
%description = {'\begin{tabular}{c} Calving \end{tabular}',...
%    '\begin{tabular}{c} Thinning \\ (ice shelf only)\end{tabular}',...
%    '\begin{tabular}{c} Thinning\\ (grounded and floating)\end{tabular}',...
%    '\begin{tabular}{c} Thinning (grounded and floating)\\ and calving \end{tabular}'};
titlestr = {'calving','ice shelf thinning','thinning (floating and grounded)','thinning and calving'};

addpath('/Volumes/Untitled/Ua/cases/PIG_IceShelf_Integrity_v2/Matlab_Tools');

CM=othercolor('RdYlBu11',220);
CM=flipdim(CM,1);
CM(end-9:end,:) = repmat([1 0 0],10,1);
CM(1:10,:) = repmat([0 0 1],10,1);
CM(110:111,:)=[1 1 1;1 1 1]*0.9;

if exist([froot,'/Matlab_Tools/FluxData.mat'])
    
    load([froot,'/Matlab_Tools/FluxData.mat']);
    
else
    
    for ii=1:length(m)

        data(ii).gs = gs{ii};
        data(ii).m = m{ii};

        load([froot,'/Inversion_MEaSUREs_2015_2016_Bedmachine_gs',gs{ii},'_m3_InverseRestartFile.mat'],'MUA','F','GF','CtrlVarInRestartFile');
        CtrlVar = CtrlVarInRestartFile; MUA_2015_2016 = MUA; GF_2015_2016=GF;
        CtrlVar.PlotXYscale = 1e3;
        speed_2015_2016 = sqrt(F.ub.^2+F.vb.^2);
        Fspeed_2015_2016 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),speed_2015_2016);
        GLgeo_2015_2016=GLgeometry(MUA.connectivity,MUA.coordinates,GF,CtrlVar);
        [Flux_2015_2016]=FluxAcrossGroundingLine_centralPIG(CtrlVar,MUA,F.ub,F.vb,0*F.ub,0*F.vb,F.h,F.rho,GLgeo_2015_2016);
        data(ii).Flux_2015_2016 = sum(Flux_2015_2016(:))/1e9;
        fprintf('Grounding line flux 2015-2106, gs=%s: %s Gt/yr \n',gs{ii},num2str(sum(Flux_2015_2016(:))/1e9));

        load([froot,'/NoPerturbationFinal_1996_gs',gs{ii},'_',m{ii},'-RestartFile.mat'],'MUA','F','GF','CtrlVarInRestartFile');
        CtrlVar = CtrlVarInRestartFile; CtrlVar.PlotXYscale = 1e3;
        MUA_1996 = MUA; GF_1996=GF;
        speed_1996 = sqrt(F.ub.^2+F.vb.^2);
        Fspeed_1996 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),speed_1996);
        [Flux_1996]=FluxAcrossGroundingLine_centralPIG(CtrlVar,MUA,F.ub,F.vb,0*F.ub,0*F.vb,F.h,F.rho,GLgeo_2015_2016);
        data(ii).Flux_1996 = sum(Flux_1996(:))/1e9;
        fprintf('Grounding line flux 1996, gs=%s: %s Gt/yr \n',gs{ii},num2str(sum(Flux_1996(:))/1e9));

        data(ii).dFlux = sum(Flux_2015_2016) - sum(Flux_1996);

        for jj=1:length(exp)

            data(ii).Perturbation(jj).exp = exp{jj};

            load([froot,'/PerturbationFinal_1996-',exp{jj},'_gs',gs{ii},'_',m{ii},'-RestartFile.mat'],'MUA','F','GF');
            data(ii).Perturbation(jj).MUA = MUA;
            
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
            data(ii).Perturbation(jj).dGLFlux = sum(Flux)-sum(Flux_1996);
            fprintf('Grounding line flux %s, gs=%s: %s Gt/yr -- Explains %s%% of the change \n',...
               exp{jj},gs{ii},num2str(sum(Flux(:))/1e9),num2str(data(ii).Perturbation(jj).dGLFlux/data(ii).dFlux*100));
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

H = fig('units','inches','width',90*12/72.27,'fontsize',16,'font','Helvetica');

%% calving, ice shelf thinning, thinning
ii=1;
for jj=1:3
    
    subplot('position',[0.1+0.3*(jj-1) 0.52 0.28 0.4]); hold on;

    PlotNodalBasedQuantities_JDR(data(ii).Perturbation(jj).MUA.connectivity,...
        data(ii).Perturbation(jj).MUA.coordinates/1e3,data(ii).Perturbation(jj).percspeedup); hold on;
%     [C,h]=tricontour(data(ii).Perturbation(jj).MUA.connectivity,data(ii).Perturbation(jj).MUA.coordinates(:,1)/1e3,...
%         data(ii).Perturbation(jj).MUA.coordinates(:,2)/1e3,data(ii).Perturbation(jj).Mask,[0.99 0.89]);
%     for ll=1:length(h)
%         h(ll).EdgeColor = [0.5 0.5 0.5];
%     end

    caxis([-110 110]);
    colormap(CM); 
    
    if strfind(data(ii).Perturbation(jj).exp,'Calving')
        plot(MUA_1996.Boundary.x/1e3,MUA_1996.Boundary.y/1e3,'--k');
    end
    plot(data(ii).Perturbation(jj).MUA.Boundary.x/1e3,data(ii).Perturbation(jj).MUA.Boundary.y/1e3,'-k');
    
    plot(xGL_central/1e3,yGL_central/1e3,'-m','linewidth',2);
    if strfind(data(ii).Perturbation(jj).exp,'Retreat')
        PlotGroundingLines(CtrlVar,MUA_1996,GF_1996,[],[],[],'--k','linewidth',1);
        PlotGroundingLines(CtrlVar,MUA_2015_2016,GF_2015_2016,[],[],[],'-k','linewidth',1);  
    else
        PlotGroundingLines(CtrlVar,MUA_1996,GF_1996,[],[],[],'-k','linewidth',1);
    end

    if jj==1
        ylabel('psy [km]','Fontsize',16);  
    else
        xlabel('psx [km]','Fontsize',16); 
        ylabel('');
        yticklabels({});
    end
    
    if ii==3
        cb=colorbar(gca,'Position',[0.92 0.58 0.015 0.35],'Location','west','AxisLocation','in',...
         'Ticks',[-100:25:100],'TickLabels',{'-100','','-50','','0','','50','','100'});
        cb.XColor='k';
        cb.YColor='k';
        cb.TickLength=0.04;
        cb.FontSize=16;
        cb.Label.String = {'100(U_{2016}-U_{1996})/U_{2016} [%]'};
        cbarrow;
    end
    
    grid on; box on; axis equal;
    xlim([xmin xmax]); ylim([ymin ymax]);

    title(titlestr{jj});

end

%% calving and thinning
subplot('position',[0.1 0.1 0.42 0.4]); hold on;

ii = 1; jj = 4;
PlotNodalBasedQuantities_JDR(data(ii).Perturbation(jj).MUA.connectivity,...
        data(ii).Perturbation(jj).MUA.coordinates/1e3,data(ii).Perturbation(jj).percspeedup); hold on;

caxis([-110 110]);
colormap(CM); 

if strfind(data(ii).Perturbation(jj).exp,'Calving')
    plot(MUA_1996.Boundary.x/1e3,MUA_1996.Boundary.y/1e3,'--k');
end
plot(data(ii).Perturbation(jj).MUA.Boundary.x/1e3,data(ii).Perturbation(jj).MUA.Boundary.y/1e3,'-k');

plot(xGL_central/1e3,yGL_central/1e3,'-m','linewidth',2);
if strfind(data(ii).Perturbation(jj).exp,'Retreat')
    PlotGroundingLines(CtrlVar,MUA_1996,GF_1996,[],[],[],'--k','linewidth',1);
    PlotGroundingLines(CtrlVar,MUA_2015_2016,GF_2015_2016,[],[],[],'-k','linewidth',1);  
else
    PlotGroundingLines(CtrlVar,MUA_1996,GF_1996,[],[],[],'-k','linewidth',1);
end

ylabel('psy [km]','Fontsize',16);
xlabel('psx [km]','Fontsize',16); 

grid on; box on; axis equal;
xlim([xmin xmax]); ylim([ymin ymax]);

title(titlestr{jj});

%% bar chart GL flux
subplot('position',[0.56 0.1 0.42 0.4]); hold on;

ii = 1;
for jj=1:4
   patch([data(ii).Flux_1996 data(ii).Flux_1996 data(ii).Perturbation(jj).GLFlux data(ii).Perturbation(jj).GLFlux],...
       [(5-jj)-0.25 (5-jj)+0.25 (5-jj)+0.25 (5-jj)-0.25],'b');
   text(data(ii).Flux_1996+data(ii).Perturbation(jj).dGLFlux+1e12,jj,...
       [num2str(data(ii).Perturbation(jj).dGLFlux/data(ii).dFlux*100,2),'%']);
end

plot([data(ii).Flux_1996 data(ii).Flux_1996],[0 length(exp)+1],'--k');
plot([data(ii).Flux_2015_2016 data(ii).Flux_2015_2016],[0 length(exp)+1],'--k');
ylim([0 length(exp)+1]);
yticks(gca,[1:length(exp)]); set(gca,'yticklabel',flipdim(titlestr,2),'TickLabelInterpreter','latex');
xticks(gca,[data(ii).Flux_1996 data(ii).Flux_2015_2016]); set(gca,'xticklabels',{'1996','2016'},'TickLabelInterpreter','latex');
xlabel('Grounding line flux','Interpreter','latex');
box on;
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
                cblabel = '(u_{perturbation}-u_{1996})/(u_{2016}-u_{1996}) [%]';
                
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


