function ASensitivity

froot = '/Volumes/mainJDeRydt/Ua/cases/PIG_IceShelf_Integrity/final/';

file(1).name='Inversion_ERS_1996_REMAthickenedfinal_gs25000_Aprior-5C_m3_AGlen-Estimate.mat';
file(2).name='Inversion_ERS_1996_REMAthickenedfinal_gs25000_m3_AGlen-Estimate.mat';
file(3).name='Inversion_ERS_1996_REMAthickenedfinal_gs25000_Aprior-25C_m3_AGlen-Estimate.mat';

file(4).name='Inversion_ERS_1996_REMAthickenedfinal_gs25000_ga0.1_m3_AGlen-Estimate.mat';
file(5).name='Inversion_ERS_1996_REMAthickenedfinal_gs25000_ga1_m3_AGlen-Estimate.mat';
file(6).name='Inversion_ERS_1996_REMAthickenedfinal_gs25000_ga100_m3_AGlen-Estimate.mat';

titlestr = {'A_{prior}=-5C','A_{prior}=-15C','A_{prior}=-25C','(g_a)_A=0.1','(g_a)_A=1','(g_a)_A=100'};

load([froot,'Inversion_ERS_1996_REMAthickenedfinal_gs25000_m3_InverseRestartFile.mat']);
CtrlVarInRestartFile.PlotXYscale = 1e3;

%% Plotting
H=fig('units','inches','width',80*12/72.27,'fontsize',16,'font','Helvetica');
%set(H,'visible','off');

for ii=1:6
    
    load([froot,file(ii).name]);
    
    row = ceil(ii/3);
    column = mod(ii-1,3)+1;

    subplot('position',[0.05+(column-1)*0.28 0.1+(row-1)*0.44 0.28 0.37]); hold on;

    PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates/1e3,log10(AGlen));

    grid on; box on;
    
    if column==1
        ylabel('psy [km]'); 
        if row==2
            set(gca,'xticklabels','');
        end
    end
    if row==1
        xlabel('psx [km]');
        if column>1
            set(gca,'yticklabels','');
        end
    end
    if row>1 && column>1
        set(gca,'xticklabels','');
        set(gca,'yticklabels','');
    end
    
    colorbar('off');
    
    title(titlestr{ii});
    
    if ii==6
        cb = colorbar(gca,'Position',[0.9 0.3 0.015 0.5],'Location','east','AxisLocation','out',...
         'Ticks',[-11:-6],'ticklabels',{'-11','-10','-9','-8','-7','-6'});
        cb.XColor='k';
        cb.YColor='k';
        cb.TickLength=0.04;
        cb.FontSize=16;
        cb.Label.String = {'log_{10}(AGlen)'}; 
    end
    caxis([-11 -6]);
 
    PlotGroundingLines(CtrlVarInRestartFile,MUA,GF,[],[],[],'k','linewidth',1);

end

pos = get(H,'Position');
set(H,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
fname = './ASensitivity';
print(H,fname,'-dpng','-r400');
