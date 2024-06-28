%% Analysis of Continuous Cultivations (Baby Bionet)
% Author: Manuel Bruch, PhD Student
% Created: 2021/04/07
% last edited: 2021/04/07

%% Tidy up
clear all
close all
clc

%% import data
FileLocation    =   '..\..\Laboratory_files\Conti_Cultures\';
FileNameCsv     =   {'Batch_BCU-M1_biomass_210325mbconti001_13.csv','Batch_BCU-M3_biomass_210413mbconti002_50.csv'};%,'ALE\Batch_BCU-M1_biomass_210513mbALE001_23.csv'};
FileNameExcel   =   {'Batch_BCU-M1_biomass_210325mbconti001_13.xlsx','Batch_BCU-M3_biomass_210413mbconti002_50.xlsx'};%,'ALE\Batch_BCU-M1_biomass_210513mbALE001_23.xlsx'};
FigureNames     =   'conti001_conti002_comparison';
FigSavePath     =   ['..\..\Laboratory_files\growth_experiments_2021\' FigureNames];
HPLCfolder      =   '..\..\Laboratory_files\HPLC_Files\20210604_Agilent_I\Results_Cmol_l\';
HPLCfileName    =   'Reactor00';

% replaces ALL ',' with '.' -- save a copy of the original file first
% the if condition prevents from re-replacing any decimal signs upon a second run
myData  =   cell(1,length(FileNameCsv));
for i = 1:length(FileNameCsv)
    if ~isfile([FileLocation FileNameExcel{i}])
        comma2point_overwrite(FileLocation,FileNameCsv{i})
        myData{i}   =   readtable([FileLocation FileNameCsv{i}],'Delimiter',';','VariableNamingRule','preserve');
        writetable(myData{i},[FileLocation FileNameExcel{i}])
    else
        myData{i} 	=   readtable([FileLocation FileNameExcel{i}],'VariableNamingRule','preserve');
    end
end

% load in offline CDW data
FileLocation    =   '..\..\Laboratory_files\C_necator_continuous.xlsx';
mySheet         =   {'Results_Conti001','Results_Conti002'};%,'Results_Conti003_ALE'};
tCDW    =   cell(1,length(FileNameCsv));
CDW     =   tCDW;
compData=   tCDW;
for i = 1:length(FileNameCsv)
    compData{i} =   readtable(FileLocation,'Sheet',mySheet{i},'Range','A:B');
    CDW{i}      =   compData{i}{:,2};
    tCDW{i}     =   compData{i}{:,1};
end

%% Data Shenanigans
% careful, some things like the pump start are set manually via hard-coding
% --> make sure to change this for a new experiment
onlineData  =   cell(1,length(FileNameCsv));
t           =   cell(1,length(FileNameCsv));
VarNames=   {'DO(%)' 'Acid_Total(ml)' 'Base_Total(ml)' 'bB_CO2(%)'....
            'bB_O2(%)' 'bB_RQ(CO2/O2)' 'bB_OUR(mmolO2/l•h)'....
            'bB_CER(mmolCO2/l•h)'};

for i = 1:length(FileNameCsv)
    t{i}    =   table2array(myData{i}(:,{'Age(Hours)'}));
    DO      =   table2array(myData{i}(:,VarNames{1}));
    acid    =   table2array(myData{i}(:,VarNames{2}));
    base    =   table2array(myData{i}(:,VarNames{3}));
    bBCO2   =   table2array(myData{i}(:,VarNames{4}));
    bBO2    =   table2array(myData{i}(:,VarNames{5}));
    RQ      =   table2array(myData{i}(:,VarNames{6}));
    OUR     =   table2array(myData{i}(:,VarNames{7}));
    CER     =   table2array(myData{i}(:,VarNames{8}));
    onlineData{i}   =   [DO acid base bBCO2];% bBO2 RQ OUR CER];
end
VarNames=   {'DO(%)' 'Acid_Total(ml)' 'Base_Total(ml)' 'bB_CO2(%)'};
VarNames=   strrep(VarNames,'(',' [');
VarNames=   strrep(VarNames,')',']');
VarNames=   strrep(VarNames,'_',' ');
VarNames=   strrep(VarNames,'O2','O_2');
VarNames=   strrep(VarNames,'•','/');
VarNames=   strrep(VarNames,'Acid Total','V(H_2SO_4)');
VarNames=   strrep(VarNames,'Base Total','V(KOH)');


myUnits     =   {' [%]',' [ml]',' [CO_2/O_2]',' [mmolO_2/l/h]',' [mmolCO_2/l/h]'};
VarNames2   =   erase(VarNames,myUnits);

HPLCData    =   cell(2,length(FileNameCsv));

%% add HPLC data to the mix
for i = 1:length(FileNameCsv)
    HPLCDataTemp    =   readtable([HPLCfolder HPLCfileName num2str(i) '.txt'],...
                        'ReadRowNames',true,'Delimiter','\t');
    
    TimePoints  =   HPLCDataTemp.Properties.RowNames;
    xHPLC  =   zeros(size(TimePoints));
    for kk = 1:length(xHPLC)
        xHPLC(kk)   =   str2double(TimePoints{kk}(5:end));
    end

    yHPLC           =   table2array(HPLCDataTemp);
    HPLCData{1,i}   =   HPLCDataTemp.Properties.VariableNames;
    HPLCData{2,i}   =   [xHPLC yHPLC];
end
%% plot these data
myFig           =   figure(1);
% myAxisNames     =   {'t [h]','CDW [g/l]'};
myLineStyles 	=   {'-x','--o',':^'};
myInterpreter   =   'tex';
myFontSize      =   [];
myTitle         =   [];
if exist('HPLCData','var')
    SqLength        =   ceil(sqrt(size(onlineData{1},2)+1+1));
else
    SqLength        =   ceil(sqrt(size(onlineData{1},2)+1));
end
tL              =   tiledlayout(SqLength,ceil((size(onlineData{1},2)+1)/SqLength));
ax              =   cell(size(onlineData{1},2)+1,1);
P               =   cell(size(onlineData{1},2),1);

ax{1}   =   nexttile;
PCDW    =   myplot(tCDW{1},CDW{1},{'',{'CDW [g/l]',''}},{'CDW'},[],myLineStyles,...
            myInterpreter,[10 10 10 12],[],ax{1},1,[0 0 0;0 0.5 0;0.5 0 0]);
PCDW    =   myplot(tCDW{2},CDW{2},{'',{'CDW [g/l]',''}},{'CDW'},[],myLineStyles,...
            myInterpreter,[10 10 10 12],[],ax{1},1,[0 0 0;0 0.5 0;0.5 0 0]);
xticks(0:96:max(max(tCDW{1}, max(tCDW{2}))+12))
xlim([0,max(max(tCDW{1}, max(tCDW{2}))+12)])
xline(12,'-r')
text(12,0.07,'\leftarrow pump activation','FontName','Arial','FontSize',10,'Color','r')

yticks(0:0.1:0.5)
ylim([0,0.5])
for i = 1:size(onlineData{1},2)
    ax{i+1}     =   nexttile;
    P           =   myplot(t{1},onlineData{1}(:,i),{'',{VarNames{i},''}},VarNames(i),...
                    [],{'-'},myInterpreter,[10 10 10 12],[],ax{i+1},1,[0 0 0;0 0.5 0;0.5 0 0]);
    P           =   myplot(t{2},onlineData{2}(:,i),{'',{VarNames{i},''}},VarNames(i),...
                    [],{'-'},myInterpreter,[10 10 10 12],[],ax{i+1},1,[0 0 0;0 0.5 0;0.5 0 0]);
    xticks(0:96:max(max(tCDW{1}, max(tCDW{2}))+12))
    xlim([0,max(max(tCDW{1}, max(tCDW{2}))+12)])
    xline(12,'-r')
    if i == 1
        yticks(0:20:100)
        ylim([0,100])
    elseif i == 2 || i == 3
        yticks(0:3:15)
        ylim([0,15])
    elseif i == 4
        yticks(0:0.1:0.3)
        ylim([0,0.3])
    elseif i == 5
        yticks(0:6:24)
        ylim([0,24])
    end
end

% plot HPLC data
ax{i+1}     =   nexttile;
P           =   myplot(HPLCData{2,1}(:,1),HPLCData{2,1}(:,2)*1000,{'',{'c_{FA} [C-mmol/l]',''}},{'c_{FA} [C-mmol/l]'},...
                [],{'-'},myInterpreter,[10 10 10 12],[],ax{i+1},1,[0 0 0;0 0.5 0;0.5 0 0]);
P           =   myplot(HPLCData{2,2}(:,1),HPLCData{2,2}(:,2)*1000,{'',{'c_{carbon} [C-mmol/l]',''}},{'c_{FA} [C-mmol/l]'},...
                [],{'-'},myInterpreter,[10 10 10 12],[],ax{i+1},1,[0 0 0;0 0.5 0;0.5 0 0]);
xticks(0:96:max(max(tCDW{1}, max(tCDW{2}))+12))
xlim([0,max(max(tCDW{1}, max(tCDW{2}))+12)])
xline(12,'-r')
yticks(0:15:90)
ylim([0,90])
    
xlabel(tL,'t [h]','FontName','Arial','FontSize',10)
lgd     =   legend({'Run 1','Run 2'}, 'FontName', 'Arial',...
                   'FontSize', 10, 'Interpreter', 'tex','NumColumns', 2);
lgd.Layout.Tile     =   'north';
set(myFig,'Position',[20 50 600 600])


%% save the figure(s)
savefig(myFig, [FigSavePath '.fig'])
exportgraphics(myFig, [FigSavePath '.emf']);
print(myFig, [FigSavePath '.svg'], '-dsvg')