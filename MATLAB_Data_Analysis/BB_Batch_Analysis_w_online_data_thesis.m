%% Analysis of Continuous Cultivations (Baby Bionet)
% Author: Manuel Bruch, PhD Student
% Created: 2021/04/07
% last edited: 2023/05/23

%% Tidy up
clear all
close all
clc

%% script parameters
FileLocationOnline  =   '..\..\Laboratory_files\Batch_Cultivations\';
FileNamesCsv        =   {'Batch_BCU-M1_biomass_Manuel020321Batch001_8.csv',...
                         'Batch_BCU-M1_biomass_Manuel021221Batch002_9.csv',...
                         'Batch_BCU-M3_biomass_Manuel021721Batch003_33.csv'};
FileNamesExcel      =   {'Batch_BCU-M1_biomass_Manuel020321Batch001_8.xlsx',...
                         'Batch_BCU-M1_biomass_Manuel021221Batch002_9.xlsx',...
                         'Batch_BCU-M3_biomass_Manuel021721Batch003_33.xlsx'};

FileLocationOffline =   ['..\..\Laboratory_files\growth_experiments_2021\',...
                        'C_necator_bioreactor_cultivation.xlsx'];
OfflineSheets       =   {'batch_triplicate_biomass',...
                         'Data_01-05_02_2021',...
                         'Data_08-13_02_2021',...
                         'Data_15-20_02_2021'};
OfflineRanges       =   {{'A2:B11', 'A21:B30'}, 'A11:B20', 'A11:B20', 'A11:B20'};
FigureLocation      =   '..\..\Laboratory_files\growth_experiments_2021\batch_fermentations\';
FigureNames         =   {'2021_WT_online_batch_data_triplicate',...
                         '2021_WT_online_batch_data_batch01',...
                         '2021_WT_online_batch_data_batch02',...
                         '2021_WT_online_batch_data_batch03',...
                         };

saveFigures         =   1;


%% load in offline data
FA_CDWconversionFactor  =   load('..\..\Laboratory_files\FA_CDWconversionFactor.mat');
FA_CDW_Factor           =   FA_CDWconversionFactor.factor;

myOfflineData   =   cell(size(OfflineSheets));
for i = 1:length(OfflineSheets)
    if iscell(OfflineRanges{i})     % load in triplicate data averages and stddevs
        avgData     =   readtable(FileLocationOffline,...
                                  'Sheet', OfflineSheets{i},...
                                  'Range', OfflineRanges{i}{1},...
                                  'ReadVariableNames', true);
        stdData     =   readtable(FileLocationOffline,...
                                  'Sheet', OfflineSheets{i},...
                                  'Range', OfflineRanges{i}{2},...
                                  'ReadVariableNames', true);
        biomassData =   cat(3, avgData{:,:}, stdData{:,:});
    else
        biomassData =   readtable(FileLocationOffline,...
                                  'Sheet', OfflineSheets{i},...
                                  'Range', OfflineRanges{i},...
                                  'ReadVariableNames', true);

        biomassData         =   biomassData{:,:};
        biomassData(:,2)    =   biomassData(:,2) * FA_CDW_Factor;
    end
    myOfflineData{i}    =   biomassData;
end


%% load online data
% replaces ALL ',' with '.' -- save a copy of the original file first
% the if condition prevents from re-replacing any decimal signs upon a second run
myData  =   cell(size(FileNamesCsv));
tOnline =   cell(size(FileNamesCsv));
for i = 1:length(FileNamesCsv)
    if ~isfile([FileLocationOnline FileNamesExcel{i}])
        comma2point_overwrite(FileLocationOnline,FileNamesCsv{i})
        myData 	    =   readtable([FileLocationOnline FileNamesExcel{i}],...
                                   'Delimiter', ';', 'VariableNamingRule',...
                                   'preserve');
        writetable(myData,[FileLocationOnline FileNamesExcel{i}])
    else
        myData{i}   =   readtable([FileLocationOnline FileNamesExcel{i}],...
                                  'VariableNamingRule','preserve');
    end
    tOnline{i}  =   table2array(myData{i}(:,{'Age(Hours)'}));
end


%% Data Shenanigans
% careful, some things like the pump start are set manually via hard-coding
% --> make sure to change this for a new experiment
VarNames=   {'DO(%)' 'Acid_Total(ml)' 'Base_Total(ml)' 'bB_CO2(%)' 'bB_O2(%)' 'bB_RQ(CO2/O2)' 'Weight1(g)' 'bB_OUR(mmolO2/l•h)' 'bB_CER(mmolCO2/l•h)' 'pH'};
DO          =   cell(1, size(FileNamesCsv, 2));     % storage cell for individual replicates
acid        =   cell(1, size(FileNamesCsv, 2));
base        =   cell(1, size(FileNamesCsv, 2));
bBCO2       =   cell(1, size(FileNamesCsv, 2));
bBO2        =   cell(1, size(FileNamesCsv, 2));
RQ          =   cell(1, size(FileNamesCsv, 2));
OUR         =   cell(1, size(FileNamesCsv, 2));
CER         =   cell(1, size(FileNamesCsv, 2));
Weight      =   cell(1, size(FileNamesCsv, 2));
pH          =   cell(1, size(FileNamesCsv, 2));
onlineData  =   cell(1, size(FileNamesCsv, 2));

for i = 1:length(FileNamesCsv)
    DO{i}           =   table2array(myData{i}(:,VarNames{1}));
    acid{i}         =   table2array(myData{i}(:,VarNames{2}));
    base{i}         =   table2array(myData{i}(:,VarNames{3}));%/1000; % converted from ml to L
    bBCO2{i}        =   table2array(myData{i}(:,VarNames{4}));
    bBO2{i}         =   table2array(myData{i}(:,VarNames{5}));
    RQ{i}           =   table2array(myData{i}(:,VarNames{6}));
    OUR{i}          =   table2array(myData{i}(:,VarNames{8}));
    CER{i}          =   table2array(myData{i}(:,VarNames{9}));
    Weight{i}       =   table2array(myData{i}(:,VarNames{7}));
    pH{i}           =   table2array(myData{i}(:,VarNames{10}));
    onlineData{i}   =   [pH{i} acid{i} base{i} DO{i} bBCO2{i}];% RQ{i} OUR{i} CER{i}];
end

VarNames2   =   {'pH' 'Acid Total' 'Base Total' 'DO' 'CO_2 In The Off-Gas'};

myUnits     =   {' [%]',' [ml]',' [CO_2/O_2]',' [g]',' [mmolO_2/l/h]',' [mmolCO_2/l/h]'};
VarNames    =   {'pH','V(H_2SO_4) [ml]','V(NH_4OH) [ml]','Dissolved O_2 [%]','c(CO_2) [%]'};


%% plot these data
for figCount = 1:length(OfflineSheets)
    figure('Renderer', 'Painters')
    myFig           =   figure(figCount);
    % myAxisNames     =   {'t [h]','CDW [g/l]'};
    myLineStyles 	=   {'-x'};
    myInterpreter   =   'tex';
    myFontSize      =   [];
    myTitle         =   [];
    if figCount == 1
        SqLength        =   ceil(sqrt(size(onlineData{1},2)+1));
        tL              =   tiledlayout(SqLength,ceil((size(onlineData{1},2)+1)/SqLength));
        ax              =   cell(size(onlineData{1},2)+1,1);
        P               =   cell(size(onlineData{1},2),1);
    else
        SqLength        =   ceil(sqrt(size(onlineData{figCount-1},2)+1));
        tL              =   tiledlayout(SqLength,ceil((size(onlineData{figCount-1},2)+1)/SqLength));
        ax              =   cell(size(onlineData{figCount-1},2)+1,1);
        P               =   cell(size(onlineData{figCount-1},2),1);
    end
    

    % Biomass
    ax{1}   =   nexttile;
    PCDW    =   myplot(myOfflineData{figCount}(:, 1, :),...
                       myOfflineData{figCount}(:, 2, :),...
                       {'',{'CDW [g/l]',''}}, {'Cell Dry Weight'}, [],...
                       myLineStyles, myInterpreter, [10 10 10 12], [],...
                       ax{1}, 1, [0 0 0]);
    % xticks(0:10:max(tCDW)+1)
    xlim([0, myOfflineData{figCount}(end, 1, 1) + 0.1])
    title('Cell Dry Weight')
    yticks(0:0.2:0.6)
    ylim([0, 0.6])

    if figCount == 1
        % truncate online data arrays to have same lengths
        tOnlineTrunc    =   cell(size(tOnline));
        onlineDataTrunc =   cell(size(tOnline));
        for lenCount = 1:length(tOnline)
            if lenCount == 1
                minLen  =   size(tOnline{lenCount}, 1);
            else
                minLen  =   min(minLen, size(tOnline{lenCount}, 1));
            end
        end
        for lenCount = 1:length(tOnline)
            tOnlineTrunc{lenCount}      =   tOnline{lenCount}(1:minLen, :);
            tOnlineTrunc{lenCount}      =   tOnline{lenCount}(1:minLen, :);
            onlineDataTrunc{lenCount}   =   onlineData{lenCount}(1:minLen, :);
        end
        xAvg    =   mean(cat(3, tOnlineTrunc{:}), 3);
        xStd    =   std(cat(3, tOnlineTrunc{:}), [], 3);
        x       =   cat(3, xAvg, xStd);
        yAvg    =   mean(cat(3, onlineDataTrunc{:}), 3);
        yStd    =   std(cat(3, onlineDataTrunc{:}), [], 3);
        y       =   cat(3, yAvg, yStd);
    else
        x       =   tOnline{figCount - 1};
        y       =   onlineData{figCount - 1};
    end

    % online data
    for i = 1:size(y,2)
        ax{i+1}     =   nexttile;
        P           =   myplot(x, y(:,i,:), {'',{VarNames{i},''}},...
                               VarNames(i), [], {'-'}, myInterpreter,...
                               [10 10 10 12], [], ax{i+1}, 1, [0 0 0]);
        %     xticks(0:10:max(tCDW)+1)
        xlim([0, myOfflineData{figCount}(end, 1, 1)+0.1])
        title(VarNames2{i}, 'FontName', 'Arial', 'FontSize', 12,...
            'Interpreter', 'tex')
        %order of online data: 'pH' 'Acid Total' 'Base Total' 'DO' 'CO_2 In The Off-Gas'
        if i == 1
            yticks(6:1:8)
            ylim([6,8])
        elseif i == 2
            yticks(0:20:60)
            ylim([0,60])
        elseif i == 3
            yticks(0:0.2:1)
            ylim([0,1])
        elseif i == 4
            yticks(0:20:100)
            ylim([0,100])
        elseif i == 5
            yticks(0:0.1:0.3)
            ylim([0,0.3])
        end


    end
    xlabel(tL, 't [d]', 'FontName', 'Arial', 'FontSize', 10, 'Interpreter',...
        'tex')
    set(myFig,'Position', [20 50 600 600])

    % save the figure(s)
    if saveFigures
        figSavePath     =   [FigureLocation FigureNames{figCount}];
        savefig(myFig,[figSavePath '.fig'])
        exportgraphics(myFig,[figSavePath '.emf'],'BackgroundColor','w');
        print(myFig, [figSavePath '.svg'], '-dsvg')
    end
end