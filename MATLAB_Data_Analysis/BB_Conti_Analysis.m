%% Analysis of Continuous Cultivations (Baby Bionet)
% Author: Manuel Bruch, PhD Student
% Created: 2021/04/07
% last edited: 2023/05/23

%% Tidy up
clear all
close all
clc

%% import data
FileLocation    =   '..\..\Laboratory_files\Conti_Cultures\ALE\';
FileNameCsv     =   {'Batch_BCU-M1_biomass_210513mbALE001_23.csv',...
                     'Batch_BCU-M1_biomass_210601mbALE002_24.csv',...
                     'Batch_BCU-M1_biomass_210802mbALE003_25.csv',...
                     'Batch_BCU-M1_biomass_mbALE004_26.csv'};
FileNameExcel   =   {'Batch_BCU-M1_biomass_210513mbALE001_23.xlsx',...
                     'Batch_BCU-M1_biomass_210601mbALE002_24.xlsx',...
                     'Batch_BCU-M1_biomass_210802mbALE003_25.xlsx',...
                     'Batch_BCU-M1_biomass_mbALE004_26.xlsx'};
FigureNames     =   {[FileLocation, 'ALE_figure_thesis']};

myYCol          =   [0.7 0 0];

% replaces ALL ',' with '.' -- save a copy of the original file first
% if condition prevents from re-replacing any decimal signs upon a second run
if ~isfile([FileLocation FileNameExcel{1}])
    comma2point_overwrite(FileLocation,FileNameCsv{1})    
    myData 	=   readtable([FileLocation FileNameCsv{1}], 'Delimiter', ';',...
                          'VariableNamingRule', 'preserve');
    writetable(myData,[FileLocation FileNameExcel{1}])
else
    myData 	=   readtable([FileLocation FileNameExcel{1}],...
                          'VariableNamingRule', 'preserve');
end

% append data from multiple files
if length(FileNameCsv) > 1
    for i = 2:length(FileNameCsv)
        if ~isfile([FileLocation FileNameExcel{i}])
            comma2point_overwrite(FileLocation,FileNameCsv{i})
            mytempData  =   readtable([FileLocation FileNameCsv{i}],...
                                      'Delimiter', ';', 'VariableNamingRule',...
                                      'preserve');
            writetable(mytempData,[FileLocation FileNameExcel{i}])
        else
            mytempData 	=   readtable([FileLocation FileNameExcel{i}],...
                                      'VariableNamingRule', 'preserve');
        end
        % get t_end from myData
        t_endTemp       =   myData{end,'Date'};
        t_startTemp     =   mytempData{1,'Date'};
        tdiff           =   hours(t_startTemp-t_endTemp);
        tadd            =   myData{end,'Age(Hours)'}+tdiff;
        mytempData{:,3} =   mytempData{:,3}+tadd;
        %%% same thing still has to be done for all the totalisers
        continuedVars   =   {'Acid_Total(ml)','Base_Total(ml)',...
                            'Bleed_Totalized(ml)','Media_Add_Totalized(ml)'};
        for j = 1:length(continuedVars)
            mytempData{:,continuedVars{j}}  = mytempData{:,...
                                                         continuedVars{j}}+...
                                                         myData{:,...
                                                         continuedVars{j}}(end);
        end
        myData  =   vertcat(myData,mytempData);
    end
end
t       =   table2array(myData(:,{'Age(Hours)'}))/24;

% load in offline CDW data
FileLocation    =   '..\..\Laboratory_files\C_necator_continuous.xlsx';
mySheet         =   'Results_Conti003_ALE';
tCDW            =   readmatrix(FileLocation,'Sheet',mySheet,'Range','A:B');
CDW             =   tCDW(:,2);
tCDW(:,2)       =   [];
tCDW            =   tCDW./24;

% load in data for pump rate changes
t0              =   readtable(FileLocation, 'Sheet', 'Data_Conti003_ALE',...
                              'Range', 'B9:B9');
ratechanges     =   readtable(FileLocation, 'Sheet',...
                              'Feed_rates_Conti003_ALE', 'Range', 'A2');
changetimes     =   datenum(ratechanges{:,1} - t0{:,:});
myrates         =   ratechanges{:,5};

%% Data Shenanigans
% careful, some things like the pump start are set manually via hard-coding
% --> make sure to change this for a new experiment
VarNames=   {'DO(%)', 'Acid_Total(ml)', 'Base_Total(ml)', 'bB_CO2(%)',...
             'bB_O2(%)', 'bB_RQ(CO2/O2)', 'Weight1(g)', 'bB_OUR(mmolO2/l•h)',...
             'bB_CER(mmolCO2/l•h)', 'pH'};
DO      =   table2array(myData(:,VarNames{1}));
acid    =   table2array(myData(:,VarNames{2}));
base    =   table2array(myData(:,VarNames{3}));
bBCO2   =   table2array(myData(:,VarNames{4}));
bBO2    =   table2array(myData(:,VarNames{5}));
RQ      =   table2array(myData(:,VarNames{6}));
OUR     =   table2array(myData(:,VarNames{8}));
CER     =   table2array(myData(:,VarNames{9}));
Weight  =   table2array(myData(:,VarNames{7}));
pH      =   table2array(myData(:,VarNames{10}));
% Weight(Weight>100)  =   NaN;
VarNames2   =   {'pH' 'Acid Total' 'Base Total' 'DO' 'CO2 In The Off-Gas'};
% VarNames=   strrep(VarNames,'Base_Total(ml)','Base Total [ml]');
% VarNames=   strrep(VarNames,'1(',' [');
% VarNames=   strrep(VarNames,'(',' [');
% VarNames=   strrep(VarNames,')',']');
% VarNames=   strrep(VarNames,'_',' ');
% VarNames=   strrep(VarNames,'O2','O_2');
% VarNames=   strrep(VarNames,'•','/');

onlineData  =   [pH acid base DO bBCO2];
myUnits     =   {' [%]',' [ml]',' [CO_2/O_2]',' [g]',' [mmolO_2/l/h]',...
                 ' [mmolCO_2/l/h]'};
% VarNames2   =   erase(VarNames,myUnits);
VarNames    =   {'pH','V(H_2SO_4) [ml]','V(KOH) [ml]','Dissolved O_2 [%]',...
                 'c(CO_2) [%]'};

%% plot these data
set(0, 'DefaultFigureRenderer', 'painters');
myFig           =   figure(1);
myLineStyles 	=   {'-x'};
myInterpreter   =   'tex';
myFontSize      =   [];
myTitle         =   [];
SqLength        =   ceil(sqrt(size(onlineData,2)+1));
tL              =   tiledlayout(SqLength,ceil((size(onlineData,2)+1)/SqLength));
ax              =   cell(size(onlineData,2)+1,1);
P               =   cell(size(onlineData,2),1);

% Biomass
ax{1}   =   nexttile;
yyaxis left
PCDW    =   myplot(tCDW, CDW, {'', 'CDW [g/l]'}, {''},...
                   [], myLineStyles, myInterpreter, [10 10 10 10], [], ax{1},...
                   1, [0 0 0]);
% xticks(0:10:max(tCDW)+1)
xlim([0,max(tCDW)+0.1])
% title('Cell Dry Weight')
yticks(0:0.1:0.5)
ylim([0, 0.5])

% alter y label position
HyLabel     =   get(gca,'ylabel');
HyLabel.VerticalAlignment   = "baseline";
yLabelPos   =   get(HyLabel,'position');
yLabelPos(1)=   yLabelPos(1) - 7;
set(HyLabel, 'position', yLabelPos)


% dilution rates
yyaxis right
set(gca, 'YColor', myYCol)
% myYCol  =   get(gca,'YColor');
stairs([changetimes;max(tCDW)], [myrates;myrates(end)], '-', 'Color', myYCol)
hold on
ylim([0, 0.162])
yticks(0:0.04:0.16)
xlim([0, 150])

% C/N ratio
if exist('CNratio', 'var')
    stairs([changetimesCN;max(tCDW)],[CNratio;CNratio(end)].*0.001,'--', 'Color', myYCol)
%     myYCol  =   get(gca,'YColor');
end

% c_FA
if exist('c_C', 'var')
    stairs([changetimesCN;max(tCDW)],[c_C;c_C(end)],'-.', 'Color', myYCol)
end
yyaxis left
set(gca,'YColor','k')


% online data
for i = 1:size(onlineData,2)
    ax{i+1}     =   nexttile;
    P           =   myplot(t, onlineData(:,i), {'', VarNames{i}},...
                           {''}, [], {'-'}, myInterpreter,...
                           [10 10 10 10], [], ax{i+1}, 1, [0 0 0]);
%     xticks(0:10:max(tCDW)+1)
    xlim([0,max(tCDW)+0.1])
%     title(VarNames2{i})
    if i == 1
        yticks(6:1:8)
        ylim([6,8])
    elseif i == 2
        yticks(0:5:20)
        ylim([0,20])
    elseif i == 3
        yticks(0:5:20)
        ylim([0,20])
    elseif i == 4
        yticks(0:20:100)
        ylim([0,100])
    elseif i == 5
        yticks(0:0.2:0.8)
        ylim([0,0.8])
    end

    % alter y label position
    HyLabel     =   get(gca,'ylabel');
    yLabelPos   =   get(HyLabel,'position');
    HyLabel.VerticalAlignment   = "baseline";
    yLabelPos(1)=   yLabelPos(1) - 7;
    set(HyLabel, 'position', yLabelPos)

    % dilution rates and carbon concentration
    yyaxis right
    stairs([changetimes;max(tCDW)],[myrates;myrates(end)],'-', 'Color', myYCol)
    hold on
    ylim([0, 0.162])
    yticks(0:0.04:0.16)
    xlim([0, 150])
    % C/N ratio
    if exist('CNratio', 'var')
        stairs([changetimesCN;max(tCDW)],[CNratio;CNratio(end)].*0.001,'--', 'Color', myYCol)
    end

    % c_FA
    if exist('c_C', 'var')
        stairs([changetimesCN;max(tCDW)],[c_C;c_C(end)],'-.', 'Color', myYCol)
    end
    set(gca,'YColor',myYCol)
    yyaxis left
    set(gca,'YColor','k')


end
xlabel(tL,'t [d]','FontName','Arial','FontSize',10,'Interpreter','tex')
ax              =   axes(myFig);
han             =   gca;
han.Visible     =   'off';
yyaxis(ax, 'right');
ylabel({'','D [h^{-1}]'}, 'FontName', 'Arial',...
    'FontSize', 10, 'Interpreter', 'tex', 'Color', myYCol)
han.YLabel.Visible = 'on';
set(myFig,'Position',[20 50 600 600])

%% save the figure(s)
savefig(myFig, [FigureNames{1} '.fig'])
exportgraphics(myFig, [FigureNames{1} '.emf']);
print(myFig, [FigureNames{1} '.svg'], '-dsvg')