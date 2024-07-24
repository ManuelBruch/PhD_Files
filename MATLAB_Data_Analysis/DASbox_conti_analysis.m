%% Analysis of Continuous Cultivations (Dasbox)
% Author: Manuel Bruch, PhD Student
% Created: 2023/12/06
% last edited: 2023/12/06

%% Tidy up
clear all
close all
clc

%% import data
% online data
runToAnalyse        =   'conti005';
curDate             =   '240722';
offlineFileLocation =   '..\..\..\Laboratory_files\Conti_Cultures\DASbox_fermentations\';
offlineFileName     =   'DASbox_conti_setup_and_data.xlsx';
switch runToAnalyse
    case 'conti003'
        FileLocation        =   '..\..\..\Laboratory_files\Conti_Cultures\DASbox_fermentations\ALE26_conti003\online_data\';
        FileNameCsv         =   'CTPCMH008571.Administrator ccd7d35f.Control.csv';
        onlineDataSaveFile  =   '20231213_online_Data.mat';
        % offline data
        sheetName           =   'Results_ALE26_conti003';
        t0Sheet             =   {'231116_DASbox_conti03_data', 'B4:B4'};
        pumpRateSheet       =   {'231116_conti_feed_rates', 'A2'};
        CNratioSheet        =   {'231116_N_concentrations', 'A7:I13'};
        fullFile            =   [offlineFileLocation, offlineFileName];
        % plotting data
        FigurePath      =   '..\..\..\Laboratory_files\Conti_Cultures\DASbox_fermentations\ALE26_conti003\figures\';
        FigureName      =   [curDate '_ALE26_conti003_DASbox_onlineData_thesis'];

    case 'conti004'
        FileLocation        =   '..\..\..\Laboratory_files\Conti_Cultures\DASbox_fermentations\ALE26_conti004\online_data\';
        FileNameCsv         =   'CTPCMH008571.CnecatorConti004.Control.csv';
        onlineDataSaveFile  =   '20240304_online_Data.mat';
        % offline data
        sheetName           =   'Results_ALE26_conti004';
        t0Sheet             =   {'240115_DASbox_conti04_data', 'B4:B4'};
        pumpRateSheet       =   {'240115_conti_feed_rates', 'A2'};
        CNratioSheet        =   {'240115_N_concentrations', 'A7:I15'};
        fullFile            =   [offlineFileLocation, offlineFileName];
        % plotting data
        FigurePath      =   '..\..\..\Laboratory_files\Conti_Cultures\DASbox_fermentations\ALE26_conti004\figures\';
        FigureName      =   [curDate '_ALE26_conti004_DASbox_onlineData_thesis'];

    case 'conti005'
        FileLocation        =   '..\..\..\Laboratory_files\Conti_Cultures\DASbox_fermentations\ALE26_conti005\online_data\';
        FileNameCsv         =   'CTPCMH008571.Administrator a9870b62.Control.csv';
        onlineDataSaveFile  =   '20240503_online_Data.mat';
        % offline data
        sheetName           =   'Results_ALE26_conti005';
        t0Sheet             =   {'240402_DASbox_conti05_data', 'B4:B4'};
        pumpRateSheet       =   {'240402_conti_feed_rates', 'A2'};
        CNratioSheet        =   {'240402_N_concentrations', 'A7:I11'};
        fullFile            =   [offlineFileLocation, offlineFileName];
        % plotting data
        FigurePath      =   '..\..\..\Laboratory_files\Conti_Cultures\DASbox_fermentations\ALE26_conti005\figures\';
        FigureName      =   [curDate '_ALE26_conti005_DASbox_onlineData_thesis'];
end

saveFigure      =   0;

if ~isfile([FileLocation onlineDataSaveFile])
    onlineData      =   readmatrix([FileLocation, FileNameCsv], 'Delimiter', ';',...
        'OutputType', 'char'); %, 'NumHeaderLines', 295);

    %% format the table
    % find the starting position of the actual table
    for i = 1:size(onlineData,1)
        try
            if strcmp(onlineData{i,2}, 'Duration') && strcmp(onlineData{i+1,2}, '0')
                tableHeaders = onlineData(i, :);
                dataTable = onlineData(i+1:end, :);
                break
            end
        end
    end

    % find the end point of the table
    for i = 1:size(dataTable,1)
        try
            if strcmp(dataTable{i,1}, '[Events]')
                dataTable = dataTable(1:i-1, :);
                break
            end
        end
    end

    % clean up column names
    for i = 1:length(tableHeaders)
        tableHeaders(i) = strrep(tableHeaders(i), 'Unit 2.', '');
    end

    dataTableT = cell2table(dataTable, "VariableNames", tableHeaders);

    % time columns to datetime
    dataTableT = convertvars(dataTableT, 'Timestamp', 'datetime');
    dataTableT = convertvars(dataTableT, 'InoculationTime []', 'datetime');

    % find inoculation begin and remove rows before that
    for t = 1:height(dataTableT)
        if ~isnat(dataTableT{t, 'InoculationTime []'})
            dataTableT = dataTableT(t:end, :);
            break
        end
    end
    % remove duration and InoculationTime columns
    dataTableT                  = removevars(dataTableT,{'Duration',...
        'InoculationTime []'});

    % convert timestamps into duration
    startTime                   = dataTableT{1, 'Timestamp'};
    dataTableT.('time [d]')     = hours(dataTableT{:, 'Timestamp'} - startTime)/24;
    dataTableT                  = removevars(dataTableT,{'Timestamp'});

    % get column headers and convert strings into numerical data
    colHeaders  = dataTableT.Properties.VariableNames;
    cols2keep   = {'time [d]', 'DO2.PV [%DO]', 'pH2.PV [pH]', 'N2.PV [rpm]',...
                   'VA2.PV [mL]', 'VB2.PV [mL]', 'FC2.PV [mL/h]'};                        % FC2.PV [mL/h]: feed pump actual value;  N2.PV [rpm]: stirrer speed;  VA2.PV [mL]: totalised acid volume; VB2.PV [mL]: totalised base volume

    myVars      = {'time [d]', 'Dissolved O_2 [%]', 'pH',...
                   'stirrer speed [rpm]',  'V(H_2SO_4) [ml]', 'V(KOH) [ml]',...
                   'pump rate [ml/h]'};

    finalTable  = dataTableT(:, cols2keep);
    finalTable.Properties.VariableNames = myVars;

    % convert strings to numbers in table
    N                   = NaN(size(finalTable));
    numberTable         = array2table(N, 'VariableNames', myVars);
    numberTable(:,1)    = finalTable(:,1);

    for col = 2:size(finalTable, 2)
        for row = 1:size(finalTable, 1)
            numberTable{row, col}   = str2double(finalTable{row, col});
        end
    end

    %% make individual tables and delete nan rows
    % artificially set starting values for pH, pump rate, stirrer speed, acid and base
    numberTable{1, 'pH'}                    = 7.2;
    numberTable{1, 'pump rate [ml/h]'}      = 0;
    numberTable{1, 'stirrer speed [rpm]'}   = 400;
    numberTable{1, 'V(H_2SO_4) [ml]'}       = 0;
    numberTable{1, 'V(KOH) [ml]'}           = 0;

    individualTables = cell(length(myVars)-1, 1);
    for i = 1:length(myVars)-1
        individualTables{i} = numberTable(:, [1, i+1]);
        individualTables{i} = rmmissing(individualTables{i});
    end

    % save the tables
    save([FileLocation, onlineDataSaveFile], 'individualTables')
else
    load([FileLocation, onlineDataSaveFile])
end

%% import offline data
tCDW            =   readmatrix(fullFile, 'Sheet', sheetName, 'Range', 'A:B');
tCDW(1, :)      =   [];
CDW             =   tCDW(:,2);
tCDW(:,2)       =   [];
tCDW            =   tCDW./24;

% load in data for pump rate changes
t0              =   readtable(fullFile, 'Sheet', t0Sheet{1}, 'Range',...
                              t0Sheet{2});
ratechanges     =   readtable(fullFile, 'Sheet', pumpRateSheet{1}, 'Range',...
                              pumpRateSheet{2});
changetimes     =   datenum(ratechanges{:,1}-t0{:, :});
changetimes     =   [0; changetimes];
myrates         =   ratechanges{:,5};
myrates         =   [0; myrates];

% load in data for the C/N ratio change
ratioChanges    =   readtable(fullFile, 'Sheet', CNratioSheet{1}, 'Range',...
                              CNratioSheet{2});
changetimesCN   =   datenum(ratioChanges{:,6}-t0{:, :});
CNratio         =   ratioChanges{:,5};
c_C             =   ratioChanges{:,1};
carbonSource    =   ratioChanges{:,8};
carbonSource    =   [num2cell(changetimesCN), carbonSource];
toDelete        =   [];
for i = 1:size(carbonSource, 1)
    if isempty(carbonSource{i,2})
        toDelete = [toDelete, i];
    end
end
carbonSource(toDelete, :) = [];


%% plotting all the data
varNames        =   {'DO' 'pH' 'stirrer speed' 'Acid Total' 'Base Total' 'pump rate'};

% get variables to plot
pltVars             =   cell(size(individualTables));
for i = 1:length(pltVars)
    pltVars{i}      =   individualTables{i}.Properties.VariableNames{2};
end
% decide which variables not to plot
% disp('Please enter a vector of the cells you do not want to plot');
% disp(varNames)
% exclude     =   input('');
exclude                     =   6; % to save this step when plotting for the thesis
individualTables(exclude)   =   [];
varNames(exclude)           =   [];
pltVars(exclude)            =   [];

figure('Renderer', 'Painters')
myFig           =   figure(1);
myLineStyles 	=   {'-x'};
myInterpreter   =   'tex';
myFontSize      =   [10 10 10 10];

SqLength        =   ceil(sqrt(size(individualTables,1) + 1));
tL              =   tiledlayout(SqLength,ceil((size(individualTables,1)+1)/SqLength));
ax              =   cell(size(individualTables,1)+1,1);
P               =   cell(size(individualTables,1),1);
stairColour     =   [0.7 0 0];
stairWidth      =   1;

% Biomass
ax{1}   =   nexttile;

% define colours manually (greyscale)
myColours   =   ones(size(carbonSource, 1), 3);
for i = 1:size(myColours, 1)
    myColours(i,:)  =   myColours(i,:) - (0.1 * (i - 1));
end
for i = 1:size(carbonSource, 1)
    % set colour
    myCol   =   myColours(i, :);
    if i < size(carbonSource, 1)
        rectangle('Position', [carbonSource{i, 1},...       % x
                               0,...                        % y
                               carbonSource{i+1, 1},...     % width
                               2000],...                    % height
                  'FaceColor', myCol, 'EdgeColor', myCol)
    else
        rectangle('Position', [carbonSource{i, 1},...       % x
                               0,...                        % y
                               1000,...                     % width
                               2000],...                    % height
                  'FaceColor', myCol, 'EdgeColor', myCol)
    end
    hold on
end

grid on
set(gca, "Layer", "top")

yyaxis left
PCDW    =   myplot(tCDW, CDW, {'', {'CDW [g/l]', ''}}, {'Cell Dry Weight'},...
                   [], myLineStyles, myInterpreter, myFontSize, [], ax{1},...
                   1, [0 0 0]);
% xticks(0:10:max(tCDW)+1)
xlim([0,max(tCDW)+0.1])
% title('Cell Dry Weight', 'FontSize', 12, 'FontName', 'Arial')
ylim([0, max(CDW) * 1.01])

% dilution rates
yyaxis right
stairs([changetimes;max(tCDW)], [myrates;myrates(end)], '-',...
       'Color', stairColour, 'LineWidth', stairWidth)
hold on

% C/N ratio
stairs([changetimesCN; max(tCDW)], [CNratio; CNratio(end)].*0.001, '--',...
       'Color', stairColour, 'LineWidth', stairWidth)

% c_FA
stairs([changetimesCN; max(tCDW)], [c_C; c_C(end)], '-.',...
       'Color', stairColour, 'LineWidth', stairWidth)

xlim([0,max(tCDW)+0.1])
ylim([0, 0.18])
yticks(0:0.03:0.18)
set(gca, 'YColor', stairColour)
yyaxis left
set(gca,'YColor','k')

% online data
for i = 1:size(individualTables,1)
    ax{i+1}     =   nexttile;

    % plot a rectangle for the carbon source
    for j = 1:size(carbonSource, 1)
        % set colour
        myCol   =   myColours(j, :);
        if j < size(carbonSource, 1)
            rectangle('Position', [carbonSource{j, 1},...       % x
                                   0,...                        % y
                                   carbonSource{j+1, 1},...     % width
                                   2000],...                    % height
                      'FaceColor', myCol, 'EdgeColor', myCol)
        else
            rectangle('Position', [carbonSource{j, 1},...       % x
                                   0,...                        % y
                                   1000,...                     % width
                                   2000],...                    % height
                      'FaceColor', myCol, 'EdgeColor', myCol)
        end
        hold on
    end

    grid on
    set(gca, "Layer", "top")

    curTable    =   individualTables{i};
    x           =   curTable{:,1};
    y           =   curTable{:,2};
    P           =   myplot(x,y,{'',{pltVars{i},''}},pltVars(i),...
                    [],{'-'},myInterpreter,myFontSize,[],ax{i+1},1,[0 0 0]);
%     xticks(0:10:max(tCDW)+1)
    xlim([0,max(tCDW)+0.1])
%     title(varNames{i}, 'FontSize', 12, 'FontName', 'Arial')

    if strcmp(varNames{i}, 'pH')
        yticks(6:1:8)
        ylim([6,8])
    else
        ylim([0,max(y)*1.01])
    end

    % dilution rates and carbon concentration
    yyaxis right
    stairs([changetimes; max(tCDW)], [myrates; myrates(end)], '-',...
       'Color', stairColour, 'LineWidth', stairWidth)
    hold on
    % C/N ratio
    stairs([changetimesCN; max(tCDW)], [CNratio; CNratio(end)].*0.001, '--',...
       'Color', stairColour, 'LineWidth', stairWidth)

    % c_FA
    stairs([changetimesCN; max(tCDW)], [c_C; c_C(end)], '-.',...
       'Color', stairColour, 'LineWidth', stairWidth)

    xlim([0,max(tCDW)+0.1])
    ylim([0, 0.18])
    yticks(0:0.03:0.18)
    set(gca, 'YColor', stairColour)
    yyaxis left
    set(gca,'YColor','k')


end
xlabel(tL,'t [d]','FontName','Arial','FontSize',10,'Interpreter','tex')
ax              =   axes(myFig);
han             =   gca;
han.Visible     =   'off';
yyaxis(ax, 'right');
ylabel({'','D [h^{-1}] (-);  C/N ratio x10^{-3} (- -);  c_{carbon} (.-) [mol/L]'}, 'FontName', 'Arial',...
    'FontSize', 10, 'Interpreter', 'tex', 'Color', stairColour)
% ylabel({'','D [h^{-1}] (-);    c_{C, feed} [mol/L] (--)'}, 'FontName', 'Arial',...
%     'FontSize', 14, 'Interpreter', 'tex', 'Color', 'r')
han.YLabel.Visible = 'on';
set(myFig,'Position',[20 50 600 600])

%% save the figure(s)
if saveFigure
    savefig(myFig,[FigurePath FigureName '.fig'])
    exportgraphics(myFig,[FigurePath FigureName '.emf'],'BackgroundColor','w');
    print(myFig, [FigurePath FigureName '.svg'], '-dsvg')
end