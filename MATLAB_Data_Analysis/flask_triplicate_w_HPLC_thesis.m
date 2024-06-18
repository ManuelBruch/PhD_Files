%% Analysis of flask growth curve triplicates with HPLC data
% Author: Manuel Bruch, PhD Student
% Created: 2020/11/22
% last edited: 2024/06/12

%% Tidy up
clear all
close all
clc

%% define script parameters
FileLocation    =   '..\..\Laboratory_files\growth_experiments_2021\';
FileName        =   'Growth_Experiments.xlsx';
FullPath        =   [FileLocation, FileName];
relevantSheet   =   {'Week_07-11_12_2020', 'Week_14-18_12_2020', 'Week_05-11_01_2021'};
relevantRange   =   'A10:M22';

% for plots:
colNum          =   3;
cSources        =   {'Gluconate', 'Na-FA', 'FA'};
pltTitles       =   {'Na-Gluconate', 'Na-Formate', 'Formic Acid'};
concentrations  =   {'20', '50', '80', '100'};
myAxisNames     =   {'t [h]', 'CDW [g/l]',...
                     'c_{carbon} [C-mol/L]'...
                     };
%                      'µ [h^{-1}]'...
FigSize         =   [20 50 600 600];
primAxCol       =   [0 0 0];
secAxCol        =   [0 0.5 0];

FigureNames     =   'flask_triplicates_w_HPLC';
FigSavePath     =   [FileLocation FigureNames];

% parameters where to save results
saveFileName    =   'Growth_Experiments_results.xlsx';
SavePath        =   [FileLocation, saveFileName];


%% calculate average CDW and µ
% import CDW conversion factors
FA_CDWconversionFactor  =   load('..\..\Laboratory_files\FA_CDWconversionFactor.mat');
FA_CDW_Factor           =   FA_CDWconversionFactor.factor;
C6_CDWconversionFactor  =   load('..\..\Laboratory_files\NaG_CDWconversionFactor.mat');
C6_CDW_Factor           =   C6_CDWconversionFactor.NaGfactor;

% load CDW data
allData     =   cell(2,length(relevantSheet));
for i = 1:length(relevantSheet)
    allData{1,i}        =   readtable(FullPath, 'Sheet', relevantSheet{i},...
                                      'Range', relevantRange);
    OriginalNames       =   allData{1, i}.Properties.VariableNames;
    % conversion from OD to CDW
    CC6                     =   find(contains(OriginalNames,'Gluconate')|...
                                     contains(OriginalNames,'Frc')| ...
                                     contains(OriginalNames,'Fructose'));
    CFormate                =   find(contains(OriginalNames,'FA'));
    allData{1,i}{:,CC6}     =   allData{1,i}{:,CC6} * C6_CDW_Factor;
    allData{1,i}{:,CFormate}=   allData{1,i}{:,CFormate} * FA_CDW_Factor;

    % growth rates
    curX        =   allData{1,i}{:,1};
    curY        =   allData{1,i}{:,2:end};
    growthRates =   slidingWindow(curX, log(curY), 2, nan);
    allMu       =   growthRates.coefficients.m;
    timeWindows =   growthRates.time_windows;
    tW          =   zeros(length(timeWindows),1);
    for j = 1:length(timeWindows)
        tW(j)       =   mean(timeWindows{j},1);
    end
    % drop the first value due to uncertainty of the real starting OD
    allMu(1,:)  =   [];
    tW(1,:)     =   [];

    % store growth rates in a matrix
    allData{2,i}    =   [tW, allMu];
    allData{1,i}    =   table2array(allData{1,i}); % assuming that OriginalNames is the same for all tables
end

% calculate averages and standard deviations across the arrays
allCDWData_mean     =   mean(cat(3,allData{1,:}), 3);
allCDWData_std      =   std(cat(3,allData{1,:}), [], 3);
allMuData_mean      =   mean(cat(3,allData{2,:}), 3);
allMuData_std       =   std(cat(3,allData{2,:}), [], 3);

tCDW    =   cat(3, allCDWData_mean(:,1), allCDWData_std(:,1));
CDW     =   cat(3, allCDWData_mean(:,2:end), allCDWData_std(:,2:end));
tMu     =   cat(3, allMuData_mean(:,1), allMuData_std(:,1));
Mu      =   cat(3, allMuData_mean(:,2:end), allMuData_std(:,2:end));

conditions  =   strrep(OriginalNames(2:end),'Na_','Na-');                       % HPLC file names are like this


%% import HPLC data
FileLocation    =   '..\..\Laboratory_files\HPLC_Files\shake_flask_averages\';
myUnits         =   {'Cmol_l'};
HPLCdata        =   cell(size(myUnits));

for i = 1:length(myUnits)
    c1  =   cell(1,length(conditions)+1);                                       % store means here
    c2  =   c1;                                                                 % store std deviations here
    for j = 1:length(conditions)
        FileName    =   [FileLocation conditions{j} '_' myUnits{i} '.txt'];     % by depending on Cond the data in the resulting matrix should be in the right order
        if isfile(FileName)
            d           =   readtable(FileName,'Delimiter','\t');
            if j == 1
                c1{j}   =   d{:,1};
                c2{j}   =   d{:,3};
            end
            if size(d{:,1},1) ~= size(c1{1},1)
                c1{j+1} =   NaN(size(c1{1}));
                c2{j+1} =   NaN(size(c1{1}));
                [sharedvals,idx]=   intersect(c1{:,1},d{:,1},'stable');
                c1{j+1}(idx)    =   d{:,2};
                c2{j+1}(idx)    =   d{:,4};
            else
                c1{j+1} =   d{:,2};
                c2{j+1} =   d{:,4};
            end
        else
            c1{j+1} =   NaN(size(c1{1}));
            c2{j+1} =   NaN(size(c1{1}));
        end
    end
    HPLCdata{i}          =   cell2mat(c1);
    HPLCdata{i}(:,:,2)   =   cell2mat(c2);
end
% since only the C-mol/L data are plotted:
tHPLC   =   HPLCdata{1}(:,1,:);
HPLC    =   HPLCdata{1}(:,2:end,:);

%% plot results
figure('Renderer', 'Painters')
myFig           =   figure(1);

tL          =   tiledlayout(length(conditions)/colNum, colNum,...
                            'TileSpacing', 'compact');%, 'Padding', 'compact');

for i = 1:length(concentrations)
    for j = 1:length(cSources)
        if j ~= length(cSources)
            curY    =  find(contains(conditions, concentrations{i})...
                            & contains(conditions, cSources{j}));
        else
            curY    =  find(contains(conditions, concentrations{i})...
                            & contains(conditions, cSources{j})...
                            & ~contains(conditions, 'Na'));
        end
        nexttile;
        % plot biomass
        yyaxis left

        myNewPlot(tCDW, CDW(:,curY,:), 'myColor', primAxCol, 'transpErr', 0.5)

        xlim([0, max(tCDW(:,:,1)) + max(tCDW(:,:,2)) + 0.1])
        ylim([0, ceil((max(CDW(:,:,1),[],'all') + max(CDW(:,:,2),[],'all')))])

        if j > 1
            set(gca,'YTickLabel',[]);
        end

        % plot HPLC and growth rate data
        yyaxis right
        set(gca, 'YColor', secAxCol)

        myNewPlot(tHPLC, HPLC(:,curY,:), 'myColor', secAxCol, 'transpErr', 0.5)
%         myNewPlot(tMu, Mu(:,curY,:), 'myColor', secAxCol,...
%                   'myLSOrder', {'--x'}, 'transpErr', 0.5)

        xlim([0, max(tHPLC(:,:,1)) + max(tHPLC(:,:,2)) + 0.1])
        ylim([0, ceil(...
                      max(...
                          max(HPLC(:,:,1),[],'all') + max(HPLC(:,:,2),[],'all')...
                      )*20 ...
                  )/20])
%                           max(Mu(:,:,1),[],'all') + max(Mu(:,:,2),[],'all')...
        if j < length(cSources)
            set(gca,'YTickLabel',[]);
        end

        yyaxis left
        set(gca, 'YColor', primAxCol)

        if i < length(concentrations)
            set(gca,'XTickLabel',[]);
        end

        if i == 1
            title(pltTitles{j}, 'FontSize', 12, 'FontName', 'Arial')
        end
    end
end

xlabel(tL, {'', myAxisNames{1}}, 'FontName', 'Arial', 'FontSize', 10,...
       'Interpreter','tex')
ylabel(tL, {myAxisNames{2}, ''}, 'FontName', 'Arial', 'FontSize', 10,...
       'Interpreter','tex')
ax              =   axes(myFig);
han             =   gca;
han.Visible     =   'off';
yyaxis(ax, 'right');
ylabel({'','', myAxisNames{3}}, 'FontName', 'Arial', 'FontSize', 10,...
       'Interpreter', 'tex', 'Color', secAxCol)
han.YLabel.Visible = 'on';

set(myFig, 'Position', FigSize)

%% save the figure(s)
savefig(myFig, [FigSavePath '.fig'])
exportgraphics(myFig, [FigSavePath '.emf']);
print(myFig, [FigSavePath '.svg'], '-dsvg')

%% write results into excel
% write biomass data
biomassAvg  =   array2table([tCDW(:,:,1) CDW(:,:,1)],...
                            'VariableNames', ['time_h' conditions]);
biomassStd  =   array2table([tCDW(:,:,2) CDW(:,:,2)],...
                            'VariableNames', ['time_h' conditions]);

[M,I]       =   max(biomassAvg{:,2:end},[],1);
maxStd      =   zeros(size(M));
for i = 1:length(I)
    maxStd(i)   =   biomassStd{I(i),i+1};
end

writecell({'average'}, SavePath, 'Sheet', 'flask_triplicate_biomass', 'Range', 'A1')
writetable(biomassAvg, SavePath, 'Sheet', 'flask_triplicate_biomass', 'Range', 'A2')
writecell({'stddev'}, SavePath, 'Sheet', 'flask_triplicate_biomass', 'Range', 'A20')
writetable(biomassStd, SavePath, 'Sheet', 'flask_triplicate_biomass', 'Range', 'A21')
writecell({'maxCDW'}, SavePath, 'Sheet', 'flask_triplicate_biomass', 'Range', 'A40')
writecell({'maxCDWstd'}, SavePath, 'Sheet', 'flask_triplicate_biomass', 'Range', 'A41')
writematrix(M, SavePath, 'Sheet', 'flask_triplicate_biomass', 'Range', 'B40')
writematrix(maxStd, SavePath, 'Sheet', 'flask_triplicate_biomass', 'Range', 'B41')
% write HPLC data
HPLCAvg     =   array2table([tHPLC(:,:,1) HPLC(:,:,1)],...
                            'VariableNames', ['time_h' conditions]);
HPLCStd     =   array2table([tHPLC(:,:,2) HPLC(:,:,2)],...
                            'VariableNames', ['time_h' conditions]);

[M,I]       =   max(HPLCAvg{:,2:end},[],1);
maxStd      =   zeros(size(M));
for i = 1:length(I)
    maxStd(i)   =   HPLCStd{I(i),i+1};
end

writecell({'average'}, SavePath, 'Sheet', 'flask_triplicate_HPLC', 'Range', 'A1')
writetable(HPLCAvg, SavePath, 'Sheet', 'flask_triplicate_HPLC', 'Range', 'A2')
writecell({'stddev'}, SavePath, 'Sheet', 'flask_triplicate_HPLC', 'Range', 'A20')
writetable(HPLCStd, SavePath, 'Sheet', 'flask_triplicate_HPLC', 'Range', 'A21')
writecell({'maxc_C'}, SavePath, 'Sheet', 'flask_triplicate_HPLC', 'Range', 'A40')
writecell({'maxc_Cstd'}, SavePath, 'Sheet', 'flask_triplicate_HPLC', 'Range', 'A41')
writematrix(M, SavePath, 'Sheet', 'flask_triplicate_HPLC', 'Range', 'B40')
writematrix(maxStd, SavePath, 'Sheet', 'flask_triplicate_HPLC', 'Range', 'B41')

% write growth rate data - calculate from µmax from individual replicates as
% µmax is reached at different time points
maxMu   =   cell(1, size(allData,2));
for i = 1:size(allData,2)
    maxMu{i}    =   max(allData{2,i}(:,2:end),[],1);
end
MuMaxAvg    =   mean(cat(3, maxMu{:}), 3);
MuMaxStd    =   std(cat(3, maxMu{:}), [], 3);

writecell({'average'}, SavePath, 'Sheet', 'flask_triplicate_mu', 'Range', 'A1')
writetable(muAvg, SavePath, 'Sheet', 'flask_triplicate_mu', 'Range', 'A2')
writecell({'stddev'}, SavePath, 'Sheet', 'flask_triplicate_mu', 'Range', 'A20')
writetable(muStd, SavePath, 'Sheet', 'flask_triplicate_mu', 'Range', 'A21')
writecell({'maxMu'}, SavePath, 'Sheet', 'flask_triplicate_mu', 'Range', 'A40')
writecell({'maxMustd'}, SavePath, 'Sheet', 'flask_triplicate_mu', 'Range', 'A41')
writematrix(MuMaxAvg, SavePath, 'Sheet', 'flask_triplicate_mu', 'Range', 'B40')
writematrix(MuMaxStd, SavePath, 'Sheet', 'flask_triplicate_mu', 'Range', 'B41')