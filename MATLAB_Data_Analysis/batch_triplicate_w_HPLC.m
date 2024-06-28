%% Analysis of batch growth curve triplicates with HPLC data
% Author: Manuel Bruch, PhD Student
% Created: 2020/11/22
% last edited: 2024/06/12

%% Tidy up
clear all
close all
clc

%% define script parameters
FileLocation    =   '..\..\Laboratory_files\growth_experiments_2021\';
FileName        =   'C_necator_bioreactor_cultivation.xlsx';
FullPath        =   [FileLocation, FileName];
relevantSheet   =   {'Data_01-05_02_2021', 'Data_08-13_02_2021', 'Data_15-20_02_2021'};
relevantRange   =   {'A11:B20', 'A11:B20', 'A11:B20'};

% for plots:
colNum          =   3;
cSources        =   {'FA'};
concentrations  =   {'80'};
myAxisNames     =   {'t [h]', 'CDW [g/l]',...
                     'c_{carbon} [C-mol/L] (-x),    µ [h^{-1}] (- -\Delta)'...
%                      {'c_{carbon} [C-mol/L] (-)',...
%                      'µ [h^{-1}] (--)'}...
                     };
FigSize         =   [20 50 600 400];
primAxCol       =   [0 0 0];
secAxCol        =   [0 0.5 0];

FigureNames     =   'batch_triplicates_w_HPLC_and_Mu';
FigSavePath     =   [FileLocation FigureNames];

% parameters where to save results
writeResults    =   1;
saveFigures     =   1;
saveFileName    =   'C_necator_bioreactor_cultivation.xlsx';
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
                                      'Range', relevantRange{i});
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
HPLCFileLocation=   '..\..\Laboratory_files\HPLC_Files\batch_cultivations\fermenter_averages\';
myUnits         =   {'Cmol_l'};
HPLCdata        =   cell(size(myUnits));

for i = 1:length(myUnits)
    c1  =   cell(1,length(conditions)+1);                                       % store means here
    c2  =   c1;                                                                 % store std deviations here
    for j = 1:length(conditions)
        FileName    =   [HPLCFileLocation 'batch_fermentation_' myUnits{i} '.txt']; % by depending on Cond the data in the resulting matrix should be in the right order
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

% plot biomass
yyaxis left

myNewPlot(tCDW, CDW, 'myColor', primAxCol, 'transpErr', 0.5)

xlim([0, max(tCDW(:,:,1)) + max(tCDW(:,:,2)) + 0.1])
ylim([0, ceil((max(CDW(:,:,1),[],'all') + max(CDW(:,:,2),[],'all')))])

% plot HPLC and growth rate data
yyaxis right
set(gca, 'YColor', secAxCol)

myNewPlot(tHPLC, HPLC, 'myColor', secAxCol, 'transpErr', 0.5)
myNewPlot(tMu, Mu, 'myColor', secAxCol,...
          'myLSOrder', {'-x', '--^'}, 'transpErr', 0.5)

xlim([0, max(tHPLC(:,:,1)) + max(tHPLC(:,:,2)) + 0.1])
% ylim([0, ceil(...
%               max(...
%                   max(HPLC(:,:,1),[],'all') + max(HPLC(:,:,2),[],'all'),...
%                   max(Mu(:,:,1),[],'all') + max(Mu(:,:,2),[],'all')...
%               )*20 ...
%           )/20])
ylim([0, 0.16])

yyaxis left
set(gca, 'YColor', primAxCol)

xlabel({'', myAxisNames{1}}, 'FontName', 'Arial', 'FontSize', 10,...
       'Interpreter','tex')
ylabel({myAxisNames{2}, ''}, 'FontName', 'Arial', 'FontSize', 10,...
       'Interpreter','tex')
% ax              =   axes(myFig);
% han             =   gca;
% han.Visible     =   'off';
% yyaxis(ax, 'right');
yyaxis right
ylabel({'', myAxisNames{3}}, 'FontName', 'Arial', 'FontSize', 10,...
       'Interpreter', 'tex', 'Color', secAxCol)
% han.YLabel.Visible = 'on';

set(myFig, 'Position', FigSize)

%% save the figure(s)
if saveFigures
    savefig(myFig, [FigSavePath '.fig'])
    exportgraphics(myFig, [FigSavePath '.emf']);
    print(myFig, [FigSavePath '.svg'], '-dsvg')
end

%% write results into excel
if writeResults
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

    writecell({'average'}, SavePath, 'Sheet', 'batch_triplicate_biomass', 'Range', 'A1')
    writetable(biomassAvg, SavePath, 'Sheet', 'batch_triplicate_biomass', 'Range', 'A2')
    writecell({'stddev'}, SavePath, 'Sheet', 'batch_triplicate_biomass', 'Range', 'A20')
    writetable(biomassStd, SavePath, 'Sheet', 'batch_triplicate_biomass', 'Range', 'A21')
    writecell({'maxCDW'}, SavePath, 'Sheet', 'batch_triplicate_biomass', 'Range', 'A40')
    writecell({'maxCDWstd'}, SavePath, 'Sheet', 'batch_triplicate_biomass', 'Range', 'A41')
    writematrix(M, SavePath, 'Sheet', 'batch_triplicate_biomass', 'Range', 'B40')
    writematrix(maxStd, SavePath, 'Sheet', 'batch_triplicate_biomass', 'Range', 'B41')
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

    writecell({'average'}, SavePath, 'Sheet', 'batch_triplicate_HPLC', 'Range', 'A1')
    writetable(HPLCAvg, SavePath, 'Sheet', 'batch_triplicate_HPLC', 'Range', 'A2')
    writecell({'stddev'}, SavePath, 'Sheet', 'batch_triplicate_HPLC', 'Range', 'A20')
    writetable(HPLCStd, SavePath, 'Sheet', 'batch_triplicate_HPLC', 'Range', 'A21')
    writecell({'maxc_C'}, SavePath, 'Sheet', 'batch_triplicate_HPLC', 'Range', 'A40')
    writecell({'maxc_Cstd'}, SavePath, 'Sheet', 'batch_triplicate_HPLC', 'Range', 'A41')
    writematrix(M, SavePath, 'Sheet', 'batch_triplicate_HPLC', 'Range', 'B40')
    writematrix(maxStd, SavePath, 'Sheet', 'batch_triplicate_HPLC', 'Range', 'B41')

    % write growth rate data - calculate from µmax from individual replicates as
    % µmax is reached at different time points
    muAvg       =   array2table([tMu(:,:,1) Mu(:,:,1)],...
                            'VariableNames', ['time_h' conditions]);
    muStd       =   array2table([tMu(:,:,2) Mu(:,:,2)],...
                                'VariableNames', ['time_h' conditions]);

    maxMu   =   cell(1, size(allData,2));
    for i = 1:size(allData,2)
        maxMu{i}    =   max(allData{2,i}(:,2:end),[],1);
    end
    MuMaxAvg    =   mean(cat(3, maxMu{:}), 3);
    MuMaxStd    =   std(cat(3, maxMu{:}), [], 3);

    writecell({'average'}, SavePath, 'Sheet', 'batch_triplicate_mu', 'Range', 'A1')
    writetable(muAvg, SavePath, 'Sheet', 'batch_triplicate_mu', 'Range', 'A2')
    writecell({'stddev'}, SavePath, 'Sheet', 'batch_triplicate_mu', 'Range', 'A20')
    writetable(muStd, SavePath, 'Sheet', 'batch_triplicate_mu', 'Range', 'A21')
    writecell({'maxMu'}, SavePath, 'Sheet', 'batch_triplicate_mu', 'Range', 'A40')
    writecell({'maxMustd'}, SavePath, 'Sheet', 'batch_triplicate_mu', 'Range', 'A41')
    writematrix(MuMaxAvg, SavePath, 'Sheet', 'batch_triplicate_mu', 'Range', 'B40')
    writematrix(MuMaxStd, SavePath, 'Sheet', 'batch_triplicate_mu', 'Range', 'B41')
end