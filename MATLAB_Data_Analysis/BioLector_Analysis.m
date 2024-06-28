%% BioLector Analysis Version 3
% newer version for analysing rfp assays, based on pure triplicate runs
% this one also calculates growth rates
% Author: Manuel Bruch, PhD student
% Creation:     29.04.2021
% Last change:  10.07.2023

%% Tidy up
clear all
close all
clc


%% define script parameters
fileLocation    =   '..\..\Laboratory_files\ALE_screening\20211124_ALE_test\';
dataFile        =   '211129_Cn_ALE_growth_test_final_TRANSFORM_CALIBIOMASS_CALIPH_CALIPO2_US.csv';
myDate          =   '211129';
ResultFile      =   [myDate '_ALE_Results.xlsx'];
sampleCodesFile =   [myDate '_Sample_Codes.csv'];
biomassFile     =   'Biomass_conv_factor.xlsx';
SignalNames     =   {'biomass_30', 'biomass_20', 'pH', 'pO2'};

% for growth rate calculation
StepWidth       =   20;

figureFolder    =   [fileLocation '\Figures\' 'Signals_' myDate '\'];
genFigName      =   '_thesis';
myAxisNames     =   {'time [h]', 'biomass [g/l]', 'biomass [g/l]', 'pH',...
                     'pO_2', 'µ [h^{-1}]'};
figSize         =   [20 50 600 650];

myLegendLocation    =   'eastoutside';
myLineStyles        =   '-';
myInterpreter       =   'tex';
myFontSize          =   [];
myTitle             =   [];

% what to plot and if to save
refPltDensity   =   5;  % select every how manieth data point is plotted to reduce figure size --> only set > 1 if it does not distort the graph
plotting        =   1;
plotIndividual  =   1;
saveFigure      =   1;

% decide whether or not to save numerical data
writingData     =   0;
overWrite       =   1;

%% data import
figureNames.Signals     =   strrep(SignalNames,' ','');
figureNames.Signals     =   strrep(figureNames.Signals,'(','_');
figureNames.Signals     =   strrep(figureNames.Signals,')','');
figureNames.Signals     =   strrep(figureNames.Signals,'=','_');
figureNames.Signals     =   [figureNames.Signals,'mu'];

% import CHANNEL names as text to filter out the Cal.Biomass signals
opts        =   detectImportOptions([fileLocation dataFile],...
                                    'NumHeaderLines',22);
opts        =   setvartype(opts,{'CHANNEL'},'char');
% import the data with the defined options
myData      =   readtable([fileLocation dataFile], opts);

SampleIDs   =   readcell([fileLocation sampleCodesFile], 'Delimiter', ',');     % requires file that maps description to each sample ID

myData.DESCRIPTION      =   [];
CalBiomassIdx           =   find(contains(myData.CHANNEL,'Cal.Biomass'));
myData(CalBiomassIdx,:) =   [];

% delete the general data at the bottom of the table:
myData(strlength(myData.CHANNEL)<1,:)   =   [];

% blank IDs based on Sample_Codes file
blankIDs    =   SampleIDs(contains(SampleIDs(:,2),'blank'),1);
blankType   =   SampleIDs(contains(SampleIDs(:,2),'blank'),2);
if find(all(contains(blankType,'w_kan')))
    kanBlank    =   find(contains(blankType,'w_kan'));
    nokanBlank  =   find(~contains(blankType,'w_kan'));
elseif find(all(contains(blankType,'lactate')))
    kanBlank    =   find(contains(blankType,'+'));
    nokanBlank  =   find(contains(blankType,'-'));
else
    kanBlank    =   [];
    nokanBlank  =   find(~contains(blankType,'w_kan'));
end

% identify samples that do not contain kan (i.e. WT)
if find(all(contains(blankType,'w_kan')))
    NonKanIDs   =   SampleIDs(contains(SampleIDs(:,2),['WT','ALE']),1);
elseif find(all(contains(blankType,'lactate')))
    NonKanIDs   =   SampleIDs(:,1);     % treat all samples with the non-lactate containing blank
else
    NonKanIDs   =   SampleIDs(~contains(SampleIDs(:,2),'blank'),1);
end

% define table headers:
myHeaders   =   strrep(SampleIDs,'= ','');
myHeaders   =   strrep(myHeaders,' ','_');
myHeaders   =   strrep(myHeaders,'{','');
myHeaders   =   strrep(myHeaders,'}','');
myHeaders   =   strrep(myHeaders,'+','');
myHeaders   =   strrep(myHeaders,'-','');

%% analyse groups of data
myTime          =   myData{2,4:end};
myTime2         =   cat(3,myTime',zeros(size(myTime')));
measuredData    =   myData(3:end,:);
% set all 'X#' to X0' if # < 10 --> G and ID will have the correct order in
% findgroups()
for i = 1:height(measuredData)
    if length(measuredData.CONTENT{i}) < 3  % test length of string
        measuredData.CONTENT{i}     =   strrep(measuredData.CONTENT{i},...
                                        'X','X0');
    end
end

% exclude certain wells
if contains(fileLocation, '20211111')
    allWells = measuredData.WELLNo_;
    measuredData(strcmp(allWells,'A02'),:)  =   [];
end

% first split by signal
[G,ID]          =   findgroups(measuredData.CHANNEL);
Signals         =   cell(length(ID)+1,1);
blankData       =   Signals;
meanSignals     =   Signals;
stdDevSignals   =   Signals;
newData         =   Signals;
TrpData         =   Signals;
SglData         =   Signals;
AllData         =   Signals;
AllTbl          =   cell(length(ID)+1,2);

SignalNames     =   [SignalNames,transpose(ID(strlength(ID)>1))];

% then calculate average for all replicates in each channel
for i = 1:length(ID)
    Signals{i}  =   measuredData(strcmp(measuredData.CHANNEL,ID{i}),:);
    % extract blanks
    blankData{i}=   Signals{i}(contains(Signals{i}.CONTENT,blankIDs),:);
    Signals{i}  = Signals{i}(~contains(Signals{i}.CONTENT,blankIDs),:);
    
    % one specific blank is far too high --> delete
    if contains(myDate,'210816') && i == 2
        blankData{i}(end,:)     =   [];
    end
    
    % subtract blanks
    [G2,ID2]    =   findgroups(blankData{i}.CONTENT);
    if strlength(ID{i}) <= 1    % only subtract blanks from biomass measurements
        meanBlanks  =   splitapply(@(x) mean(x,1),...
                        table2array(blankData{i}(:,4:end)),G2);
        
        % substract blank without kan
        Signals{i}{contains(Signals{i}.CONTENT,NonKanIDs),4:end}  =   ...
            bsxfun(@minus, Signals{i}{contains(Signals{i}.CONTENT,...
            NonKanIDs),4:end}, meanBlanks(nokanBlank,:));

        % substract blank with kan
        if ~isempty(NonKanIDs)
            Signals{i}{~contains(Signals{i}.CONTENT,NonKanIDs),4:end}  =...
                bsxfun(@minus, Signals{i}{~contains(Signals{i}.CONTENT,...
                NonKanIDs),4:end}, meanBlanks(kanBlank,:));
        end
    end
    
    [G2,ID2]    =   findgroups(Signals{i}.CONTENT);
    meanSignals{i}      =   splitapply(@(x) mean(x,1),table2array(...
                            Signals{i}(:,4:end)),G2);
    meanSignals{i}      =   transpose(meanSignals{i});
    % applies the function in the first argument to all groups in the
    % second argument. The function is defined as an "anonymous function"
    % (@(InputVar1,InputVar2,...)) followed by what the function does with
    % the inputs
    stdDevSignals{i}    =   splitapply(@(x) std(x,[],1),table2array(...
                            Signals{i}(:,4:end)),G2);
    stdDevSignals{i}    =   transpose(stdDevSignals{i});
    newData{i}          =   cat(3,meanSignals{i},stdDevSignals{i});
    
    % find all measurements that are not done in triplicates: 
    % (https://uk.mathworks.com/matlabcentral/answers/354020-find-values-in-list-that-appear-exactly-once)
    tmp                 =   sort(G2);       
        % sorts array in ascending numerical order
    idp                 =   diff(tmp)>0;    
        % tests if neighbouring elements in an array have a difference > 0
    single_runsIDX      =   tmp([true;idp]&[idp;true]);
    
    % triplicate data:
    TrpData{i}  =   newData{i}(:,setdiff(1:end,single_runsIDX),:);
    % single run data:
    SglData{i}  =   newData{i}(:,single_runsIDX,:);
    
    % combine both
    AllData{i}  =   [TrpData{i},SglData{i}];
    
    % also write triplicate data into table format:
    for j = 1:2
        TrpTbl      =   array2table([myTime2(:,:,j) TrpData{i}(:,:,j)],...
                        'VariableNames',['Time_h';myHeaders(~contains(...
                        myHeaders(:,2),'blank') & ~contains(myHeaders(:,...
                        2),'ALE'),2)]);
        % add on any sglData
        SglTable    =   array2table([myTime2(:,:,j) SglData{i}(:,:,j)],...
                        'VariableNames',['Time_h';myHeaders(contains(...
                        myHeaders(:,2),'ALE'),2)]);
        try
            AllTbl{i,j} =   join(TrpTbl,SglTable,'Keys','Time_h');
        catch
            AllTbl{i,j} =   [TrpTbl,SglTable(:,2:end)];
        end
    end
end

% remove "raw" pO2 and pH data
removeThese     =   [find(strcmp(SignalNames,'pO2')),find(strcmp(...
                    SignalNames,'pH'))];
TrpData(removeThese)    =   [];
SglData(removeThese)    =   [];
ID(removeThese)         =   [];
SignalNames(removeThese)=   [];
AllTbl(removeThese,:)   =   [];

%% calculate growth rates
% convert the biomass to g/l
if strcmp(myDate,'210816')
    myODs   =   readmatrix(['..\..\Laboratory_files\BioLector\',myFolder,...
                '\Biomass_conv_factor.xlsx'],'Range','B6:P6');
    Fs      =   myODs./TrpData{1}(1,:,1);
    avgF    =   mean(Fs);
    stdF    =   std(Fs);
    FA_CDWconversionFactor  =   ...
        load('..\..\Laboratory_files\FA_CDWconversionFactor.mat');
    FA_CDW_Factor           =   FA_CDWconversionFactor.FAfactor;
    ItsTheFinalFactor       =   avgF * FA_CDW_Factor;
    save(['..\..\Laboratory_files\BioLector\' myFolder '\conv_factor.mat'],...
        'ItsTheFinalFactor')
elseif strcmp(myDate,'211129')
    myODs   =   readmatrix([fileLocation biomassFile],'Range','B30:AR30');
    Fs      =   myODs./AllTbl{1}{end,2:end}; % for a gain of 30
    avgF    =   mean(Fs);
    stdF    =   std(Fs);
    FA_CDWconversionFactor  =   ...
        load('..\..\Laboratory_files\FA_CDWconversionFactor.mat');
    FA_CDW_Factor           =   FA_CDWconversionFactor.factor;
    ItsTheFinalFactor       =   avgF * FA_CDW_Factor;
    save([fileLocation 'conv_factor.mat'], 'ItsTheFinalFactor')
else
    ItsTheFinalFactor       =   load('..\..\Laboratory_files\ALE_screening\20211124_ALE_test\conv_factor.mat');
    ItsTheFinalFactor       =   ItsTheFinalFactor.ItsTheFinalFactor;
end
for i = 1:size(AllTbl,1)-1
    if ~contains(ID{i}, 'Cal.')
        try
            AllTbl{i,1}{:,2:end}  =   AllTbl{i,1}{:,2:end} * ItsTheFinalFactor;
            AllTbl{i,2}{:,2:end}  =   AllTbl{i,2}{:,2:end} * ItsTheFinalFactor;
        end
    end
end

% check which signal is biomass_30:
biomassPos                          =   find(contains(SignalNames,'_30'));
allMyBiomass                        =   Signals{biomassPos}{:,4:end};
% allMyBiomass(allMyBiomass <= 0)     =   nan;    % otherwise log will result in imaginary numbers
% rather: add a large number as that should not change the slope but shift the
% entire curve upwards
allMyBiomass                        =   allMyBiomass + 10;

% store data for all in the last cell of the "Signals" array
SWrates     =   slidingWindow(transpose(myTime),log(transpose(...
                allMyBiomass) * ItsTheFinalFactor),StepWidth,nan);
SWtimes     =   zeros(length(SWrates.time_windows),1);
Signals{end}=   SWrates.coefficients.m;

for i = 1:length(SWtimes)
   SWtimes(i)   =   mean(SWrates.time_windows{i}); 
end
SWtimes2    =   cat(3,SWtimes,zeros(size(SWtimes)));

[G3,ID3]            =   findgroups(Signals{biomassPos}.CONTENT);
meanSignals{end}    =   splitapply(@(x) mean(x,2),Signals{end},G3');
stdDevSignals{end}  =   splitapply(@(x) std(x,[],2),Signals{end},G3');
newData{end}        =   cat(3,meanSignals{end},stdDevSignals{end});
AllData{end}        =   newData{end};
for j = 1:2
    AllTbl{end,j} =   array2table([SWtimes2(:,:,j) AllData{end}(:,:,j)],...
        'VariableNames',['Time_h';myHeaders(1+length...
        (blankIDs):end,2)]);
end

SignalNames         =   [SignalNames,'mu'];

disp('These are the µ(max)')
disp(max(AllTbl{end,1}{:,2:end}))


%% plot data
if plotting
    if plotIndividual
        for i = 1:length(SignalNames)
            % define data for current figure
            x       =   cat(3, AllTbl{i,1}{:, 1},...
                               AllTbl{i,2}{:, 1});
            allY    =   cat(3, AllTbl{i,1}{:, 2:end},...
                               AllTbl{i,2}{:, 2:end});
            myYmax  =   max(allY(:,:,1),[],'all');
            myYmin  =   min(allY(:,:,1),[],'all');

            % define figure parameters
            myFig   =   figure(i);
            set(gcf,'renderer','painters')

            colNum  =   6;
            t       =   tiledlayout(ceil((size(allY, 2) -1) / colNum), colNum,...
                                    'TileSpacing', 'none');
            % get vector with plot numbers that will be at the left edge
            edgePlots   =   1:colNum:ceil((size(allY, 2) -1) / colNum) * colNum;
            for j = 2:size(allY,2)
                ax      =   nexttile;
                if all(allY(:,2,2)) == 0 % for single replicates, all stddevs will be 0 in the table
                    y       =   allY(:,j,1);
                    x       =   x(:,1,1);
                else
                    y       =   allY(:,j,:);
                end
                % plot the WT triplicate
                myNewPlot(x(1:refPltDensity:end),...
                          allY(1:refPltDensity:end,1,:),...
                          'myyLim', [myYmin, myYmax],...
                          'myLSOrder', {'-'}, 'myColor', [0 0 0],...
                          'transpErr', 0.1)
                % plot the current ALE colony
                myNewPlot(x, y, 'myyLim', [myYmin, myYmax],...
                          'myLSOrder', {'-'}, 'myColor', [0.7 0 0])
                ylim([myYmin myYmax])

                if ~ismember(j-1, edgePlots)
                    set(gca,'YTickLabel',[])
                end
                if j < edgePlots(end)
                    set(gca,'XTickLabel',[])
                end
            end
            %define axis labels
            xlabel(t, {'', myAxisNames{1}}, 'FontName', 'Arial', 'FontSize', 10)
            ylabel(t, {myAxisNames{i+1}, ''}, 'FontName', 'Arial',...
                   'FontSize', 10)
            lgd             =   legend({'H16', 'ALE colony'}, 'FontName',...
                                       'Arial', 'FontSize', 10,...
                                       'Interpreter', 'tex', 'NumColumns', 2);
            lgd.Layout.Tile =   'north';

            set(myFig, 'Position', figSize)
            curFigName  =   figureNames.Signals{i};
            % save the figure
            if saveFigure
                FigSavePath     =   [figureFolder curFigName genFigName];
                savefig(myFig, [FigSavePath '.fig'])
                exportgraphics(myFig, [FigSavePath '.emf']);
                print(myFig, [FigSavePath '.svg'], '-dsvg')
            end
        end
    end
end
%% save data in excel sheet
if writingData
    myWriteFile     =   [fileLocation, ResultFile];
    TableTypes      =   {'mean','std'};
    figureNames.Signals     =   [figureNames.Signals 'mu'];

    if overWrite
        recycle on
        delete(myWriteFile);
    end

    for i = 1:size(AllTbl,1)
        for j = 1:size(AllTbl,2)
                    writetable(AllTbl{i,j},myWriteFile,'Sheet',[...
                        figureNames.Signals{i} '_' TableTypes{j}])
        end
    end
end