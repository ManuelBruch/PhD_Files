function [ResultTable,SampleName] = HPLC_Sample_Analysis(FileLocation,myDate,SampleNr,MySamples,StdRetTimes,myStandards,SampleAmount,timePoints,GeneralFileName)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% determine internal sample number
ISN         =   SampleNr;
SampleName  =   MySamples{1};
Results     =   zeros(size(MySamples,1),size(myStandards,1));               % Rows: Timepoints; Columns: Standards

%loop through all lines in MySamples (corresponding to all time points)
for i = 1:length(MySamples)
    % load in data
    SampleFile  =   [FileLocation myDate '\peak_tables\' GeneralFileName num2str(ISN) '.txt']; % changed this so it may be controlled from the main script (12.02.2021)
    warning('off','all')
    RawData     =   readcell(SampleFile,'Delimiter','\t');
    warning('on','all')
    
    % find peak table
    PKTstart    =   find(cellfun(@(c) ischar(c) && ~isempty(strfind(c, '[Peak Table(Detector A-Ch1)]')), RawData));
    PKTend      =   find(cellfun(@(c) ischar(c) && ~isempty(strfind(c, '[Peak Table(PDA-Ch1)]')), RawData));
    PeakTable   =   cell2table(RawData(PKTstart+3:PKTend-1,:),'VariableNames',RawData(PKTstart+2,:));
    
    % get retention times and areas of all peaks
    SampleRetTimes  =   table2array(PeakTable(:,2));
    SamplePeakAreas =   table2array(PeakTable(:,5));
    
    % loop through all standards and check RetTimes
    for j = 1:length(StdRetTimes)
        RetTimeDiff     =   SampleRetTimes-StdRetTimes(j);
        PotPeakIdx      =   find(RetTimeDiff >= -0.05 & RetTimeDiff <= 0.05);
        if ~isempty(PotPeakIdx)
            Concentration   = (SamplePeakAreas(PotPeakIdx) - ...
                myStandards{j,2}.StandardCurve(2,1)) / ...
                myStandards{j,2}.StandardCurve(1,1);
            Results(i,j)    = Concentration;                                % may run into trouble if multiple Peaks are found
        end
    end
    ISN     =   ISN+SampleAmount;
end

ResultTable     =   array2table(Results,'VariableNames',myStandards(:,1));

%% plot the results
SubPlotRows     =   ceil(sqrt(length(StdRetTimes)));
SubPlotCols     =   floor(sqrt(length(StdRetTimes)));
myAxisNames     =   {'t [h]','concentration [g/l]'};

myFig   =   figure(SampleNr);
for i = 1:length(StdRetTimes)
    ax       =   subplot(SubPlotRows,SubPlotCols,i);
    myLineStyles 	=   {'-x','--^',':o','-.d'};
%     P{i,j}   =  myplot(timePoints,Results(:,i),myAxisNames,...
%                 MySamples,[],myLineStyles,[],[],[myStandards{i,1}],ax,0);
    P{i,j}   =  myplot(timePoints,Results(:,i),myAxisNames,...
                myStandards(i,1),[],myLineStyles,[],[],[myStandards{i,1}],ax,0);
    ylim([0 Inf])
    sgtitle(['Analysis of Growth on ' SampleName],'FontName','Arial','FontSize',14,'Interpreter','none')
end
set(myFig,'Position',[20 50 1200 600])
% save the figure(s)
savefig(myFig,[FileLocation myDate '\Condition_' SampleName '.fig'])
exportgraphics(myFig,[FileLocation myDate '\Condition_' SampleName '.emf']);
