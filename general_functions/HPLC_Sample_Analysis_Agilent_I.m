function [ResultTable,SampleName] = HPLC_Sample_Analysis_Agilent_I(FileLocation,myDate,CurFileNrs,MySamples,StdRetTimes,myStandards,SampleAmount,timePoints,GeneralFolderName)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% determine internal sample number
SampleName  =   MySamples{end};
Results     =   zeros(size(MySamples,1),size(myStandards,1));               % Rows: Timepoints; Columns: Standards

%loop through all lines in MySamples (corresponding to all time points)
for i = 1:length(MySamples)
    % load in data
    if ~isnan(CurFileNrs(i))
        if CurFileNrs(i) < 10
            SampleFile  =   [FileLocation myDate '_Agilent_I\' GeneralFolderName '0' num2str(CurFileNrs(i)) '.D\Report.txt']; % changed this so it may be controlled from the main script (12.02.2021)
        else
            SampleFile  =   [FileLocation myDate '_Agilent_I\' GeneralFolderName num2str(CurFileNrs(i)) '.D\Report.txt'];
        end
        warning('off','all')
        RawData     =   importfile_Agilent_I(SampleFile);
        warning('on','all')
        
        % delete all rows containing NaN values to get peak tables for DAD and
        % RI
        RawData(isnan(RawData{:,1}),:)    =   [];
        % check where "1" appears for the second time in the first column,
        % indicating the start of the RI peak table
        idx                 =   find(RawData{:,1}==1,1,'last');
        PeakTable           =   RawData(idx:end,:);
        SampleRetTimes      =   PeakTable{:,2};
        SamplePeakAreas     =   str2double(PeakTable{:,6});
        
        % get retention times and areas of all peaks
%         SampleRetTimes  =   table2array(PeakTable(:,2));
%         SamplePeakAreas =   table2array(PeakTable(:,5));
        
        % loop through all standards and check RetTimes
        for j = 1:length(StdRetTimes)
            RetTimeDiff     =   SampleRetTimes-StdRetTimes(j);
            PotPeakIdx      =   find(RetTimeDiff >= -0.15 & RetTimeDiff <= 0.15);
            if ~isempty(PotPeakIdx)
                Concentration   = (SamplePeakAreas(PotPeakIdx) - ...
                    myStandards{j,2}.StandardCurve(2,1)) / ...
                    myStandards{j,2}.StandardCurve(1,1);
                if length(Concentration) > 1
                    Concentration   =   max(Concentration);
                end
                Results(i,j)    = Concentration;                                % may run into trouble if multiple Peaks are found
            end
        end
    else
        Results(i,:)    = NaN;
    end
end

ResultTable     =   array2table(Results,'VariableNames',myStandards(:,1));

%% plot the results
SubPlotRows     =   ceil(sqrt(length(StdRetTimes)));
SubPlotCols     =   floor(sqrt(length(StdRetTimes)));
myAxisNames     =   {'t [h]','concentration [g/l]'};

myFig   =   figure(CurFileNrs(end));
for i = 1:length(StdRetTimes)
    ax       =   subplot(SubPlotRows,SubPlotCols,i);
    myLineStyles 	=   {'-x','--^',':o','-.d'};
%     P{i,j}   =  myplot(timePoints,Results(:,i),myAxisNames,...
%                 MySamples,[],myLineStyles,[],[],[myStandards{i,1}],ax,0);
    P        =  myplot(timePoints,Results(:,i),myAxisNames,...
                myStandards(i,1),[],myLineStyles,[],[],[myStandards{i,1}],ax,0);
    ylim([0 Inf])
    sgtitle(['Analysis of Growth on ' SampleName],'FontName','Arial','FontSize',14,'Interpreter','none')
end
set(myFig,'Position',[20 50 1200 600])
% save the figure(s)
savefig(myFig,[FileLocation myDate '_Agilent_I\Condition_' SampleName '.fig'])
exportgraphics(myFig,[FileLocation myDate '_Agilent_I\Condition_' SampleName '.emf']);
