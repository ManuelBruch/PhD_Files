%% Analysis of HPLC Data from the Shimadzu HPLC
% Author: Manuel Bruch, PhD Student
% Created: 2021/06/17
% last edited: 2021/11/24

%% Tidy up
clear all
close all
clc

%% import sample descriptions
FileLocation    =   '..\..\Laboratory_files\HPLC_Files\';
GeneralFile     =   'HPLC_samples.xlsx';
relevantSheet   =   '20_05_2024';
FigureNames     =   {''};

% change the following four lines manually
myDate              =   '20240520';
GeneralFileName     =   '240520_MB0';
% relevantSheet   =   '09-11_12_2020';
% relevantRange   =   {'C3:O9'};
relevantRange       =   {'B3:E12'};

Samples             =   readcell([FileLocation GeneralFile],'Sheet',...
                        relevantSheet,'Range',relevantRange{1});
SampleNr            =   1;

%% analyse standards if desired
% Standards       =   find(contains(Samples(:,end),'standard')); % does not work
% for some reason
stdIdx  =   [];
for i = 1:size(Samples,1)
    try
        if contains(Samples{i,end},'standard')
            stdIdx      =   [stdIdx, i];
        end
    end
end
StdSamples      =   Samples(stdIdx,:);

% assumption: 
% - samples are always run first
% - the last cell in a row contains the name of the standard
MySamplesIdx    =   setdiff(1:size(Samples,1),stdIdx);
numStdSamples   =   numel(StdSamples(:,1:end-1));

if isempty(MySamplesIdx)
    MySamples   =   [];
else
    MySamples   =   Samples(MySamplesIdx,:);
end
% determine smaple indexes for standard samples:
stdFileNr   =   reshape(1:numStdSamples, size(StdSamples(:,1:end-1),2), ...
                size(StdSamples(:,1:end-1),1)).';


%% extract time point containing columns
if ~isempty(MySamples)
    % to circumvent troubles with missing cells change the last row
    % temporarily
    fristRow    =   MySamples(1,:);
    fristRow(cellfun(@(x) isa(x,'missing'), fristRow))    =   {'I am empty'};
    timeCols    =   contains(fristRow,'=');       % Assumes that samples end in the last row
    timeColsPos =   find(timeCols);
    times       =   MySamples(:,timeCols);
    MySamples   =   MySamples(:,~timeCols);
    timePoints  =   zeros(size(MySamples,1),length(timeColsPos));
    for j = 1:size(timePoints,2)
        for i = 1:size(timePoints,1)
            eqPosition      =   strfind(times{i,j},'=');
            TimeStr         =   times{i,j}(eqPosition+2:end);
            timePoints(i,j) =   str2double(TimeStr);
        end
    end
    if isempty(timePoints)  % assumption: if no time is given, a concentration is listed
        ConcCols    =   contains(fristRow,'mM');
        ConcColsPos =   find(ConcCols);
        if ~isempty(ConcColsPos)
            Concs       =   MySamples(:,ConcCols);
            MySamples   =   MySamples(:,~ConcCols);
            
            myConcs     =   zeros(size(MySamples,1),length(ConcColsPos));
            for j = 1:size(myConcs,2)
                for i = 1:size(myConcs,1)
                    myConcs(i,j)    =   sscanf(Concs{i,j},'%f');
                end
            end
        else
            % some arbitrary denomnation can be found in the last column
            runName     =   MySamples(:,end);
            MySamples   =   MySamples(:,1:end-1);
        end
    end
end

%% find empty fields in the sample table and create an array of sample numbers
if ~isempty(MySamples)
    mask        =   cellfun(@ismissing,MySamples(:,:),'UniformOutput',false);
    FileNr      =   zeros(size(mask));
    counter     =   numStdSamples + 1;
    for i = 1:size(mask,1)
        for j = 1:size(mask,2)
            if ~mask{i,j}
                FileNr(i,j) =   counter;
                counter     =   counter+1;
            else
                FileNr(i,j) =   nan;
            end
        end
    end
    SampleFileNrs   =   FileNr(MySamplesIdx - size(StdSamples,1),:);
else
    SampleFileNrs   =   [];
end

%% Standard Analysis
disp('Would you like to analyse any standards?')
disp('Yes: type 1')
disp('No: type 0')
StandardAnalysis    =   input('');
allRetTimes     =   zeros(length(stdIdx),1);

% Analyse all standards
if StandardAnalysis
    for i = 1:length(stdIdx)
        Standard    =   Samples(stdIdx(i),:);
        relFileNrs  =   stdFileNr(stdIdx(i),:);
        
        % determine retention time and variance
%         disp('Would you like to enter the target retention time manually or shall I try to figure it out myself?')
%         disp('Yes: type 1')
%         disp('No: type 0')
%         a   =   input('');
%         
%         if a
%             RetTime     =   input('Please enter your desired retention time [min]: ');
%         else
%             RetTime     =   [];
%         end
        RetTime     =   [];
        
        [SampleNr,PeakAreas,Concentrations,StandardCurve,PeakTable,RetTime,myFig,Rsq]   =   ...
            HPLC_Standard_Analysis_Shimadzu(FileLocation,myDate,Standard,RetTime,i,StdSamples{i,end},GeneralFileName,relFileNrs);
        
        allRetTimes(i)  =   RetTime;
        
        save([FileLocation '\Standard_Parameters\' myDate '_' Standard{end} '.mat'],...
            'StandardCurve','RetTime','Concentrations','PeakAreas','Rsq')
        set(myFig,'Position',[20 50 1200 600])
        % save the figure(s)
        savefig(myFig,['..\..\Matlab_figures\Original_Figures\HPLC_standards\' myDate '_' Standard{end} '.fig'])
        exportgraphics(myFig,['..\..\Matlab_figures\EMF_Files\HPLC_standards\' myDate '_' Standard{end} '.emf'],'BackgroundColor','w');
    end
end

%% Prepare sample analysis
% import standard data
if ~isempty(SampleFileNrs)
    disp('Are all standards run at the same date as the samples?')
    disp('Yes: type 1')
    disp('No: type 0')
    a   =   input('');
    if ~a
        myDate2  =  input('Please enter the date on which the standards were run (YYYYMMDD): ','s');
    else
        myDate2  =  myDate;
    end
    
    % check which data are available
    myFiles     =   what([FileLocation 'Standard_Parameters\']);
    FileNames   =   myFiles.mat;
    % have to select all filenames belonging to the set date
    FileNames   =   FileNames(contains(FileNames,myDate2));
    
    %%% ONLY FOR MY CASE (don't want to stop each time for the prompt)
%     FileNames   =   FileNames(~contains(FileNames,'sodium_formate'));
    
    myStandards =   cell(size(FileNames,1),2);
    StdRetTimes =   zeros(size(FileNames));
    
    for i = 1:length(FileNames)
        myStandards{i,1}    =   FileNames{i}(10:end-4);
        myStandards{i,2}    =   load([FileLocation 'Standard_Parameters\' FileNames{i}]);
        StdRetTimes(i)      =   myStandards{i,2}.RetTime;
    end
    
    %% adjust Standards to be analysed
    if any(contains(FileNames,'avg'))
        disp('The chosen range contains averages for existing standards.')
        disp('Would you like to only load in these averages?')
        disp('Yes: type 1')
        disp('No: type 0')
        a   =   input('');
        if a
            averagedStds    =   contains(FileNames,'avg');
            myStandards     =   myStandards(averagedStds,:);
            StdRetTimes     =   StdRetTimes(averagedStds);
        end
    end
    allStandards 	=   myStandards;
    allRetTimes     =   StdRetTimes;
    % myStandards     =   myStandards([1,4,5,6],:);
    % StdRetTimes     =   StdRetTimes([1,4,5,6]);
    %
    % myStandards{contains(myStandards(:,1),'formic'),1}  =   'formic_acid';
    % myStandards(:,1)  =   strrep(myStandards(:,1),'_',' ');
    
    disp('You have successfully loaded in the following standards:')
    disp(myStandards(:,1))
    
    %% Sample Analysis
    % find all unique strings in the MySamples array:
    groupNames  =   unique(MySamples(cellfun('isclass',MySamples,'char')));
    for i = 1:numel(MySamples)
        if ismissing(MySamples{i})
            MySamples{i} = 'IamNotASample';
        end
    end
    
    for i = 1:length(groupNames)
        SampleName  =   groupNames(i);
        samplePos   =   strfind(MySamples, groupNames{i});
        [row, col]  =   find(cellfun(@(x)isequal(x,1),samplePos));
        CurFileNrs  =   zeros(length(row),1);
        for j = 1:length(row)
            CurFileNrs(j)   =   SampleFileNrs(row(j),col(j));
        end
        if ~isempty(timePoints)
            CurTime     =   timePoints(row,i);    % for the plots
            CurTimeStr  =   times(row,i);         % for naming the rows of the resulting table
            NotATime    =   find(isnan(CurTime));
            CurTime(NotATime)   =   [];
            CurTimeStr(NotATime)=   [];
            %         CurFileNrs          =   FileNr(MySamplesIdx,i);
            [Results,SampleName]=   HPLC_Sample_Analysis_Shimadzu(FileLocation,myDate,CurFileNrs,SampleName,StdRetTimes,myStandards,size(MySamples,2),CurTime,GeneralFileName);
            SampleNr    =   SampleNr+1;
            Results.Properties.RowNames     =   CurTimeStr;
            writeFile   =   [FileLocation myDate '\Results\' SampleName '.txt'];
            writetable(Results,writeFile,'WriteRowNames',true,'Delimiter','\t')
        end
        if exist('myConcs','var') && ~exist('runName','var')
            CurFileNrs          =   FileNr(MySamplesIdx,i);
            [Results,SampleName]=   HPLC_Sample_Analysis_Shimadzu_noT(FileLocation,myDate,CurFileNrs,SampleName,StdRetTimes,myStandards,size(MySamples,2),myConcs,GeneralFileName);
            SampleNr    =   SampleNr+1;
            Results.Properties.RowNames     =   Concs;
            writeFile   =   [FileLocation myDate '\Results\' SampleName '.txt'];
            writetable(Results,writeFile,'WriteRowNames',true,'Delimiter','\t')
        end
        if exist('runName','var')
            [Results,SampleName]=   HPLC_Sample_Analysis_Shimadzu_noT(FileLocation,myDate,CurFileNrs,SampleName,StdRetTimes,myStandards,size(MySamples,2),runName,GeneralFileName);
            SampleNr    =   SampleNr+1;
            Results.Properties.RowNames     =   runName;
            writeFile   =   [FileLocation myDate '\Results\' SampleName '.txt'];
            writetable(Results,writeFile,'WriteRowNames',true,'Delimiter','\t')
        end
    end
end
