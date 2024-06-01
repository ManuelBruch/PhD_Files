function [SampleNr,PeakAreas,Concentrations,StandardCurve,PeakTable,RetTime,myFig,Rsq] = HPLC_Standard_Analysis_Agilent_I(FileLocation,myDate,Standard,RetTime,SampleNumber,SampleName,GeneralFolderName,FileNr)
%HPLC_Standard_Analysis Summary of this function goes here
%   Detailed explanation goes here

%% loop through all files belonging to one sample series
mask            =   cellfun(@ismissing,Standard,'UniformOutput',false);
Concentrations  =   cell2mat(Standard(~cell2mat(mask(1:end-1))));
PeakAreas       =   zeros(1,length(Concentrations));
SampleName      =   strrep(SampleName,'_standard','');
SampleName      =   strrep(SampleName,'_',' ');

for j = 1:length(Concentrations)
    if FileNr(j) < 10
        SampleFile  =   [FileLocation myDate '_Agilent_I\' GeneralFolderName '0' num2str(FileNr(j)) '.D\Report.txt'];
    else
        SampleFile  =   [FileLocation myDate '_Agilent_I\' GeneralFolderName num2str(FileNr(j)) '.D\Report.txt'];
    end
    warning('off','all')
    
    %%%%%%% copied code from auto-generated importfile function
%     %% Read columns of data as text:
%     % For more information, see the TEXTSCAN documentation.
%     formatSpec = '%4s%8s%3s%2s%8s%11s%11s%s%[^\n\r]';
%     startRow = 1;
%     endRow = inf;
%     %% Open the text file.
%     fileID = fopen(SampleFile,'r','n','UTF16-LE');
%     % Skip the BOM (Byte Order Mark).
%     fseek(fileID, 2, 'bof');
%     
%     %% Read columns of data according to the format.
%     % This call is based on the structure of the file used to generate this code. If an error occurs for a different file, try regenerating the code from the Import Tool.
%     dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
%     for block=2:length(startRow)
%         frewind(fileID);
%         dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
%         for col=1:length(dataArray)
%             dataArray{col} = [dataArray{col};dataArrayBlock{col}];
%         end
%     end
    %%%%%%%
    
    RawData     =   importfile_Agilent_I(SampleFile);
    warning('on','all')
    
    % delete all rows containing NaN values to get peak tables for DAD and
    % RI
    RawData(isnan(RawData{:,1}),:)    =   [];
    % check where "1" appears for the second time in the first column,
    % indicating the start of the RI peak table
    idx     =   find(RawData{:,1}==1,1,'last');
    PeakTable   =   RawData(idx:end,:);
    RetTimes    =   PeakTable{:,2};
    PeakAreasTemp   =   str2double(PeakTable{:,6});
    
    %find or set peak
    if j == 1
        if isempty(RetTime)
            RetTime     =   RetTimes(PeakAreasTemp==max(PeakAreasTemp));
        end
    end
    % find peak closest to retention time in first file
    C       =   unique(abs(RetTimes-RetTime),'stable');
    [~,b]   =   min(C);
    PeakAreas(j)    =   PeakAreasTemp(b);
end
SampleNr    =   FileNr(j);

% present results as table
myIDs   =   transpose(1:length(Concentrations));
ResultTable     = table(myIDs,transpose(Concentrations),transpose(PeakAreas),'VariableNames',{'ID','Concentration [g/l]','Peak Area [a.u.]'});
disp(SampleName)
disp(ResultTable)

% Ask if some concentrations shall be omitted and delete them from the
% array
% disp('Would you like to exclude certain concentrations and respective Peak Areas?')
% disp('Yes: type 1')
% disp('No: type 0')
% a   =   input('');
a = 0;
if a
   delC     =   input('Enter the IDs of excluded data as a vector: ');
   Concentrations(delC)     =   [];
   PeakAreas(delC)          =   [];
end

% Plot linear fit to data
myFig   =   figure(SampleNumber);
FS    	=   [11 11 11 11];
myTitle =   SampleName;

% calculate R^2
% as forced through 0 --> no y-intercept neccessary
% https://uk.mathworks.com/help/matlab/data_analysis/linear-regression.html

b1      =   Concentrations'\PeakAreas';
yCalc   =   Concentrations*b1;
Rsq     =   1 - sum((PeakAreas-yCalc).^2)/sum((PeakAreas-mean(PeakAreas)).^2);
disp(['R^2 = ' num2str(Rsq)])

%try forcing the curve through 0
Concentrations  =   [0 Concentrations];
PeakAreas       =   [0 PeakAreas];
xerr    =   zeros(1,length(Concentrations));
yerr    =   ones(1,length(PeakAreas));                                      % emulating excel behaviour with minimal error in Y data
% yerr    =   ones(1,length(PeakAreas)).*PeakAreas*0.01;                      % assuming an error of 1%
yerr(1) =   0;
[P,SP]  =   linfitxy(Concentrations,PeakAreas,xerr,yerr);                   % Copyright (c) 2017, Julien Browaeys AND Tristan Beau - All rights reserved.

set(gca,'FontName','Arial')
grid on
title(myTitle,'FontName','Arial','FontSize',FS(4))
xlabel('Concentration [g/l]','FontName','Arial','FontSize',FS(2))
ylabel('Peak Area [a. u.]','FontName','Arial','FontSize',FS(2))
xl  =   xlim;
yl  =   ylim;
mystr       =   {'PA = m*c + b',...
                ['m = ' num2str(P(1)) ' \pm ' num2str(SP(1))], ...
                ['b = ' num2str(P(2)) ' \pm ' num2str(SP(2))], ...
                ['R^2 = ' num2str(Rsq)]};
text(xl(1)+0.1,yl(2)-yl(2)*0.1,mystr,'FontSize',11,'FontName', 'Arial')

StandardCurve   = [P;SP];

% delete concentration 0
Concentrations(1)   =   [];
PeakAreas(1)        =   [];
end

