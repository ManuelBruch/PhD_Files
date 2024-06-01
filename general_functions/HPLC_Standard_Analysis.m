function [SampleNr,PeakAreas,Concentrations,StandardCurve,PeakTable,RetTime] = HPLC_Standard_Analysis(FileLocation,myDate,SampleNr,Standard,RetTime,SampleNumber,SampleName)
%HPLC_Standard_Analysis Summary of this function goes here
%   Detailed explanation goes here

%% loop through all files belonging to one sample series
mask            =   cellfun(@ismissing,Standard,'UniformOutput',false);
Concentrations  =   cell2mat(Standard(~cell2mat(mask(1:end-1))));
PeakAreas       =   zeros(1,length(Concentrations));
SampleName      =   strrep(SampleName,'_standard','');
SampleName      =   strrep(SampleName,'_',' ');

for j = SampleNr:SampleNr+length(Concentrations)-1
    SampleFile  =   [FileLocation myDate '\full_reports\MANUELBRUCH09122020_1292020_' num2str(j) '.txt'];   % try to use the just_peak_tables folder and compare if it works
    warning('off','all')
    RawData     =   readcell(SampleFile,'Delimiter','\t');
    warning('on','all')
    
    % find peak table
    PKTstart    =   find(cellfun(@(c) ischar(c) && ~isempty(strfind(c, '[Peak Table(Detector A-Ch1)]')), RawData));
    PKTend      =   find(cellfun(@(c) ischar(c) && ~isempty(strfind(c, '[Peak Table(PDA-Ch1)]')), RawData));
%     PKTstart    =   find(contains(RawData(:,1),'[Peak Table(Detector A-Ch1)]'));
%     PKTend      =   find(contains(RawData(:,1),'[Peak Table(PDA-Ch1)]'));
    PeakTable   =   cell2table(RawData(PKTstart+3:PKTend-1,:),'VariableNames',RawData(PKTstart+2,:));
    
    %find or set peak
    if j == SampleNr
        if isempty(RetTime)
            RelevantTime    =   find(table2array(PeakTable(:,2))>2.1);      % manually set minimum Peak time to 2.1 min, due to a huge peak in my files before that
            PeakNumber      =   find(table2array(PeakTable(:,5))==max(table2array(PeakTable(RelevantTime,5))));
            RetTime         =   table2array(PeakTable(PeakNumber,2));
        end
    end
    % find peak closest to retention time in first file
    C       =   unique(abs(table2array(PeakTable(:,2))-RetTime),'stable');
    [~,b]   =   min(C);
    PeakAreas(j-SampleNr+1)    =   table2array(PeakTable(b,5));
end
SampleNr    =   SampleNr+length(Concentrations);

% present results as table
myIDs   =   transpose(1:length(Concentrations));
ResultTable     = table(myIDs,transpose(Concentrations),transpose(PeakAreas),'VariableNames',{'ID','Concentration [g/l]','Peak Area [a.u.]'});
disp(ResultTable)

% Ask if some concentrations shall be omitted and delete them from the
% array
disp('Would you like to exclude certain concentrations and respective Peak Areas?')
disp('Yes: type 1')
disp('No: type 0')
a   =   input('');

if a
   delC     =   input('Enter the IDs of excluded data as a vector: ');
   Concentrations(delC)     =   [];
   PeakAreas(delC)          =   [];
end

% Plot linear fit to data
figure(SampleNumber)
FS    	=   [11 11 11 11];
myTitle =   SampleName;

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
                ['b = ' num2str(P(2)) ' \pm ' num2str(SP(2))]} ...
                ['R^2 = ' num2str(Rsq)];
text(xl(1)+0.1,yl(2)-yl(2)*0.1,mystr,'FontSize',11,'FontName', 'Arial')

StandardCurve   = [P;SP];

% delete concentration 0
Concentrations(1)   =   [];
PeakAreas(1)        =   [];
end

