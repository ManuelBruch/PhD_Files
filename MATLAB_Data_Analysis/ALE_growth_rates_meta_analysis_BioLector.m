%% BioLector Processed Data Analysis
% newer version for analysing rfp assays, based on pure triplicate runs
% this one also calculates growth rates
% Author: Manuel Bruch, PhD student
% Creation:     04.12.2021
% Last change:  04.12.2021

%% Tidy up
clear all
close all
clc

%% data import
myFolder    =   '20211124_ALE_test';
myDate2     =   '211129';
ResultFile  =   ['..\..\Laboratory_files\ALE_screening\',myFolder,'\',...
                myDate2 '_ALE_Results.xlsx'];
figureFolder=   ['..\..\Laboratory_files\ALE_screening\' myFolder '\Figures\'];
figureName  =   'biolector_growth_rate_analysis';
overWrite   =   1;

mu          =   readtable(ResultFile,'Sheet','mu_mean');
mu_error    =   readtable(ResultFile,'Sheet','mu_std');

%% extract mu max from each column
[M,I]       =   max(mu{:,2:end});
ALEmu       =   M(2:end);
% get error of µ(max, WT)
WTerror     =   mu_error{I(1),2};
WTmuUpper   =   M(1) + WTerror;
WTmuLower   =   M(1) - WTerror;

%% plot the mu max
x       =   [1:size(ALEmu,2)];
figure('Renderer', 'Painters')
myFig   =   figure(1);
l1      =   yline(M(1), '-k');
hold on
l2      =   yline(WTmuUpper,'--k');
hold on
l3      =   yline(WTmuLower,'--k', 'H16', 'LabelVerticalAlignment', 'bottom');
hold on
set(l3,'FontSize',8,'FontName','arial')
% l3      =   yline(WTmuUpper*2,'--r','2x WT_{max}');
% set(l2,'FontSize',11,'FontName','arial')
hold on
% box on
grid on
s       =   scatter(x,ALEmu,15,'filled');

ylim([0,ceil(max(ALEmu,[],'all')*100)/100])
ytickformat('%.2f')
yticks(linspace(0,ceil(max(ALEmu,[],'all')*100)/100,6))

ylabel({'\mu_{max} [h^{-1}]';''},'FontSize',10,'FontName','arial')
xlabel('ALE colony','FontSize',10,'FontName','arial')

% plot the five highest values in a different colour and with labels
[B, I]  =   maxk(ALEmu,5);
s2      =   scatter(I, B, 20, [0.7 0 0], 'filled');
dx      =   0.3;
dy      =   0;
% create labels
myLabels=   cell(size(B));
for i = 1:length(B)
    myLabels{i}     =   [' ALE',num2str(I(i))];
end
pltLbl  =   text(I + dx,B + dy,myLabels,'FontSize',8,'FontName','arial');

% plot a line that looks like it connects upper and lower error of the WT to
% make it look like an errorbar
xLimits     =   xlim;
xErrBar     =   xLimits(2)/2;
plot([xErrBar xErrBar],[WTmuLower, WTmuUpper],'-k')
ax = gca;
ax.FontSize = 9;
ax.FontName = 'Arial';


%% save the figure
set(myFig,'Position',[20 50 600 400])

figureFolder1   =   figureFolder;
%     figureFolder1   =   ['..\..\Matlab_figures\Original_Figures\BioLector\'....
%                         ,myDate];
if ~exist(figureFolder1, 'dir')
    mkdir(figureFolder1)
end
figureFolder2   =   figureFolder;
%     figureFolder2   =   ['..\..\Matlab_figures\EMF_Files\BioLector\' myDate];
if ~exist(figureFolder2, 'dir')
    mkdir(figureFolder2)
end

if ~isfile([figureFolder1 figureName '.fig']) || overWrite
    savefig(myFig,[figureFolder1, figureName, '.fig'])
    exportgraphics(myFig,[figureFolder2, figureName, '.emf'])
    print(myFig, [figureFolder2, figureName, '.svg'], '-dsvg')
end