%% Calculating a cell dry weight conversion factor
% Author: Manuel Bruch, PhD Student
% Created: 2021/03/12
% last edited: 2021/03/21

%% Tidy up
clear all
close all
clc

%% import data
FileName    =   '..\..\Laboratory_files\C_necator_continuous.xlsx';
mySheet     =   'CDW_FA_est_03_06_22';                                          % Excel sheet in which data are found
myRange     =   {'B8:D8','B10:D10','B11:D11'};                                  % cells in which data are found

finalODs    =   readmatrix(FileName,'Sheet',mySheet,'Range',myRange{1});        
Volume      =   readmatrix(FileName,'Sheet',mySheet,'Range',myRange{2})./1000;  % conversion of sample volume from ml to L
mPellet     =   readmatrix(FileName,'Sheet',mySheet,'Range',myRange{3});        % in g

%% calculate the factor

factors     =   mPellet./(Volume.*finalODs);                                    % factor will have the unit g/L/OD --> OD measurements can be multiplied with this factor to yield g/L
factor      =   mean(factors);
stdDev      =   std(factors);

%% save the results
save('..\..\Laboratory_files\FA_CDWconversionFactor.mat','factor','stdDev')