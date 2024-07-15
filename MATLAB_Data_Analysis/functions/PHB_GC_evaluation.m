function outTable = PHB_GC_evaluation(inTable)
%PHB_GC_evaluation This function takes the offline parameters of a
%cultivation in which PHB was produced and calculates C/N & C/P ratios, PHB
%percentage, titre and STY and returns an extended table with these
%parameters.
%   required table headers in the input table are:
%   cultivation_type
%   c_C [mmol/L]
%   c_N [mmol/L]
%   c_P [mmol/L]
%   D_h-1_conti
%   sampling_time_h_batch
%   V_harvested_ml
%   m_pellet_mg
%   m_pellet_used_mg
%   peak_area_IS
%   peak_area_C4


% PHB quantification equation parameters
methanolysis_conv_factor    =    0.85;
betaIS                      =    7;     % beta of the internal standard (methyl benzoate) (beta = (number of C atoms) - 0.5*(number of O atoms)) --> C8 O2
betaC4                      =    3.5;   % beta of the C4 monomer (Methyl 3-hydroxybutyric acid) --> C5 O3
molarConcIS                 =    0.000735;  % molar concentration of the internal standard (0.1 mg/ml, MW = 136.1 g/mol)
dehydrogenatedC4molarMass   =   86.09;  % molar mass of 3-hydroxybutyric acid minus H2O (lost in polymerisation to PHB)
sampleVolume                =    2;     % sample volume

%% calculate additional sample stats
% C/N ratio
fullData    =   [inTable array2table(inTable.c_C_mmol_L_ ./ inTable.c_N_mmol_L_,...
                 'VariableNames',{'C_N_ratio'})];
% C/P ratio
fullData    =   [fullData array2table(fullData.c_C_mmol_L_ ./ fullData.c_P_mmol_L_,...
                 'VariableNames',{'C_P_ratio'})];

% CDW [g/L]
fullData    =   [fullData array2table(fullData.m_pellet_mg ./ fullData.V_harvested_ml,...
                 'VariableNames',{'CDW_g_L'})];

% calculate PHB content based on formula created in the lab
monomerMass =   (fullData.peak_area_C4 ./ fullData.peak_area_IS) * ...
                (betaIS/betaC4) * ...
                molarConcIS * dehydrogenatedC4molarMass * sampleVolume;
fullData    =   [fullData array2table(monomerMass,...
                 'VariableNames',{'monomer_mass_mg'})];

PHBpercent  =   (fullData.monomer_mass_mg ./ fullData.m_pellet_used_mg) * ...
                100 / methanolysis_conv_factor;
fullData    =   [fullData array2table(PHBpercent,...
                 'VariableNames',{'PHB_percent'})];

% PHB titre
fullData    =   [fullData array2table(fullData.CDW_g_L .* (fullData.PHB_percent / 100),...
                 'VariableNames',{'c_PHB_g_L'})];

% PHB STY
PHB_STY     =   zeros(height(fullData), 1);
for i = 1:height(fullData)
    if contains(fullData.cultivation_type{i}, 'batch')
        PHB_STY(i)  =   fullData.c_PHB_g_L(i) / fullData.sampling_time_h_batch(i);
    elseif contains(fullData.cultivation_type{i}, 'conti')
        PHB_STY(i)  =   fullData.c_PHB_g_L(i) * fullData.D_h_1_conti(i);
    end
end
outTable    =   [fullData array2table(PHB_STY,...
                 'VariableNames',{'PHB_STY_g_L_h'})];
end