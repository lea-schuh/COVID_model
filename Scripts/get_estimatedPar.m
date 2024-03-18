%get estimated parameters 

clearvars;
clc;

%get data from Ke et al. 2022
data_Ke = readtable('Data_Ke2022.xlsx');

%get all patient IDs from data table
ID = unique(data_Ke.('Ind'));

load('sol');

icount = 1;
for ID_opt = 1:length(sol)
    if ismember(ID_opt,[24,35,41,48])
    else
        %determine individual-specific parameter from estimations for our model
        pB(ID_opt) = sol{ID_opt}.P(1);
        pV(ID_opt) = sol{ID_opt}.P(2);
        dB(ID_opt) = sol{ID_opt}.P(3);
        sigma(ID_opt) = sol{ID_opt}.P(4);
        I(ID_opt) = ID(ID_opt);
    end
end
