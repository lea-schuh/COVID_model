%data fitting to Ke et al. 2022
clearvars;
clc;

%do plot of short-and long-term fit | 1 - yes | 0 - no
do_plot = 1;

%% prep data and retrieve fixed parameters
%get data from Ke et al. 2022
data_Ke = readtable('Data_Ke2022.xlsx');

%get all patient IDs from data table
ID = unique(data_Ke.('Ind'));

%create data structure - get specific patient information (time of relative
%measurement and CN values per individual)
icount = 1;
for i_ID = ID'

    %find all rows corresponding to patient ID
    row_ID{icount} = find(data_Ke.('Ind') == i_ID);

    %exclude NaN values of patient data
    a = data_Ke.('Nasal_CN')(row_ID{icount});
    row_ID1{icount} = row_ID{icount}(~isnan(a));

    %ensure that we only look at time points >=0 (after shifting data
    %by +6 days)
    b = data_Ke.('Time')(row_ID1{icount})+6;
    row_ID2{icount} = row_ID1{icount}(b>=0);

    %summarize data in data_ID structure
    data_ID{icount} = -data_Ke.('Nasal_CN')(row_ID2{icount});

    %set the CN values to -42 (detection threshold) if lower than -42
    data_ID{icount}(data_ID{icount}<-42)=-42;

    %shift data by +6 days (in raw data, peak viral load is centered at 0)
    time_ID{icount} = data_Ke.('Time')(row_ID2{icount})+6;
    icount = icount+1;

end

%determine fixed parameter values
S0 = 8*10^7; %total number of epithelial cells in nose at t=0, Ke et al., 2022
dN = 1/11; %death rate of all target cells, Tomasetti et al., 2017
pN = S0*dN; %production of new epithelial cells
b0 = 4.92*10^(-9); %infectivity rate, Ke et al., 2022
dI = 2.45; %death of infected cells, Ke et al., 2022
dV = 10; %deactivation virus, Ke et al., 2022

%% fitting of model to single-individual CN values from Ke et al.
%fit over every patient seperately
for ID_opt = 1:length(ID)
%for ID_opt = 1

    clearvars -except data_Ke ID ID_opt sol data_ID time_ID S0 dN pN b0 dI dV do_plot

    %get current individual
    ID_opt

    %number of single optimzation runs
    n_start = 10;

    %determine lower (lb) and upper (ub) parameter boundaries
    lb = [-10,2,-4,-2]; %pB, pV, dB, sigma
    ub = [-4,3,0,1];

    %use latin hypercube sampling to sample n_start initial values
    par0_all = lhsdesign_modified(n_start,lb,ub); %latinhypercube sampled startign parameters

    %optimize for each n_start
    for i = 1:n_start
        [par_opt(i,:),fval(i)] = fmincon(@(par)logLikelihood(data_ID{ID_opt},time_ID{ID_opt},par),par0_all(i,:),[],[],[],[],lb,ub);
    end

    %sort all n_starts according to best runs (highest fval)
    [fval_sorted,ind] = sort(fval);
    sol{ID_opt}.fval_sorted = fval_sorted;

    %set free parameters to best estimates from optimzation and
    %re-trasnform from log10 parameter space
    pB = 10^par_opt(ind(1),1);
    pV = 10^par_opt(ind(1),2);
    dB = 10^par_opt(ind(1),3);

    %collect the final optmimized parameter values
    sol{ID_opt}.P = 10.^par_opt(ind(1),:); %best optimized parameters according to logL
    sol{ID_opt}.Pall = 10.^par_opt(ind,:); %optimized parameters of all 10 runs

    %save sol structure
    % save('.\Results\sol','sol');

    if do_plot == 1

        %plot the results for short- and long-term dynamics
        for j = 1:2

            if j == 1
                tspan = 0:0.1:20; %short term
            else
                tspan = 0:1:360; %long term
            end

            %determine the initial values and B_thres for simulation
            y0 = [S0, 1, 0, 0];
            B_thres = 1-dI*dV/(b0*S0*(pV-dI));
            options = odeset('NonNegative',[1,2,3,4]); %specify non-negative values
            [t,y] = ode45(@(t,y) odefcn_SARSCoV2_infection(t,y,b0,dI,pV,dV,pN,dN,pB,dB,B_thres), tspan, y0,options);

            figure
            subplot(size(y,2)+1,1,1)
            plot(t,y(:,1),'Color','k','Linewidth',2) %S

            subplot(size(y,2)+1,1,2)
            plot(t,y(:,2),'Color','k','Linewidth',2) %I

            subplot(size(y,2)+1,1,3)
            plot(t,y(:,3),'Color','k','Linewidth',2) %V

            subplot(size(y,2)+1,1,4)
            if j == 1
                plot(time_ID{ID_opt}, data_ID{ID_opt}, '.', 'Color', 'k','Markersize',10) %data
                hold on
            end
            y_short = y(:,3);
            %if values too small, fix at 1 (numerical problems)
            y_short(y_short<1)=1;
            plot(t,-(log10(y_short)-11.35)/(-0.25),'Color','k','Linewidth',2) %Ke 2022 nasal

            subplot(size(y,2)+1,1,5)
            plot(t,y(:,4), '-', 'Color', 'k','Linewidth',2) %B
            hold on
            plot(t,1-(dI*dV/(b0*S0*(pV-dI)))*ones(length(t),1),'--','Color', 'k','Linewidth',2) %B_thres
            ylim([0,1])
            hold on

            set(gcf, 'Position',  [100, 100, 250, 600])
        end
    end
end

%% get average number of viral particles produced by an infected cell across all individuals

load('sol')

j = 1;
for i_ind = 1:length(sol)
    if ismember(i_ind,[24,35,41,48]) %remove the 4 individuals as Ke et al.
    else
        pV(j) = sol{i_ind}.P(2);
        j =j+1;
    end
end

dI = 2.45;
mean(pV/dI)
std(pV/dI)

figure
boxplot(pV/2.45)

median(pV/dI)
prctile(pV/dI,25)
prctile(pV/dI,75)
