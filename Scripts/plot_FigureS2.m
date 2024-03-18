%plot Figure S2
%% plot Figure S2- plot the sensitivity analysis
clearvars;
clc;

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

    ind_t_diag(icount) = find(data_Ke.('Index')(row_ID{icount})==1,1,'first');
    t_diag(icount) = data_Ke.('Time')(row_ID{icount}(ind_t_diag(icount)));

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

load('sol')

i_ind = 1;

%determine individual-specific parameter from estimations for our model
pB = sol{i_ind}.P(1);
pV = sol{i_ind}.P(2);
dB = sol{i_ind}.P(3);
sigma = sol{i_ind}.P(4);

P_all = [dN*0.9:dN*0.01:dN*1.1;
    b0*0.9:b0*0.01:b0*1.1;
    dI*0.9:dI*0.01:dI*1.1;
    pV*0.9:pV*0.01:pV*1.1;
    dV*0.9:dV*0.01:dV*1.1;
    pB*0.9:pB*0.01:pB*1.1;
    dB*0.9:dB*0.01:dB*1.1];

% % figure('visible','off');
figure
% %compute the fits for our model and Ke et al. refractory  model and mean
% %sqaured errors (MSEs)
% 
% %specify the layout of the figure - only plot specific fits and according
% %to random or specific order
t = tiledlayout(1,7,'TileSpacing','Compact','Padding','Compact');
% tile = 1:7;

icount = 1;
for ID_opt = 1:7
    %for ID_opt = 1

    nexttile
    % nexttile(tile(icount))
    % icount = icount+1;

    %determine time span for evaluation
    tspan = 0:1:20; %short term
    % tspan = 0:0.1:90; %short term

    for i_par = 1:size(P_all,2)

        %i_par

        dN = P_all(1,11);
        b0 = P_all(2,11);
        dI = P_all(3,11);
        pV = P_all(4,11);
        dV = P_all(5,11);
        pB = P_all(6,11);
        dB = P_all(7,11);

        if ID_opt == 1
            dN = P_all(ID_opt,i_par);
        elseif ID_opt == 2
            b0 = P_all(ID_opt,i_par);
        elseif ID_opt == 3
            dI = P_all(ID_opt,i_par);
        elseif ID_opt == 4
            pV = P_all(ID_opt,i_par);
        elseif ID_opt == 5
            dV = P_all(ID_opt,i_par);
        elseif ID_opt == 6
            pB = P_all(ID_opt,i_par);
        elseif ID_opt == 7
            dB = P_all(ID_opt,i_par);
        end

        sigma = sol{1}.P(4);

        y0 = [S0, 1, 0, 0];
        B_thres = 1-dI*dV/(b0*S0*(pV-dI));
        options = odeset('NonNegative',[1,2,3,4]); %specify non-negative values
        [t,y] = ode45(@(t,y) odefcn_SARSCoV2_infection(t,y,b0,dI,pV,dV,pN,dN,pB,dB,B_thres), tspan, y0, options);

        %for each day retrieve the measured CN value (if measured at all)
        y_short = [];
        ind = [];

        %map the simulated time points to the observed time points of data
        for i_time = 1:length(time_ID{i_ind})
            ind = find(tspan == time_ID{i_ind}(i_time));
            y_short(i_time) = y(ind,3);
        end

        %if values too small, fix at 1 (otherwise numerical problems)
        y_short((y_short<1))=1;
        %y_new = log10(y_short);

        %calculate CN/CT values to compare viral load with data given the
        %conversion by Ke 2022
        y_new = -(log10(y_short)-11.35)/(-0.25); %Ke 2022

        %number of data points for this individual
        n = length(data_ID{i_ind});

        %determine value of negative logL
        logL = 0.5*n*log(2*pi*sigma^2)+0.5*sum((data_ID{i_ind}-y_new').^2./sigma^2);

        logL_all(ID_opt,i_par) = logL;

    end

    plot(P_all(ID_opt,:),-logL_all(ID_opt,:),'-','LineWidth',1,'Color','k')
    hold on
    plot(P_all(ID_opt,11),-logL_all(ID_opt,11),'*','Color','k')

    % %specify axes and plot
    % if ismember(ID_opt, 1:12)
    %     xlim([0,20])
    ylim([-42,-36])
    %     %ylim([0,8])
    %     if ismember(ID_opt, [1,2,3,4,5,6,7,8])
    %         xticks([0,10,20])
    %         xticklabels({})
    %     else
    %         xticks([0,10,20])
    %     end
    %     if ismember(ID_opt,[1,5,9])
    %         yticks([-40,-30,-20,-10])
    %         yticklabels({40,30,20,10})
    %         %yticks([0,2,4,6,8])
    %     else
    %         yticks([-40,-30,-20,-10])
    %         yticklabels({})
    %     end
    % 
    % end

end

clear figure_property;
figure_property.units = 'inches';
figure_property.format = 'pdf';
figure_property.Preview= 'none';
figure_property.Width= '8'; % Figure width on canvas
figure_property.Height= '2'; % Figure height on canvas
figure_property.Units= 'inches';
figure_property.Color= 'rgb';
figure_property.Background= 'w';
figure_property.FixedfontSize= '10';
figure_property.ScaledfontSize= 'auto';
figure_property.FontMode= 'scaled';
figure_property.FontSizeMin= '10';
figure_property.FixedLineWidth= '1';
figure_property.ScaledLineWidth= 'auto';
figure_property.LineMode= 'none';
figure_property.LineWidthMin= '1';
figure_property.FontName= 'Arial';% Might want to change this to something that is available
figure_property.FontWeight= 'auto';
figure_property.FontAngle= 'auto';
figure_property.FontEncoding= 'latin1';
figure_property.PSLevel= '3';
figure_property.Renderer= 'painters';
figure_property.Resolution= '600';
figure_property.LineStyleMap= 'none';
figure_property.ApplyStyle= '0';
figure_property.Bounds= 'tight';
figure_property.LockAxes= 'off';
figure_property.LockAxesTicks= 'off';
figure_property.ShowUI= 'off';
figure_property.SeparateText= 'off';
chosen_figure=gcf;
set(chosen_figure,'PaperUnits','inches');
set(chosen_figure,'PaperPositionMode','auto');
set(chosen_figure,'PaperSize',[str2num(figure_property.Width) str2num(figure_property.Height)]); % Canvas Size
set(chosen_figure,'Units','inches');
hgexport(gcf,'.\Figures\FigureS2_raw.pdf',figure_property); %Set desired file name


