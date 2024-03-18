%plot all fits
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

%retrieve the individual-specific parameters from Ke et al.
par_Ke = readtable('Parameters_Ke2022.txt');

icount = 1;
for i_ID = ID'
    row_ID{icount} = find(par_Ke.('ID') == i_ID);
    t0_all(icount) = table2array(par_Ke(row_ID{icount},2));
    beta_all(icount) = table2array(par_Ke(row_ID{icount},3));
    delta_all(icount) = table2array(par_Ke(row_ID{icount},4));
    pi_all(icount) = table2array(par_Ke(row_ID{icount},5));
    Phi_all(icount) = table2array(par_Ke(row_ID{icount},6));
    roh_all(icount) = table2array(par_Ke(row_ID{icount},7));
    icount = icount+1;
end

%determine fixed parameter values
S0 = 8*10^7; %total number of epithelial cells in nose at t=0, Ke et al., 2022
dN = 1/11; %death rate of all target cells, Tomasetti et al., 2017
pN = S0*dN; %production of new epithelial cells
b0 = 4.92*10^(-9); %infectivity rate, Ke et al., 2022
dI = 2.45; %death of infected cells, Ke et al., 2022
dV = 10; %deactivation virus, Ke et al., 2022

load('sol');

%determine fixed parameter values for Ke et al. refractory model
%(saliva)
c  = 10;
k = 4;

% figure('visible','off');
figure
%compute the fits for our model and Ke et al. refractory  model and mean
%sqaured errors (MSEs)

%specify the layout of the figure - only plot specific fits and according
%to random or specific order
t = tiledlayout(10,6,'TileSpacing','Compact','Padding','Compact');
% tile = 1:12;

icount = 1;
for ID_opt = 1:length(sol)

    %determine time span for evaluation
    tspan = 0:0.1:20; %short term
    % tspan = 0:0.1:90; %short term

    %determine individual-specific parameters for Ke et al. refractory
    %model
    t0 = t0_all(ID_opt);
    beta = beta_all(ID_opt)*10^-8;
    delta = delta_all(ID_opt);
    pi = pi_all(ID_opt);
    Phi= Phi_all(ID_opt)*10^-6;
    roh = roh_all(ID_opt);

    %determine individual-specific parameter from estimations for our model
    pB = sol{ID_opt}.P(1);
    pV = sol{ID_opt}.P(2);
    dB = sol{ID_opt}.P(3);

    if ~isnan(t0)
        nexttile

        options = odeset('NonNegative',[1,2,3,4]); %specify non-negative values

        %get fits of Ke et al.
        if ~isnan(t0)
            y0 = [S0, 0, 1, 0, 0]; %specify the initial values
            [t,y_Ke] = ode45(@(t,y) odefcn_Ke_refractory(t,y,beta,delta,pi,Phi,roh,k,c), tspan, y0, options);
        end

        %get fits of our model
        y0 = [S0, 1, 0, 0]; %specify the initial values
        B_thres = 1-dI*dV/(b0*S0*(pV-dI)); %specify B_thres
        [t,y] = ode45(@(t,y) odefcn_SARSCoV2_infection(t,y,b0,dI,pV,dV,pN,dN,pB,dB,B_thres), tspan, y0,options);

        %plot detection threshold
        plot([0,20],[-42,-42],'--','Color', [0,0,0]./255,'LineWidth',1)
        hold on

        %get Ke et al. values of viral load
        y_short_Ke = y_Ke(:,5);
        %for numerics, limit lowest viral load to 1
        y_short_Ke(y_short_Ke<1)=1;
        %plot Ke et al. fit
        plot(tspan-(t0-t_diag(ID_opt))+6,-(log10(y_short_Ke)-11.35)/(-0.25),'Color',[150,150,150]./255,'Linewidth',1) %Ke 2022 nasal


        hold on

        %get values of viral load
        y_short = y(:,3);
        %for numerics, limit lowest viral load to 1
        y_short(y_short<1)=1;
        %plot fit
        plot(tspan,-(log10(y_short)-11.35)/(-0.25),'Color', [255,102,0]./255,'Linewidth',1) %Ke 2022 nasal
        hold on

        %plot CN values (fitted data)
        plot(time_ID{ID_opt}, data_ID{ID_opt}, 'o', 'Color', [0,0,0]./255,'Markersize',2,'MarkerFaceColor',[0,0,0]./255) %data
        hold on

        %specify axes and plot
        xlim([0,20])
        ylim([-45,-10])

        xticks([0,10,20])
        xticklabels({})

        yticks([-40,-30,-20,-10])
        yticklabels({})

        title(sprintf('#%d',ID(ID_opt)),'FontWeight','normal')
        %t.TickLength = [0.02 0.035];
        box off
    end
end

clear figure_property;
figure_property.units = 'inches';
figure_property.format = 'pdf';
figure_property.Preview= 'none';
figure_property.Width= '8'; % Figure width on canvas
figure_property.Height= '11'; % Figure height on canvas
figure_property.Units= 'inches';
figure_property.Color= 'rgb';
figure_property.Background= 'w';
figure_property.FixedfontSize= '8';
figure_property.ScaledfontSize= 'auto';
figure_property.FontMode= 'fixed';
figure_property.FontSizeMin= '8';
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
set(chosen_figure,'PaperSize',[8.27,11.69]); % Canvas Size
set(chosen_figure,'Units','inches');
hgexport(gcf,'.\Figures\FigureS1_raw.pdf',figure_property); %Set desired file name