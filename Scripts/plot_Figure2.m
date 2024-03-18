%plot Figure 2
%% plot Figure 2A - Ke et al. and our model fits vs data
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

%load('sol1');
load('sol')

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
t = tiledlayout(3,4,'TileSpacing','Compact','Padding','Compact');
tile = 1:12;

icount = 1;
for ID_opt = 1:length(sol)
%for ID_opt = 1

    %if individual should be plotted in Figure 2
    if ismember(ID_opt, 1:12)
        %create new figure panel
        nexttile(tile(icount))
        icount = icount+1;
    end

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
    if ismember(ID_opt, 1:12)
        plot([0,20],[-42,-42],'--','Color', [0,0,0]./255,'LineWidth',1) 
        hold on
    end

    if ~isnan(t0)
        %get Ke et al. values of viral load
        y_short_Ke = y_Ke(:,5);
        %for numerics, limit lowest viral load to 1
        y_short_Ke(y_short_Ke<1)=1;
        %plot Ke et al. fit
        if ismember(ID_opt, 1:12)
           plot(tspan-(t0-t_diag(ID_opt))+6,-(log10(y_short_Ke)-11.35)/(-0.25),'Color',[150,150,150]./255,'Linewidth',1) %Ke 2022 nasal
        end
    end
   
    hold on

    %get values of viral load
    y_short = y(:,3);
    %for numerics, limit lowest viral load to 1
    y_short(y_short<1)=1;
    %plot fit
    if ismember(ID_opt, 1:12)
        plot(tspan,-(log10(y_short)-11.35)/(-0.25),'Color', [255,102,0]./255,'Linewidth',1) %Ke 2022 nasal
        %plot(tspan,log10(y_short),'Color', [255,102,0]./255,'Linewidth',1) %Ke 2022 nasal
    end
    hold on

    %plot CN values (fitted data)
    if ismember(ID_opt, 1:12)
        plot(time_ID{ID_opt}, data_ID{ID_opt}, 'o', 'Color', [0,0,0]./255,'Markersize',2,'MarkerFaceColor',[0,0,0]./255) %data
        %plot(time_ID{ID_opt}, 0.25*data_ID{ID_opt}+11.35, 'o', 'Color', [0,0,0]./255,'Markersize',2,'MarkerFaceColor',[0,0,0]./255) %data
        hold on
    end

    %specify axes and plot
    if ismember(ID_opt, 1:12)
        xlim([0,20])
        ylim([-45,-10])
        %ylim([0,8])
        if ismember(ID_opt, [1,2,3,4,5,6,7,8])
            xticks([0,10,20])
            xticklabels({})
        else
            xticks([0,10,20])
        end
        if ismember(ID_opt,[1,5,9])
            yticks([-40,-30,-20,-10])
            yticklabels({40,30,20,10})
            %yticks([0,2,4,6,8])
        else
            yticks([-40,-30,-20,-10])
            yticklabels({})
        end
        title(sprintf('#%d',ID(ID_opt)),'FontWeight','normal')
        %t.TickLength = [0.02 0.035];
        box off
        if ID_opt == 5
            ylabel('Nasal CN values')
        end
        if ID_opt == 10
            xlabel('                                       Days since infection')
        end
    end

    %determine the MSEs
    if ~isnan(t0)
        y_short_Ke = [];
        y_short_Schuh = [];
        y_new_Ke = [];
        y_new_Schuh = [];
        ind1 = [];
        ind2 = [];

        %map the simulated time points to the observed time points of data
        for i_time = 1:length(time_ID{ID_opt})

            if -(t0-t_diag(ID_opt)) < time_ID{ID_opt}(i_time)-6
                ind1 = find(round(tspan-(t0-t_diag(ID_opt)),2) == time_ID{ID_opt}(i_time)-6);
                ind2 = find(round(tspan-6,2) == time_ID{ID_opt}(i_time)-6);
                y_short_Ke(i_time) = y_Ke(ind1,5);
                y_short_Schuh(i_time) = y(ind2,3);
            else
                data_ID{ID_opt}(i_time) = [];
                y_short_Ke(i_time) = NaN;
                y_short_Schuh(i_time) = NaN;
            end

        end

        %if values too small, fix at 1 (numerical problems)
        y_short_Ke((y_short_Ke<1))=1;
        y_short_Schuh((y_short_Schuh<1))=1;

        y_new_Ke = -(log10(y_short_Ke)-11.35)/(-0.25); %Ke 2022
        y_new_Schuh = -(log10(y_short_Schuh)-11.35)/(-0.25); %our model

        n = length(data_ID{ID_opt}); %number of data points

        %determine value of mean squared error for both fits
        MSE_Ke(ID_opt) = 1/length(data_ID{ID_opt})*sum((y_new_Ke(~isnan(y_new_Ke))'-data_ID{ID_opt}).^2);
        MSE_Schuh(ID_opt) = 1/length(data_ID{ID_opt})*sum((y_new_Schuh(~isnan(y_new_Schuh))'-data_ID{ID_opt}).^2);

    else
        
        %if Ke et al. did not fit this individual, set everything to zero
        MSE_Ke(ID_opt) = 0;
        MSE_Schuh(ID_opt) = 0;

    end

end
clear figure_property;
figure_property.units = 'inches';
figure_property.format = 'pdf';
figure_property.Preview= 'none';
figure_property.Width= '5'; % Figure width on canvas
figure_property.Height= '4'; % Figure height on canvas
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
hgexport(gcf,'.\Figures\Figure2A_raw.pdf',figure_property); %Set desired file name

%% Figure 2B and C- plot MSE comparisons

figure
% figure('visible','off');
t = tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');

nexttile
plot([0,max(MSE_Schuh)+10],[0,max(MSE_Schuh)+10],'--','Color','k','Linewidth',1)
hold on
scatter(MSE_Ke(MSE_Ke>0),MSE_Schuh(MSE_Schuh>0),'o','MarkerEdgeColor', 'none','SizeData',10,'MarkerFaceColor',[0,0,0]./255,'MarkerFaceAlpha',0.5)
% hold on
% scatter(MSE_Ke(1:12),MSE_Schuh(1:12),'o','MarkerEdgeColor', 'none','SizeData',10,'MarkerFaceColor','r','MarkerFaceAlpha',0.5)

xlim([0,50])
ylim([0,50])
xticks([0,10,20,30,40,50])
yticks([0,10,20,30,40,50])
xlabel('MSE_{Ke}')
ylabel('MSE_{model}')
box off

nexttile
%compute the log10 values of the MSEs for both fits
MSE_Ke_new_n = log10(MSE_Ke(MSE_Ke>0));
MSE_Schuh_new_n = log10(MSE_Schuh(MSE_Schuh>0));

plot([0,max(MSE_Schuh_new_n)+1],ones(2,1)*mean(MSE_Schuh_new_n-MSE_Ke_new_n),'--','Color','k','Linewidth',1)
hold on
plot([0,max(MSE_Schuh_new_n)+1],ones(2,1)*(mean(MSE_Schuh_new_n-MSE_Ke_new_n)+2*std(MSE_Schuh_new_n-MSE_Ke_new_n)),':','Color','k','Linewidth',1)
hold on
plot([0,max(MSE_Schuh_new_n)+1],ones(2,1)*(mean(MSE_Schuh_new_n-MSE_Ke_new_n)-2*std(MSE_Schuh_new_n-MSE_Ke_new_n)),':','Color','k','Linewidth',1)
hold on
scatter((MSE_Ke_new_n+MSE_Schuh_new_n)./2,MSE_Schuh_new_n-MSE_Ke_new_n,'o','MarkerEdgeColor', 'none','SizeData',10,'MarkerFaceColor',[0,0,0]./255,'MarkerFaceAlpha',0.5)
% hold on
% scatter((MSE_Ke_new_n(11)+MSE_Schuh_new_n(11))./2,MSE_Schuh_new_n(11)-MSE_Ke_new_n(11),'o','MarkerEdgeColor', 'none','SizeData',10,'MarkerFaceColor','r','MarkerFaceAlpha',0.5)
xlim([0,2])
ylim([-0.45,0.45])
xticks([0,1,2])
yticks([-0.4,-0.2,0,0.2,0.4])
xlabel('mean(log_{10} MSE_{Ke}, log_{10} MSE_{model})')
ylabel('MSE_{model} - MSE_{Ke}')
box off

upper_val = 10.^(mean(MSE_Schuh_new_n-MSE_Ke_new_n)+2*std(MSE_Schuh_new_n-MSE_Ke_new_n))
lower_val = 10.^(mean(MSE_Schuh_new_n-MSE_Ke_new_n)-2*std(MSE_Schuh_new_n-MSE_Ke_new_n))

se_bias = sqrt(std(MSE_Schuh_new_n-MSE_Ke_new_n)^2/length(MSE_Ke_new_n));

t = tinv(1-0.025,length(MSE_Ke_new_n)-1);

bias = mean(MSE_Schuh_new_n-MSE_Ke_new_n)

bias2std_plus = mean(MSE_Schuh_new_n-MSE_Ke_new_n)+2*std(MSE_Schuh_new_n-MSE_Ke_new_n)
bias2std_minus = mean(MSE_Schuh_new_n-MSE_Ke_new_n)-2*std(MSE_Schuh_new_n-MSE_Ke_new_n)

bias_lower = mean(MSE_Schuh_new_n-MSE_Ke_new_n) - (t*se_bias);
bias_upper = mean(MSE_Schuh_new_n-MSE_Ke_new_n) + (t*se_bias);

clear figure_property;
figure_property.units = 'inches';
figure_property.format = 'pdf';
figure_property.Preview= 'none';
figure_property.Width= '3'; % Figure width on canvas
figure_property.Height= '4'; % Figure height on canvas
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
hgexport(gcf,'.\Figures\Figure2BC_raw.pdf',figure_property); %Set desired file name

