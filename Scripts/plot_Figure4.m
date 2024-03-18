%% plot Figure 4A - plot of fitted long-term dynamics

clearvars;
clc;

%get data (time and CN values per individual)
data_Ke = readtable('Data_Ke2022.xlsx');

%get all patient IDs
ID = unique(data_Ke.('Ind'));

%determine fixed parameter values
S0 = 8*10^7; %total number of epithelial cells in nose at t=0, Ke et al., 2022
dN = 1/11; %death rate of all target cells, Tomasetti et al., 2017
pN = S0*dN; %production of new epithelial cells
b0 = 4.92*10^(-9); %infectivity rate, Ke et al., 2022
dI = 2.45; %death of infected cells, Ke et al., 2022
dV = 10; %deactivation virus, Ke et al., 2022

load('sol');

tspan = 0:0.001:90; %long-term infection dynamics

%specify colors
C = [{[200,200,200]./255}, {[255,102,0]./255}, {[0,0,0]./255}, {[95,141,211]./255}, {[255,85,153]./255}, {[113,55,200]./255}];
LW = 1;
FA = 0.1;

j = 1;

for ID_opt = 1:length(sol)

    if ismember(ID_opt,[24,35,41,48]) %remove the 4 individuals as Ke et al.

    else
        %specify the individual-specific parameters from fit
        pB = sol{j}.P(1);
        pV = sol{j}.P(2);
        dB = sol{j}.P(3);


        %determine the initial values and B_thres for simulation
        y0_n = [S0, 1, 0, 0];
        B_thres_n = 1-dI*dV/(b0*S0*(pV-dI));
        options = odeset('NonNegative',[1,2,3,4]); %specify non-negative values
        [t,y] = ode45(@(t,y) odefcn_SARSCoV2_infection(t,y,b0,dI,pV,dV,pN,dN,pB,dB,B_thres_n), tspan, y0_n,options);

        y_short = y(:,3);
        y_short(y_short<1)=1;

        y(:,3) = y_short;

        Y1(:,j) = y(:,1);
        Y2(:,j) = y(:,2);
        Y3(:,j) = y(:,3);
        Y4(:,j) = y(:,4);

        j = j+1;
    end
end

figure
t_reinf = tiledlayout(1,4,'TileSpacing','Compact','Padding','Compact');
%nt = [1,6,11,16,2,7,12,17,3,8,13,18,4,9,14,19,5,10,15,20];
%nt = [1,5,9,13,2,6,10,14,3,7,11,15,4,8,12,16];


nexttile
for i_ind = 1:size(Y1,2)
    plot(tspan,Y1(:,i_ind)./S0,'Color',[C{2},FA],'Linewidth',LW) %long-term dynamics
    hold on
end
xlabel('Days since primary infection')
ylabel('{\it S}')
xticks([0,45,90])
ylim([0,1])
yticks([0,0.5,1])
xlim([0,90])
box off

nexttile
for i_ind = 1:size(Y2,2)
    plot(tspan,Y2(:,i_ind)./max(Y2(:,i_ind)),'Color',[C{2},FA],'Linewidth',LW) %long-term dynamics
    hold on
end
xlabel('Days since primary infection')
ylabel('{\it I}')
xticks([0,45,90])
ylim([0,1])
yticks([0,0.5,1])
xlim([0,90])
box off

nexttile
for i_ind = 1:size(Y3,2)
    plot(tspan,Y3(:,i_ind)./max(Y3(:,i_ind)),'Color',[C{2},FA],'Linewidth',LW) %long-term dynamics
    hold on
end
xlabel('Days since primary infection')
ylabel('{\it V}')
xticks([0,45,90])
ylim([0,1])
yticks([0,0.5,1])
xlim([0,90])
box off

nexttile
for i_ind = 1:size(Y4,2)
    plot(tspan,Y4(:,i_ind),'Color',[C{2},FA],'Linewidth',LW) %long-term dynamics
    hold on
    % plot([0,max(tspan)],[1-dI*dV/(b0*S0*(sol{i_ind}.P(2)-dI)),1-dI*dV/(b0*S0*(sol{i_ind}.P(2)-dI))],'--','Color',[C{3},FA],'LineWidth',LW)
    % hold on
end


xlabel('Days since primary infection')
ylabel('{\it B}')
xticks([0,45,90])
ylim([0,1])
yticks([0,0.5,1])
xlim([0,90])
box off

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
hgexport(gcf,'.\Figures\Figure4_raw.pdf',figure_property); %Set desired file name
