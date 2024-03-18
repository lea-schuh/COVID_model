%% plot Figure 6 - reproduction number over immune response B and over course of reinfection
clearvars;
clc;

%determine fixed parameter values 
S0 = 8*10^7; %total number of epithelial cells in nose at t=0, Ke et al., 2022
dN = 1/11; %death rate of all target cells, Tomasetti et al., 2017
pN = S0*dN; %production of new epithelial cells
b0 = 4.92*10^(-9); %infectivity rate, Ke et al., 2022
dI = 2.45; %death of infected cells, Ke et al., 2022
dV = 10; %deactivation virus, Ke et al., 2022

%determine individual-specific parameter from estimations for our model
pV = 200;
pB = 10^(-8);
dB = 0.01;

B_thres = 1-dI*dV/(b0*S0*(pV-dI)); %specify B_thres

% figure('visible','off');
figure

subplot(1,2,1)
B = 0:0.01:1;
icount = 1;
for i = 0:0.01:1
    R0(icount) = (pV*b0*(1-i)*S0)/(dI*(dV+b0*(1-i)*S0));
    icount = icount+1;
end
plot(B,R0,'Color',[0,0,0]./255,'Linewidth',1)
hold on
yline(1,'--','Color',[0,0,0]./255,'Linewidth',1)
hold on
xline(B_thres,':','Color',[0,0,0]./255,'Linewidth',1)

%specify axes and plot
xlim([0,1])
ylim([0,Inf])

xticks([0,0.2,0.4,0.6,0.8,1])
% yticks([0,0.2,0.4,0.6,0.8,1])

box off
ylabel('R')
xlabel('B')

subplot(1,2,2)

%get fits of our model
y0 = [S0, 1, 0, 0]; %specify the initial values
options = odeset('NonNegative',[1,2,3,4]); %specify non-negative values
tspan = 0:0.1:20;
% tspan = 0:0.1:90;
[t,y] = ode45(@(t,y) odefcn_SARSCoV2_infection(t,y,b0,dI,pV,dV,pN,dN,pB,dB,B_thres), tspan, y0,options);

R = (pV*b0*(1-y(:,4))*S0)./(dI*(dV+b0*(1-y(:,4))*S0));

plot(tspan,R,'Color',[0,0,0]./255,'Linewidth',1)
hold on
yline(1,'--','Color',[0,0,0]./255,'Linewidth',1)

xlim([0,20])
ylim([0,Inf])

xticks([0,10,20])
% yticks([0,0.2,0.4,0.6,0.8,1])

box off
ylabel('R')
xlabel('Days since infection')

clear figure_property;
figure_property.units = 'inches';
figure_property.format = 'pdf';
figure_property.Preview= 'none';
figure_property.Width= '5'; % Figure width on canvas
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
hgexport(gcf,'.\Figures\Figure6_raw.pdf',figure_property); %Set desired file name
