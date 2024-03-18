function logL = logLikelihood(data,time,par)

%input:
%data - CN values per individual
%time - measurement time points, with peak viral load assumed to be at 6
%days post infection
%par - parameters to be optimized to receive best fit to Ke et al. data

%output
%logL - negative logLikelihood value between fit with estimated parameters
%and Ke et al. data for optimization

%specify the free parameters and re-transform from log10 space
pB = 10^par(1); %rate at which B increases
pV = 10^par(2); %viral production rate
dB = 10^par(3); %rate at which B approaches B_thres
sigma = 10^par(4); %Gaussian noise parameter

%determine fixed parameter values (same as defined in main_opt_Ke2022)
S0 = 8*10^7; %total number of epithelial cells in nose at t=0, Ke et al., 2022
dN = 1/11; %death rate of all target cells, Tomasetti et al., 2017
pN = S0*dN; %production of new epithelial cells
b0 = 4.92*10^(-9); %infectivity rate, Ke et al., 2022
dI = 2.45; %death of infected cells, Ke et al., 2022
dV = 10; %deactivation virus, Ke et al., 2022

%determine the initial values and B_thres for simulation
% if i == 1
y0 = [S0, 1, 0, 0]; %S, I, V, B - infection initiated with a single infected cell
% elseif i == 2
%     y0 = [S0, 10, 0, 0];
% elseif i == 3
%     y0 = [S0, 100, 0, 0];
% elseif i == 4
%     y0 = [S0, 1000, 0, 0];
% elseif i == 5
%     y0 = [S0, 10000, 0, 0];
% end

B_thres = 1-dI*dV/(b0*S0*(pV-dI));
tspan = [0 20]; %time span of solving ODE
options = odeset('NonNegative',[1,2,3,4]); %specify non-negative values
sol = ode45(@(t,y) odefcn_SARSCoV2_infection(t,y,b0,dI,pV,dV,pN,dN,pB,dB,B_thres), tspan, y0,options);

%evaluate the ODE at daily time points 0 to 20
y = deval(sol,0:1:20);
yT = y';

%for each day retrieve the measured CN value (if measured at all)
tspan = 0:1:20;
y_short = [];
ind = [];

%map the simulated time points to the observed time points of data
for i_time = 1:length(time)
    ind = find(tspan == time(i_time));
    y_short(i_time) = yT(ind,3);
end

%if values too small, fix at 1 (otherwise numerical problems)
y_short((y_short<1))=1;
%y_new = log10(y_short);

%calculate CN/CT values to compare viral load with data given the
%conversion by Ke 2022
y_new = -(log10(y_short)-11.35)/(-0.25); %Ke 2022
%data = 0.25*data+11.35;

%number of data points for this individual
n = length(data);

%determine value of negative logL
logL = 0.5*n*log(2*pi*sigma^2)+0.5*sum((data-y_new').^2./sigma^2);

end