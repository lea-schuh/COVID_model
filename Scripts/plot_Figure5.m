%% plot Figure 5 - plot reinfection dynamics
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
C = [{[255,127,42]./255}, {[255,102,0]./255}, {[0,0,0]./255}, {[33,68,120]./255}, {[255,204,0]./255}, {[113,55,200]./255}];
LW = 1;
FA = 0.1;

ID_opt = 1;

%determine a vertical tile numeration for plot
nt = 1:4;

% figure('visible','off');
figure
t_reinf = tiledlayout(1,4,'TileSpacing','Compact','Padding','Compact');

for i = 1:2
    tspan = 0:0.001:90; %long-term infection dynamics

    %specify the individual-specific parameters from fit
    pB = sol{ID_opt}.P(1);
    pV = sol{ID_opt}.P(2);
    dB = sol{ID_opt}.P(3);

    if ismember(i,[1,3])
        b02 = 2*b0;
    else
        b02 = b0;
    end

    if ismember(i,[2,3])
        CI = 0.5;
    else
        CI = 1;
    end

    %determine the initial values and B_thres for simulation
    y0_n = [S0, 1, 0, 0];
    B_thres_n = 1-dI*dV/(b0*S0*(pV-dI));
    options = odeset('NonNegative',[1,2,3,4]); %specify non-negative values
    [t,y] = ode45(@(t,y) odefcn_SARSCoV2_infection(t,y,b0,dI,pV,dV,pN,dN,pB,dB,B_thres_n), tspan, y0_n,options);

    y02_n = [y(end,1), 1, 0, y(end,4)*CI]; %determien the initial values, take S0 and B0 from long-term fit
    B_thres2 = 1-dI*dV/(b02*S0*(pV-dI));
    tspan2 = 0:0.001:60; %reinfection time span for simulations
    [t_reinf,y_reinf] = ode45(@(t,y) odefcn_SARSCoV2_infection(t,y,b02,dI,pV,dV,pN,dN,pB,dB,B_thres2), tspan2, y02_n,options);

    y_short = y(:,3);
    y_short(y_short<1)=1;

    y_short_reinf = y_reinf(:,3);
    y_short_reinf(y_short_reinf<1)=1;

    y(:,3) = y_short;
    y_reinf(:,3) = y_short_reinf;

    %for each of the 4 populations, S, I, V, and B
    for j_comp = 3:4

        %next panel
        nexttile

        %plot the line highlighting the introduction of the new virus
        %variant
        plot([max(tspan),max(tspan)],[0,1],':','Color',C{3},'Linewidth',LW)
        hold on

        %for S, I, and V
        if j_comp < 4
            if j_comp == 1 %for S, normalization according to S0
                plot(tspan,y(:,j_comp)./S0,'Color',C{1},'Linewidth',LW) %long-term dynamics
                hold on
                plot(max(tspan)+tspan2,y_reinf(:,j_comp)./S0,'Color',C{2},'Linewidth',LW) %reinfection
                hold on
            else %for I and V, normalization according to max I or max V
                plot(tspan,y(:,j_comp)./max(y(:,j_comp)),'Color',C{1},'Linewidth',LW) %long-term dynamics
                hold on
                plot(max(tspan)+tspan2,y_reinf(:,j_comp)./max(y(:,j_comp)),'Color',C{i+3},'Linewidth',LW) %reinfection
                hold on
            end
        else %for B
            plot(tspan,y(:,j_comp),'Color',C{1},'Linewidth',LW) %long-term dynamics
            hold on
            plot(max(tspan)+tspan2,y_reinf(:,j_comp),'Color',C{i+3},'Linewidth',LW) %reinfection
            hold on
        end

        %add B_thres to B plot
        if j_comp == 4
            plot([0,max(tspan)],[B_thres_n,B_thres_n],'--','Color',C{1},'LineWidth',LW)
            hold on
            plot([max(tspan),max(tspan)+max(tspan2)],[B_thres2,B_thres2],'--','Color',C{i+3},'LineWidth',LW)
            hold on
        end

        xlim([0,max(tspan)+max(tspan2)])
        if j_comp == 4
            xticks([0,90,150])
        else
            xticks([0,90,150])
            % xticklabels({})
            ylim([0,1])
        end
        if i == 1
            % if j_comp == 3
            %     yticks([0,1,2])
            % else
            yticks([0,0.5,1])
            % end

            if j_comp == 1
                ylabel('Norm. {\it S}')
            elseif j_comp  == 2
                ylabel('Norm. {\it I}')
            elseif j_comp == 3
                ylabel('Norm. {\it V}')
            else
                ylabel('{\it B}')
            end
        else
            % if j_comp ~= 3
            yticks([0,0.5,1])
            % yticklabels({})
            % else
            %     yticks([0,1,2])
            %     yticklabels({})
            % end
        end
        if i == 2
            if j_comp == 4
                xlabel('                                       Days since primary infection')
            end
        end
        % if j_comp == 1
        %      title(sprintf('#%d',ID(ID_opt)),'FontWeight','normal')
        % end
        box off
    end

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
hgexport(gcf,'.\Figures\Figure5_raw.pdf',figure_property); %Set desired file name
