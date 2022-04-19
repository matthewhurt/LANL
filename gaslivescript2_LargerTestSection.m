%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% _____________ Gas Compression Code, vfinal 04.19.22 _______________%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Constants
% Median: AIR
Ma = [3, 5, 7, 10];
mDot_WT = [7267 12545 11752 7332];                  % [kg/s]
altitude = 3350;                                        % [m], 11000 [ft]
ro_alt = 0.86435;                                       % [kg/m3], 11000 [ft]
P_alt = 7.012*10^4;                                     % [Pa], 11000 [ft]
T_alt = 266.38;                                         % [K], 11000 [ft]
%altitude = 2225;                                        % [m], 7300 [ft]
%ro_alt = 0.98502;                                       % [kg/m3], 7300 [ft]
%P_alt = 7.739*10^4;                                     % [Pa], 7300 [ft]
%T_alt = 273.54;                                         % [K], 7300 [ft]

R = 8.314/0.02897;                                      % [J/(kg-K)]
gammaa = 1.4;
Cp = (gammaa*R)/(gammaa-1);                               % calorically perfect gas [J/(kg-K)]
Cv = Cp-R;
T_liq = 273.15-196;                                     % air liquefaction temperature
T_convert = (1+((gammaa-1)/2).*Ma);
P_convert = (1+((gammaa-1)/2).*Ma).^(-gammaa/(gammaa-1));
time = 30;                                              % [sec]
mass = mDot_WT.*time;
P_tank_max = 2.7579e7;                                 % [Pa]
T_tank_max = 300;                                       % K  atmospheric temperature in Los Alamos NM for storage tank calcs




%% State 4 - Wind Tunnel 

% initial
V_4i = 11.17;                               % volume of Nozzle
A_4i = 1.98*180;                            % cross Sectional Area: Nozzle Entrance
mDot_4i = mDot_WT;                          % mass flow rate: Nozzle Entrance
P_4i_tot = [128.56, 1587.27, 6209.75, 21219.62]*10^3;                     % total pressure: Nozzle Entrance
P_4i = P_4i_tot.*P_convert;
T_4i_tot = T_liq.*(1+(((gammaa-1)/2).*Ma.^2));  % temperature: Nozzle Entrance, stdAtm(Pressure)
T_4i_stat = T_liq;  % static temperature [K]
ro_4i = P_4i_tot./(R.*T_4i_tot);

% exit
M_4e = 0.95;                                % estimated Mach at WT exit
V_4e = 7.66;                                % volume of Diffuser
A_4e = 1.98^2;                              % cross Sectional Area: Diffuser Exit
mDot_4e = mDot_4i;                          % mass flow rate: Diffuser Exit
P_4e_tot = [42.21, 97.96, 95.33, 64.61]*10^3;                      % total pressure: Diffuser Exit
P_4e = P_4e_tot.*P_convert;
T_4e_tot = T_4i_tot;                        % temperature: Diffuser Exit, stdAtm(Pressure)
ro_4e = P_4e_tot./(R.*T_4e_tot);
T_4e_stat = T_4e_tot/(1 + (M_4e^2)*((gammaa-1)/2));  % static temperature [K]
h_4 = Cp*(T_4i_tot-T_4e_tot);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ____________________ AFTER WT, WORKING FORWARDS ___________________%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% State 5 - Heat Exchanger (After WT)
% Assumptions: Cp = const, Cv = const, ideal gas (P=ro*RT)
P_5 = P_4e_tot;
T_5e = T_alt*0.5;
T_5e_stat = T_5e;
Wdot_5 = 0;                     % Power generated/used is zero
mDot_5i = mDot_4e;
h_5i = T_4e_tot*Cp;
u_5i = T_4e_tot*Cv;
h_5e = T_5e*Cp;
ro_5 = P_5/(R*T_5e);
vel_5 = mDot_5i/(A_4e*ro_5);
h_5 = h_5e-h_5i;
Qdot_5 = -h_5.*mDot_5i;          % [J/s]




%% State 6 - Vacuum Chamber
% Assumptions: Cp = const, Cv = const, ideal gas (P=ro*RT)

Qdot_6 = 0;                                     % Heat Transfer Rate is zero
Wdot_6 = 0;
h_6i = h_5e;
P_6i = 1000;                                     % Initial pressure before run
P_6f = P_5;                                      % Final Pressure after run
ro_6 = ((P_6i./P_5).*(ro_5.^gammaa)).^(1/gammaa);
T_6 = P_6i./(R.*ro_6);
T_6stat = T_6;
h_6e = Cp*T_6;
h_6 = h_6e-h_6i;
Vol_6 = mDot_WT * time *R*T_5e./(P_6f-P_6i);
fprintf('The volume of the Vacuum Chamber is : %d meters^3.\n', max(Vol_6));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% __________________ BEFORE WT, WORKING BACKWARDS ___________________%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% State 1 - Compressor
% Assumptions: Cp = const, Cv = const, ideal gas (P=ro*RT)

Qdot_1 = 0;                     % Power generated/used is zero
AV_1 = 27000/((3.281^3)*60);    % [m^3/s]
ro_1i = ro_alt;
mDot_1 = AV_1*ro_1i;            % [kg/s]
T_1i = T_alt;
P_1e = P_4i_tot;                % Inlet of the Nozzle
P_1i = P_alt;
ro_1e = ((P_1e./P_1i).*(ro_1i.^gammaa)).^(1/gammaa);
h_1i = Cp.*T_1i;
T_1e = P_1e/(R*ro_1e);
T_1e_stat = T_1e;
h_1e = Cp*T_1e;
h_1 = h_1i-h_1e;
u_1 = Cv*(T_1i-T_1e);

Wdot_1 = mDot_1*h_1;
Wdot_1_1run = mDot_WT*h_1;  % [J/s]



%% State 2 - Storage Tank
% Assumptions: Cp = const, Cv = const, ideal gas (P=ro*RT)

Wdot_2 = 0;             % Power generated/used is zero
P_2 = P_4i_tot;         % Needs to be inlet of Nozzle, assume volume is large enough for constant 
T_2 = T_alt;            % Overtime returns back to atmospheric
T_2stat = T_2;
Qdot_2 = 0;
ro_2 = P_2/(R*T_2);
h_2 = Cp*T_2;
u_2 = Cv*T_2;

V_small = (mDot_WT*time)./ (((P_2./(R*T_tank_max)).*((P_tank_max./P_2).^(1/gammaa) - 1)));  % storage tank volume for max pressure
fprintf('Selecting a higher tank pressure, tank volume in cubic meeters are :\n\t\tMach = 3\tMach = 5\tMach = 7\tMach = 10\n');
disp(V_small);
% calculating pressure after run from max pressure.
rho_2 = (1/V_small(4)).*(P_tank_max*V_small(4)./(R*T_tank_max) - mDot_WT.*time);
P_final = P_tank_max.*((1 + ((mDot_WT*time)./(rho_2.*V_small(4)))).^-gammaa);
fprintf('The final pressure in the storage tank after a run starting at max pressure for each mach number is\n')
disp(P_final);

%% State 3 - Heat Exchanger (Prior to WT)
% Assumptions: Cp = const, Cv = const, ideal gas (P=ro*RT)

Wdot_3 = 0;                     % Power generated/used is zero

% exit
mDot_3e = mDot_4i;              % mass flow rate: Diffuser Exit
h_3e = T_4i_tot*Cp;
u_3e = T_4i_tot*Cv;
M_3e = 0.95;                                % estimated Mach at WT inlet
T_3e_stat = T_4i_tot/(1+(((gammaa-1)/2)*M_3e^2));  % static temperature [K]
% inlet
T_3 = T_2;
h_3i = T_3*Cp;
h_3 = h_3i-h_3e;  % I believe the order of this subtraction needs to be swapped!!!
u_3 = Cv*(T_3-T_4i_tot);

Qdot_3 = -h_3.*mDot_3e;
Total_pwr = (Qdot_5 + Qdot_3)*30;  % [Joules]
Total_pwr_kwhr = (Qdot_5 + abs(Qdot_3) + abs(max(Wdot_1_1run)))*30/(3.6e6);  % [kw-hr]
Los_alamos_powr = 154708.07e3;                 % kilowatt hours annually
annual_power_percent = (Total_pwr_kwhr*5*52)/Los_alamos_powr*100;
%% I changed P_3i to P_4i_tot
P_3i = P_4i_tot;                    % Inlet same as storage, exit depends on liquefaction temperature 
ro_3 = ro_2;                    % Assuming incompressible state


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ______________________ LET'S PLOT THIS THING ______________________%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i = 1;
Mach3.ro = [ro_1i, ro_1e(i), ro_2(i), ro_3(i), ro_4i(i), ro_4e(i), ro_5(i), ro_6(i)];
Mach3.T  = [T_1i, T_1e(i), T_2, T_3, T_4i_tot(i), T_4e_tot(i), T_5e(i), T_6(i)];
Mach3.T_stat  = [T_1i/T_1i, T_1e_stat/T_1e, T_2stat/T_2, T_3e_stat(i)/T_4i_tot(i), T_4i_stat/T_4i_tot(i), T_4e_stat(i)/T_4e_tot(i), T_5e/T_5e, T_6stat(i)/T_6(i)];
Mach3.P  = [P_1i, P_1e(i), P_2(i), P_3i(i), P_4i_tot(i), P_4e_tot(i), P_5(i), P_6f(i)];
Mach3.h  = [h_1i, h_1e, h_2, h_3(i), Cp*T_4i_tot(i), Cp*T_4e_tot(i), h_5(i), h_6(i)];
Mach3.v  = 1./Mach3.ro;

i = 2;
Mach5.ro = [ro_1i, ro_1e(i), ro_2(i), ro_3(i), ro_4i(i), ro_4e(i), ro_5(i), ro_6(i)];
Mach5.T  = [T_1i, T_1e, T_2, T_3, T_4i_tot(i), T_4e_tot(i), T_5e, T_6(i)];
Mach5.T_stat  = [T_1i/T_1i, T_1e_stat/T_1e, T_2stat/T_2, T_3e_stat(i)/T_4i_tot(i), T_4i_stat/T_4i_tot(i), T_4e_stat(i)/T_4e_tot(i), T_5e/T_5e, T_6stat(i)/T_6(i)];
Mach5.P  = [P_1i, P_1e(i), P_2(i), P_3i(i), P_4i_tot(i), P_4e_tot(i), P_5(i), P_6f(i)];
Mach5.h  = [h_1i, h_1e, h_2, h_3(i), Cp*T_4i_tot(i), Cp*T_4e_tot(i), h_5(i), h_6(i)];
Mach5.v  = 1./Mach5.ro;

i = 3;
Mach7.ro = [ro_1i, ro_1e(i), ro_2(i), ro_3(i),ro_4i(i), ro_4e(i), ro_5(i), ro_6(i)];
Mach7.T  = [T_1i, T_1e, T_2, T_3, T_4i_tot(i), T_4e_tot(i), T_5e, T_6(i)];
Mach7.T_stat  = [T_1i/T_1i, T_1e_stat/T_1e, T_2stat/T_2, T_3e_stat(i)/T_4i_tot(i), T_4i_stat/T_4i_tot(i), T_4e_stat(i)/T_4e_tot(i), T_5e/T_5e, T_6stat(i)/T_6(i)];
Mach7.P  = [P_1i, P_1e(i), P_2(i), P_3i(i), P_4i_tot(i), P_4e_tot(i), P_5(i), P_6f(i)];
Mach7.h  = [h_1i, h_1e, h_2, h_3(i), Cp*T_4i_tot(i), Cp*T_4e_tot(i), h_5(i), h_6(i)];
Mach7.v  = 1./Mach7.ro;

i = 4;
Mach10.ro = [ro_1i, ro_1e(i), ro_2(i), ro_3(i),ro_4i(i), ro_4e(i), ro_5(i), ro_6(i)];
Mach10.T  = [T_1i, T_1e, T_2, T_3, T_4i_tot(i), T_4e_tot(i), T_5e, T_6(i)];
Mach10.T_stat  = [T_1i/T_1i, T_1e_stat/T_1e, T_2stat/T_2, T_3e_stat(i)/T_4i_tot(i), T_4i_stat/T_4i_tot(i), T_4e_stat(i)/T_4e_tot(i), T_5e/T_5e, T_6stat(i)/T_6(i)];
Mach10.P  = [P_1i, P_1e(i), P_2(i), P_3i(i), P_4i_tot(i), P_4e_tot(i), P_5(i), P_6f(i)];
% Mach10.h  = [h_1i, h_1e, h_2, h_3(i), Cp*T_4i_tot(i), Cp*T_4e_tot(i), h_5(i), h_6(i)];
Mach10.h  = [h_1i, h_1e, h_2, h_3e(i), Cp*T_4i_tot(i), Cp*T_4e_tot(i), h_5e, h_6e(i)];
Mach10.v  = 1./Mach10.ro;


states = [1, 2, 3, 4, 5, 6, 7, 8];
Stateslabels = [cellstr('Compressor Inlet'), cellstr(sprintf('Compressor Exit\nTank Inlet')), cellstr(sprintf('Tank Exit\nHeat Exchanger 1 Inlet')), cellstr(sprintf('Heat Exchanger 1 Exit\nSettling Chamber Inlet')), cellstr(sprintf('Settling Chamber Exit\nWind Tunnel Inlet')), cellstr(sprintf('Diffuser Exit\nHeat Exchanger II Inlet')), cellstr(sprintf('Heat Exchanger II Exit\nVacuum Chamber Inlet')), cellstr(sprintf('Vacuum Chamber Exit'))];
labels = [cellstr('1i'), cellstr('1e'), cellstr('2'), cellstr('3'), cellstr('4i'), cellstr('4e'), cellstr('5'), cellstr('6')];
dx = 0.15;
dy = 0.15;

close all

% Pressure vs. Specific Volume
figure
scatter(Mach3.v, Mach3.P)
hold on
scatter(Mach5.v, Mach5.P)
hold on 
scatter(Mach7.v, Mach7.P)
hold on 
scatter(Mach10.v, Mach10.P)
grid on
xlabel('Specific Volume, v [m^3/kg]','FontSize', 12)
ylabel('Pressure, P [Pa]','FontSize', 12)
title('Pressure vs. Specific Volume','FontSize', 15)
legend('Mach 3', 'Mach 5', 'Mach 7', 'Mach 10')
% saveas(gcf,'Graphs_11000ft/PressureVolume_AllMa.png')

% Pressure vs. Temperature
figure
scatter(Mach3.T, Mach3.P)
hold on
scatter(Mach5.T, Mach5.P)
hold on 
scatter(Mach7.T, Mach7.P)
hold on 
scatter(Mach10.T, Mach10.P)
grid on
xlabel('Temperature, T [K]','FontSize', 12)
ylabel('Pressure, P [Pa]','FontSize', 12)
title('Pressure vs. Temperature','FontSize', 15)
legend('Mach 3', 'Mach 5', 'Mach 7', 'Mach 10')
% saveas(gcf,'Graphs_11000ft/PressureTemperature_AllMa.png')

% Temperature vs. Specific Volume
figure
scatter(Mach3.v, Mach3.T)
hold on
scatter(Mach5.v, Mach5.T)
hold on 
scatter(Mach7.v, Mach7.T)
hold on 
scatter(Mach10.v, Mach10.T)
grid on
xlabel('Specific Volume, v [m^3/kg]','FontSize', 12)
ylabel('Temperature, T [K]','FontSize', 12)
title('Temperature vs. Specific Volume','FontSize', 15)
legend('Mach 3', 'Mach 5', 'Mach 7', 'Mach 10')
% saveas(gcf,'Graphs_11000ft/TemperatureVolume_AllMa.png')

% Enthalpy vs. Specific Volume
figure
scatter(Mach3.v, Mach3.h)
hold on
scatter(Mach5.v, Mach5.h)
hold on 
scatter(Mach7.v, Mach7.h)
hold on 
scatter(Mach10.v, Mach10.h)
grid on
xlabel('Specific Volume, v [m^3/kg]','FontSize', 12)
ylabel('Enthalpy, h [kJ/kg]','FontSize', 12)
title('Enthalpy vs. Specific Volume','FontSize', 15)
legend('Mach 3', 'Mach 5', 'Mach 7', 'Mach 10')
% saveas(gcf,'Graphs_11000ft/EnthalpyVolume_AllMa.png')

% Enthalpy vs. Specific Volume Mach 3
figure
scatter(Mach3.v, Mach3.h)
grid on
xlabel('Specific Volume, v [m^3/kg]','FontSize', 12)
ylabel('Enthalpy, h [kJ/kg]','FontSize', 12)
title('Enthalpy vs. Specific Volume for Mach 3','FontSize', 15)
text(Mach3.v + dx, Mach3.h + dy, labels, 'FontSize', 12)
% enthalpy, pressure at states
fig = figure();
ax = axes(fig);
plot(ax, states, Mach10.h, '--o')
hold(ax, 'on')
grid on
xlabel('States', 'FontSize', 12)
ylabel('Enthalpy, h [kJ/kg]','FontSize', 12)
title('Enthalpy vs. State','FontSize', 15)
text(ax, states, Mach10.h, Stateslabels, 'FontSize', 12)
fig = figure();
ax = axes(fig);
hold(ax, 'on')
plot(ax, states, Mach3.P, '--o')
text(ax, states, Mach3.P, labels, 'FontSize', 12)
plot(ax, states, Mach5.P, '--o')
text(ax, states, Mach5.P, labels, 'FontSize', 12)
plot(ax, states, Mach7.P, '--o')
text(ax, states, Mach7.P, labels, 'FontSize', 12)
plot(ax, states, Mach10.P, '--o')
text(ax, states, Mach10.P, labels, 'FontSize', 12)
hold(ax, 'on')
legend('Mach 3', 'Mach 5', 'Mach 7', 'Mach 10')
grid on
xlabel('States', 'FontSize', 12)
ylabel('Pressure, P [Pa]','FontSize', 12)
title('Pressure vs. State','FontSize', 15)
% saveas(gcf,'Graphs_11000ft/EnthalpyVolume_Ma3.png')
%%%%%%%%%%%%%%%%%%%%%%%%%
% Temperature at states %
%%%%%%%%%%%%%%%%%%%%%%%%%
fig = figure();
ax = axes(fig);
hold(ax, 'on')
plot(ax, states, Mach3.T, '--o')
text(ax, states, Mach3.T, labels, 'FontSize', 12)
plot(ax, states, Mach5.T, '--o')
text(ax, states, Mach5.T, labels, 'FontSize', 12)
plot(ax, states, Mach7.T, '--o')
text(ax, states, Mach7.T, labels, 'FontSize', 12)
plot(ax, states, Mach10.T, '--o')
text(ax, states, Mach10.T, labels, 'FontSize', 12)
hold(ax, 'on')
legend('Mach 3', 'Mach 5', 'Mach 7', 'Mach 10')
grid on
xlabel('States', 'FontSize', 12)
ylabel('Total Temperature, T [K]','FontSize', 12)
title('Total Temperature vs. State','FontSize', 15)

fig = figure();
ax = axes(fig);
hold(ax, 'on')
plot(ax, states, Mach3.T_stat, '--o')
text(ax, states, Mach3.T_stat, labels, 'FontSize', 12)
plot(ax, states, Mach5.T_stat, '--o')
text(ax, states, Mach5.T_stat, labels, 'FontSize', 12)
plot(ax, states, Mach7.T_stat, '--o')
text(ax, states, Mach7.T_stat, labels, 'FontSize', 12)
plot(ax, states, Mach10.T_stat, '--o')
text(ax, states, Mach10.T_stat, labels, 'FontSize', 12)
hold(ax, 'on')
legend('Mach 3', 'Mach 5', 'Mach 7', 'Mach 10')
grid on
xlabel('States', 'FontSize', 12)
ylabel('T/T_o','FontSize', 12)
title('T/T_o vs. State','FontSize', 15)


% Enthalpy vs. Specific Volume Mach 5
figure
scatter(Mach5.v, Mach5.h)
grid on
xlabel('Specific Volume, v [m^3/kg]','FontSize', 12)
ylabel('Enthalpy, h [kJ/kg]','FontSize', 12)
title('Enthalpy vs. Specific Volume for Mach 5','FontSize', 15)
text(Mach5.v + dx, Mach5.h + dy, labels, 'FontSize', 12)
% saveas(gcf,'Graphs_11000ft/EnthalpyVolume_Ma5.png')

% Enthalpy vs. Specific Volume Mach 7
figure
scatter(Mach7.v, Mach7.h)
grid on
xlabel('Specific Volume, v [m^3/kg]','FontSize', 12)
ylabel('Enthalpy, h [kJ/kg]','FontSize', 12)
title('Enthalpy vs. Specific Volume for Mach 7','FontSize', 15)
text(Mach7.v + dx, Mach7.h + dy, labels, 'FontSize', 12)
% saveas(gcf,'Graphs_11000ft/EnthalpyVolume_Ma7.png')

% Enthalpy vs. Specific Volume Mach 10
figure
scatter(Mach10.v, Mach10.h)
grid on
xlabel('Specific Volume, v [m^3/kg]','FontSize', 12)
ylabel('Enthalpy, h [kJ/kg]','FontSize', 12)
title('Enthalpy vs. Specific Volume for Mach 10','FontSize', 15)
text(Mach10.v + dx, Mach10.h + dy, labels, 'FontSize', 12)
% saveas(gcf,'Graphs_11000ft/EnthalpyVolume_Ma10.png')



% Making Data Tables
format shortG
rowNames = [cellstr('Ma 3'), cellstr('Ma 5'), cellstr('Ma 7'), cellstr('Ma 10')];
pressure    = array2table([round(Mach3.P,2); round(Mach5.P,2); round(Mach7.P,2); round(Mach10.P,2)], "VariableNames",labels, 'RowNames',rowNames)
temperature = array2table([round(Mach3.T,2); round(Mach5.T,2); round(Mach7.T,2); round(Mach10.T,2)], "VariableNames",labels, 'RowNames',rowNames)
density     = array2table([round(Mach3.ro,2); round(Mach5.ro,2); round(Mach7.ro,2); round(Mach10.ro,2)], "VariableNames",labels, 'RowNames',rowNames)
enthalpy    = array2table([round(Mach3.h,2); round(Mach5.h,2); round(Mach7.h,2); round(Mach10.h,2)], "VariableNames",labels, 'RowNames',rowNames)

heatTrans = array2table([round(Qdot_3, 2); round(Qdot_5, 2)], "VariableNames", rowNames, 'RowNames', {'Qdot 3' 'Qdot 5'})
