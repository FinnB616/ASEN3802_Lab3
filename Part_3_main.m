clear;
clc;
close all;

b = 33 + (1/3);%in ft
c_root = 5 + (1/3);%in ft
c_tip = 3 + (8.5/12);%in ft

V=100;

root_foil = "0012";
tip_foil = "2412";

Num_panels = 72; %total number of panels in the airfoil
Num_odd_terms_pllt = 100;
geo_r = 1;
geo_t = 0;

%Deconstruct input strings into 3 integrers to use in the NACA_Airfoils
%function
m_root = str2double(extractBetween(root_foil,1,1))/100;
p_root = str2double(extractBetween(root_foil,2,2))/10;
t_root = str2double(extractBetween(root_foil,3,4))/100;
m_tip = str2double(extractBetween(tip_foil,1,1))/100;
p_tip = str2double(extractBetween(tip_foil,2,2))/10;
t_tip = str2double(extractBetween(tip_foil,3,4))/100;

%Create naca airfoil coordinates for the root and tip (in ft)
[r_x_b, r_y_b, r_x_camber, r_y_camber] = NACA_Airfoil(m_root,p_root,t_root,c_root,Num_panels);
[t_x_b, t_y_b, t_x_camber, t_y_camber] = NACA_Airfoil(m_tip,p_tip,t_tip,c_tip,Num_panels);

for i=1:15
    CL_root(i)=Vortex_Panel(r_x_b, r_y_b, V, i-6);
    CL_tip(i)=Vortex_Panel(t_x_b, t_y_b, V, i-6);
    alpha(i) = i-6;
end

%use a linear fit to find the lift slopes and convert to 1/rad
root_fit = polyfit(alpha,CL_root,1);
a0_r = root_fit(1)*(180/pi);
tip_fit = polyfit(alpha,CL_tip,1);
a0_t = tip_fit(1)*(180/pi);

%find aerodynamic twist (same thing as zero lift aoa) at root and tip using linear interpolation of CL
%and alpha
aero_r = interp1(CL_root, alpha, 0);
aero_t = interp1(CL_tip, alpha, 0);

%run pllt code for 0 effective aoa and 1 degree to get 3d lift slope. Pllt
%always makes linear lift slope so only need to run 2 anlges of attack.
[e,CL,CDi] = PLLT(b,a0_t,a0_r,c_tip,c_root,aero_t,aero_r,geo_t,geo_r,Num_odd_terms_pllt);
[e_2,CL_2,CDi_2] = PLLT(b,a0_t,a0_r,c_tip,c_root,aero_t,aero_r,geo_t+1,geo_r+1,Num_odd_terms_pllt);

%find 3D lift slope in 1/deg
Lift_slope_3D = CL_2-CL;


%% Deliverable 1: Table of convergence (CL and CDi vs N)
alpha_test = 4; % AoA for convergence study
N_values = 1:1:40; % odd terms

CL_vals = zeros(size(N_values));
CDi_vals = zeros(size(N_values));

for i = 1:length(N_values)
    N = N_values(i);
    
    [~, CL_temp, CDi_temp] = PLLT(b,a0_t,a0_r,c_tip,c_root,aero_t,aero_r,geo_t+alpha_test,geo_r+alpha_test,N);
    
    CL_vals(i) = CL_temp;
    CDi_vals(i) = CDi_temp;
end

% Reference values (largest N)
CL_ref = CL_vals(end);
CDi_ref = CDi_vals(end);

% Relative errors
CL_error = abs((CL_vals - CL_ref)/CL_ref);
CDi_error = abs((CDi_vals - CDi_ref)/CDi_ref);

% Determine N for required tolerances
tolerances = [0.1, 0.01, 0.001];

fprintf('\nDeliverable 1 Table:\n')

for t = 1:length(tolerances)
    tol = tolerances(t);
    
    iCL = find(CL_error < tol,1);
    iCDi = find(CDi_error < tol,1);
    
    fprintf('\nTolerance = %.3f\n', tol);
    fprintf('CL: N = %d, CL = %.4f\n', N_values(iCL), CL_vals(iCL));
    fprintf('CDi: N = %d, CDi = %.5f\n', N_values(iCDi), CDi_vals(iCDi));
end

%% Deliverable 2: Plot CL and CDi vs N with error markers
figure
plot(N_values, CL_vals, 'LineWidth',2)
hold on

% Mark tolerance locations
labels = {'10%','1%','0.1%'};
for t = 1:length(tolerances)
    tol = tolerances(t);
    iCL = find(CL_error < tol,1); % lines at N=2,4,6
    xline(N_values(iCL),'--',labels{t});
end

xlim([1 10]) % to zoom in and see differences
ylim([.4 0.6])
xlabel('Number of odd terms, N')
ylabel('Lift Coefficient')
title('Convergence of C_L')
xline(NaN,'--'); % dummy for legend
legend('C_L','Relative error thresholds','Location','best')
grid on



figure
hold on
plot(N_values, CDi_vals, 'LineWidth',2)

% Mark tolerance locations
labels = {'10%','1%','0.1%'};
for t = 1:length(tolerances)
    tol = tolerances(t);
    iCDi = find(CDi_error < tol,1); % lines at N=2,4,6
    xline(N_values(iCDi),'--',labels{t});
end

xlim([1 10]) % to zoom in and see differences
ylim([0 0.05])
xlabel('Number of odd terms, N')
ylabel('Induced Drag Coefficinet')
title('Convergence of C_{D,i}')
xline(NaN,'--'); % dummy for legend
legend('C_{D,i}','Relative error thresholds','Location','best')
grid on

%% Deliverable 3: Table of forces and efficiency
rho = 0.001756; % slug/ft^3
V_ft = V * 1.68781; % ft/s

S = b * (c_root + c_tip)/2;

N_final = N_values(find(CL_error < 0.001,1));

[~, CL_final, CDi_final] = PLLT(b,a0_t,a0_r,c_tip,c_root,aero_t,aero_r,geo_t+alpha_test,geo_r+alpha_test,N_final);

L = 0.5 * rho * V_ft^2 * S * CL_final;
Di = 0.5 * rho * V_ft^2 * S * CDi_final;

cd_profile = 0.01;
CD_total = CDi_final + cd_profile;

LD = CL_final / CD_total;

fprintf('\nDeliverable 3 Table:\n')
fprintf('Lift = %.2f lb\n', L);
fprintf('Induced Drag = %.2f lb\n', Di);
fprintf('L/D = %.2f\n', LD);

%% Deliverable 4: Total drag vs angle of attack
alpha_sweep = -2:1:10;

CDi_array = zeros(size(alpha_sweep));
CD_total_array = zeros(size(alpha_sweep));

for i = 1:length(alpha_sweep)
    
    [~, CL_temp, CDi_temp] = PLLT(b,a0_t,a0_r,c_tip,c_root,aero_t,aero_r,geo_t+alpha_sweep(i),geo_r+alpha_sweep(i),N_final);
    
    CDi_array(i) = CDi_temp;
    CD_total_array(i) = CDi_temp + cd_profile;
end

figure
plot(alpha_sweep, CD_total_array, 'LineWidth',2)
hold on
plot(alpha_sweep, CDi_array, '--', 'LineWidth',2)
xlabel('Angle of attack (deg)')
ylabel('Drag coefficient')
title('Total and Induced Drag vs Angle of Attack')
legend('Total drag','Induced drag')
grid on

%% Deliverable 5: L/D vs angle of attack
LD_array = zeros(size(alpha_sweep));

for i = 1:length(alpha_sweep)
    
    [~, CL_temp, CDi_temp] = PLLT(b,a0_t,a0_r,c_tip,c_root,aero_t,aero_r,geo_t+alpha_sweep(i),geo_r+alpha_sweep(i),N_final);
    
    CD_total_temp = CDi_temp + cd_profile;
    LD_array(i) = CL_temp / CD_total_temp;
end

figure
plot(alpha_sweep, LD_array, 'LineWidth',2)
xlabel('Angle of attack (deg)')
ylabel('Lift to drag ratio (L/D)')
title('Aerodynamic Efficiency vs Angle of Attack')
grid on
