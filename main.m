%% Aero Lab Part 1 MAIN (Tasks 1-4)
clear;
clc;
close all;

%% Task 1: Airfoil Geometry

% Parameters
c = 1;
N = 100; % 50 panels per surface (100 total)

%% NACA 0021 (symmetric)
m1 = 0.0;
p1 = 0.0;
t1 = 0.21;
[x1, y1, xc1, yc1] = NACA_Airfoil(m1, p1, t1, c, N);

%% NACA 2421 (cambered) 
m2 = 0.02;
p2 = 0.4;
t2 = 0.21;
[x2, y2, xc2, yc2] = NACA_Airfoil(m2, p2, t2, c, N);

%% Plot
figure
hold on

% NACA 0021
plot(x1, y1, 'b-', 'LineWidth', 1.5)
plot(x1, y1, 'bo', 'MarkerSize', 4)

% NACA 2421 
plot(x2, y2, 'r-', 'LineWidth', 1.5)
plot(x2, y2, 'ro', 'MarkerSize', 4)

% Camber line (only for 2421)
plot(xc2, yc2, 'k--', 'LineWidth', 2)
axis equal
grid on
xlabel('x/c')
ylabel('y/c')
title('NACA 0021 and NACA 2421 Airfoils')
legend('0021 Surface','0021 Panel Points', ... 
  '2421 Surface','2421 Panel Points','2421 Camber Line', 'Location','best')

%% Task 2: Convergence Study

%% Task 3: Effect of Airfoil Thickness on Lift

%given minimum panels needed in task 2 = 70  
N_min = 70;

% random freestream velo
VINF = 100; %m/s?

% aoa alpha ranges from -5 to 15 degrees
alpha = -5:1:15;

% NACA 0006
[x31, y31, xc31, yc31] = NACA_Airfoil(0.0, 0.0, 0.06, c, N_min);

% NACA 0012
[x32, y32, xc32, yc32] = NACA_Airfoil(0.0, 0.0, 0.12, c, N_min);

% NACA 0018
[x33, y33, xc33, yc33] = NACA_Airfoil(0.0, 0.0, 0.18, c, N_min);

% lift coeff arrays
cl_06 = zeros(size(alpha));
cl_12 = zeros(size(alpha));
cl_18 = zeros(size(alpha));

% call VPM for each AF at each aoa
for i = 1:length(alpha)
    cl_06(i)=Vortex_Panel(x31,y31,VINF,alpha(i));
    cl_12(i)=Vortex_Panel(x32,y32,VINF,alpha(i));
    cl_18(i)=Vortex_Panel(x33,y33,VINF,alpha(i));
end

% symmetric AF -> thin airfoil theory
cl_TAT = 2*pi*deg2rad(alpha);

% use fit to find lift slope and L=0 aoa
idx = (alpha >= -4) & (alpha <= 8);

p06 = polyfit(alpha(idx), cl_06(idx), 1);
p12 = polyfit(alpha(idx), cl_12(idx), 1);
p18 = polyfit(alpha(idx), cl_18(idx), 1);

a0_06 = p06(1);
a0_12 = p12(1);
a0_18 = p18(1);

alphaL0_06 = -p06(2)/p06(1);
alphaL0_12 = -p12(2)/p12(1);
alphaL0_18 = -p18(2)/p18(1);

a0_TAT = 2*pi*(pi/180);
alphaL0_TAT = 0;

%Plotting
figure
hold on
grid on
box on


plot(alpha, cl_06, 'b-', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'VPM NACA 0006')
plot(alpha, cl_12, 'r-', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'VPM NACA 0012')
plot(alpha, cl_18, 'g-', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'VPM NACA 0018')
plot(alpha, cl_TAT, 'k--', 'LineWidth', 2, 'DisplayName', 'Thin Airfoil Theory')

xlabel('\alpha (deg)')
ylabel('c_l')
title('Effect of Airfoil Thickness on Lift')
legend('Location','best')

%% Task 4: Effect of Airfoil Camber on Lift

% Function inputs
VINF = 100;
ALPHA = -5:15;


% From Task 2: optimal number of panels is 70
N = 70;

% NACA 0012
m1 = 0;
p1 = 0;
t1 = 0.12;
[x1, y1, xc1, yc1] = NACA_Airfoil(m1, p1, t1, c, N);

% NACA 2412 (Moderately Cambered Airfoil)

m1 = 0.02;
p1 = 0.4;
t1 = 0.12;
[x2, y2, xc2, yc2] = NACA_Airfoil(m2, p2, t2, c, N);

% NACA 4412 (Significantly Cambered Airfoil)

m3 = 0.04;
p3 = 0.4;
t3 = 0.12;
[x3, y3, xc3, yc3] = NACA_Airfoil(m3, p3, t3, c, N);


% Vortex Panel Function for varying alpha

for i = 1:length(ALPHA)
    CL1(i) = Vortex_Panel(x1, y1, VINF, ALPHA(i)); % Symmetric
    CL2(i) = Vortex_Panel(x2, y2, VINF, ALPHA(i)); % Moderately Cambered
    CL3(i) = Vortex_Panel(x3, y3, VINF, ALPHA(i)); % Significantly Cambered
end

% interpolating alpha to find exactly when CL is predicted to be zero

alpha0_0012 = interp1(CL1, ALPHA, 0);
alpha0_2412 = interp1(CL2, ALPHA, 0);
alpha0_4412 = interp1(CL3, ALPHA, 0);

% Thin Airfoil Theory, assuming x-inter to be 0, -2, -4 respectfully

% TAT uses radians
alpha_rad = deg2rad(ALPHA);

% CL / alpha is always 2pi with different x intercepts 

CL_TAT_0012 = 2*pi*(alpha_rad - deg2rad(0));
CL_TAT_2412 = 2*pi*(alpha_rad - deg2rad(-2));
CL_TAT_4412 = 2*pi*(alpha_rad - deg2rad(-4));

% Use poly fit function to find slope of each airfoil 
slope_0012 = polyfit(ALPHA, CL1, 1);
slope_2412 = polyfit(ALPHA, CL2, 1);
slope_4412 = polyfit(ALPHA, CL3, 1);
figure

% Vortex Panel
plot(ALPHA, CL1, 'b-', 'LineWidth', 1)
hold on
plot(ALPHA, CL2, 'r-', 'LineWidth', 1)
plot(ALPHA, CL3, 'g-', 'LineWidth', 1)

% Thin Airfoil Theory 
plot(ALPHA, CL_TAT_0012, 'b--')
plot(ALPHA, CL_TAT_2412, 'r--')
plot(ALPHA, CL_TAT_4412, 'g--')

xlabel('\alpha (deg)')
ylabel('Cl')
title('Effect of Airfoil Camber on Lift')
legend('NACA 0012','NACA 2412', 'NACA 4412', 'Thin Airfoil 0012', 'Thin Airfoil 2412', 'Thin Airfoil 4412' )
grid on
hold off
