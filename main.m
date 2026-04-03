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



