function [e,CL,CDi] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)
%% Input:
% b = wingspan (ft)
% a0_t/r = lift slope at tip/root (per rad)
% c_t/r = chord at tip/root (ft)
% aero_t/r = zero-lift AoA at tip/root (deg)
% geo_t/r = geometric AoA at tip/root (deg)
% N = number of Odd terms in series
%% Output:
% e = span efficiency factor
% CL = lift coefficient
% CDi = induced drag coefficient

%% Linear variation across span 
% cos(0) at tip and cos(pi/2) at root
% theta = ipi/2N
theta = (1:N) * pi / (2*N);
% chord
c = c_r + (c_t - c_r).*cos(theta);
% lift slope
a0 = a0_r + (a0_t - a0_r).*cos(theta);
% geometric AoA (convert to radians)
alpha_geo = deg2rad(geo_r + (geo_t - geo_r).*cos(theta));
% zero-lift AoA (convert to radians)
alpha_L0 = deg2rad(aero_r + (aero_t - aero_r).*cos(theta));

%% [M][A] = [b]
M = zeros(N,N);
alpha_eff = alpha_geo - alpha_L0; % [b]

for i = 1:N
    for j = 1:N
        n = 2*j - 1; % Only odd terms
        % First term (geometric AoA contribution)
        term1 = (4*b)/(a0(i)*c(i)) * sin(n*theta(i)); %a0=2pi in main
        % Second term (induced AoA contribution)
        term2 = n * sin(n*theta(i)) / sin(theta(i));
        % Fill matrix
        M(i,j) = term1 + term2;
    end
end

A = M \ alpha_eff';

% AR and wing area calc
S = b * (c_r + c_t) / 2;
AR = b^2 / S; 

% Lift Coefficient
CL = A(1) * pi * AR;

% Induced Drag
sum_term = 0;
for j = 1:N
    n = 2*j - 1;
    sum_term = sum_term + n * A(j)^2;
end
CDi = pi * AR * sum_term;

% Span Efficiency
delta = 0;
for j = 2:N % start at second term
    n = 2*j - 1;
    delta = delta + n * A(j)^2;
end
delta = delta / (A(1)^2);

e = 1 / (1 + delta);
end