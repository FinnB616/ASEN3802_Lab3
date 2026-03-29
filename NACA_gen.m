clear;
clc;
close all;

%[x_b,y_b] = function NACA_Airfoils(m,p,t,c,N)

m=.04;
p=.4;
t=.12;
c=1;
N=1000;

%create n/2 different x positions
x_positions = linspace(0,c,round(N/2));

%find the thickness at each x position
y_t = ((t./0.2).*c).* ((.2969.*sqrt(x_positions./c)) - (.1260.*(x_positions./c)) - (.3516.*((x_positions./c).^2)) + (.2843.*((x_positions./c).^3)) - (.1036.*((x_positions./c).^4)));

%split the positions into before the max thickness and after
front_positions = x_positions(x_positions < (p*c));
rear_positions = x_positions(x_positions >= (p*c));

%create the mean camber line (different eqautions infront and behind the
%max thickness)
y_c = (m.*front_positions./(p.^2)).*((2.*p) - (front_positions./c));
y_c = [y_c (m.*(c-rear_positions)./((1-p).^2)).*(1+(rear_positions./c)-(2.*p))];

%find the slope of the camber line using the derivitives of the expression
%to find the mean camber line
dyc_dx = (2*m/p) - (2.*m.*front_positions./(p^2 * c));
dyc_dx = [dyc_dx (((2*m)/((1-p).^2)).*(p-(rear_positions./c)))];
%convert the slope to an angle
zeta = atan(dyc_dx);

%create vectors for the upepr and lower surfaces of the wing
x_U = x_positions - (y_t.*sin(zeta));
y_U = y_c + y_t.*cos(zeta); 
x_L = x_positions + (y_t.*sin(zeta));
y_L = y_c - y_t.*cos(zeta);

%combine the upper and lower lines into a single back to front to back line
x_b = [flip(x_L) x_U(2:end)];
y_b = [flip(y_L) y_U(2:end)]; 

%% debugging plots
plot(x_b,y_b);
axis equal
