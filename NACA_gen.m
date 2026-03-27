clear;
clc;
close all;

%[x_b,y_b] = function NACA_Airfoils(m,p,t,c,N)

m=.04;
p=.4;
t=.12;
c=1;
N=1000;

%create n/2 x positions
x_positions = linspace(0,c,N/2);

y_t = ((t./0.2).*c).* ((.2969.*sqrt(x_positions./c)) - (.1260.*(x_positions./c)) - (.3516.*((x_positions./c).^2)) + (.2843.*((x_positions./c).^3)) - (.1036.*((x_positions./c).^4)));


front_positions = x_positions(x_positions < (p*c));
rear_positions = x_positions(x_positions >= (p*c));

y_c = (m.*front_positions./(p.^2)).*((2.*p) - (front_positions./c));
y_c = [y_c (m.*(c-rear_positions)./((1-p).^2)).*(1+(rear_positions./c)-(2.*p))];

for i=1:(length(x_positions)-1)
    zeta(1)=0;
    zeta(i+1) = atan((y_c(i+1)-y_c(i))/(x_positions(i+1)-x_positions(i)));
end
x_U = x_positions - (y_t.*sin(zeta));
y_U = y_c + y_t.*cos(zeta);
x_L = x_positions + (y_t.*sin(zeta));
y_L = y_c - y_t.*cos(zeta);

x_b = [flip(x_L) x_U(2:end)];
y_b = [flip(y_L) y_U(2:end)]; 

%% debugging plots
% figure;
% plot(x_positions,y_t);
% 
% figure;
% plot(x_positions,y_c)
% 
% plot(x_U,y_U);
% hold on
% plot(x_L,y_L);
% axis equal
plot(x_b,y_b);
axis equal

