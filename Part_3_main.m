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
Num_odd_terms_pllt = 10;
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
a0_r = root_fit(1)*(pi/180);
tip_fit = polyfit(alpha,CL_tip,1);
a0_t = tip_fit(1)*(pi/180);

%find aerodynamic twist (same thing as zero lift aoa) at root and tip using linear interpolation of CL
%and alpha
aero_r = interp1(CL_root, alpha, 0);
aero_t = interp1(CL_tip, alpha, 0);


[e,CL,CDi] = PLLT(b,a0_t,a0_r,c_tip,c_root,aero_t,aero_r,geo_t,geo_r,Num_odd_terms_pllt);