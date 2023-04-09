
L_tot = 22; %in, total length
L = 20; %in, length between supports

Oak_Tens_Avg = 17873; %psi
Oak_Tens_Min = 12283; %psi
Oak_Shear_Avg = 2873; %psi
Oak_Shear_Min = 2230; %psi
OakGlue_Shear_Avg = 1391; %psi
OakGlue_Shear_Min = 1012; %psi

Pine_Tens_Avg = 14327; %psi
Pine_Tens_Min = 5533; %psi
Pine_Shear_Avg = 1492; %psi
Pine_Shear_Min = 1425; %psi
PineGlue_Shear_Avg = 989; %psi
PineGlue_Shear_Min = 525; %psi

E_Oak = 1.8*10^6; %psi
E_Pine = 1.5*10^6; %psi

dens_Oak = .024; %lb/in^3
dens_Pine = .014; %lb/in^3

a = 12; %in
b = 8; %in

F = 2000; %lb

flange_width = 2; %in
flange_height = 3/4; %in
web_width = 3/4; %in
web_height = 4-2*flange_height; %in

flange_n = E_Oak/E_Oak; %E chosen/ E material
web_n = E_Oak/E_Oak; %E chosen/ E material

%converted I cross section
I = 1/12*(flange_width*flange_n*(2*flange_height+web_height)^3-(flange_width-web_width)*web_n*web_height^3);


syms x;
%using lower E value, E pine
deflectLeft = -F*b*x/(6*E_Pine*I*L)*(L^2-b^2-x^2);
deflectRight = -F*a*(L-x)/(6*E_Pine*I*L)*(L^2-a^2-(L-x)^2);
slopeLeft = diff(deflectLeft);
slopeRight = diff(deflectRight);

MLeft = diff(slopeLeft)*E_Pine*I;
MRight = diff(slopeRight)*E_Pine*I;

VLeft = diff(MLeft);
VRight = diff(MRight);

%{
%This won't work because there's a discontinuity at the point I'm trying
to get
% I will hardcode in the value instead
maxM_x_left = max(solve(VLeft));
maxM_x_right = min(solve(VRight));
%}

maxM_x = a;
x = maxM_x;
maxM_left = subs(MLeft);
x = maxM_x;
maxM_right = subs(MRight);

maxMArr = [maxM_left;maxM_right];
[~,X] = max(abs(maxMArr));
M_max =  maxMArr(X);

if abs(VLeft)>abs(VRight)
    V_max = VLeft;
else
    V_max = VRight;
end

maxDeflect_x = solve(slopeLeft);
x = max(maxDeflect_x);
maxDeflect = vpa(subs(deflectLeft));

%Now do max stress calculations for tensile, shear, and glue shear


