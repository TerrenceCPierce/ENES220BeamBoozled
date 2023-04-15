clear;
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

flange_width = 1.125; %in
flange_height = 0.1875; %in
web_width = 0.5625; %in
web_height = 1.8125; %in

flange_n = E_Oak/E_Oak; %E chosen/ E material
web_n = E_Oak/E_Pine; %E chosen/ E material

%Converted cross sections
flange_width_con = flange_width/flange_n; %in
web_width_con = web_width/web_n; %in

%converted I cross section
I = 1/12*(flange_width_con*(2*flange_height+web_height)^3-(flange_width_con-web_width_con)*web_height^3);

syms x; syms F;
%using lower E value, E pine
deflectLeft = -F*b*x/(6*E_Oak*I*L)*(L^2-b^2-x^2);
deflectRight = -F*a*(L-x)/(6*E_Oak*I*L)*(L^2-a^2-(L-x)^2);
slopeLeft = diff(deflectLeft);
slopeRight = diff(deflectRight);

MLeft = diff(slopeLeft)*E_Oak*I;
MRight = diff(slopeRight)*E_Oak*I;

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
M_max =  maxMArr(1);

if abs(diff(VLeft))>abs(diff(VRight))
V_max = VLeft;
else
V_max = VRight;
end



%Now do max stress calculations for tensile, shear, and glue shear

%This assumes strength is the same in compression and tension
%this assumption may not hold for pine
bendStress_flange = -M_max*(flange_height+web_height/2)/(I*flange_n);
bendStress_web = -M_max*(web_height/2)/(I*web_n);

%shear and glue
Q_joint = flange_width_con*flange_height*(flange_height/2+web_height/2);
Q_max = Q_joint + web_width_con*web_height/2*web_height/4;

shearStress_max = V_max*Q_max/(web_n*I*web_width_con); %middle of I
shearStress_flange = V_max*Q_joint/(flange_n*I*flange_width_con); %connection between flange and web, on flange side
shearStress_glueWeb = V_max*Q_joint/(web_n*I*web_width_con); %connection between flange and web, on web side
%not sure if I should multiple last equation by web_n or not

%calculate safety factors
%need to take derivatives and compare
SF_bendStress = min(abs(Oak_Tens_Avg*F/bendStress_flange),abs(Pine_Tens_Avg*F/bendStress_web))/F;
SF_shearStress = min(abs(Oak_Shear_Avg*F/shearStress_flange),abs(Pine_Shear_Avg*F/shearStress_max))/F;
SF_shearStressGlue = min(abs(OakGlue_Shear_Avg*F/shearStress_flange),abs(PineGlue_Shear_Avg*F/shearStress_glueWeb))/F;

F_max_bendStress = vpa(solve(SF_bendStress==1,F));
F_max_shearStress = vpa(solve(SF_shearStress==1,F));
F_max_shearStressGlue = vpa(solve(SF_shearStressGlue==1,F));

F_max = double(min([F_max_bendStress,F_max_shearStress,F_max_shearStressGlue]));

flangeVolume = L_tot*flange_width*flange_height;
FlangeWeight = flangeVolume*dens_Oak;

webVolume = L_tot*web_width*web_height;
webWeight = webVolume*dens_Pine;

weight = 2*FlangeWeight + webWeight;

str2Weight = F_max/weight;

syms F; syms x;
slopeLeftmax = subs(slopeLeft,F,double(F_max));
maxDeflect_x = solve(slopeLeftmax,x);
x_val_max = double(max(maxDeflect_x));
syms F; syms x;
maxDeflect = subs(deflectLeft,[x, F],[x_val_max, F_max]);