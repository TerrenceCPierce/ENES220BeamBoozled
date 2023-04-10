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

bestStr2Weight = 0;
bestConfig = [0,0,0,0];
for flange_width = 3/16:1/16:2 %29 times
    for flange_height = 3/16:1/16:3/4 %9 times
        for web_width = 3/16:1/16:min(flange_width,3/4) %9 times or less
            for web_height= 3/16:1/16:2*flange_width-2*flange_height %58 times or less

                flange_n = E_Pine/E_Pine; %E chosen/ E material
                web_n = E_Pine/E_Pine; %E chosen/ E material
                
                %converted I cross section
                I = 1/12*(flange_width*flange_n*(2*flange_height+web_height)^3-(flange_width-web_width)*web_n*web_height^3);
                
                
                syms x; syms F;
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
                M_max =  maxMArr(1);
                
                if abs(diff(VLeft))>abs(diff(VRight))
                    V_max = VLeft;
                else
                    V_max = VRight;
                end
                

                
                %Now do max stress calculations for tensile, shear, and glue shear
                
                %This assumes strength is the same in compression and tension
                %this assumption may not hold for pine
                bendStress_flange = -M_max/I*(flange_height+web_height/2);
                bendStress_web = -M_max/I*(web_height/2);
                
                %shear and glue
                Q_joint = flange_width*flange_n*flange_height*(flange_height/2+web_height/2);
                Q_max = Q_joint + web_width*web_n*web_height*web_height/4;
                
                shearStress_max = V_max*Q_max/(web_n*I*web_width); %middle of I
                shearStress_flange = V_max*Q_joint/(flange_n*I*flange_width); %connection between flange and web
                shearStress_glue = V_max*Q_joint/(web_n*I*web_width); %connection between flange and web
                %not sure if I should multiple last equation by web_n or not
                
                %calculate safety factors
                SF_bendStress = min(abs(Pine_Tens_Avg*F/bendStress_flange),abs(Pine_Tens_Avg*F/bendStress_web))/F;
                SF_shearStress = min(abs(Pine_Shear_Avg*F/shearStress_max),abs(Pine_Shear_Avg*F/shearStress_flange))/F;
                SF_shearStressGlue = PineGlue_Shear_Avg/abs(shearStress_glue);

                F_max_bendStress = vpa(solve(SF_bendStress==1,F));
                F_max_shearStress = vpa(solve(SF_shearStress==1,F));
                F_max_shearStressGlue = vpa(solve(SF_shearStressGlue==1,F));

                F_max = min([F_max_bendStress,F_max_shearStress,F_max_shearStressGlue]);

                flangeVolume = L_tot*flange_width*flange_height;
                FlangeWeight = flangeVolume*dens_Pine;
                
                webVolume = L_tot*web_width*web_height;
                webWeight = webVolume*dens_Pine;
                
                weight = 2*FlangeWeight + webWeight;
                
                %This is not accurate because this F is not the F that the beam breaks
                str2Weight = F_max/weight;
                
                if(str2Weight> bestStr2Weight && F_max > 1000 && F_max < 2500)
                    bestStr2Weight = str2Weight;
                    bestConfig = [flange_width,flange_height,web_width,web_height];
                end
            end
        end
    end
end

%{
                maxDeflect_x = solve(slopeLeft,x);
                x = max(maxDeflect_x);
                maxDeflect = subs(deflectLeft);
%}