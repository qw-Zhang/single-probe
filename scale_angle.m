function [ang1,ang2] = scale_angle(phi_real_1,phi_real_2)
ang1 = round(rad2deg(phi_real_1));
ang2 = round(rad2deg(phi_real_2));
if ang1 <= 0
    ang1 = ang1 + 360;
elseif ang1 > 360
    ang1 = ang1 - 360;
end
if ang2 <= 0
    ang2 = ang2 + 360;
elseif ang2 > 360
    ang2 = ang2 - 360;
end
end