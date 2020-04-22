function [ang_P_1,ang_P_2] = generate_ang_pattern(dis_1,dis_2,r,phi_1,phi_2)
[ang_1, ang_2] = scale_angle(phi_1,phi_2);
sin_ang_1 = r*sin(deg2rad(ang_1))/dis_1;
sin_ang_2 = r*sin(deg2rad(ang_2))/dis_2;
ang_P_1 = real(asin(sin_ang_1));
ang_P_2 = real(asin(sin_ang_2));
[ang_P_1, ang_P_2] = scale_angle((pi - ang_P_1),(ang_P_2));
% ang_P_1 = ang_P_1+1;ang_P_2 = ang_P_2+2;
end