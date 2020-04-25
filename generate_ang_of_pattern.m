function [ang_P_1,ang_P_2] = generate_ang_of_pattern(dis_1,dis_2,d,r,phi_1_rad,phi_2_rad)
phi = (phi_2_rad - phi_1_rad)/2;
%four part within [-pi,pi] of phi 
if (phi < pi/2)&&(phi > -pi/2)
    sin_top_1_rad = d*sin(phi_1_rad)/dis_1;
    ang_P_1 = real(asin(sin_top_1_rad)) + phi_1_rad;
    if phi>0
        ang_P_2 = real(asin(r*sin(phi_2_rad)/dis_2));
    else
        ang_P_2 = pi - real(asin(r*sin(phi_2_rad)/dis_2));
    end
else
    if phi>0
        phi_1_rad = -phi_1_rad;
        ang_P_1 = 2*pi - (real(asin(d*sin(phi_1_rad)/dis_1)) + phi_1_rad);
        phi_2_rad = 2*pi - (phi_2_rad + pi/2);
        ang_P_2 = 2*pi - (pi - phi_2_rad - (pi/2 + real(asin(d*sin(phi_2_rad)/dis_2))));
    else
        phi_1_rad = 2*pi - phi_1_rad;
        ang_P_1 = 2*pi - (real(asin(d*sin(phi_1_rad)/dis_1)) + phi_1_rad);
        phi_2_rad = abs(phi_2_rad);
        ang_P_2 = 2*pi - (pi - phi_2_rad - real(asin(d*sin(phi_2_rad)/dis_2)));
    end
end
[ang_P_1, ang_P_2] = scale_angle(ang_P_1, ang_P_2);
end