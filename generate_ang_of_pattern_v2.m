function [ang_P_1,ang_P_2] = generate_ang_of_pattern_v2(top,pos_ant_1,pos_ant_2)
%generate angle with porobe location and antennas
ang(1,1) = top(1) - pos_ant_1(1);
ang(1,2) = top(2) - pos_ant_1(2);
ang(2,1) = top(1) - pos_ant_2(1);
ang(2,2) = top(2) - pos_ant_2(2);

ang_P_1 = atan2(ang(1,2), ang(1,1));
ang_P_2 = atan2(ang(2,2), ang(2,1));
[ang_P_1, ang_P_2] = scale_angle(ang_P_1, ang_P_2);
end