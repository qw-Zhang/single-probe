function P_az = generate_pattern(fc)
%     % azimuth for dipole
%     P_az = patternAzimuth(ant,plotFrequency,0,'Azimuth',-180:1:180);
%     

% Create a dipole antenna
%% Antenna Properties 
% Design antenna at frequency 2535000000Hz
antennaObject = design(dipole,fc);
% Properties changed 
antennaObject.Width = 0.03;
antennaObject.Length = 0.2;
% Update load properties 
antennaObject.Load.Impedance = 50;
antennaObject.Load.Frequency = fc;

%% Antenna Analysis 
% Define plot frequency 
plotFrequency = fc;
% Define frequency range 
% freqRange = (2281.5:25.35:2788.5) * 1e6;
% % show for dipole
% figure;
% show(antennaObject) 
% % pattern for dipole
% figure;
% pattern(antennaObject, plotFrequency) 
% % azimuth for dipole
% figure;
P_az = patternAzimuth(antennaObject, plotFrequency,0,'Azimuth',-180:1:180);

