fc = 2.535e9;
ant = design(dipole ,fc);
%generate the pattern gain of the antenna.
%parameter:antenna,frequency,al,ez.
%G(dBi)=10lgGi G(dBd)=10lgG dBi=dBd+2.15
[D,al,ez] = pattern(ant,fc,0:1:360,-90:1:90);
[pat_az] = patternAzimuth(ant,fc);

%patternAzimuth(object,frequency,elevation)
%patternElevation(object,frequency,azimuth) 

%%
pcolor(al,ez,D)
shading interp;
colorbar
%%
contourf(al,ez,D)
colorbar