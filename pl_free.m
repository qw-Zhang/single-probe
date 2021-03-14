function [mag,db] = pl_free(dis_m,wavelength)
c = 3e8;
fre = c / wavelength;
dis_km = dis_m / 1e3;

db = 32.5 + 20*log10(fre) + 20*log10(dis_km);
mag = db2mag(db);
end
