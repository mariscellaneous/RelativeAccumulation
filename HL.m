function [depth, rho, age] = HL(acc,temp,rho0)

% uses the Herron and Langway FDM to determine depth-density profile of
% dry snow

R = 8.314;   % gas constant
rhoi = 0.917;   % density of ice

depth = (0:0.01:100)';

% rate constants
k0 = 11 * exp((-10160 / (R * temp)));
k1 = 575 * exp((-21400 / (R * temp)));

w1 = log(rho0 / (rhoi - rho0));
w2 = log(0.55 / (rhoi - 0.55));

% depth to 2nd stage of densification
h55 = (1 / (rhoi * k0)) * (w2 - w1);
% age of firn at transition to 2nd stage
t55 = (1 / (k0 * acc)) * log((rhoi - rho0) ./ (rhoi - 0.55));


i1 = find(depth <= h55);
i2 = find(depth > h55);
d1 = depth(i1);
d2 = depth(i2);

% determine densities and ages for the first stage of densification
z0 = exp(rhoi * k0 * d1 + w1);
r1 = (rhoi * z0) ./ (1 + z0);
a1 = (1 / (k0 * acc)) * log((rhoi - rho0) ./ (rhoi - r1));
% determine densities and ages for the second stage of densification
z1 = exp((rhoi * k1 * (d2 - h55)) / sqrt(acc) + w2);
r2 = (rhoi * z1) ./ (1 + z1);
a2 = (1 / (k1 * sqrt(acc))) * log((rhoi - 0.55) ./ (rhoi - r2)) + t55;

% combine!
rho = [r1; r2];
age = [a1; a2];