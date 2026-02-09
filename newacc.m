function [nacc,nD,naccA] = newacc(rho0,temp,iacc,twt,age)

% given an initial guess of accumulation rate, as well as rho0 and temp
% which do not vary, determine a new accumulation rate that is consistent
% with measured twtt to a marker at depth

% first get depth-density initial guess
[depth, rho] = HL(iacc,temp,rho0);  % depth in meters, rho in g/cm^3
% determine the permittivity
[e] = perm(rho);
% then convert modeled depths to twtt based on the permittivity
tt = d2twtt(e);   % seconds
% get the cumulative load (mass) based on the density
cl = cumsum(rho * 1);   % results are in g/cm^2

% now we have a depth to twtt conversion, so we can calculate a new
% accumulation rate based on the radar twtt
nD = interp1(tt,depth,twt);       % depth of horizon (meters)
% determine the cumulative mass above that depth
nC = interp1(tt,cl,twt);       % load in g/cm^2

% calculate the new long-term accumulation rate
nacc = (nC ./ age) / 100;   % accumulation rate in m w.e. per year

% calculate the new annual accumulation rate
naccA = (diff([0;nC]) ./ diff([0;age])) / 100;   % accumulation rate in m w.e. per year