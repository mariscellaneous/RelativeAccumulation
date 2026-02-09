function [facc,nD,naccA] = convacc(rho0,temp,iacc,twt,age)

% iterate based on a guessed accumulation to get density and accumulation
% consistent with a density model and the depth to radar horizon

thres = 0.005;      % threshold for ending iteration
v = 1000;
acc1 = iacc;
% iterate until you converge on an accumulation rate that is consistent
% with the HL FDM and a radar derived twtt
while v > thres
    [nacc,nD,naccA] = newacc(rho0,temp,acc1,twt,age);
    v = abs(acc1 - nacc(end)); % if this is within threshold, stop
    acc1 = nacc(end);
end
facc = nacc;

