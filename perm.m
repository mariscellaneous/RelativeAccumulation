function [e] = perm(rho)

% determine the permittivity of the firn based on its density
n = 1; % if 0, Looyenga.  If 1, Kovacs. 
rhoi = 0.917;   % density of ice
eice = 3.15;    % permittivity of air
eair = 1;       % permittivity of ice

if n == 0
    e = ((eice^(1/3) - eair^(1/3)) * (rho / rhoi) + eair^(1/3)).^3; 
else
    e = (1 + 0.845 * rho).^2;
end

