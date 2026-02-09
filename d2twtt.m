function tt = d2twtt(e)

% use the depths and permittivity to convert to twtt

c = 299792458;  % m/s

deltt = (2 * sqrt(e) .* 0.01) / c;  % HL model is at 0.01 m vertical rez.
tt = cumsum(deltt);