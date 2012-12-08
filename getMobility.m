function mu = getMobility(C, loc, tsteps)

nsteps = size(loc,2);
% if nsteps = 1000 and tsteps = 4: loc(5:1000) - loc(1:996)
d1 = loc((1+tsteps):nsteps)-loc(1:(nsteps-tsteps));
time1 = C.tstep * tsteps;
avgx2 = mean(d1.^2);

mu = (1.6e-19 * (4e-8).^2 /2/1.38e-23/298/48.8e-15) * avgx2/time1;


