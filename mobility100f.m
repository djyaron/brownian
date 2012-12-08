%% Common
clear classes;
Clib.type    = 'perES';
Clib.nangles = 100;
Clib.Vgs     = 1;
Clib.betaES  = -20;
Clib.beta1   = 1;
Clib.tstep   = 1;
Clib.temp    = 298;
Clib.ESoverlap = true;
Clib.ESeps     = 0.1;

nangles = [50 100];
Vgs = [0.3 1];
betaES = [-20];
beta1 = [1];
tstep = [1 10 0.2 0.05];
nsteps = [200000 20000 600000 1000000];
%ESeps = [0.1 0.2 0.3 0.4 0.5 0.7 0.9 1.0];
ESeps = [0.1 0.2 0.3 0.4 0.5];


C = TrajConfig;
C.ESoverlap     = true;
C.ESforces      = 1;
C.optWidth      = 0;
C.periodic      = true;
nsave = [0 1000 0 0 1];

lib = Library('data3\mobility100f','perES');

%% Generate data

for i1 = 1:size(nangles,2)
for i2 = 1:size(Vgs,2)
for i3 = 1:size(betaES,2)
for i4 = 1:size(beta1,2)
for i5 = 1:size(tstep,2)
for i6 = 1:size(ESeps,2)

Clib.nangles       = nangles(i1);
Clib.Vgs           = Vgs(i2);
Clib.betaES        = betaES(i3);
Clib.beta1         = beta1(i4);
Clib.tstep         = tstep(i5);
Clib.ESeps         = ESeps(i6);

C.nangles          = Clib.nangles;
C.Vgs              = Clib.Vgs;
C.betaES           = Clib.betaES;
C.beta1            = Clib.beta1;
C.tstep            = Clib.tstep;
C.temp             = Clib.temp;
C.ESeps            = Clib.ESeps;

found = lib.find(Clib);
nfound = size(found,1);
if (nfound > 0)   
disp(['angles ',num2str(C.nangles), ...
   ' Vgs ', num2str(C.Vgs), ...
   ' betaEs ', num2str(C.betaES), ...
   ' beta1 ', num2str(C.beta1), ...
   ' tstep ', num2str(C.tstep), ...
   ' found'
   ]);
else
disp(['angles ',num2str(C.nangles), ...
   ' Vgs ', num2str(C.Vgs), ...
   ' betaEs ', num2str(C.betaES), ...
   ' beta1 ', num2str(C.beta1), ...
   ' tstep ', num2str(C.tstep), ...
   ' not found'
   ]);   
t1 = TrajSegment(C,nsave,0,0);
tic;
olaps = t1.runTraj( nsteps(i5) );
toc;
lib.store(Clib, t1);
end

end
end
end
end
end
end
%% load data

for i1 = 1:size(nangles,2)
for i2 = 1:size(Vgs,2)
for i3 = 1:size(betaES,2);
for i4 = 1:size(beta1,2);
for i5 = 1:size(tstep,2);

Clib.nangles       = nangles(i1);
Clib.Vgs           = Vgs(i2);
Clib.betaES        = betaES(i3);
Clib.beta1         = beta1(i4);
Clib.tstep         = tstep(i5);


disp(['angles ',num2str(Clib.nangles), ...
   ' Vgs ', num2str(Clib.Vgs), ...
   ' betaEs ', num2str(Clib.betaES), ...
   ' beta1 ', num2str(Clib.beta1), ...
   ' tstep ', num2str(Clib.tstep), ...
   ]);
ifound = lib.find(Clib,false);

t1 = lib.retrieve(ifound);
c1 = t1.data('cent');
d1 = diff(c1);
dtest = [d1+2*pi; d1; d1-2*pi];
d2 = min(abs(dtest),[],1);
hist(d2,200);
input junk;
end
end
end
end
end

%% analysis

for i1 = 1:size(nangles,2)
for i2 = 1:size(Vgs,2)
for i3 = 1:size(betaES,2)
for i4 = 1:size(beta1,2)
for i5 = 1:size(tstep,2)
for i6 = 1:size(ESeps,2)

Clib.nangles       = nangles(i1);
Clib.Vgs           = Vgs(i2);
Clib.betaES        = betaES(i3);
Clib.beta1         = beta1(i4);
Clib.tstep         = tstep(i5);
Clib.ESeps         = ESeps(i6);

ifound = lib.find(Clib,false);

% disp(['angles ',num2str(Clib.nangles), ...
%    ' Vgs ', num2str(Clib.Vgs), ...
%    ' betaEs ', num2str(Clib.betaES), ...
%    ' beta1 ', num2str(Clib.beta1), ...
%    ' tstep ', num2str(Clib.tstep), ...
%    ' ESeps ', num2str(Clib.ESeps), ...
%    ' found ', num2str(ifound), ...
%    ]);
if (~isempty(ifound))
t1 = lib.retrieve(ifound);
c1 = t1.data('cent');
c1 = c1(1,:);
d1 = diff(c1);
dtest = [d1+2*pi; d1; d1-2*pi];
d2 = min(abs(dtest),[],1);
% d2 is in radians, need to convert to unit cells
d2 = d2 .* t1.C.nangles/(2*pi);

if (Clib.tstep < 1)
   disp('>');
   nsteps = 1/Clib.tstep;
   d3 = [];
   currentPoint = 1;
   while((currentPoint+nsteps) < length(d2))
      d3(end+1) = sum(d2(currentPoint:(currentPoint+nsteps-1)));
      currentPoint = currentPoint+nsteps;
   end
   avgx2 = sum(d3.^2)/size(d3,2);
   mob = (1.6e-19 * (4e-8).^2 /2/1.38e-23/298/48.8e-15) * avgx2/...
      (Clib.tstep * nsteps);   
else
   avgx2 = sum(d2.^2)/size(d2,2);
   mob = (1.6e-19 * (4e-8).^2 /2/1.38e-23/298/48.8e-15) * avgx2/Clib.tstep;
end
%disp([' mobility = ',num2str(mob)]);
disp(['angles ',num2str(Clib.nangles), ...
   ' Vgs ', num2str(Clib.Vgs), ...
   ' betaEs ', num2str(Clib.betaES), ...
   ' beta1 ', num2str(Clib.beta1), ...
   ' tstep ', num2str(Clib.tstep), ...
   ' ESeps ', num2str(Clib.ESeps), ...
   ' found ', num2str(ifound), ...
   ' avgx2 ', num2str(avgx2), ...
   ' mobility ',num2str(mob), ...
   ]);
end
%hist(d2,200);

%input junk;
end
end
end
end
end
end

%% How does x2 change with summation of steps ?
%ESeps = [0.1 0.2 0.3 0.4 0.5];
col = {'r','y','g','b','k'};
%tstep = [1 10 0.2 0.05];
sym = {'x','.','^','o'};
for i1 = 1:size(ESeps,2)
   for i2 = 1:size(tstep,2)
   
Clib.nangles       = 100;
Clib.Vgs           = 0.3;
Clib.betaES        = -20;
Clib.beta1         = 1;
Clib.tstep         = tstep(i2);
Clib.ESeps         = ESeps(i1);

ifound = lib.find(Clib,false);
t1 = lib.retrieve(ifound);
angles = t1.data('cent');
%y1 = c1(2,:);
%state = floor(y1);
%olap  = y1 - state;
%olap = olap(2:size(olap,2));
loc = perToLinear(t1.C.nangles, angles);

plotx = zeros(100,1);
ploty = zeros(100,1);
ic = 1;
for time = 1:100
   bin = floor(time/t1.C.tstep);
   plotx(ic,1) = Clib.tstep * bin; % * 48.8/1000;
   ploty(ic,1) = getMobility(t1.C, loc, bin);
   ic = ic+1;
end
figure(400 + i1);
hold on;
plot(plotx,ploty,[col{i1},sym{i2}]);
end
end

%% Since results are converged by limit of above plot (time = 100), look at
% this point only
nangles = [50 100];
Vgs = [0.3 1];
betaES = [-20];
beta1 = [1];
tstep = [0.05 0.2 1 10];
ESeps = [0.1 0.2 0.3 0.4 0.5];

mob = zeros(size(nangles,2), size(ESeps,2), size(tstep,2));

for i1 = 1:size(nangles,2)
lib = Library(['mobility',num2str(nangles(i1)),'f'],'perES');
   
for i2 = 1:size(ESeps,2)
for i3 = 1:size(tstep,2)
   
Clib.nangles       = nangles(i1);
Clib.Vgs           = 0.3;
Clib.betaES        = -20;
Clib.beta1         = 1;
Clib.tstep         = tstep(i3);
Clib.ESeps         = ESeps(i2);

ifound = lib.find(Clib,false);
t1 = lib.retrieve(ifound);
angles = t1.data('cent');
loc = perToLinear(t1.C.nangles, angles);
time = 100;
bin = floor(time/t1.C.tstep);
mob(i1,i2,i3) = getMobility(t1.C, loc, bin);
end
end
end
%% plot mob 
nangles = [50 100];
sym = {'x','o'};
tstep = [0.05 0.2 1 10];
col = {'r','k','b','g'};
ESeps = [0.1 0.2 0.3 0.4 0.5];

for i1 = 1:size(nangles,2)
for i2 = 1:size(ESeps,2)
for i3 = 1:size(tstep,2)
   figure(10);
   hold on;
   plot(ESeps(i2), mob(i1,i2,i3),[col{i3},sym{i1}]);
end
end
end

%% What is happening with the cutoffs?
for nangles = 50 % 50:50:100
Clib.nangles       = nangles;
Clib.Vgs           = 1.0;
Clib.betaES        = -20;
Clib.beta1         = 1;
Clib.tstep         = 0.05;

lib = Library(['data3\mobility',num2str(Clib.nangles),'f'],'perES');

if Clib.nangles == 100
   sym = '^';
   ESeps = [0.1 0.2 0.3 0.4 0.5];
else
   sym = '>';
   ESeps = [0.1 0.2 0.3 0.4 0.5];
end

col = {'b','r','k','b','g','c'};

for i1 = 1 %:size(ESeps,2)
Clib.ESeps         = ESeps(i1);

ifound = lib.find(Clib,false);
if (ifound == [])
   disp('data not found');
end
t1 = lib.retrieve(ifound);
c1 = t1.data('cent');
angles = c1(1,:);
% flag is stored in center(2,:)
y1 = c1(2,:);
% integer component is the state selected
state = floor(y1);
% remaining portion is the overlap
olap  = y1 - state;
% dist(n) holds loc(n+1) - loc(n)
% so state(n+1) and olap(n+1) are appropriate
olap = olap(2:size(olap,2));
state = state(2:size(state,2));
loc = perToLinear(t1.C.nangles, angles);
dist = loc(1,2:size(loc,2)) - loc(1,1:(size(loc,2)-1));

%[~,ist1] = find(state < 1.5);
%olap = olap(ist1);
%dist = dist(ist1);

% examine distances as a function of overlap
step = 0.025;
ols = Clib.ESeps:step:(1.0-step);  % available data starts at ESeps
nplot = size(ols,2);
xplot = zeros(nplot,1);
yplot = zeros(nplot,1);
yplot2 = zeros(nplot,1);
mu = zeros(nplot,1);
i=0;
for ol = ols
   i = i+1;
   xplot(i) = ol+(step/2);
   % select points within a range of overlap
   [~,ind1] = find((ol<olap) & (olap<(ol+step)));
   % find average x^2 for steps with this overlap
   yplot(i) = mean(dist(ind1).^2);
   yplot2(i) = size(ind1,2);
end
% get fraction of steps with this overlap
yplot2n = yplot2/sum(yplot2);
for i2 = 1:nplot
   % fraction of steps, assuming we keep states with overlaps > ols(i2)
   ytemp = yplot2(i2:nplot)/sum(yplot2(i2:nplot));
   % average x2 for this cut off
   avgx2 = sum(ytemp.*yplot(i2:nplot));
   time1 = Clib.tstep;
   mu(i2) = (1.6e-19 * (4e-8).^2 /2/1.38e-23/298/48.8e-15) * avgx2/time1;
end
figure(10);
hold on;
plot(xplot,yplot,[col{i1},sym,'-']);
figure(11);
hold on;
plot(xplot,log10(yplot2n),[col{i1},sym,'-']);
figure(12);
hold on;
plot(xplot,yplot.*yplot2n,[col{i1},sym,'-']);
figure(13);
hold on;
plot(xplot,mu,[col{i1},sym,'-']);
end
end

figure(10)
title('avg x2 versus overlap');
figure(11)
title('log(# steps) versus overlap');
figure(12);
title('contribution to x2 versus overlap');
figure(13);
title('mobility as function of overlap cutoff');