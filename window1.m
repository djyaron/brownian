%% First attempt at the window approach to mobility 
%clear classes;
dataroot = 'c:\matdl\brownian';
% Create default library structure for window mobility calcs
Clib.type    = 'wind1';
Clib.nangles = 100;
Clib.Vgs     = 1;
Clib.betaES  = -20;
Clib.beta1   = 1;
Clib.tstep   = 1;
Clib.temp    = 298;
Clib.wsize   = 9;

% Variables to be looped over
nangles = 100; % [50 100];
Vgs = 1; %[0.3 1];
betaES = [-20];
beta1 = [1];
tstep = 1;% [1 10 0.2 0.05];
nsteps = 2e5; %[200000 20000 600000 1000000];
wsize = 8;
nruns = 2;

C = TrajConfig;
C.ESoverlap     = true;
C.ESforces      = 1;
C.optWidth      = 0;
C.periodic      = true;
nsave = [0 1000 1 0 1];
nener = 3; % save ground state and first two excited states
nwf = 0;

lib = Library(dataroot,'wind1');

%% Generate data

for i1 = 1:length(nangles)
for i2 = 1:length(Vgs)
for i3 = 1:length(betaES)
for i4 = 1:length(beta1)
for i5 = 1:length(tstep)
for i6 = 1:length(wsize)

% Structure used to label the data for library storage/retreival   
Clib.nangles       = nangles(i1);
Clib.Vgs           = Vgs(i2);
Clib.betaES        = betaES(i3);
Clib.beta1         = beta1(i4);
Clib.tstep         = tstep(i5);
Clib.wsize         = wsize(i6);

% Copy the configuration variables into the structure used to configure the
% simulation
C.nangles          = Clib.nangles;
C.Vgs              = Clib.Vgs;
C.betaES           = Clib.betaES;
C.beta1            = Clib.beta1;
C.tstep            = Clib.tstep;
C.temp             = Clib.temp;
C.wsize            = Clib.wsize;

% See how many steps are in the library
found = lib.find(Clib);
nfound = size(found,1);
if (nfound >= nruns)   
disp(['angles ',num2str(C.nangles), ...
   ' Vgs ', num2str(C.Vgs), ...
   ' betaEs ', num2str(C.betaES), ...
   ' beta1 ', num2str(C.beta1), ...
   ' tstep ', num2str(C.tstep), ...
   ' wsize ', num2str(C.wsize), ...
   ' found'
   ]);
else
   disp(['angles ',num2str(C.nangles), ...
      ' Vgs ', num2str(C.Vgs), ...
      ' betaEs ', num2str(C.betaES), ...
      ' beta1 ', num2str(C.beta1), ...
      ' tstep ', num2str(C.tstep), ...
      ' wsize ', num2str(C.wsize), ...
      ' not found'
      ]);
   t1 = TrajSegment(C,nsave,0,0);
   t1.nener = nener;
   t1.nwf = nwf;
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

%%
i = 1;
nplot = 1000;
x = Analysis.removePeriodic(t1.cent(1,:),t1.nangles);
while 1
   plot(x(i:(i+nplot)),'r.');
   i = i+nplot;
   input junk;
end
%%  Plots of a trajectory (assuming nener = nwf = 2)
time1 = t1.time('ener','fs');
figure(100)
clf;
plot(time1,t1.ener(1,:),'k.');
hold on;
plot(time1,t1.ener(2,:),'r.');
figure(200)
plot(time1,t1.cent(1,:),'b.');
title('center');
figure(205)
plot(time1,t1.cent(2,:),'g.');
title('wf width');
disp(['average wf width = ',num2str(mean(t1.cent(2,:)))]);
% for i = 1:size(t1.wf,3)
%    figure(300)
%    plot(t1.wf(:,1,i));
%    figure(301)
%    plot(t1.wf(:,2,i));
%    input junk;
% end


%% load data

for i1 = 1:length(nangles)
for i2 = 1:length(Vgs)
for i3 = 1:length(betaES)
for i4 = 1:length(beta1)
for i5 = 1:length(tstep)
for i6 = 1:length(wsize)

% Structure used to label the data for library storage/retreival   
Clib.nangles       = nangles(i1);
Clib.Vgs           = Vgs(i2);
Clib.betaES        = betaES(i3);
Clib.beta1         = beta1(i4);
Clib.tstep         = tstep(i5);
Clib.wsize         = wsize(i6);

disp(['angles ',num2str(Clib.nangles), ...
   ' Vgs ', num2str(Clib.Vgs), ...
   ' betaEs ', num2str(Clib.betaES), ...
   ' beta1 ', num2str(Clib.beta1), ...
   ' tstep ', num2str(Clib.tstep), ...
   ]);

ifound = lib.find(Clib);
figure(100)
clf;
for ifile = ifound(:)'
   t1 = lib.retrieve(ifile);
   sumLengths = [1:1000:100000];
   [mu, t] = Analysis.muFromC(t1,sumLengths);
   t = t * 48.8/1000;
   figure(100);
   hold on;
   plot(t,mu,'r.')
end
xlabel('time (ps)');
ylabel('mobility');

end
end
end
end
end
end



