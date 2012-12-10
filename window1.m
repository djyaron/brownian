%% First attempt at the window approach to mobility 
clear classes;
dataroot = 'c:\matdl\brownian';
% Create default library structure for window mobility calcs

% Variables to be looped over
nangles = 100; % [50 100];
Vgs = 1; %[0.3 1];
betaES = [-20];
beta1 = [1];
tstep = 1;% [1 10 0.2 0.05];
nsteps = 1e6; %[200000 20000 600000 1000000];
wsize = 3:20; %[3 4 5 6 7 8 9 10 11 12 15 20];
nruns = 15;

nsave = [0 1000 10 0 10];
nener = 3; % save ground state and first two excited states
nwf = 0;

%% Generate data

for i1 = 1:length(nangles)
for i2 = 1:length(Vgs)
for i3 = 1:length(betaES)
for i4 = 1:length(beta1)
for i5 = 1:length(tstep)
parfor i6 = 1:length(wsize)

% Structure used to label the data for library storage/retreival   
Clib = [];
Clib.type = 'wind1';
Clib.nangles       = nangles(i1);
Clib.Vgs           = Vgs(i2);
Clib.betaES        = betaES(i3);
Clib.beta1         = beta1(i4);
Clib.tstep         = tstep(i5);
Clib.temp          = 298;
Clib.wsize         = wsize(i6);

lib = Library([dataroot,'\c',num2str(Clib.nangles),'\w',num2str(Clib.wsize)] ...
   ,'wind');

% Copy the configuration variables into the structure used to configure the
% simulation
C = TrajConfig;
C.ESoverlap     = false;
C.ESforces      = 1;
C.optWidth      = 0;
C.periodic      = true;
C.nangles          = Clib.nangles;
C.Vgs              = Clib.Vgs;
C.betaES           = Clib.betaES;
C.beta1            = Clib.beta1;
C.tstep            = Clib.tstep;
C.temp             = Clib.temp;
C.wsize            = Clib.wsize;

% See how many files are in the library
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
      disp([...
         ' angles ',num2str(C.nangles), ...
         ' Vgs ', num2str(C.Vgs), ...
         ' betaEs ', num2str(C.betaES), ...
         ' beta1 ', num2str(C.beta1), ...
         ' tstep ', num2str(C.tstep), ...
         ' wsize ', num2str(C.wsize), ...
         ' found ', num2str(nfound) ...
         ]);
   for irun = (nfound+1):nruns
      t1 = TrajSegment(C,nsave,nener,nwf);
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
end

%%
% i = 1;
% nplot = 1000;
% x = Analysis.removePeriodic(t1.cent(1,:),t1.nangles);
% while 1
%    plot(x(i:(i+nplot)),'r.');
%    i = i+nplot;
%    input junk;
% end
% %%  Plots of a trajectory (assuming nener = nwf = 2)
% time1 = t1.time('ener','fs');
% figure(100)
% clf;
% plot(time1,t1.ener(1,:),'k.');
% hold on;
% plot(time1,t1.ener(2,:),'r.');
% figure(200)
% plot(time1,t1.cent(1,:),'b.');
% title('center');
% figure(205)
% plot(time1,t1.cent(2,:),'g.');
% title('wf width');
% disp(['average wf width = ',num2str(mean(t1.cent(2,:)))]);
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
   xw = []; yw=[]; sdw=[];
   figure(500); clf;
   figure(100); clf;
   figure(600); clf;
for i6 = 1:length(wsize)

% Structure used to label the data for library storage/retreival  
Clib = [];
Clib.type = 'wind1';
Clib.nangles       = nangles(i1);
Clib.Vgs           = Vgs(i2);
Clib.betaES        = betaES(i3);
Clib.beta1         = beta1(i4);
Clib.tstep         = tstep(i5);
Clib.temp          = 298;
Clib.wsize         = wsize(i6);

lib = Library([dataroot,'\c',num2str(Clib.nangles),'\w',num2str(Clib.wsize)] ...
   ,'wind');

ifound = lib.find(Clib);
disp(['angles ',num2str(Clib.nangles), ...
   ' Vgs ', num2str(Clib.Vgs), ...
   ' betaEs ', num2str(Clib.betaES), ...
   ' beta1 ', num2str(Clib.beta1), ...
   ' tstep ', num2str(Clib.tstep), ...
   ' wsize ', num2str(Clib.wsize),...
   ' found ',num2str(length(ifound)) ...
   ]);

%sumLengths = [1:1000:100000];
sumLengths = [1:5:200];
% mu(files, timewindow)
mu = zeros(length(ifound),length(sumLengths));
nhistBins = 400;
ignore = 5000;
egapx = zeros(length(wsize),nhistBins);
egapy = zeros(length(wsize),nhistBins);
ic = 0;
gaps = [];
wfwidth = [];
for ifile = ifound(:)'
   ic = ic+1;
   t1 = lib.retrieve(ifile);
   gaps = [gaps, t1.ener(3,ignore:end) - t1.ener(2,ignore:end)];
   wfwidth = [wfwidth,t1.cent(2,ignore:end)];
   [mu(ifile,:), t] = Analysis.muFromC(t1,sumLengths,ignore);
   figure(100);
   hold on;
   t = t * 48.8/1000;
   plot(t,mu(ifile,:),'r.')
end
% Plot gap between E1 and E2
[egapy(i6,:),egapx(i6,:)] = hist(gaps,nhistBins);
egapy(i6,:) = egapy(i6,:)/max(egapy(i6,:));
figure(500);
hold on;
plot(egapy(i6,:)+Clib.wsize,egapx(i6,:),'b-');
[wfwy,wfwx] = hist(wfwidth,nhistBins);
wfwy = wfwy/max(wfwy);
figure(600);
hold on;
plot(wfwy+Clib.wsize,wfwx,'c-');

avgMu = mean(mu);
sdMu = std(mu)/sqrt(length(ifound));

xw(end+1) = Clib.wsize;
yw(end+1) = avgMu(1,end);
sdw(end+1) = sdMu(1,end);

figure(200);
hold on;
errorbar(t,avgMu,sdMu,'k')
   
xlabel('time (ps)');
ylabel('mobility');


end
figure(300);
errorbar(xw,yw,sdw);
end
end
end
end
end

%%
% Test trajectory, to make sure we know what we are doing.
trs = cell(0,0);
for nangles = [50 100]
   for stepSize = [2.0 4.0]
%nangles = 100;
%stepSize = 2.0; % gaussian random steps of this magnitude each time step
nstep = 1e6;

tr.nangles = nangles;
tr.C.tstep = 1;
tr.t = 1:nstep;
c = zeros(1,nstep);
rnum = stepSize * randn(nstep,1);

c(1) = 0;
for i = 1:(nstep-1)
   c(i+1) = c(i) + rnum(i);
   if (c(i+1) > nangles)
      c(i+1) = c(i+1) - nangles;
   end
   if (c(i+1) < 0)
      c(i+1) = c(i+1) + nangles;
   end
end
tr.cent = c;
trs{end+1} = tr;
   end
end
%%
for it = 1:length(trs)
sumlengths = 1:100;
ignore = 200;
[mu, t] = Analysis.muFromC(trs{it},sumLengths,ignore);
figure(55)
hold on;
plot(t,mu);
end

%% What is the depenence if we just move randomly within the window
% Test trajectory, to make sure we know what we are doing.
trs = cell(0,0);
winds = [];
for wind = 3:50
%nangles = 100;
%stepSize = 2.0; % gaussian random steps of this magnitude each time step
nstep = 1e6;
winds(end+1) = wind;

tr.nangles = nangles;
tr.C.tstep = 1;
tr.t = 1:nstep;
c = zeros(1,nstep);
rnum = wind * (-1 + 2 * rand(nstep,1));

c(1) = 0;
for i = 1:(nstep-1)
   c(i+1) = c(i) + rnum(i);
   if (c(i+1) > nangles)
      c(i+1) = c(i+1) - nangles;
   end
   if (c(i+1) < 0)
      c(i+1) = c(i+1) + nangles;
   end
end
tr.cent = c;
trs{end+1} = tr;

end
%%
muWind = [];
for it = 1:length(trs)
sumlengths = 1:100;
ignore = 200;
[mu, t] = Analysis.muFromC(trs{it},sumLengths,ignore);
figure(55)
hold on;
plot(t,mu);
muWind(end+1) = mu(end);
end
figure(56)
plot(winds,muWind);
