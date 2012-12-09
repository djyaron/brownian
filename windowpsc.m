function windowpsc

% First attempt at the window approach to mobility 
%clear classes;
dataroot = '/brashear/yaron/matdl/browniantest';
% Create default library structure for window mobility calcs

% Variables to be looped over
nangles = 100; % [50 100];
Vgs = 0.3; %[0.3 1];
betaES = [-10];
beta1 = [1];
tstep = [1 0.2];% [1 10 0.2 0.05];
nsteps = [50 100]; %[1e6 5e6]; %[200000 20000 600000 1000000];
wsize = 3:20; %[3 4 5 6 7 8 9 10 11 12 15 20];
nruns = 15;

nsave = [0 1000 10 0 10];
nener = 3; % save ground state and first two excited states
nwf = 0;



% Save type of calcs so we can 
calcsTodo = cell(0,0);
for i1 = 1:length(nangles)
for i2 = 1:length(Vgs)
for i3 = 1:length(betaES)
for i4 = 1:length(beta1)
for i5 = 1:length(tstep)
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

a.clib = Clib;
if (i5 == 1)
   a.nsteps = nsteps(i5);
   a.nsave = [0 1000 10 0 10];
else
   a.nsteps = nsteps(i5);
   a.nsave = [0 5000 50 0 50];
end
calcsTodo{end+1} = a;

end
end
end
end
end
end
%%
parfor icalc = 1:length(calcsTodo)

Clib = calcsTodo{icalc}.clib;
   

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
% found = lib.find(Clib);
% nfound = size(found,1);
% if (nfound >= nruns)   
% disp(['angles ',num2str(C.nangles), ...
%    ' Vgs ', num2str(C.Vgs), ...
%    ' betaEs ', num2str(C.betaES), ...
%    ' beta1 ', num2str(C.beta1), ...
%    ' tstep ', num2str(C.tstep), ...
%    ' wsize ', num2str(C.wsize), ...
%    ' found'
%    ]);
% else
%       disp([...
%          ' angles ',num2str(C.nangles), ...
%          ' Vgs ', num2str(C.Vgs), ...
%          ' betaEs ', num2str(C.betaES), ...
%          ' beta1 ', num2str(C.beta1), ...
%          ' tstep ', num2str(C.tstep), ...
%          ' wsize ', num2str(C.wsize), ...
%          ' found ', num2str(nfound) ...
%          ]);

lib = Library([dataroot,'/gs',num2str(Clib.Vgs),...
   '/es',num2str(abs(Clib.betaES)),'/w',num2str(Clib.wsize)] ...
   ,'wind');
nfound = 0;
   for irun = (nfound+1):nruns
      t1 = TrajSegment(C,calcsTodo{icalc}.nsave,nener,nwf);
      tic;
      olaps = t1.runTraj( calcsTodo{icalc}.nsteps );
      toc;
      lib.store(Clib, t1);
   end
end


end