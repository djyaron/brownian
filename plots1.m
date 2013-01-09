function plots1()
dataroot = 'c:\matdl\brownian';
Clib.type = 'wind1';
Clib.nangles = 100; 
Clib.Vgs = 0; 
Clib.betaES = -70; 
Clib.beta1 = 1; 
Clib.tstep = 1;
Clib.temp = 298;
Clib.wsize = 10;

wsize = [10 20];
sumLengths = [200];
ignore = 5000;
nruns = 7;
genData = 1;

nsave = [0 1000 10 0 10];
nener = 3; % save ground state and first two excited states
nwf = 0;
nsteps = 1e6; 

xtype = 'betaES';
xvals = [-30:-10:-60]; %[-30 -40 -50 -60];
linetype = 'Vgs';
lvals = [0.0 0.1 0.3];


%lt = {'ro','go','bo','ko','co','r^','g^','b^','k^','c^'};
lt = {'r','g','b','k','c','r','g','b','k','c'};
legnd = cell(0,0);
close all;

for iline = 1:length(lvals)
   xp = []; yp=[]; ypw=[];
   Clib = setfield(Clib,linetype,lvals(iline));
   for ix = 1:length(xvals)
      xp(end+1) = xvals(ix);
      Clib = setfield(Clib,xtype,xvals(ix));
      ysave = zeros(1,length(wsize));
      for iw = 1:length(wsize)
         Clib.wsize = wsize(iw);
         lib = Library([dataroot,'\gs',num2str(Clib.Vgs),...
            '\es',num2str(abs(Clib.betaES)),'\w',num2str(Clib.wsize)] ...
            ,'wind');
         
         ifound = lib.find(Clib);

         nfound = size(ifound,1);
         disp(['angles ',num2str(Clib.nangles), ...
            ' Vgs ', num2str(Clib.Vgs), ...
            ' betaEs ', num2str(Clib.betaES), ...
            ' beta1 ', num2str(Clib.beta1), ...
            ' tstep ', num2str(Clib.tstep), ...
            ' wsize ', num2str(Clib.wsize), ...
            ' found ', num2str(nfound)
            ]);
         
         % mu(files, timewindow)
         mu = zeros(length(ifound),length(sumLengths));
         ic = 0;
         trajs = cell(0,0);
         for ifile = ifound(:)'
            ic = ic+1;
            t1 = lib.retrieve(ifile);
            [mu(ic,:), t] = Analysis.muFromC(t1,sumLengths,ignore);
         end
         
         ysave(iw) = mean(mu);
         %sdMu = std(mu)/sqrt(length(ifound));
      end
      yp(end+1) = mean(ysave);
      ypw(end+1) = (max(ysave)-min(ysave))/2;
   end
   figure(300);
   hold on;
   errorbar(xp,yp,ypw,lt{iline});
   legnd{end+1} = [linetype,' = ',num2str(lvals(iline))];
end
figure(300);
xlabel(xtype);
ylabel('mobility');
legend(legnd);


end

%          if ((nfound < nruns) && genData)
%             fprintf(1,'%s','generating data ');
%             C = TrajConfig;
%             C.ESoverlap     = false;
%             C.ESforces      = 1;
%             C.optWidth      = 0;
%             C.periodic      = true;
%             C.nangles          = Clib.nangles;
%             C.Vgs              = Clib.Vgs;
%             C.betaES           = Clib.betaES;
%             C.beta1            = Clib.beta1;
%             C.tstep            = Clib.tstep;
%             C.temp             = Clib.temp;
%             C.wsize            = Clib.wsize;
%             
%             todo = nruns-nfound;
%             tsave = cell(todo,1);
%             parfor irun = 1:todo
%                fprintf(1,'%i ',irun);
%                nsave1 = nsave/Clib.tstep;
%                t1 = TrajSegment(C,nsave1,nener,nwf);
%                olaps = t1.runTraj( nsteps );
%                tsave(irun) = t1;
%             end
%             for irun = 1:todo
%                lib.store(Clib, t1);
%             end
% 
%             ifound = lib.find(Clib);
%             
%             nfound = size(ifound,1);
%             disp(['try again: angles ',num2str(Clib.nangles), ...
%                ' Vgs ', num2str(Clib.Vgs), ...
%                ' betaEs ', num2str(Clib.betaES), ...
%                ' beta1 ', num2str(Clib.beta1), ...
%                ' tstep ', num2str(Clib.tstep), ...
%                ' wsize ', num2str(Clib.wsize), ...
%                ' found ', num2str(nfound)
%                ]);
%          end


