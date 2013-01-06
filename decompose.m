load('mobilitySave.mat');
ndata = length(resSave);
fitme = FitFunc(resSave{ndata}.w,resSave{ndata}.mob);

% pguess = [1 5 0.1];
% for idata = 1% :(ndata-1)
%    d1 = resSave{idata};
%    x = fitme.doFit(d1.mob, pguess);
%    mreal(idata) = x(1);
%    plot(d1.w,d1.mob,'ro-');
%    hold on;
%    plot(d1.w,fitme.pred(x),'bo-');
% end

close all;
pguess = [5 0.07 0.0005];
for idata = 1% 1:ndata
   d1 = resSave{idata};
   x = fitme.doFit(d1.mob, pguess);
   parab = fitme.parabToPlot(x);
   rest = d1.mob - parab;
   figure(1)
   hold on;
   plot(d1.w,d1.mob,'ko-');
   figure(1)
   hold on;
   plot(d1.w,rest,'bo-');
end
