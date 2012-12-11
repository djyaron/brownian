classdef Analysis
   methods (Static)
      function res = takeDiff(x,ndiff)
         % for n = 1000 and ndiff = 1
         %  x(2:1000) - x(1:999)
         % for n = 1000 and ndiff = 2
         %  x(3:1000) - x(1:998)
         res = x((ndiff+1):end) - x(1:(end-ndiff));
      end
      function res = removePeriodic(x,nangles)
         % takes a list of centers, in terms of unit cell numbers
         % that lie between 0 and nangles, and returns a list of 
         % real displacements
         % 1) calculate raw differences
         d1 = Analysis.takeDiff(x,1);
         % 2) for differences <-nangles or >nangles, 
         %    set correct difference
         i1 = find(abs(d1+nangles) < abs(d1));
         d1(i1) = d1(i1) + nangles;
         i1 = find(abs(d1-nangles) < abs(d1));
         d1(i1) = d1(i1) - nangles;
         % 3) reconstruct trajectory from differences 
         res = cumsum(d1);
      end
      function [mu, t] = muFromC(t1,sumLengths,ignore)
         % calculates mobility from a trajSegment t1
         % it ignores the first "ignore" time steps
         % sumLengths is a list of times, used to get the 
         % rms(x2) where x2 = (x(t+sumLength)-x(t)).
         x = Analysis.removePeriodic(t1.cent(1,ignore:end),t1.nangles);
         mu = zeros(size(sumLengths));
         t = zeros(size(sumLengths));
         ic = 1;
         for iwind = sumLengths
            d1 = Analysis.takeDiff(x,iwind);
            time1 = t1.C.tstep * iwind;
            avgx2 = mean(d1.^2);         
            mu(ic) = (1.6e-19 * (4e-8).^2 /2/1.38e-23/298/48.8e-15) ...
               * avgx2/time1;
            t(ic) = time1;
            ic = ic+1;
         end
      end
      function [xp,yp] = jumpDistribution(trajs,sumlength,ignore)
         jumps = [];
         for i = 1:length(trajs)
            t1 = trajs{i};
            x = Analysis.removePeriodic(t1.cent(1,ignore:end),t1.nangles);
            d1 = Analysis.takeDiff(x,sumlength);
            jumps = [jumps,d1(:)'];
         end
         [yp,xp] = hist(jumps,50);
         yp = yp/max(yp);
      end
   end
   
   
end

