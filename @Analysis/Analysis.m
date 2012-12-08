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
         % take differences
         d1 = Analysis.takeDiff(x,1);
         i1 = find(abs(d1+nangles) < abs(d1));
         d1(i1) = d1(i1) + nangles;
         i1 = find(abs(d1-nangles) < abs(d1));
         d1(i1) = d1(i1) - nangles;
         % apply periodic conditions
         res = cumsum(d1);
      end
      function [mu, t] = muFromC(t1,sumLengths)
         x = Analysis.removePeriodic(t1.cent(1,:),t1.nangles);
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
   end
end

