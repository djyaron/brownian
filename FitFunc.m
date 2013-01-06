classdef FitFunc < handle
   properties
      mobility;
      w     % list of window sizes
      m100  % mobility versus w for betaES = 100
   end
   methods
      function res = FitFunc(w,m100)
         res.w = w;
         res.m100 = m100;
      end
      function res = doFit(obj,mobVersusW,pguess)
         obj.mobility = mobVersusW;
         [x,resnorm,residual,exitflag,output] =...
            lsqnonlin(@obj.toMinimize2,pguess);
         res = x;
      end
      function res = pred1(obj,pars)
         w100Multiplier = pars(1);
         x0 = pars(2);
         curvature = pars(3);
         res = w100Multiplier .* obj.m100 + ...
            curvature .* heaviside(obj.w - x0) .* ...
            (obj.w - x0).^2;
      end
      function res = toMinimize(obj,pars)
         res = obj.mobility - obj.pred(pars);
      end
      function res = toMinimize2(obj,par)
         x0 = par(1);
         y0 = par(2);
         curvature = par(3);
         shiftDown = obj.mobility - y0;
         shiftDown(shiftDown<0.0) = 0.0;
         res = shiftDown - obj.parab(par);
      end
      function res = parab(obj,par)
         x0 = par(1);
         curvature = par(3);
         res = curvature .* heaviside(obj.w - x0) .* ...
            (obj.w - x0).^2;
      end
      function res = parabToPlot(obj,par)
         x0 = par(1);
         y0 = par(2);
         curvature = par(3);
         res =y0 + curvature .* heaviside(obj.w - x0) .* ...
            (obj.w - x0).^2;
      end
   end
end