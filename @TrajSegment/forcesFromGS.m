function [forcesGS, E] = forcesFromGS(V,angles,periodic)

%      [F, E]    = forcesFromGS(V, angles)

% INPUT:  V      = Rotational barrier in ground state.
%         angles = List of all angles of the rings.
%         periodic = true for periodic boundary conditions [default is
%           false]
% OUTPUT: F      = Forces on each ring.
%         E      = Total energy of ground state, relative to planar.
% DEBUG CONDITION: if periodic == 2, then use harmonic potential

% from angles, we know the length of the chain (i.e # angles)
if (nargin < 3)
   periodic = false;
end

[nangles,junk]= size(angles);

forcesGS      = zeros(nangles,1);
E=0;

if (periodic == 2)
   for iangle = 2:nangles
      E = E + (4*V/pi^2)*(angles(iangle)-angles(iangle-1))^2;
   end
   
   forcesGS(1,1) = +(8*V/pi^2)*(angles(2)-angles(1));
   for iangle = 2:(nangles-1)
      forcesGS(iangle,1) = -(8*V/pi^2)*(angles(iangle)-angles(iangle-1)) ...
         +(8*V/pi^2)*(angles(iangle+1)-angles(iangle));
   end
   forcesGS(nangles,1) = -(8*V/pi^2)*(angles(nangles)-angles(nangles-1));
   
else % regular cosine potential
if (V ~= 0)
   if (periodic)
      forcesGS(1,1) = -V* sin(2*(angles(1)-angles(nangles))) ...
         + V*sin(2*(angles(2)-angles(1)));
   else
      forcesGS(1,1) = V*sin(2*(angles(2)-angles(1)));
   end
   for iangle = 2:(nangles-1)
      forcesGS(iangle,1) = -V*sin(2*(angles(iangle)-angles(iangle-1))) ...
         +V*sin(2*(angles(iangle+1)-angles(iangle)));
   end
   if (periodic)
      forcesGS(nangles,1) = -V*sin(2*(angles(nangles)-angles(nangles-1))) ...
         + V*sin(2*(angles(1)-angles(nangles)));
   else
      forcesGS(nangles,1) = -V*sin(2*(angles(nangles)-angles(nangles-1)));
   end
   for iangle = 2:nangles
      E = E + (V/2.0)*(1-cos(2*(angles(iangle)-angles(iangle-1))));
   end
   if (periodic)
      E = E + (V/2.0)*(1-cos(2*(angles(nangles)-angles(1))));
   end
end

end