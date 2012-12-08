function [forces,E,c,wf,flag] = forcesFromESwind(beta,angles, ...
   cen,wsize)
% Window version of the ES calculation
% INPUT: beta = coupling between adjacent rings when planar (should be <0)
%       angles = list of all angles of the rings
%    The following two are optimization variables. We first diagonalize
%    the subblock around the cen position, and then see if charge is
%    localized here.
%         cen = previous center of charge
%         wsize = diagonalize matrix from cen-window...cen+window
% Output: F = force on each ring
%         E = excitation energies (lowest is the state used for forces)
%         c = wavefunction of state used for the force calculation
%         wf = excited state wavefunctions (ordered as in E)
%         flag = width of the wavefunction

% from angles, we know the length of the chain (i.e # angles)
flag=0;
nangles= length(angles);
% Set up the entire Hamiltonian
hamil= zeros(nangles,nangles);
for i=1:nangles-1
   hamil(i,i+1) = beta*cos(angles(i+1)-angles(i));
   hamil(i+1,i) = hamil(i,i+1);
end
% because periodic, we couple ring 1 and N
hamil(1,nangles) = beta*cos(angles(nangles)-angles(1));
hamil(nangles,1) = hamil(1,nangles);

% pull out the range of cells around the current location
r1 = getRange(round(cen),nangles,wsize);
Hsub = hamil(r1,r1);
% diagonalize, and sort eigenvalues (sort may not be needed, but matlab
% just states they take what lapack returns, and not sure lapack will 
% always call dsyev
[V,D] = eig(Hsub);
[E order] = sort(diag(D));
Vsub = V(:,order);
% I now have the wavefunctions in the range r1, stored in Vsub
% I'll copy this into the correct location of V(nangles, :)
% V will then hold the wavefunctions on the range 1..nangles
% but I only have a subset of states calculated
wf = zeros(nangles,length(E));
wf(r1,:) = Vsub;
   
% Analytic derivatives of the charge energy
c = wf(:,1);  % use lowest state wavefunction
flag = 1/sum(c.^4);
forces = zeros(nangles,1);
forces(1) = 2*beta*sin(angles(1)-angles(nangles)) *c(1)*c(nangles) ...
   -2*beta*sin(angles(2)-angles(1))*c(1)*c(2);
for i=2:nangles-1
   forces(i) = 2*beta*sin(angles(i)-angles(i-1))*c(i)*c(i-1) ...
      -2*beta*sin(angles(i+1)-angles(i))*c(i)*c(i+1);
end
forces(nangles) = 2*beta*sin(angles(nangles)-angles(nangles-1))...
   *c(nangles)*c(nangles-1) ...
   -2*beta*sin(angles(1)-angles(nangles))*c(1)*c(nangles);

end

function res = getRange(cen,nangles,wsize)
% pull out the sub-block corresponding to range cen-window... cen+window
r1 = (cen-wsize):(cen+wsize);
% displace values that are > nangles, down by nangles
tooHigh = find(r1>nangles);
r1(tooHigh) = r1(tooHigh) - nangles;
tooLow = find(r1<1);
r1(tooLow) = r1(tooLow) + nangles;

res = r1;
end
