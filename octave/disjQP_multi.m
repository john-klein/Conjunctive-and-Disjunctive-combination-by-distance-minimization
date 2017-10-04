function [mout, OBJ, INFO, LAMBDA]=disjQP_multi(mass,dist)
	
% This function computes a disjunctive combination of two input mass functions by minimizing a distance.
% Distances are computed between the total conflict mass function and all functions belonging to those less comitted than both input mass functions.
% The minimization is guaranteed to converge to a unique solution.
% inputs:
%  - mass : a M x N matrix containing M mass functions. Each line of the matrix is an input mass functions. N is the size of the power set.
%  - dist : a string indicating which distance to use for approximation. Possible choices are 'pl' and 'q'.
%    Note that here k=2 since we will be using quadratic programming.
%
% outputs:
%  - mout : the result of the combination of m1 with m2. mout > m1 and m2 for some partial order linked with the chosen distance.
%  - OBJ, INFO, LAMBDA : some convergence information of the quadratic programming, see function qp for more information on these parameters.

if (nargin~=2)
  error('missing input parameters.')
end

[M,N]=size(mass);
n = log(N)/log(2);

if ((n-floor(n))>0)
  error('the size of mass functions must be a power of 2.')
end

if (abs(sum(sum(mass')==ones(1,M)-1))>1e-10)
  error('matrix mass does not contain mass functions.')
end

%Building the incidence matrix M.
M=[1 1 ; 0 1];
if (n>=2)
  for i=2:n
    M=kron([1 1 ; 0 1],M);
  end
end 

%Initial mass function in the Conditional Subspace.
m0=zeros(N,1);
m0(N)=1;

%Total conflict mass function
mconf=zeros(N,1);
mconf(1)=1;

%Vaccuous mass function
migno=zeros(N,1);
migno(N)=1;


%Bounds for mass functions
mass_lb=zeros(N,1);
mass_ub=ones(N,1);

%Selecting a distance
if ((strcmp(dist,'pl')==1))
  %Let us compute the optimum in the plausibility space
  %The  matrix 1-J*M' maps mass functions to plausibility functions.
  J=fliplr(eye(N));
  A=1-J*(M');
  %There is no reverse transfer matrix for plausibilities
  %Instead revA will compute mass function values from a plausibility for non empty subsets.
  R=repmat(migno',N,1);
  revA=inv(M')*J*(R-eye(N));
  %Boundary conditions for plausibilities (open world assumption)
  lower_bound=max((A*mass')');
  upper_bound=ones(N,1);
  upper_bound(1)=0;
  %Vectors are valid plausibility functions if applying revA yields a valid mass function for non empty sets: m(A) in [0;1] for A non empty.
  %We also need that m(emptyset) is in [0;1]. This condition is ensured as long as pl(Omega) is in [0;1].
elseif (strcmp(dist,'q')==1)
  %The matrix M maps mass functions to commonality functions.
  A=M;
  %Reverse passage matrix is revA
  revA=inv(M);
  %Boundary conditions for commonalities
  lower_bound=max((A*mass')');
  upper_bound=ones(N,1);
  %Vectors are valid commonality functions if applying revA yields a valid mass function

else
  error('undefined distance.')
end

%Calling solver
[mout, OBJ, INFO, LAMBDA] =qp (m0, A'*A, -A'*A*mconf, mass_ub', 1, mass_lb, mass_ub, lower_bound, A, upper_bound);
%Erasing very small masses dued to computation noise
mout=mout.*(mout>1e-10);

