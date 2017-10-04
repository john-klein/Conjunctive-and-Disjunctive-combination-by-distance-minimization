function [mout]=conjQP_multi(mass,dist)
	
% This functions computes a conjunctive combination of two input mass functions by minimizing a distance.
% Distances are computed between the vaccuous mass function and all functions belonging to those more comitted than both input mass functions.
% The minimization is guaranteed to converge to a unique solution.
% inputs:
%  - m1 : a mass function to be combined 
%  - m2 : a mass function to be combined
%  - dist : a string indicating which distance to use for approximation. Possible choices are 'pl', 'b' and 'q'.
%    Note that here k=2 since we will be using quadratic programming.
%
% outputs:
%  - mout : the result of the combination of m1 with m2. mout < m1 and m2 for some partial order linked with the chosen distance.


if (license('test', 'optimization_toolbox')==0)
    error('missing optimization toolbox. The function quadprog belongs to the optimization toolbox.')
end

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
m0(1)=1;

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
  lower_bound=zeros(N,1);
  upper_bound=min((A*mass')');
  upper_bound(1)=0;
  %Vectors are valid plausibility functions if applying revA yields a valid mass function for non empty sets: m(A) in [0;1] for A non empty.
  %We also need that m(emptyset) is in [0;1]. This condition is ensured as long as pl(Omega) is in [0;1].
elseif ((strcmp(dist,'b')==1))
  %Let us compute this optimum in the implicability space
  %The transpose of matrix M maps mass functions to implicability functions.
  A=M';
  %Reverse passage matrix is revA
  revA=inv(M');
  %Boundary conditions for implicabilities
  lower_bound= max((A*mass')');
  lower_bound(N)=1;
  upper_bound=ones(N,1);
  %Calling solver
  %Vectors are valid implicability functions if applying revA yields a valid mass function
elseif (strcmp(dist,'q')==1)
  %The matrix M maps mass functions to commonality functions.
  A=M;
  %Reverse passage matrix is revA
  revA=inv(M);
  %Boundary conditions for commonalities
  lower_bound=zeros(N,1);
  lower_bound(1)=1;
  upper_bound=min((A*mass')');
  %Vectors are valid commonality functions if applying revA yields a valid mass function
else
  error('undefined distance.')
end

%Calling solver
options = optimoptions('quadprog',...
    'Algorithm','interior-point-convex','MaxIterations',1000,'OptimalityTolerance',1e-14,'ConstraintTolerance',1e-14,'Display','off');
[mout] =quadprog ( A'*A, -A'*A*migno, A, upper_bound, mass_ub', 1, mass_lb, mass_ub, m0,options);

%Erasing very small masses dued to computation noise
mout=mout.*(mout>1e-10);

