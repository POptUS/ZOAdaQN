
function y = calfun_sample(x,sigma,seed)
%     This is a modified version of the subroutine calfun.m
%     available at
%     https://github.com/POptUS/BenDFO
%
% Inputs:
%       x 	array of length n
%	sigma  	(optional) scalar defines the standard deviation of noise
% 	seed	(optional) scalar the defines the rand/randn seed used
% Outputs:
%       f 	scalar function value at x
%
%
%     Additional problem descriptors are passed through the global
%     variables:
%       m a positive integer (length of output from dfovec).
%          m must not exceed n.
%       nprob is a positive integer that defines the number of the problem.
%          nprob must not exceed 22.
%       probtype is a string specifying the type of problem desired:
%           'smooth' corresponds to smooth (noise-free) problems
%           'absuniform' corresponds to stochastic uniform absolute noise
%           'absnormal' corresponds to stochastic Gaussian absolute noise
%           'reluniform' corresponds to stochastic uniform relative noise
%           'relnormal' corresponds to stochastic Gaussian relative noise
%	**Note: the noise is applied independently to each component before
%		the components are squared and summed, additional variance
%		control will necessarily need to account for the value m	
%
%
%     To store the evaluation history, additional variables are passed 
%     through global variables. These may be commented out if a user 
%     desires. They are:
%       nfev is a non-negative integer containing the number of function 
%          evaluations done so far (nfev=0 is a good default).
%          after calling calfun, nfev will be incremented by one.
%       np is a counter for the test problem number. np=1 is a good
%          default if only a single problem/run will be done.
%       fvals is a matrix containing the history of function
%          values, the entry fvals(nfev+1,np) being updated here.
%

global m nprob probtype fvals nfev np

n = size(x,1); % Problem dimension

% Generate the vector
fvec = dfovec(m,n,x,nprob); 

% Calculate the function value
switch probtype
    case 'reluniform'
	if nargin<2
	        sigma=10^-3;
	elseif nargin >=3
		rand('seed',seed)
	end
        z = sigma*sqrt(3)*(-ones(m,1)+2*rand(m,1));
        fvec = fvec.*(1+z);
        y = sum(fvec.^2);
    case 'relnormal'
	if nargin<2
	        sigma=10^-3;
	elseif nargin >=3
		randn('seed',seed)
	end
        z = sigma*randn(m,1);
        fvec = fvec.*(1+z);
        y = sum(fvec.^2);
    case 'absuniform'
	if nargin<2
	        sigma=10^-3;
	elseif nargin >=3
		rand('seed',seed)
	end
        z = sigma*sqrt(3)*(-ones(m,1)+2*rand(m,1));
        fvec = fvec+z;
        y = sum(fvec.^2);
    case 'absnormal'
	if nargin<2
	        sigma=10^-3;
	elseif nargin >=3
		randn('seed',seed)
	end
        z = sigma*randn(m,1);
        fvec = fvec+z;
        y = sum(fvec.^2);
    case 'smooth'
        y = sum(fvec.^2);
end

% Update the function value history
nfev = nfev +1;
fvals(nfev,np) = y;

% Optional truncation:
%if y>1e64
%  display('Function value exceeds 10^64')
%end
