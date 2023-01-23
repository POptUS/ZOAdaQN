function [y, ysmooth] = calfun_batch(x,sigma,batch,problem)
%     This is a modified version of the subroutine calfun.m
%     available at
%     https://github.com/POptUS/BenDFO
%
% Inputs:
%       x 	array of length n
%	sigma  	scalar defines the standard deviation of noise
% 	batch	consists of set of indices
%           
% Outputs:
%       y 	```: stochastic function value
%       ysmooth: smooth deterministic function value
%
%     Additional problem descriptors are passed through the global
%     variables:
%       problem.m a positive integer (length of output from dfovec).
%          m must not exceed n.
%       problem.nprob is a positive integer that defines the number of the problem.
%          nprob must not exceed 22.
%       problem.probtype is a string specifying the type of problem desired:
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

%global m nprob probtype fvals nfev np

m = problem.m;
nprob = problem.nprob;
probtype = problem.probtype;
n = size(x,1); % Problem dimension
S = length(batch); % batch size
if nprob<100
    fvec = dfovec(m,n,x,nprob);
else
    fvec = mghvec(m,n,x,nprob);
    %fvec = mghvec(m,n,x,nprob-100); 
end

% Calculate the function value
switch probtype
    case 'reluniform' % Didn't finish coding uniform
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
        end
        if nargout>1
            ysmooth = sum(fvec.^2);
        end
        z=zeros(m,S);
        for i=1:S
            randn('seed',batch(i));
            z(:,i) = sigma*randn(m,1);
        end
        fvec = fvec.*(1+z);
        y = sum(fvec.^2);
        y=(1/(1 + sigma^2))*(sum(y)/length(y));
    case 'absuniform' % Didn't finish coding uniform
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
        end
        if nargout>1
            ysmooth = sum(fvec.^2);
        end
        z=zeros(m,S);
        for i=1:S
            randn('seed',batch(i));
            z(:,i) = sigma*randn(m,1);
        end
	%elseif nargin >=3
	%	randn('seed',seed)
	%end
    %    z = sigma*randn(m,1);
        fvec = fvec+z;
        y = sum(fvec.^2);
        y=(sum(y)/length(y)) - m*sigma^2;
    case 'smooth'
        y = sum(fvec.^2);
end

% Update the function value history
%nfev = nfev +1;
%fvals(nfev,np) = y;

% Optional truncation:
%if y>1e64
%  display('Function value exceeds 10^64')
%end
