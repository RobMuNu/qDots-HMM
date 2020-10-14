function logprob=fa_log(o, trans, ipi)
% Forward algorithm in log probabilities to avoid underflow issues
% 
% INPUTS:
% o 			= 	Given emitted sequence labelled in numerics e.g. [1 2 2 1 2 1]
% ipi 			= 	initial probability of states in a row vector
% trans(i,j,k) 	= 	(N x N x |alphabet|) sized transition matrix 
% 					(N x N x each possible emitted symbol). 
% 					trans(:,:,1) contain transitions emitting a '1', etc.
% 					make sure symbols used in 'o' and 'trans' are consistent.
%
% OUTPUTS:
% logprob		= 	log(P(o|trans))
%
% NOTES:
% The only 'assumption' this code makes is that the hmm's fed into it is ergodic
% This means that the stationary distribution will have no asymptotic zero elements
%
% Optimisation of this code is left as an exercise to the reader 


% ---------------------------
%            Setup
% ---------------------------
% Length of observed sequence
T=length(o); 
% logPi row vector containing all log pi elements
% eMs are ergodic so pi should never have states that are never visited
logPi = log(ipi);
% Log elements of transition tensor
logTrans = log(trans); % !! Will incur -Inf elements !!

% ---------------------------
%       Initialisation
% ---------------------------
% Define log-forward variable
b(1,:) = logPi;

% ---------------------------
%         Recursion
% ---------------------------
for t = 2:(T+1)

	for j = 1:length(trans(1,:,1))
		% make vector of logsum exponent arguments
		d = b(t-1,:)' + logTrans(:,j,o(t-1));  
		% Remove -Inf elements
		d = d(find(d~=(-Inf)));

		% If no logsumexp args exist then no valid paths exist which end on state j
		% hence write -Inf to the forward element which will be removed in the next pass
		if numel(d)==0
			b(t,j) = -Inf;
		else
			% Update forward variable with logsumexp trick
			b(t,j) = max(d) + log(sum(exp(d - max(d))));
		end
	
	end
end 

% ---------------------------
%          Summation
% ---------------------------
% Do a quick check to see if path is impossible or not
% Path is impossible if all the final forward elements are -Inf
if unique(b(T+1,:)) == -Inf
	disp(['Path given by data is impossible under model ' inputname(2)])
	logprob = -Inf;
else
	logprob = max(b(T+1,:)) + log(sum(exp(b(T+1,:) - max(b(T+1,:)))));
end




