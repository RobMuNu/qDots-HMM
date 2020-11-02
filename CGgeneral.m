function output=CGgeneral(dat,x,method,overflow)
% Description:
% Converts raw photon trajectory dataset into a compressed
% (interleaved, leading ON) waiting-time dataset.
% Compression sends raw waiting times into intervals defined
% by exponential, or powerlaw rules:
% EXP:			[x^k, x^(k+1)) (k the variable, x the chosen base)
% PL:			[x^k, (x+1)^k) (x the variable, k the chosen power)


% INPUTS:		dat - numerical, binary (0,1) photon trajectory. 1 line.
%				x - base or exponent
%				method - 'exp' or 'pl' compression rules
%				overflow - which bin number to make the overflow

% OUTPUT:		Compressed (interleaved, leading ON) dataset


% Check inputs ok
if method ~= 'exp' & method ~= 'pl'
	error("CG method must be either 'exp' or 'pl'");
	return;

elseif overflow < 0 
	error('Overflow must be 0 or positive integer');
	return;

elseif method == 'exp' & x <= 1
	error("Exponent base 'x' must be positive integer > 1 for 'exp' method");
	return;

elseif method == 'pl' & x <= 0 
	error("Power 'x' must be positive integer for 'pl' method");
	return;

end


% Perform compression methods
if method == 'exp'
	if overflow == 0
		
		% EXP bin index rule for matlab
		% k(w) = floor(log_x(w))+1

	elseif overflow > 0
		% check whether its greater than the number of bins that would have 
		% been made automatically (hard to do?)

	end

elseif method == 'pl'
	if overflow == 0

		% PL bin index rule for matlab
		% x(w) = floor(w^(1/k))

	elseif overflow > 0 
		% check whether its greater than the number of bins that would have 
		% been made automatically (hard to do?)

	end

end



	% Function to 
  indexRule = @(w) floor(log2(w))+1;
  zipfPMF = @(x,a) x.^(-a) ./ genHarm(waitMax,a);


	% Function to calculate log of any base
	 function logbase = logb(x,base)
	 	logbase = log(x)./log(base);
	 end



	 % Func. to read binary trajectory and convert to wait times
	 function waitDat = traj2wait(traj)
	 	% traj must already be a numerical, binary dataset! 1 line!
        % No idea what these do
        findOne = find(diff([0,traj,0]==1)); % For one
        findZero = find(diff([1,traj,1]==0)); % For zero
        % Find starting indices of *blocks* of ones and zeros
        idxStartOne = findOne(1:2:end-1);  % Starting indices of 1's blocks
        idxStartZero = findZero(1:2:end-1);  % Start indices
        % The generalised waiting times
        waitOne = findOne(2:2:end)-idxStartOne;  % Consecutive ones? counts
        waitZero = findZero(2:2:end)-idxStartZero;  % Consecutive zeros? counts

        % Now we want to interleave on /off waiting times together
        % Initialise
        waitDat = zeros(1,length(waitOne)+length(waitZero));
        % First check what the starting traj symbol is
        if traj(1) == 1
        	% Then interleave leading with waitOne
        	waitDat((1:length(waitOne))*2 - 1) = waitOne;
        	waitDat((1:length(waitZero))*2) = waitZero;

        else
        	% Interleave leading with waitZero and remove the first waiting time
        	% That way we know the waiting times began in the ON state
        	waitDat((1:length(waitZero))*2 - 1) = waitZero;
        	waitDat((1:length(waitOne))*2) = waitOne;
        	waitDat(1) = [];
        end
	 end






