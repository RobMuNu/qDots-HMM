function output = CGgeneral(dat,x,method,overflow)
% Description:
% Converts raw photon trajectory dataset into a compressed
% (interleaved, leading ON) waiting-time dataset.
% Compression sends raw waiting times into intervals defined
% by exponential, or powerlaw rules:
% EXP:			[x^k, x^(k+1)) (k the variable, x the chosen base)
% PL:			[a^x, (a+1)^x) (a the variable, x the chosen power)


% INPUTS:		dat - numerical, binary (0,1) photon trajectory. 1 line.
%				x - base or exponent
%				method - 'ex' or 'pl' compression rules
%				overflow - which bin number to make the overflow

% OUTPUT:		Compressed (interleaved, leading ON) dataset 
%				alphanumeric string of compressed alphabet symbols


	% Check inputs ok
	if ~strcmpi(method,'ex') & ~strcmpi(method,'pl')
		error("CG method must be either 'ex' or 'pl'");
		return;

	elseif overflow < 0 
		error('Overflow must be 0 or positive integer');
		return;

	elseif strcmpi(method,'ex') & x <= 1
		error("Exponent base 'x' must be positive integer > 1 for 'ex' method");
		return;

	elseif strcmpi(method,'pl') & x <= 0 
		error("Power 'x' must be positive integer for 'pl' method");
		return;

	end

	% Perform compression methods
	alphaFull = 'a':'z';
	if overflow > 1
		% Brainless implementation of overflow bin
		% Will work for powerlaw or exponential bin intervals
		% as long as you don't have absurdly long waiting times
		% If the overflow bin is set to a waiting time greater than
		% what is present in the data, then this will do nothing to the output
		alphaFull(overflow+1:end) = alphaFull(overflow);
	end

	if method == 'ex'	
		% Waiting time to EXP bin index rule for matlab
		% k(w) = floor(log_x(w))+1
		[waitDat,~] = traj2wait(dat);
		CGexpBins = floor(logbase(waitDat,x))+1;
		output{1,1} = sprintf('Compressed data under %s method\n',method);
		output{1,2} = alphaFull(CGexpBins);


	elseif method == 'pl'
		% Waiting time to PL bin index rule for matlab
		% x(w) = floor(w^(1/k))
		[waitDat,~] = traj2wait(dat);
		CGplBins = floor(waitDat.^(1/x));
		output{1,1} = sprintf('Compressed data under %s method\n',method);
		output{1,2} = alphaFull(CGplBins);

	end

	% Nested functions needed for this to work
	% Function to calculate log of any base
	 function logb = logbase(x,base)
	 	logb = log10(x)./log10(base);
	 end

	 % Func. to read binary trajectory and convert to wait times
	 function [waitDat,waitMax] = traj2wait(traj)
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

        % Find the longest waiting time
        waitMax = max(max(waitOne),max(waitZero));
	 end

end







