function nLLg=nllgHMM(datH0,tH0,tH1,piH0,piH1,a)
% Finds the LRT threshold g for two HMM hypotheses H0 and H1 
% given significance level a. Uses log-forward algorithm to calculate LRT.
%
% INPUTS:
% datH0 = an NxM array of H0 data. M trials of length N
% tH0 = transition tensor of H0
% tH1 = transition tensor of H1
% a = the desired significance level
%
% OUTPUTS:
% nLLg = the nLLRT threshold (-2*log(gamma)) for the HMMs at sig level a
% 
% NOTES:
% The computation of the negative log likelihood will take some time so sit tight


if (a*100) < ((1/length(datH0(1,:)))*100)
	disp('Significance level likely to be too small given number of indep. data sets...')
%	disp('Press any key to continue...')
%	pause;
end


% Compute negative Log-likelihood nLL = -2*log(LRT)
% Vpa used for arbitrary precision
nLL = zeros(1,length(datH0(1,:)));
wBar = waitbar(0,'Computing all log-likelihoods...');
for i = 1:length(datH0(1,:))
	nLL(i) = double(-2*log(exp(vpa(fa_log(datH0(:,i),tH1,piH1)))/...
			exp(vpa(fa_log(datH0(:,i),tH0,piH0)))));
	waitbar(i/length(datH0(1,:)),wBar)
end
close(wBar);


% Check edge cases for full impossible paths
% Impossible H1 ==> nLL = +Inf
% I.e. if H1 is always impossible then the null should be
% infinitely more likely. The converse should never happen for thresholding
if length(unique(nLL)) == 1 & unique(nLL) == Inf
	nLLg = Inf; % Strictly speaking the threshold is undefined


else
	% Tabulate nLL results
	nLLTab = tabulate(nLL);

	% Check to see if all the log-likelihood results are unique
	if length(unique(nLLTab(:,3)))==1

		% Check to see if the significance level is too big
		if (a*100) < unique(nLLTab(:,3))
			disp(sprintf('Not enough trials to achieve sig level %g',a))
			disp(sprintf('Tightest significance available: %g',unique(nLLTab(:,3))/100))
			disp('Running Log-likelihood at new sig level...')
			% Update significance level
			a = unique(nLLTab(:,3))/100;
		end

		% Since it's a flat distribution then all we need to do is index
		disp('[ Method ] Direct Indexing')
		nLLg = nLLTab(floor(a*length(datH0(1,:))),1);


	% If results are not unique then we'll have to do a sum loop
	else
		
		% THIS IS SO BROKEN
		% FOR A PROBABILITY SUM LOOP THE IF STATEMENT BELOW MAKES NO SENSE
		% First check to see if the significance level is too big
		if (a*100) < unique(nLLTab(:,3))
			disp(sprintf('Not enough trials to achieve sig level %g',a))
			disp(sprintf('Tightest significance available: %g',unique(nLLTab(:,3))/100))
			disp('Running Log-likelihood at new sig level...')
			% Update significance level
			a = unique(nLLTab(:,3))/100;
		end

		disp('[ Method ] Probability Sum Loop')
		sumProb = 0;
		iProb = 0;
		while 1
			iProb = iProb + 1;
			sumProb = sumProb + (nLLTab(iProb,3)/100);
			if sumProb > a
				% Ugly fix to stop overcounting
				% I think this breaks if it overcounts on the first shot
				% and when all results are in 1 bin (isnt this rare af?)
				% apparently it happens :S
				% This can happen during stress testing where the length of
				% data is very short. 
				sumProb = sumProb - (nLLTab(iProb,3)/100); 
				iProb = iProb - 1;
				break;
			end
		end
		nLLg = nLLTab(iProb,1);

	end
end

disp(sprintf('Results:\n\tnLL Threshold: %g\n\tSignificance level: %g',nLLg,a))









