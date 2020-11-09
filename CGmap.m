% Script that makes a heatmap of the fraction of datasets
% where CSSR is preferred over ARP 
% as a function of CG bins and max memory length



% Run in the directory that has all the gviz files as well as the CG data


for x = 2:10
	% For all quantum dots @ CG bin size
	strBaseExtension = sprintf('*-EXb%02d',x);

	% 29 datafiles coarse grained under this CG bin size
	datFilePattern = fullfile(pwd,strBaseExtension);
	datList = dir(datFilePattern);


	for dSet = 1:length(datList)

		% For a specific quantum dot @ CG bin size
		strDotName = datList(dSet).name;
		gvFilePattern = fullfile(pwd, horzcat(strDotName,'*.dot'));
		gvList = dir(gvFilePattern);


		% Construct the ARP for the dot @ this CG base



		for lam = 1:length(gvList)
			% For a specific quantum dot @ CG bin size & CSSR lambda
			strCSSRName = gvList(lam).name;

			% Get CSSR model for the dot @ this CG base & lambda





		end


	end 

end


% Locally defined functions 
function datNum = importDat(datName)
% IMPORTDAT Local function that imports a coarse-grained datafile
% and translates it into a numerical rather than alphabetical sequence
	dat = fileread(datName);

	% Grab interleaved data emission statistics
    datTable = tabulate(dat'); % transpose for row vector data
    datFreqs = cell2mat(datTable(:,2)); % get letter frequencies
    datProbs = datFreqs./length(dat); % get letter probabilities
    datTable(:,4) = num2cell(datProbs);
    [~, datSrt] = sort([datTable{:,1}],'ascend');
    datSrtTab = datTable(datSrt,:);
    datAlpha = [datSrtTab{:,1}]; % This ensures that the emission alphabet is sorted
    datProbs = [datSrtTab{:,4}];
   
    % The forward algorithm won't recognize alphanumeric sequences
    % We will need to convert the path to numerical emissions
    datRange = 1:length(datAlpha); 
    datNum = zeros(1,length(dat));
    for i = 1:length(datAlpha)
        datNum(strfind(dat,datAlpha(i))) = datRange(i);
    end   
end


function datAlpha = getAlphabet(datName)
	dat = fileread(datName);

	% Grab interleaved data emission statistics
    datTable = tabulate(dat'); % transpose for row vector data
    datFreqs = cell2mat(datTable(:,2)); % get letter frequencies
    datProbs = datFreqs./length(dat); % get letter probabilities
    datTable(:,4) = num2cell(datProbs);
    [~, datSrt] = sort([datTable{:,1}],'ascend');
    datSrtTab = datTable(datSrt,:);
    datAlpha = [datSrtTab{:,1}]; % This ensures that the emission alphabet is sorted
end


function [ttDOT, piDOT, tmDOT] = dot_to_transition2(filename,alphabet,doFuzz)
	% Exctracts a transition tensor, and the asymptotics of the flattened matrix
	% from a GraphViz file http://www.research.att.com/sw/tools/graphviz
	%
	% INPUT:  'filename' - the file in DOT format containing the graph layout
	%                       where graph is a unifilar HMM
	%                       Labels should satisfy <symbol:probability> format
	%           'alphabet' - string vector of all the CG2k alphabet
	%                       symbols expressed in the HMM (same as alphabet
	%                       file used by CSSR to generate DOT file
	%           'doFuzz' - flag which asks the function to make a "fuzzy"
	%                       transition tensor. This writes a very small
	%                       probability into zero elements of the transition
	%                       tensor which ensures full sequence support while
	%                       maintaining some semblance of structure.
	%
	% OUTPUT:   'ttDOT' - the transition tensor of the unifilar HMM 
	%           'piDOT' - the vector of asymptotics simulated from the HMM
	%           'tmDOT' - flattened transition tensor (the transition matrix)

	%
	% NOTES: Only tested on directed, ergodic & unifilar graphs.
	%       Only tested on DOT output generated from CSSR algorithm
	%       Assumes some metadata in the DOT file is always structured
	%        the same way 
	%       The rows of the transition matrix pulled from DOT don't always
	%       Sum to 1. This is usually due to rounding errors, but is rectified
	%       via the dtmc function. 

	% Read file into cell array of lines ignoring c style comments
	dotFile = textread(filename,'%s','delimiter','\n','commentstyle','c');  
	dot_lines = strvcat(dotFile);    

	 % Check for valid DOT file
	if isempty(strfind(dot_lines(1,:), 'graph '))
	   error('* * * File does not appear to be in valid DOT format. * * *');
	   return;
	end;

	% Cleanup
	% Remove metadata present from file header and footer
	dot_lines = dot_lines(7:end-1,:);

	% Get number of states
	nStates = split(dot_lines(end,:),' ');
	nStates = str2num(nStates{1})+1;
	nAlpha = length(alphabet);

	% Initialise transition tensor and asymptotic vector
	ttDOT = zeros(nStates,nStates,nAlpha);
	piDOT = zeros(1,nStates);


	% Pull labels and fill matrix
	for lidx = 1:length(dot_lines(:,1))
	    % Split each line into a 12x1 cell array of space delimited values 
	    dotline = split(dot_lines(lidx,:),' '); 

	    % Write in the transition tensor based on the following predictable fields
	    % {1} FROM
	    % {3} TO
	    % {6} Symbol emitted
	    % {7} Transition probability

	    aTemp = erase(dotline{6},':');
	    aTemp = erase(aTemp,'"');
	    if isempty(strfind(alphabet,aTemp))
	        error('No matching DOT symbol to alphabet. Check supplied alphabet same as the one which generated the dot file.');
	        return;
	    end
	    ttDOT(str2double(dotline{1})+1,...
	        str2double(dotline{3})+1,...
	        strfind(alphabet,aTemp)) = str2double(dotline{7});
	    
	end


	% Do fuzzy version?
	if doFuzz == 0
	    % Evaluate state asymptotics from transition matrix
	    % Flatten transition tensor across alphabet dimension
	    tmDOT = zeros(nStates,nStates);
	    tmDOT = sum(ttDOT,3);
	    dtmcDOT = dtmc(tmDOT);
	    piDOT = asymptotics(dtmcDOT);

	elseif doFuzz == 1

	    % Fill in zero transitions with a small error
	    % Let's just make it 1/100th of the smallest existing nonzero transition
	    ttDOT(find(ttDOT==0)) = min(ttDOT(ttDOT>0))/100;
	 

	    % Evaluate state asymptotics from transition matrix
	    % Flatten transition tensor across alphabet dimension
	    tmDOT = zeros(nStates,nStates);
	    tmDOT = sum(ttDOT,3);
	    dtmcDOT = dtmc(tmDOT);
	    piDOT = asymptotics(dtmcDOT);


	else 
	    error('doFuzz should be 0 or 1')
	    return;
	end



	disp('dot_to_transition2: Check there are no CG2k alphabet letters skipped!')
end


function [ttPLARP, piPLARP] = makeARP(zeros(2,2,max(length(datAlphaOn),length(datAlphaOff))));
    % remember symbol 1 = off, symbol 2 = on
    % Only need to fill in the off-diagonals
    % For convention set A->B = On (matrix element 1,2)


    % Note for waiting time HMMs, only the rows of the emission-flattened matrix sum to 1
    ttPLARP(1,2,datRangeOn) = datProbsOn; % Fill A->B (On waiting times)
    ttPLARP(2,1,datRangeOff) = datProbsOff;  % Fill B->A (Off waiting times)
    % !!! This fails if alphabet letters are skipped cf. line 63 !!!


    % Check if on and off aphabets are asymmetrical to make fuzzy ARP
    % Fill shorter state trailing zeros with a micro transition
    % Let's just make it 1/100th of the smallest existing nonzero transition
    if length(datAlphaOn) < length(datAlphaOff)
        % ON shorter, need to fill in tt(1,2) with small micro transition
        ttPLARP(1,2,find(ttPLARP(1,2,:)==0)) = min(ttPLARP(ttPLARP>0))/100;
        % Renormalise
        ttPLARP(1,2,:) = ttPLARP(1,2,:)./sum(ttPLARP(1,2,:));
        
    elseif length(datAlphaOn) > length(datAlphaOff)
        % OFF shorter, need to fill in tt(2,1)
        ttPLARP(2,1,find(ttPLARP(2,1,:)==0)) = min(ttPLARP(ttPLARP>0))/100;
        % Renormalise
        ttPLARP(2,1,:) = ttPLARP(2,1,:)./sum(ttPLARP(2,1,:));

    end


    % Remember interleaved data begins at ON by construction, 
    % so pi should be (1,0) for the forward algorithm
    piPLARP = [1,0]; % Since A->B implies On transition


    % Generate PL-ARP data
    % Doing this as two alternating biased dice
    % The fast way is to index every ODD element right away and sample
    % then index every EVEN element right away and sample
    simPLARP = zeros(length(dat),nDatSets);

    distPLARPon = [ttPLARP(1,2,:)];
    distPLARPoff = [ttPLARP(2,1,:)];
    % Shave off trailing zeros on the distribution if they exist
    % (Shouldn't be the case if we do fuzzy ARP)
    distPLARPon = distPLARPon(1:find(distPLARPon,1,'last')); 
    distPLARPoff = distPLARPoff(1:find(distPLARPoff,1,'last'));
    % Need to define these for the PLARP echo 
    distCGPLon = distPLARPon;
    distCGPLoff = distPLARPoff;

    for n = 1:nDatSets
        % Sample On wait times
        simPLARP(1:2:end,n) = randsample(...
            1:max(max(datRangeOn),max(datRangeOff)),...
            length(dat(1:2:end)),...
            true,...
            distPLARPon);

        % Sample Off wait times
        simPLARP(2:2:end,n) = randsample(...
            1:max(max(datRangeOn),max(datRangeOff)),...
            length(dat(2:2:end)),...
            true,...
            distPLARPoff);
    end
    clear n

end


function nLL=nllHMM(dat,tH0,tH1,piH0,piH1)

	nLL = double(-2*log(exp(vpa(fa_log(dat,tH1,piH1)))/...
			exp(vpa(fa_log(dat,tH0,piH0)))));

end


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
		% Strictly speaking the threshold is undefined but this assignment will
		% produce the right results
		nLLg = Inf;


	else
		% Tabulate nLL results
		nLLTab = tabulate(nLL);

		% Check to see if all the log-likelihood results are unique
		if length(unique(nLLTab(:,3)))==1

			% Check to see if the significance level is too tight for number of trials
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
			
			% First check to see if the significance level is too big
			% In this instance we only really need the significance to not overshoot the
			% first bin

			if (a*100) < unique(nLLTab(1,3)) 

				
				disp(sprintf('Not enough trials to achieve sig level %g',a))
				disp(sprintf('Tightest significance available: %g',unique(nLLTab(1,3))/100))
				disp('Running Log-likelihood at new sig level...')
				% Update significance level
				a = unique(nLLTab(1,3))/100; 
			end

			disp('[ Method ] Probability Sum Loop')
			sumProb = 0;
			iProb = 0;
			while 1
				iProb = iProb + 1;
				sumProb = sumProb + (nLLTab(iProb,3)/100);
				if sumProb > a
					% Ugly fix to stop overcounting
					sumProb = sumProb - (nLLTab(iProb,3)/100); 
					iProb = iProb - 1;
					break;
				end
			end
			nLLg = nLLTab(iProb,1);

		end
	end

	disp(sprintf('Results:\n\tnLL Threshold: %g\n\tSignificance level: %g',nLLg,a))
end















