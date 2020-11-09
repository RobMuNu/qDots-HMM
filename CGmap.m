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
		alphaFilePattern = fullfile(pwd,...
							horzcat(erase(strDotName,...
							erase(strBaseExtension,'*')),...
							'-ALPHA'));
		alphaVec = dir(alphaFilePattern); % Should only ever index 1 file 
		gvFilePattern = fullfile(pwd, horzcat(strDotName,'*.dot'));
		gvList = dir(gvFilePattern);


		% Construct and simulate the ARP for this data
		[ttARP, piARP, simARP] = makeARP(strDotName,1000);


		for lam = 1:length(gvList)
			% For a specific quantum dot @ CG bin size & CSSR lambda
			strCSSRName = gvList(lam).name;

			% Read the correct line from the alphabet vector
			fid = fopen(alphaVec.name);
			cellAline = textscan(fid,'%s',1,'delimiter','\n','headerlines',lam-1);
			cellAvec = cellAline{1};
			if isempty(cellAvec)
				error('Tried to index a lambda longer than the available alphabet vector')
				return;
			end
			alphaString = cellAvec{1};
			fclose(fid);

			% Get CSSR model for the dot @ this CG base & lambda
			[ttCSSR, piCSSR, tmCSSR] = dot_to_transition2(strCSSRName,alphaString,1)

			% Simulate CSSR
			simCSSR = simulateCSSR(!!!!!!!!!!,1000,ttCSSR,tmCSSR)


		end


	end 

end


% -------------------------------------
% -------------------------------------

% Locally defined functions 

% -------------------------------------
% -------------------------------------

function [datNum, datProbs, datProbsOn, datProbsOff] = importDat(datName)
	% IMPORTDAT Local function that imports a coarse-grained datafile
	% and translates it into a numerical rather than alphabetical sequence
	%
	% Also grabs probabilities
	dat = fileread(datName);
	datOn = dat(1:2:end);
	datOff - dat(2:2:end);

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

	% Grab split data emission statistics (same procedure as above)
	datTableOn = tabulate(datOn');
	datTableOff = tabulate(datOff');
	datFreqsOn = cell2mat(datTableOn(:,2)); % get letter frequencies
	datFreqsOff = cell2mat(datTableOff(:,2)); 
	datProbsOn = datFreqsOn./length(datOn); % get letter probabilities
	datProbsOff = datFreqsOff./length(datOff); 
	datTableOn(:,4) = num2cell(datProbsOn);
	datTableOff(:,4) = num2cell(datProbsOff);
	[~, datSrtOn] = sort([datTableOn{:,1}],'ascend');
	[~, datSrtOff] = sort([datTableOff{:,1}],'ascend');
	datSrtTabOn = datTableOn(datSrtOn,:);
	datSrtTabOff = datTableOff(datSrtOff,:);
	datAlphaOn = [datSrtTabOn{:,1}]; % This ensures that the emission alphabet is sorted
	datAlphaOff = [datSrtTabOff{:,1}]; 
	datProbsOn = [datSrtTabOn{:,4}];
	datProbsOff = [datSrtTabOff{:,4}];

	datRangeOn = 1:length(datAlphaOn); 
	datNumOn = zeros(1,length(datOn));
	for i = 1:length(datAlphaOn)
	    datNumOn(strfind(datOn,datAlphaOn(i))) = datRangeOn(i);
	end

	datRangeOff = 1:length(datAlphaOff); 
	datNumOff = zeros(1,length(datOff));
	for i = 1:length(datAlphaOff)
	    datNumOff(strfind(datOff,datAlphaOff(i))) = datRangeOff(i);
	end
end


function [datAlpha, datAlphaOn, datAlphaOff] = getAlphabet(datName)
	dat = fileread(datName);	
	datOn = dat(1:2:end);
	datOff = dat(2:2:end);

	% Grab interleaved data emission statistics
    datTable = tabulate(dat'); % transpose for row vector data
    datFreqs = cell2mat(datTable(:,2)); % get letter frequencies
    datProbs = datFreqs./length(dat); % get letter probabilities
    datTable(:,4) = num2cell(datProbs);
    [~, datSrt] = sort([datTable{:,1}],'ascend');
    datSrtTab = datTable(datSrt,:);
    datAlpha = [datSrtTab{:,1}]; % This ensures that the emission alphabet is sorted

    datTableOn = tabulate(datOn');
    datTableOff = tabulate(datOff');
    datFreqsOn = cell2mat(datTableOn(:,2)); % get letter frequencies
    datFreqsOff = cell2mat(datTableOff(:,2)); 
    datProbsOn = datFreqsOn./length(datOn); % get letter probabilities
    datProbsOff = datFreqsOff./length(datOff); 
    datTableOn(:,4) = num2cell(datProbsOn);
    datTableOff(:,4) = num2cell(datProbsOff);
    [~, datSrtOn] = sort([datTableOn{:,1}],'ascend');
    [~, datSrtOff] = sort([datTableOff{:,1}],'ascend');
    datSrtTabOn = datTableOn(datSrtOn,:);
    datSrtTabOff = datTableOff(datSrtOff,:);
    datAlphaOn = [datSrtTabOn{:,1}]; % This ensures that the emission alphabet is sorted
    datAlphaOff = [datSrtTabOff{:,1}]; 
end


function [ttDOT, piDOT, tmDOT] = dot_to_transition2(filename,alphabetVec,doFuzz)
	% Exctracts a transition tensor, and the asymptotics of the flattened matrix
	% from a GraphViz file http://www.research.att.com/sw/tools/graphviz
	%
	% INPUT:  'filename' - the file in DOT format containing the graph layout
	%                       where graph is a unifilar HMM
	%                       Labels should satisfy <symbol:probability> format
	%           'alphabetVec' - string vector of all the CG2k alphabet
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
	nAlpha = length(alphabetVec);

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
	    if isempty(strfind(alphabetVec,aTemp))
	        error('No matching DOT symbol to alphabet. Check supplied alphabet same as the one which generated the dot file.');
	        return;
	    end
	    ttDOT(str2double(dotline{1})+1,...
	        str2double(dotline{3})+1,...
	        strfind(alphabetVec,aTemp)) = str2double(dotline{7});   
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
end


function simCSSR = simulateCSSR(datNum,nDatSets,ttDOT,tmDOT)

	% Simulation step
	dtmcCSSR = dtmc(tmDOT);
	% Now we need to get symbol emissions
    simCSSR = zeros(length(datNum),nDatSets);
    for j = 1:nDatSets 
        s = simulate(dtmcCSSR,length(datNum)+1);
        for t = 1:length(datNum)

            % Find the nonzero probability symbols that can be emitted
            z1 = find(ttDOT(s(t), s(t+1),:));
            z2 = reshape(ttDOT(s(t),s(t+1),find(ttDOT(s(t),s(t+1),:))),[1,length(z1)]);
            eProbs = [z1, z2'];

            % eProbs is now a probability distribution of possible emissions given the
            % state transition at the current step. Draw from it
            if length(eProbs(:,1)) == 1
                simCSSR(t,j) = eProbs(1,1);
            else
                simCSSR(t,j) = randsample(eProbs(:,1),1,true,eProbs(:,2));
            end
        end
    end
end


function [ttARP, piARP, simARP] = makeARP(datName,nDatSets)
	% INPUTS:	filename of compressed data
	%
	% OUTPUTS:	ARP transition tensor
	%			ARP initial seeding vector
	%			Simulated ARP data


    % remember symbol 1 = off, symbol 2 = on
    % Only need to fill in the off-diagonals
    % For convention set A->B = On (matrix element 1,2)

    % Get split alphabet from data
 	[datNum, ~, datProbsOn, datProbsOff] = importDat(datName);
    [~, datAlphaOn, datAlphaOff] = getAlphabet(datName);

    % Initialise memory
    ttARP = zeros(2,2,max(length(datAlphaOn),length(datAlphaOff)));

    % Note for waiting time HMMs, only the rows of the emission-flattened matrix sum to 1
    ttARP(1,2,1:length(datAlphaOn)) = datProbsOn; % Fill A->B (On waiting times)
    ttARP(2,1,1:length(datAlphaOff)) = datProbsOff;  % Fill B->A (Off waiting times)
    % !!! This fails if alphabet letters are skipped cf. line 63 !!!


    % Check if on and off aphabets are asymmetrical to make fuzzy ARP
    % Fill shorter state trailing zeros with a micro transition
    % Let's just make it 1/100th of the smallest existing nonzero transition
    if length(datAlphaOn) < length(datAlphaOff)
        % ON shorter, need to fill in tt(1,2) with small micro transition
        ttARP(1,2,find(ttARP(1,2,:)==0)) = min(ttARP(ttARP>0))/100;
        % Renormalise
        ttARP(1,2,:) = ttARP(1,2,:)./sum(ttARP(1,2,:));
        
    elseif length(datAlphaOn) > length(datAlphaOff)
        % OFF shorter, need to fill in tt(2,1)
        ttARP(2,1,find(ttARP(2,1,:)==0)) = min(ttARP(ttARP>0))/100;
        % Renormalise
        ttARP(2,1,:) = ttARP(2,1,:)./sum(ttARP(2,1,:));

    end

    % Remember interleaved data begins at ON by construction, 
    % so pi should be (1,0) for the forward algorithm
    piARP = [1,0]; % Since A->B implies On transition

    % Generate PL-ARP data
    % Doing this as two alternating biased dice
    % The fast way is to index every ODD element right away and sample
    % then index every EVEN element right away and sample
    simARP = zeros(length(datNum),nDatSets);

    distPLARPon = [ttARP(1,2,:)];
    distPLARPoff = [ttARP(2,1,:)];
    % Shave off trailing zeros on the distribution if they exist
    % (Shouldn't be the case if we do fuzzy ARP)
    distPLARPon = distPLARPon(1:find(distPLARPon,1,'last')); 
    distPLARPoff = distPLARPoff(1:find(distPLARPoff,1,'last'));

    for n = 1:nDatSets
        % Sample On wait times
        simARP(1:2:end,n) = randsample(...
            1:max(max(1:length(datAlphaOn)),max(1:length(datAlphaOff))),...
            length(datNum(1:2:end)),...
            true,...
            distPLARPon);

        % Sample Off wait times
        simARP(2:2:end,n) = randsample(...
            1:max(max(1:length(datAlphaOn)),max(1:length(datAlphaOff))),...
            length(datNum(2:2:end)),...
            true,...
            distPLARPoff);
    end
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















