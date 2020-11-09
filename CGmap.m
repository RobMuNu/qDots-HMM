% Script that makes a heatmap of the fraction of datasets
% where CSSR is preferred over ARP 
% as a function of CG bins and max memory length



% Run in the directory that has all the gviz files as well as the CG data

nTrials = 1000;
sig = 0.01;

% Initialize the heatmap matrix
matAcceptMemory = zeros(10,10,2);

wBar = waitbar(0,'Doing a button of comparisons...');
wCount = 0;

for x = 2:10
	% For all quantum dots @ CG bin size
	strBaseExtension = sprintf('*-EXb%02d',x);

	% 29 datafiles coarse grained under this CG bin size
	datFilePattern = fullfile(pwd,strBaseExtension);
	datList = dir(datFilePattern);


	for dSet = 1:length(datList)

		% For a specific quantum dot @ CG bin size
		strDotName = datList(dSet).name;
		[datNum,~,~,~] = importDat(strDotName);
		alphaFilePattern = fullfile(pwd,...
							horzcat(erase(strDotName,...
							erase(strBaseExtension,'*')),...
							'-ALPHA'));
		alphaVec = dir(alphaFilePattern); % Should only ever index 1 file 
		gvFilePattern = fullfile(pwd, horzcat(strDotName,'*.dot'));
		gvList = dir(gvFilePattern);

		% Read the correct line from the alphabet vector
		fid = fopen(alphaVec.name);
		cellAline = textscan(fid,'%s',1,'delimiter','\n','headerlines',(x-1)-1);
		cellAvec = cellAline{1};
		if isempty(cellAvec)
			error('Tried to index a lambda longer than the available alphabet vector')
			return;
		end
		alphaString = cellAvec{1};
		fclose(fid);


		% Construct and simulate the ARP for this data
		[ttARP, piARP, simARP] = makeARP(strDotName,nTrials);


		for lam = 1:length(gvList)
			% For a specific quantum dot @ CG bin size & CSSR lambda
			strCSSRName = gvList(lam).name;

			

			% Get CSSR model for the dot @ this CG base & lambda
			[ttCSSR, piCSSR, tmCSSR] = dot_to_transition2(strCSSRName,alphaString,1);
			% Simulate CSSR
			simCSSR = simulateCSSR(datNum,nTrials,ttCSSR,tmCSSR);

			% Compute the log-likelihood threshold between 
			% ARP and this memory-length CSSR
			thresh = nllgHMM(simARP,ttARP,ttCSSR,piARP,piCSSR,sig);
			% Compute the actual likelihood for data
			nLL = nllHMM(datNum,ttARP,ttCSSR,piARP,piCSSR);

			% Write into the proportion matrix
			if nLL < thresh 
				% Accept CSSR
				matAcceptMemory(x,lam,1) = matAcceptMemory(x,lam,1) + 1;
			else
				% Reject CSSR
				matAcceptMemory(x,lam,2) = matAcceptMemory(x,lam,2) + 1;
			end
			wCount = wCount + 1;
			waitbar(wCount/783,wBar)
		end
	end 
end
close(wBar);

