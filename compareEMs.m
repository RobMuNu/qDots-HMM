% Compare two CSSR HMMs at different max memory lengths to the data
% to see whether longer memory is a better fit or not
% This script uses local functions which need at least R2016b MATLAB


% Get compressed blinking quantum dot data
datFilePattern = fullfile(pwd, '*-SML'); % Change to whatever pattern you need.
datList = dir(datFilePattern);

% Make folder to put results in
resName = horzcat('resCompareEMs-',datestr(now,30));
mkdir(resName);


for dSet = 1:length(datList)
	datName = datList(dSet).name;
	datPath = fullfile(datList(dSet).folder, datName);
	fprintf(1, 'Now reading %s\n', datName);

	% Get graphviz dot files inferred from CSSR
	gvizFilePattern = fullfile(pwd, horzcat(datName,'*.dot'));
	gvizList = dir(gvizFilePattern);

	% Pick the memory lengths you want to test between
	% Minimum memory length 3
	dotNameMin = gvizList(3).name;
	dotNameMax = gvizList(length(gvizList)).name;

	% Now can we call the nLL AND dot_to_transition2 at the same time?
	[ttMin, piMin, ~] = dot_to_transition2(dotNameMin,getAlphabet(datName),1);
	[ttMax, piMax, ~] = dot_to_transition2(dotNameMax,getAlphabet(datName),1);
	
    nLL = nllHMM(importDat(datName),ttMin,ttMax,piMin,piMax);

    cd(fullfile(pwd,resName))
    fidResults = fopen('eMcomparison.txt', 'a+');
	% nLL(importDat(datName),ttH0,ttH1,piH0,piH1)
    fprintf(fidResults, 'Comparison:\n');
    fprintf(fidResults, 'H0:\t%s\n', dotNameMin);
    fprintf(fidResults, 'H1:\t%s\n', dotNameMax);
    fprintf(fidResults, 'nLL:\t%1.2f\n', nLL);
    if nLL > 0
        fprintf(fidResults, '%s More likely given data', dotNameMin);
    elseif nLL < 0
        fprintf(fidResults, '%s More likely given data', dotNameMax);
    elseif nLL == 1
        fprintf(fidResults, 'Both eMs equally likely');
    else
        fprintf(fidResults, '[ERROR] nLL calculation failed');
    end

    cd ..

end

fclose(fidResults);
clear dotNameMax dotNameMin ttMin ttMax piMin piMax nLL fidResults resName
clear datFilePattern datList


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








