% Perform hypothesis between HMM architectures based off a data set
% INPUT:        A CG-compressed dataset (interleaved, leading ON)
%               The dot file for the CSSR graph obtained from the CG qd data
%               Desired significance level
%               Number of independent trials to run for thresholding            
%
% OUTPUT:       Confusion matrices for all Hypothesised test pairs
%               Rejection decisions of all test pairs based off real CG data
%               (Optional) Statistical power vs. Data length plots 

% TEST PAIRS:   PL-ARP vs FAIR COIN
%               PL-ARP vs BIASED COIN
%               PL-ARP vs DUPLICATED PL-ARP
%               PL-ARP vs CSSR
    
% NOTES:        
%               Req. CG datafile(s) (joint and state split) in dir
%               (Opt) Req. Binarised QD trace datafile in running dir 
%
%               Req. Log forward algorithm fa_log.m in running dir
%               Req. Random n-sided dice roller roll.m in running dir
%               Req. Log-likelihood calculator nllHMM.m in running dir
%               Req. nLL threshold calculator nllgHMM.m in running dir
%               Req. CSSR dot file to transition matrix dot_to_transition2.m
%               Req. Generalised harmonic number function genHarm.m
%
%               All the data simulated is in the CG compression!!



% ---------------------------
%          Inputs
% ---------------------------
% Get compressed blinking quantum dot data
datFilePattern = fullfile(pwd, '*-SML'); % Change to whatever pattern you need.
datList = dir(datFilePattern);

%promptTrials = 'Enter amount of desired trials: ';
%nDatSets = input(promptTrials);
%promptSig = 'Enter desired significance level: ';
%siglevel = input(promptSig);
nDatSets = 1000;
siglevel = 0.01;

% Make folder to put results in
resName = horzcat('results-',datestr(now,30));
mkdir(resName);


% ---------------------------
%       Dataset Loop
% ---------------------------
for dSet = 1:length(datList)
    datName = datList(dSet).name;
    datPath = fullfile(datList(dSet).folder, datName);
    fprintf(1, 'Now reading %s\n', datName);

    % ---------------------------
    %      Data Importing
    % ---------------------------
    % Pull in CG (interleaved) datafile
    dat = fileread(datName);
    % Make CG (split) data (req. for PL-ARP construction)
    datOn = dat(1:2:end);
    datOff = dat(2:2:end);
    warning('Make sure all CG2k data has no empty line at EOF')


    % Grab interleaved data emission statistics
    datTable = tabulate(dat'); % transpose for row vector data
    datFreqs = cell2mat(datTable(:,2)); % get letter frequencies
    datProbs = datFreqs./length(dat); % get letter probabilities
    datTable(:,4) = num2cell(datProbs);
    [~, datSrt] = sort([datTable{:,1}],'ascend');
    datSrtTab = datTable(datSrt,:);
    datAlpha = [datSrtTab{:,1}]; % This ensures that the emission alphabet is sorted
    datProbs = [datSrtTab{:,4}];
    clear datTable datFreqs datSrt % cleanup

    % The forward algorithm won't recognize alphanumeric sequences
    % We will need to convert the path to numerical emissions
    datRange = 1:length(datAlpha); 
    datNum = zeros(1,length(dat));
    for i = 1:length(datAlpha)
        datNum(strfind(dat,datAlpha(i))) = datRange(i);
    end
    clear i

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
    clear datTableOn datTableOff datFreqsOn datFreqsOff datSrtOn datSrtOff % cleanup

    % We will need to convert the path to numerical emissions for forward algorithm
    % !!! Warning: This method fails if an alphabet letter is skipped !!!
    %               e.g. Alphabet a,b,d will produce wait times 1,2,3 NOT 1,2,4 

    % Maybe this doesn't actually matter unless your trying to match the
    % wait times under the same binning between two HMMs which we aren't
    datRangeOn = 1:length(datAlphaOn); 
    datNumOn = zeros(1,length(datOn));
    for i = 1:length(datAlphaOn)
        datNumOn(strfind(datOn,datAlphaOn(i))) = datRangeOn(i);
    end
    clear i

    datRangeOff = 1:length(datAlphaOff); 
    datNumOff = zeros(1,length(datOff));
    for i = 1:length(datAlphaOff)
        datNumOff(strfind(datOff,datAlphaOff(i))) = datRangeOff(i);
    end
    clear i
    warnText = sprintf('Check that no CG2k alphabet letters are skipped!\nIf there exist skipped letters, <datProbs> variables need zeros inserted!');
    warning(warnText)
    %alphaFull = 'a':'z';
    % Maybe we can use diff(alphabet). If we somehow use the correct letter positions
    % then unique(diff(alphabet)) will be 1
    % Say our alphabet is [a,b,d] --> [1,2,4] --> diff([1,2,4]) = [1,2]
    % therefore we can check unique(diff(alphabet)) == 1;
    clear warnText

    % ---------------------------
    %     Fair Coin Creation
    % ---------------------------
    % SUPPORT:  Full
    % A 1 state HMM with flat emission distribution
    % Equivalent to rolling a length(datAlpha)-sided dice


    ttFCoin = zeros(1,1,length(datAlpha));
    ttFCoin(1,1,:) = deal(1/length(datAlpha));
    piFCoin = 1; % Trivial asymptotics

    % Generate fair coin data
    simFCoin = roll(nDatSets,length(datAlpha),length(dat),1);


    % ---------------------------
    %     Biased Coin Creation
    % ---------------------------
    % SUPPORT:  Full
    % 1 state HMM with data emission probs as distribution 

    ttBCoin = zeros(1,1,length(datAlpha));
    for symbol = 1:length(datAlpha)
        ttBCoin(1,1,symbol) = datProbs(symbol);
    end
    clear symbol
    piBCoin = 1; % Trivial asymptotics

    % Generate biased coin data
    % nDatSets indep datasets, each at same length of blinking data
    for n = 1:nDatSets
        simBCoin(:,n) = randsample(...
        datRange,...
        length(dat),...
        true,...
        [ttBCoin(1,1,:)]);
    end
    clear n

    % -----------------------------------------
    %           PL-ARP Creation
    % -----------------------------------------
    % SUPPORT:  Full (for case 2)
    % 2 State alternating PL-HMM

    % For the emission / transition distributions...
    % CASE 1: Pull distributions off the split data & make "fuzzy" ARP
    %           This helps when the powerlaw fits from case 2 don't perform well.
    % CASE 2: Fit truncated power law to waiting times first

    plCase = 1;
    if plCase == 1
        % ---- Case 1 ----
        ttPLARP = zeros(2,2,max(length(datAlphaOn),length(datAlphaOff)));
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


    elseif plCase == 2
        % ---- Case 2 ----

        % Create waiting time trajectory from the binary quantum dot intensity trace
        % This is only really used for method 2 ARP creation
        datInfo=dir('*-cen');
        datNames={datInfo.name};
        traj = fileread(datNames{1});
        % Quickly convert the string into a numerical array
        traj = traj-'0';
        % Convert the binary trace to state-split waiting times
        findOn = find(diff([0,traj,0]==1)); % For one
        findOff = find(diff([1,traj,1]==0)); % For zero
        % Find starting indices of *blocks* of ones and zeros
        idxStartOn = findOn(1:2:end-1);  % Starting indices of 1's blocks
        idxStartOff = findOff(1:2:end-1);  % Start indices
        % The generalised waiting times
        waitOn = findOn(2:2:end)-idxStartOn;  % Consecutive ones counts
        waitOff = findOff(2:2:end)-idxStartOff;  % Consecutive zeros counts
        clear findOn findOff idxStartOff idxStartOn traj datInfo datNames


        waitMax = max([waitOn,waitOff]);

        waitTabOn = tabulate(waitOn);
        waitTabOff = tabulate(waitOff);
        waitFreqsOn = waitTabOn(:,2);
        waitFreqsOff = waitTabOff(:,2);
        waitTabOn(:,4) = waitFreqsOn./length(waitOn);
        waitTabOff(:,4) = waitFreqsOff./length(waitOff);
        waitTimeOn = waitTabOn(:,1);
        waitTimeOff = waitTabOff(:,1);
        waitProbsOn = [waitTabOn(:,4)];
        waitProbsOff = [waitTabOff(:,4)];
        clear waitFreqsOn waitFreqsOn waitTabOn waitTabOff


        % Take log wait times and probabilities for PL fitting
        % Remove -Inf log-probabilities and corresp. waiting times from influencing fit
        wLogOn = log(waitTimeOn);
        wLogProbOn = log(waitProbsOn);
        wLogOn(isinf(wLogProbOn)) = [];
        wLogProbOn = wLogProbOn(find(wLogProbOn~=(-Inf)));

        wLogOff = log(waitTimeOff);
        wLogProbOff = log(waitProbsOff);
        wLogOff(isinf(wLogProbOff)) = [];
        wLogProbOff = wLogProbOff(find(wLogProbOff~=(-Inf)));

        % Define anonymous Zipf distribution with a fixed support which matches
        % the size of the largest alphabet between Wait-on or Wait-off
        zipfPMF = @(x,a) x.^(-a) ./ genHarm(waitMax,a);
        % Run MLE distribution fit for a pinned maximum wait time
        [plParamOn, plCIOn] = mle(waitOn,'pdf',zipfPMF,'start',5,'lower',1);
        [plParamOff, plCIOff] = mle(waitOff,'pdf',zipfPMF,'start',5,'lower',1);

            % Optional plotting to see how well the MLE fits
            % close all
            % subplot(1,2,1)
            % plot(wLogOn,wLogProbOn,'+g');
            % hold on
            % fplot(@(x) -plParamOn*x - log(genHarm(waitMax,plParamOn)),...
            %     [0,max(wLogOff)]);

            % subplot(1,2,2)
            % plot(wLogOff,wLogProbOff,'+r');
            % hold on
            % fplot(@(x) -plParamOff*x - log(genHarm(waitMax,plParamOff)),...
            %     [0,max(wLogOff)]);
        clear wLogOn wLogProbOn wLogOff wLogProbOff waitTimeOn waitTimeOff
        clear waitProbsOn waitProbsOff


        % Create waiting time distributions based off MLE fits
        distPLWaitOn = zipfPMF(1:waitMax,plParamOn);
        distPLWaitOff = zipfPMF(1:waitMax,plParamOff);
        % Create the full supp CG2k-dist directly from full supp PL-dist
        % A surjective map f: PL(w)-->CG(k)
        indexRule = @(w) floor(log2(w))+1;
        distCGPLon = zeros(1,indexRule(waitMax));
        distCGPLoff = zeros(1,indexRule(waitMax));
        for w = 1:waitMax
            distCGPLon(indexRule(w)) = distCGPLon(indexRule(w))...
                                     + distPLWaitOn(w);
            distCGPLoff(indexRule(w)) = distCGPLoff(indexRule(w))...
                                     + distPLWaitOff(w);
        end
        clear w indexRule distPLWaitOn distPLWaitOff zipfPMF waitMax

        % Simulate data
        simPLARP = zeros(length(datNum),nDatSets);
        for n = 1:nDatSets
            % Sample On wait times
            simPLARP(1:2:end,n) = randsample(...
                1:length(distCGPLon),...
                length(dat(1:2:end)),...
                true,...
                distCGPLon);

            % Sample Off wait times
            simPLARP(2:2:end,n) = randsample(...
                1:length(distCGPLoff),...
                length(dat(2:2:end)),...
                true,...
                distCGPLoff);
        end
        clear n

        % Create PL-ARP structure (CG2k)
        if length(distCGPLon) ~= length(distCGPLoff)
            error('Powerlaw-fit distributions must share same CG2k support!');
            return;
        else
            ttPLARP = zeros(2,2,max(length(distCGPLon),length(distCGPLoff)));
            % Only need to fill in the off-diagonals
            % For convention set A->B = On 
            % Note for waiting time HMMs, only the rows of the emission-flattened matrix
            % sum to 1
            ttPLARP(1,2,1:length(distCGPLon)) = distCGPLon; % Fill A->B (On waiting times)
            ttPLARP(2,1,1:length(distCGPLoff)) = distCGPLoff;  % Fill B->A (Off waiting times)
        
            % Get initial seed distribution assuming leading symbol is ON
            piPLARP = [1,0];

        end
        warning('Check: that the powerlaw fits for ON and OFF make sense!')
       
    else
        error('PL-ARP generation must be method 1 or 2');
        return;
    end
    clear plCase 
    clear plCIOff plCIOn plParamOff plParamOn


    % -----------------------------------------
    %        Duplicated PL-ARP Creation
    % -----------------------------------------

    doMethod = 1;
    if doMethod == 1
        % HMM will draw from ON and OFF independently, then "echo" that
        % pair of outcomes before drawing independently again.

        % Let's not assume that the ON and OFF alphabet are the same size
        CGAsizeOn = length(distCGPLon); % Size of ON alphabet |A_on|
        CGAsizeOff = length(distCGPLoff); % Size of OFF alphabet |A_off|

        % The number of states follows this rule:
        % 2*|A_on|*|A_off| + |A_on| + 1
        nEchostates = 2*CGAsizeOn*CGAsizeOff + CGAsizeOn + 1;
        ttPLARPecho = zeros(nEchostates,nEchostates,max(CGAsizeOn,CGAsizeOff));

        % Fill the available transitions
        % There's not really a nice and symmetric way of doing this as far as I can tell
        % so here we go:

            % Sample ON elements
            for b = 1:CGAsizeOn
                ttPLARPecho(1,b+1,b) = distCGPLon(b);
            end

            % Sample OFF elements
            for b = 1:CGAsizeOn
                for c = 1:CGAsizeOff
                    ttPLARPecho(1+b,CGAsizeOn+1+(CGAsizeOff*(b-1))+c,c) = distCGPLoff(c);
                end
            end

            % Deterministic Echo ON elements
            d1 = 1+CGAsizeOn+(CGAsizeOn*CGAsizeOff); 
            for b = 1:CGAsizeOn
                for c = 1:CGAsizeOff
                    ttPLARPecho(CGAsizeOn+1+c+(CGAsizeOff*(b-1)),...
                        d1+c+(CGAsizeOff*(b-1)),b) = 1;
                end
            end

            % Deterministic Echo OFF elements
            for c = 1:CGAsizeOff
                ttPLARPecho(d1+c:CGAsizeOff:end,1,c) = 1;
            end

        piPLARPecho = zeros(1,nEchostates);
        piPLARPecho(1) = 1;

        % Simulate
        % Quick and dirty method
        % Check whether original dataset is even or odd (fixes assignment issues)
        if rem(length(datNum),2) == 0
            % Even assignment
            % Allocate memory
            simPLARPecho = zeros(2*length(simPLARP(:,1)),nDatSets);
            for n = 1:nDatSets
                % Reshape the data vector into on/off pairs (col1 on, col2 off)
                % Then repeat the on/off pairs
                % Reshape back into a vector
                simPLARPecho(:,n) = reshape(repelem(reshape(simPLARP(:,n),2,[])',2,1)',1,[])';
            end
        else
            % Odd assignment
            % Allocate memory
            simPLARPecho = zeros(2*length(simPLARP(1:end-1,1)),nDatSets);
            % Just shave off the last datapoint to make it even
            for n = 1:nDatSets
                simPLARPecho(:,n) = reshape(repelem(reshape(simPLARP(1:end-1,n),2,[])',2,1)',1,[])';
            end
        end
        % Cut to fit original dataset
        simPLARPecho = simPLARPecho(1:length(datNum),:);

        clear nEchostates distCGPLoff distCGPLon b c d1 n


    elseif doMethod == 2
        % Every waiting time is repeated. If ON was 3, then OFF will be 3
        % So if the interleaved CG2k sequence is 
        % a1,b1,a2,b2,...aN,bN
        % a = on waiting time
        % b = off waiting time
        % then we can make a model where bk = ak forall k
        % that is both unifilar and has full support. It's a bit lame though since
        % all OFF transitions are fully deterministic but hey.
        % effectively what we end up if we split the on/off sequences and find their
        % distributions is that the OFF sequence is an *exact* copy of the ON sequence
        % You'd think that this model could be immediately ruled out if we observed
        % a OFF distribution != ON distribution.
        % Thats a feature. However, what if we took the off transitions and rotated or permuted
        % them about the ON HMM state? (Draw it if that helps)
        % You'd end up with an OFF distribution != the on distribution
        % the probabilities would be the same but with the bins permuted about.
        % Anyway this is just a strawman so it doesn't really matter too much


        % Make it share the same alphabet size as PL-ARP by
        % making the ON transition probabilities the same as the fitted powerlaw
        % from ttPLARP

        maxCGPLalpha = max(length(distCGPLon),length(distCGPLoff)); % args are probs the same
        ttPLARPecho = zeros(maxCGPLalpha+1,maxCGPLalpha+1,maxCGPLalpha);
        for b = 1:maxCGPLalpha
            ttPLARPecho(1,b+1,b) = distCGPLon(b);
            ttPLARPecho(b+1,1,b) = 1;
        end
        % Asymptotics are trivial if we assume leaing ON always
        piPLARPecho = zeros(1,maxCGPLalpha+1);
        piPLARPecho(1) = 1;

        % Simulate
        % Quick and dirty method. Take simulated PLARP on bins and double them
        simPLARPecho = zeros(length(datNum),nDatSets);

        for n = 1:nDatSets
            % Copy every second element to every fourth one until half the seq length
            simPLARPecho(1:4:end,n) = simPLARP(1:2:end/2,n);
            simPLARPecho(2:4:end,n) = simPLARP(2:2:end/2,n);
        end

        % Check whether original dataset is even or odd (fixes assignment issues)
        if rem(length(datNum),2) == 0
            % Even assignment
            for n = 1:nDatSets
                simPLARPecho(1:2:end,n) = simPLARP(1:2:end,n);
                % Do the echo
                simPLARPecho(2:2:end,n) = simPLARP(1:2:end,n);
            end
        else
            % Odd assignment
            for n = 1:nDatSets
                simPLARPecho(1:2:end,n) = simPLARP(1:2:end,n);
                % Do the echo
                simPLARPecho(2:2:end,n) = simPLARP(1:2:end-1,n);
            end
        end

        clear maxCGPLalpha distCGPLoff distCGPLon b n
    end


    % ---------------------------------
    %     Statistical power loop
    % ---------------------------------
    % Run the thresholding and confusion matrix as a function of data length

    %dLenRange = horzcat([5:5:100],[200,300,400,500]);
    dLenRange = length(datNum);
    powVec = zeros(length(dLenRange),4);
    idxDLen = 1;
    if length(dLenRange) > 1
        wBar2 = waitbar(0,'Data length power scanning...');
    end
    for dLen = dLenRange

        % ---------------------------------
        %     Hypothesis pair thresholds
        % ---------------------------------
        % Collate all datasets into a big array for ease of indexing
        simSets = zeros(dLen,nDatSets,5);
        simSets(:,:,1) = simFCoin(1:dLen,:);
        simSets(:,:,2) = simBCoin(1:dLen,:);
        simSets(:,:,3) = simPLARPecho(1:dLen,:); 
        % CSSR defined after CSSR creation
        simSets(:,:,5) = simPLARP(1:dLen,:);

        % Calculate the thresholds for all the hypothesis pairs
        threshPair = zeros(1,4);
        nLLPair = zeros(1,4);
        % Pair 1    H0: FCoin vs. H1: PLARP
        threshPair(1) = nllgHMM(simSets(:,:,1),...
                                ttFCoin,ttPLARP,...
                                piFCoin,piPLARP,...
                                siglevel);
        nLLPair(1) = nllHMM(datNum, ttFCoin, ttPLARP, piFCoin, piPLARP);

        % Pair 2    H0: BCoin vs. H1: PLARP
        threshPair(2) = nllgHMM(simSets(:,:,2),...
                                ttBCoin,ttPLARP,...
                                piBCoin,piPLARP,...
                                siglevel);
        nLLPair(2) = nllHMM(datNum, ttBCoin, ttPLARP, piBCoin, piPLARP);

        % Pair 3    H0: PLARPecho vs. H1: PLARP
        threshPair(3) = nllgHMM(simSets(:,:,3),...
                                ttPLARPecho,ttPLARP,...
                                piPLARPecho,piPLARP,...
                                siglevel);
        nLLPair(3) = nllHMM(datNum, ttPLARPecho, ttPLARP, piPLARPecho, piPLARP);


        % ------------------------------
        %          GraphViz Import
        % ------------------------------
        % Get graphviz dot files inferred from CSSR
        gvizFilePattern = fullfile(pwd, horzcat(datName,'*.dot'));
        gvizList = dir(gvizFilePattern);
        % Loop over each maximum memory length CSSR model you have
        for lam = 1:length(gvizList)
            % Name of the CSSR-generated dot file
            dotName = gvizList(lam).name;
            fprintf(1, 'CSSR at max length: %d\n', lam);

            % ------------------------------
            %          CSSR Creation
            % ------------------------------
            doFuzz = 1;
            if doFuzz ~= 0 & doFuzz ~= 1
                error('doFuzz should be 0 or 1 only')
                return;
            end
            [ttCSSR, piCSSR, tmCSSR] = dot_to_transition2(dotName,datAlpha,doFuzz);

            % Simulate (this is gonna be a pain to do)
            dtmcCSSR = dtmc(tmCSSR);
            % Now we need to get symbol emissions
            simCSSR = zeros(length(datNum),nDatSets);
            for j = 1:nDatSets 
                s = simulate(dtmcCSSR,length(datNum)+1);
                for t = 1:length(datNum)

                    % Find the nonzero probability symbols that can be emitted
                    z1 = find(ttCSSR(s(t), s(t+1),:));
                    z2 = reshape(ttCSSR(s(t),s(t+1),find(ttCSSR(s(t),s(t+1),:))),[1,length(z1)]);
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
            clear j t dtmcCSSR s z1 z2 eProbs doFuzz

            % ----------------------------------
            % Pair 4    H0: CSSR vs. H1: PLARP
            % ----------------------------------
            simSets(:,:,4) = simCSSR(1:dLen,:);
            threshPair(4) = nllgHMM(simSets(:,:,4),...
                                    ttCSSR,ttPLARP,...
                                    piCSSR,piPLARP,...
                                    siglevel);
            nLLPair(4) = nllHMM(datNum, ttCSSR, ttPLARP, piCSSR, piPLARP);


            % TESTING / DEBUGGING
            % simSets(:,:,1) = simFCoin(1:dLen,:);
            % simSets(:,:,2) = simBCoin(1:dLen,:);
            % simSets(:,:,3) = simPLARPecho(1:dLen,:); 
            % simSets(:,:,4) = simPLARP(1:dLen,:);
            % simSets(:,:,5) = simCSSR(1:dLen,:);

            % Calculate the thresholds for all the hypothesis pairs
            % !!! TESTING / DEBUGGING VERSION BELOW !!!
            % threshPair = zeros(1,4);
            % % Pair 1    CSSR vs FCoin
            % threshPair(1) = nllgHMM(simSets(:,:,5),...
            %                         ttCSSR,ttFCoin,...
            %                         piCSSR,piFCoin,...
            %                         siglevel);
            % % Pair 2    H0: CSSR vs Bcoin
            % threshPair(2) = nllgHMM(simSets(:,:,5),...
            %                         ttCSSR,ttBCoin,...
            %                         piCSSR,piBCoin,...
            %                         siglevel);
            % % Pair 3    H0: CSSR vs PLARPecho
            % threshPair(3) = nllgHMM(simSets(:,:,5),...
            %                         ttCSSR,ttPLARPecho,...
            %                         piCSSR,piPLARPecho,...
            %                         siglevel);
            % % Pair 4    H0: CSSR vs. H1: PLARP
            % threshPair(4) = nllgHMM(simSets(:,:,5),...
            %                         ttCSSR,ttPLARP,...
            %                         piCSSR,piPLARP,...
            %                         siglevel);

            % TESTING / DEBUGGING
            % ttAll = {ttFCoin, ttBCoin, ttPLARPecho, ttPLARP, ttCSSR};
            % piAll = {piFCoin, piBCoin, piPLARPecho, piPLARP, piCSSR};


            % ---------------------------
            %     Confusion Matrix
            % ---------------------------
            % Do this based on the number of trials
            % The full training set consists of ALL the simulated sets from H0 and H1
            % for a full number of trials as 2*nDatSets (cause why not)
            % Need all thresholds, a set of labeled datasets, nllHMM to test them with

            % Make a cell array of all the transition tensors
            ttAll = {ttFCoin, ttBCoin, ttPLARPecho, ttCSSR, ttPLARP};
            piAll = {piFCoin, piBCoin, piPLARPecho, piCSSR, piPLARP};

            % Initialise the confusion matrix
            doConfusion = 0;
            if doConfusion == 1
                matConf = zeros(2,2,4);
                for pair = 1:4 % default 1:4 (4 hypothesis pairs)
                    wBar = waitbar(0,'Counting confusion matrix...');

                    for j = 1:nDatSets

                        % run the nll test (for BOTH hypotheses! and fill the confusion matrixs)
                        confTrueH0 = nllHMM(simSets(:,j,pair),...
                                ttAll{pair},ttAll{5},...
                                piAll{pair},piAll{5});

                        confTrueH1 = nllHMM(simSets(:,j,5),...
                                ttAll{pair},ttAll{5},...
                                piAll{pair},piAll{5});

                        % % ----------TESTING / DEBUGGING ---------
                        % confTrueH0 = nllHMM(simSets(:,j,5),...
                        %         ttAll{5},ttAll{pair},...
                        %         piAll{5},piAll{pair});

                        % confTrueH1 = nllHMM(simSets(:,j,pair),...
                        %         ttAll{5},ttAll{pair},...
                        %         piAll{5},piAll{pair});
                        % % ---------------------------------------

                        % Fill in confusion matrix for this round
                        % True H1
                        if confTrueH1 < threshPair(pair)
                            % detect H1 | H1 true (TP)
                            matConf(1,1,pair) = matConf(1,1,pair)+1;
                        else
                            % detect H0 | H1 true (FN)
                            matConf(1,2,pair) = matConf(1,2,pair)+1;
                        end

                        % True H0
                        if confTrueH0 < threshPair(pair)
                            % detect H1 | H0 true (FP)
                            matConf(2,1,pair) = matConf(2,1,pair)+1;
                        else
                            % detect H0 | H0 true (TN)
                            matConf(2,2,pair) = matConf(2,2,pair)+1;
                        end
                        waitbar(j/nDatSets,wBar)

                    end
                    close(wBar);

                    powVec(idxDLen,pair) = matConf(1,1,pair);
                    
                end

                waitbar(idxDLen/length(dLenRange),wBar2)
                idxDLen = idxDLen + 1;
                % Normalise to get probabilities
                matConf = matConf./nDatSets;

                clear j pair nllTemp wBar

            end % Confusion

            % This snipped might not belong here since moving loops around.
            if length(dLenRange) > 1
                close(wBar2);
            end

            % Normalise TP counts to get power in probability
            powVec = powVec./nDatSets;
            


            % Write results to file
            % Make sure you write it to the results directory!
            cd(fullfile(pwd,resName))
            fidResults = fopen('results.txt', 'a+');

            fprintf(fidResults, 'Data:\t%s\n', dotName);
            fprintf(fidResults, 'H1 list:\n');
            fprintf(fidResults, '\tFDice\n\tBDice\n\tPLARPecho\n\tCSSR\n');
            fprintf(fidResults, 'nLL Thresholds:\n');
            fprintf(fidResults, '\t%1.2f\n', threshPair)
            fprintf(fidResults, 'nLL:\n');
            fprintf(fidResults, '\t%1.2f\n', nLLPair);
            fprintf(fidResults, 'Accept PL-ARP:\n');
            fprintf(fidResults, '\t%d \t%d \t%d \t%d\n',...
                                nLLPair(1)<threshPair(1),...
                                nLLPair(2)<threshPair(2),...
                                nLLPair(3)<threshPair(3),...
                                nLLPair(4)<threshPair(4));
            fprintf(fidResults, '---------------------------------------\n\n\n')
            fclose(fidResults);
            cd ..


            % ------------------------------------------
            %     (Optional) Statistical power plots
            % ------------------------------------------
            close all
            % Plot only if we're doing a stat power scan
            if length(dLenRange) > 1
                % --- 
                subplot(2,2,1)
                plot(dLenRange,powVec(:,1), ':or',...
                        'LineWidth',2,...
                        'MarkerFaceColor','r');
                hold on
                f(1) = yline(0.8);
                f(1).LineWidth = 2.2;
                f(1).Color = [0.5 0.5 0.5];
                f(1).LineStyle = '--';
                f(1).Alpha = 0.7;
                % Formatting
                xlabel('CG Data length')
                ylabel('Statistical power $(1-FN)$','Interpreter','latex')
                legend('Fair dice vs. PL-ARP',...
                        'Location','Southeast');
                set(gca,'FontSize',14);
                set(gca,'YLim',[0,1.1]);

                % ---
                subplot(2,2,2)
                plot(dLenRange,powVec(:,2), ':o',...
                        'LineWidth',2,...
                        'Color',[1 0.651 0],...
                        'MarkerFaceColor',[1 0.651 0]); 
                hold on
                f(1) = yline(0.8);
                f(1).LineWidth = 2.2;
                f(1).Color = [0.5 0.5 0.5];
                f(1).LineStyle = '--';
                f(1).Alpha = 0.7;
                % Formatting
                xlabel('CG Data length')
                ylabel('Statistical power $(1-FN)$','Interpreter','latex')
                legend('Biased dice vs. PL-ARP',...
                        'Location','Southeast');
                set(gca,'FontSize',14);
                set(gca,'YLim',[0,1.1]);

                % ---
                subplot(2,2,3)
                plot(dLenRange,powVec(:,3), ':o',...
                        'LineWidth',2,...
                        'Color',[0.13 0.54 0.13],...
                        'MarkerFaceColor',[0.13 0.54 0.13]);
                hold on
                f(1) = yline(0.8);
                f(1).LineWidth = 2.2;
                f(1).Color = [0.5 0.5 0.5];
                f(1).LineStyle = '--';
                f(1).Alpha = 0.7;
                % Formatting
                xlabel('CG Data length')
                ylabel('Statistical power $(1-FN)$','Interpreter','latex')
                legend('Double PL-ARP vs. PL-ARP',...
                        'Location','Southeast');
                set(gca,'FontSize',14);
                set(gca,'YLim',[0,1.1]);

                % ---
                subplot(2,2,4)
                plot(dLenRange,powVec(:,4), ':ob',...
                        'LineWidth',2,...
                        'MarkerFaceColor','b');
                hold on
                f(1) = yline(0.8);
                f(1).LineWidth = 2.2;
                f(1).Color = [0.5 0.5 0.5];
                f(1).LineStyle = '--';
                f(1).Alpha = 0.7;
                % Formatting
                xlabel('CG Data length')
                ylabel('Statistical power $(1-FN)$','Interpreter','latex')
                legend('CSSR vs. PL-ARP',...
                        'Location','Southeast');
                set(gca,'FontSize',14);
                set(gca,'YLim',[0,1.1]);
            end 


            clear promptSig promptTrials n datRange datRangeOn datRangeOff 

        end % Lambda

    end % Stat power

    clear simFCoin simBCoin simPLARP simCSSR simPLARPecho
    clear idxDLen dLen

end % Data loop

clear resName
 