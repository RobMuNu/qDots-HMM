% Generate 2k encoded speeches from photon trajectory
% INPUT: Binarised photon trajectories with some common fname pattern
% OUTPUT: ONE 2k encoded speech (on and off together) sharing a common alphabet

% NOTES: We discard leading zeros of any trajectory to make them all start at one.
% That way, the CG2k speech symbols alternate between ON and OFF, where the first
% Symbol is always ON.



% Read all datafiles in the directory and store the names in an cell
datInfo=dir('*-cen');
datNames={datInfo.name};


% Define the full alphabet (covers wait times up to 62 million steps lol)
alphaFull = 'a':'z';

% Make a results folder
%[yr,mo,da]=ymd(datetime);
%[hr,mn,se]=hms(datetime);
%mkdir dat-CG2k-nosplit

close all

% For each GY photon traj:
for dat = 1:numel(datNames)
%for dat = 10:10 % for testing

	traj = fileread(datNames{dat}); 

    % Quickly convert the char array into numerical one
    traj = traj-'0';

	% Convert the entire sequence to a trajectory of durations / waiting times
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
    waitTraj = zeros(1,length(waitOne)+length(waitZero));
    % First check what the starting traj symbol is
    if traj(1) == 1
    	% Then interleave leading with waitOne
    	waitTraj((1:length(waitOne))*2 - 1) = waitOne;
    	waitTraj((1:length(waitZero))*2) = waitZero;

    else
    	% Interleave leading with waitZero and remove the first waiting time
    	% That way we know the waiting times began in the ON state
    	waitTraj((1:length(waitZero))*2 - 1) = waitZero;
    	waitTraj((1:length(waitOne))*2) = waitOne;
    	waitTraj(1) = [];
    end

	% Find the longest waiting time
	waitMax = max(max(waitOne),max(waitZero));

	% Coarse grain the durations vector, by defining alphabet symbols 
	% a0,a1… under the ak = [2^k,2^k+1) rule
	% The LARGEST alphabet interval is defined by [2^aLower, 2^aUpper)
	aLower = floor(log2(waitMax));
	aUpper = aLower +1; % this is also max|A|

	% Now replace wait times with the alphabet bin indices e.g. 5 --> a(3)
	CG2kbins = floor(log2(waitTraj))+1;

	% Take the alphabet at those bin indices and write in the symbols eg. a(3)-->'c'
	% ------------------------------------
	CG2kall = alphaFull(CG2kbins);       % <-- This is the encoded speech
	% ------------------------------------

	% Get the alphabet array of the compressed speech
	CG2kalpha = unique(CG2kall);

	% Get the *actual* value of |A| 
	sizeCG2kalpha = numel(unique(CG2kall));
	% Calculate the length of the coarse grained data speech
	nCGdat = numel(CG2kall);

	% Evaluate maximum Lambda = log2(Ncg)/log2(|A|)
	LambdaAll = floor(log2(nCGdat)/log2(sizeCG2kalpha));


    % ----------ON / OFF SPLIT------------
    % ------------------------------------
    % Find the longest of those durations
    waitMaxOne = max(waitOne);
    waitMaxZero = max(waitZero);
    % Coarse grain the durations vector, by defining alphabet symbols 
    % a0,a1… under the ak = [2^k,2^k+1) rule
    % The LARGEST alphabet interval is defined by [2^aLower, 2^aUpper)
    aLowerOne = floor(log2(waitMaxOne));
    aUpperOne = aLowerOne + 1; % this is also max|A_on|
    aLowerZero = floor(log2(waitMaxZero));
    aUpperZero = aLowerZero + 1; % this is also max|A_off|
    % Now replace wait times with the alphabet bin indices e.g. 5 --> a(3)
    CG2kbinsOne = floor(log2(waitOne))+1;
    CG2kbinsZero = floor(log2(waitZero))+1;
    % ------------------------------------
    % Take the alphabet at those bin indices and write in the symbols eg. a(3)-->'c'
    CG2kOne = alphaFull(CG2kbinsOne);    % <- These are the compressed speeches
    CG2kZero = alphaFull(CG2kbinsZero);  %
    % ------------------------------------
    % Store the compressed alphabets for ON and OFF 
    CG2kalphaOne = unique(CG2kOne);
    CG2kalphaZero = unique(CG2kZero);
    % Get their *actual* |A| size
    sizeCG2kalphaOne = numel(unique(CG2kOne));
    sizeCG2kalphaZero = numel(unique(CG2kZero));
    % Calculate the length of the coarse grained data speech
    nCGdatOne = numel(CG2kOne);
    nCGdatZero = numel(CG2kZero);

    % Evaluate maximum Lambda = log2(Ncg)/log2(|A|)
    LambdaOne = log2(nCGdatOne)/log2(sizeCG2kalphaOne);
    LambdaZero = log2(nCGdatZero)/log2(sizeCG2kalphaZero);
    % ------------------------------------
    % -------------END SPLIT--------------


    % -----------PLOTMAKING---------------
    % Compressed alphabet histograms
    CG2ktable = tabulate(CG2kall'); % need to transpose if CG2k is a row vector
    CG2kfreqs = cell2mat(CG2ktable(:,2)); % get letter frequencies
    CG2kprobs = CG2kfreqs./length(CG2kall); % get letter probabilities
    CG2ktable(:,4) = num2cell(CG2kprobs);
    [~, srtJoin] = sort([CG2ktable{:,1}],'ascend');
    srtTabJoin = CG2ktable(srtJoin,:);
    
    CG2ktableOne = tabulate(CG2kOne');
    CG2kfreqsOne = cell2mat(CG2ktableOne(:,2));
    CG2kprobsOne = CG2kfreqsOne./length(CG2kOne);
    CG2ktableOne(:,4) = num2cell(CG2kprobsOne);
    [~, srtOne] = sort([CG2ktableOne{:,1}],'ascend');
    srtTabOne = CG2ktableOne(srtOne,:);

    CG2ktableZero = tabulate(CG2kZero');
    CG2kfreqsZero = cell2mat(CG2ktableZero(:,2));
    CG2kprobsZero = CG2kfreqsZero./length(CG2kZero);
    CG2ktableZero(:,4) = num2cell(CG2kprobsZero);
    [~, srtZero] = sort([CG2ktableZero{:,1}],'ascend');
    srtTabZero = CG2ktableZero(srtZero,:);

    figure('Position',[157 157 806 549]);
    

    subplot(2,3,4); % Bottom left (ON)
    figBarOne = bar(categorical(CG2ktableOne(:,1)),CG2kprobsOne);
    %
    text(1:sizeCG2kalphaOne,[srtTabOne{:,4}],CG2kalphaOne(:),...
        'vert','bottom','horiz','center'); 
    %
    set(gca, 'XTickLabel', cellstr(num2str([unique(CG2kbinsOne)]')))
    xlabel('On intvl. $k\in[2^{k-1},2^k)$','Interpreter','latex')



    subplot(2,3,5); % Bottom mid (OFF)
    figBarZero = bar(categorical(CG2ktableZero(:,1)),CG2kprobsZero);
    %
    text(1:sizeCG2kalphaZero,[srtTabZero{:,4}],CG2kalphaZero(:),...
        'vert','bottom','horiz','center'); 
    %
    set(gca, 'XTickLabel', cellstr(num2str([unique(CG2kbinsZero)]')))
    xlabel('Off intvl. $k\in[2^{k-1},2^k)$','Interpreter','latex')
 

    subplot(2,3,6); % Bottom right (JOINT)
    figBar = bar(categorical(CG2ktable(:,1)),CG2kprobs);
    %
    text(1:sizeCG2kalpha,[srtTabJoin{:,4}],CG2kalpha(:),...
        'vert','bottom','horiz','center'); 
    %
    set(gca, 'XTickLabel', cellstr(num2str([unique(CG2kbins)]')))
    xlabel('Joint intvl. $k\in[2^{k-1},2^k)$','Interpreter','latex')


    % -----------------------------------
    % On/Off duration & probability plots
    subplot(2,3,1); % Top left (ON)
    waitTableOne = tabulate(waitOne);
    waitTableOne(:,4) = waitTableOne(:,3)./100; % to get probabilities
    figDurOne = scatter(waitTableOne(:,1), waitTableOne(:,4));
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlabel('On state duration (ms)','Interpreter','latex')
    ylabel('Normalised event count','Interpreter','latex')

    subplot(2,3,2); % Top mid (OFF)
    waitTableZero = tabulate(waitZero);
    waitTableZero(:,4) = waitTableZero(:,3)./100; % to get probabilities
    figDurZero = scatter(waitTableZero(:,1), waitTableZero(:,4));
    set(gca,'xscale','log');
    set(gca,'yscale','log');  
    xlabel('Off state duration (ms)','Interpreter','latex')

    % Intensity peaks plot
    cd ..
    cd 'traces and intensity histos'
    figAxPeaks = subplot(2,3,3); % Specify axes we want to paste into
    figPeakLoad = hgload(erase(convertCharsToStrings(datNames{dat}),"-cen"));
    copyobj(allchild(get(figPeakLoad,'CurrentAxes')),figAxPeaks);
    cd ..
    cd 'thresh-central'
    figPeakLoad.Name = 'killme';
    close killme
    xlabel('Photon counts / 10ms time bin','Interpreter','latex')
    ylabel('Event count','Interpreter','latex')

    qName = convertCharsToStrings(datNames{dat});
    sgtitle(qName);

    % Display Cmu values on subplots
    figAxes = findall(gcf,'type','axes'); 
    figPos = get(figAxes,'position'); % cell array of the positions of subplots

    % Printing Cmu ON values onto plot
    dirCmuOne = dir('**/Cmu*-s0');
    fnameCmuOne = {dirCmuOne.name};
    foldCmuOne = {dirCmuOne.folder};
    CmuOne = fileread(strcat(foldCmuOne{dat},'/',fnameCmuOne{dat}));
    annotation(gcf,'textbox',...
        [0.580772287553494,0.395607235142119,0.037841191066998,0.040983606557377],...
        'String',CmuOne,'FitBoxToText','on','HorizontalAlignment','right')

    % Print Cmu OFF value onto plot
    dirCmuZero = dir('**/Cmu*-s1');
    fnameCmuZero = {dirCmuZero.name};
    foldCmuZero = {dirCmuZero.folder};
    CmuZero = fileread(strcat(foldCmuZero{dat},'/',fnameCmuZero{dat}));
    annotation(gcf,'textbox',...
        [0.299975186104218,0.395607235142119,0.037841191066998,0.040983606557377],...
        'String',CmuZero,'FitBoxToText','on','HorizontalAlignment','right')

    % Print Cmu JOINT value onto plot 
    dirCmuJoin = dir('**/Cmu*-CG2k');
    fnameCmuJoin = {dirCmuJoin.name};
    foldCmuJoin = {dirCmuJoin.folder};
    CmuJoin = fileread(strcat(foldCmuJoin{dat},'/',fnameCmuJoin{dat}));
    annotation(gcf,'textbox',...
        [0.862810083791851,0.393785741517347,0.037841191066998,0.040983606557377],...
        'String',CmuJoin,'FitBoxToText','on','HorizontalAlignment','right')

    % ------------------------------------    

    % Write into a text file to produce a vector of LamdaMax's
    %fidAlpha = fopen('CGlambdaAll.txt', 'a+');
    %fprintf(fidAlpha, '%d\n', LambdaAll);
    %fclose(fidAlpha);

 	% Write compressed speeches to datafiles	
    %qName = convertCharsToStrings(datNames{dat});
 	%dlmwrite(qName+'-CG2k', CG2kall, ''); % dot CG2k speech
 	%dlmwrite(qName+'-alpha', CG2kalpha, ''); % dot CG2k alphabet
    cd 'CG2k-summary'
    print(gcf,qName+'-fig','-dpdf','-bestfit')
    close all
    cd ..

end

% Clean up some output
clear qName idxStartOne idxStartZero findOne findZero



