% Generate 2k encoded speeches from photon trajectory
% INPUT: Binarised photon trajectories with some common fname pattern
% OUTPUT: two 2k encoded speeches (on/off each) plus their alphabets

% NOTES: We either double the alphabet size, or reduce our data
% by some factor. The latter is far less detrimental for Lmax so that
% is what we are doing.



% Read all datafiles in the directory and store the names in an cell
datInfo=dir('*-cen');
datNames={datInfo.name};

% Define the full alphabet (covers seq lengths up to 62 million lol)
alphaFull = 'a':'z';

% Make a results folder
[yr,mo,da]=ymd(datetime);
[hr,mn,se]=hms(datetime);
%mkdir CG2kData

% For each GY photon traj:
for dat = 1:numel(datNames)
	traj = fileread(datNames{dat}); 

    % Quickly convert the string into numerical array
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
	% Take the alphabet at those bin indices and write in the symbols eg. a(3)-->'c'
	% ------------------------------------
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

	fid = fopen('CGcheckAllGY.txt', 'w+');
    fid0 = fopen('CGlambda0.txt', 'a+');
    fid1 = fopen('CGlambda1.txt', 'a+');
	% If Lambda >= |A| then mark data as “2^k CG safe”
	if (LambdaOne >= sizeCG2kalphaOne) & (LambdaZero >= sizeCG2kalphaZero)
		% Then dataset is fully safe for ON
		safeState = convertCharsToStrings(datNames{dat})+'   ----- 2 safe';
		fprintf(fid, '%s\n', safeState);

	elseif (LambdaOne < sizeCG2kalphaOne) & (LambdaZero >= sizeCG2kalphaZero)
		safeState = convertCharsToStrings(datNames{dat})+'   ----- Off safe';
		fprintf(fid, '%s\n', safeState);

	elseif (LambdaOne >= sizeCG2kalphaOne) & (LambdaZero < sizeCG2kalphaZero)
		safeState = convertCharsToStrings(datNames{dat})+'   ----- On safe';
		fprintf(fid, '%s\n', safeState);

	elseif (LambdaOne < sizeCG2kalphaOne) & (LambdaZero < sizeCG2kalphaZero)
		safeState = convertCharsToStrings(datNames{dat})+'   ----- unsafe';
		%fprintf(fid, '%s\n', safeState);
		fprintf(fid, '%s\n', convertCharsToStrings(datNames{dat}));
		fprintf(fid,...
		 'On state: \tData length=%d\n\t\tMax Lambda=%.4f\n\t\t|A|=%d\n',...
		 nCGdatOne, LambdaOne, sizeCG2kalphaOne);
		fprintf(fid,...
		 'Off state: \tData length=%d\n\t\tMax Lambda=%.4f\n\t\t|A|=%d\n\n',...
		 nCGdatZero, LambdaZero, sizeCG2kalphaZero);
     
  
        fprintf(fid0, '%d\n',floor(LambdaZero));
        fprintf(fid1, '%d\n',floor(LambdaOne));
	end
	fclose(fid);
    fclose(fid0);
    fclose(fid1);
 	
 	% Write compressed speeches to datafiles	
%     qName = convertCharsToStrings(datNames{dat});
% 	dlmwrite(qName+'-s1', CG2kOne,''); % on
% 	dlmwrite(qName+'-s0', CG2kZero,''); % off
% 	dlmwrite(qName+'-alpha1', CG2kalphaOne,''); % alphabet for 1
% 	dlmwrite(qName+'-alpha0', CG2kalphaZero,''); % alphabet for 0

end

% Clean up some output
clear qName fid idxStartOne idxStartZero safeState findOne findZero



