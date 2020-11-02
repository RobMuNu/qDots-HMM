% Script that will generate all EXP-CG datasets 
% iterated over some exponential base (no overflow)

% Output will be sent to a results file
% That will contain labelled exp-CG data and the respective alphabet
% so that we can readily run CSSR over all of them


datInfo = dir('*-cen');
datNames = {datInfo.name};

% Make results directory
resName = horzcat('resCGgen-',datestr(now,30));
mkdir(resName);

for dat = 1:numel(datNames)
	qName = convertCharsToStrings(datNames{dat});
	fprintf(1, 'Now reading %s\n',datNames{dat});
	traj = fileread(datNames{dat});
	traj = traj-'0';


	for b = 2:10
		fprintf(1, 'Exponential coarse graining at base %d\n', b);

		outCG  = CGgeneral(traj,b,'ex',0);
		datCG = outCG{1,2};
		alphaCG = unique(datCG);
		lamCG = floor(log2(numel(datCG))/log2(numel(alphaCG)));


		cd(fullfile(pwd,resName))
		fidAlpha = fopen(qName+'-ALPHA', 'a+');
		fidLambda = fopen(qName+'-LAMBDA', 'a+');
		fprintf(fidAlpha, '%s\n', alphaCG);
		fprintf(fidLambda, '%d\n', lamCG);

		nameCG = sprintf(qName+'-EXb%d',b);
		dlmwrite(nameCG, datCG, '');

		fclose(fidAlpha);
		fclose(fidLambda);
		cd ..
		
	end % exp base


end % datafile in dir


clear dat b traj outCG nameCG