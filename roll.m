function out=roll(num,side,roll,trial)
% Roll a dice with the following parameters
% num number of dice to roll simultaneously
% side number of sides
% roll number of rolls
% trial number of independent trials


out = randi([1 side],[roll num trial]);