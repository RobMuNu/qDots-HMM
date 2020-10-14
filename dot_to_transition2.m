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
