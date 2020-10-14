function p = pr_heavy(o,a,pi)
% INPUT
% o = vector of emission symbols
% a = transition tensor (rows sum to 1, third index marks emission transitions)
% pi = initial distribution over states
% OUTPUT
% p = P(emission sequence | model)

% get number of states
nState = length(a(1,:,1));
% get length of observation seq
L = length(o);


% Initialise 
for i = 1:nState
    fwd(1,i) = pi(i);
end

% Do the recursion
for t = 2:(L+1)
    
    for j = 1:nState
        pathSum = 0;
        
        for i = 1:nState
            pathSum = pathSum + fwd(t-1,i)*a(i,j,o(t-1));
        end
        
        fwd(t,j) = pathSum;
    end
     
end

p = 0;
for i = 1:nState
    p = p+fwd(L+1,i);
end


end