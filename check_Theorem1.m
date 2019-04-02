function [id, is_complete, is_3items, is_distinct] = check_Theorem1(Q, varargin)

% This function checks if the necessary and sufficient identifiability 
% conditions of Q under DINA model are satisfied

% Q = [eye(3); 1 1 0; 1 0 1; 0 1 1; 1 1 1];

% Q = [eye(8);...
%      0 0 1 1 1 0 1 1;...
%      0 1 0 1 0 1 1 1;...
%      1 0 0 0 1 1 1 1;...
%      1 1 1 1 1 1 0 1];

% time0 = cputime;


% set defaults for optional inputs
optargs = {eps 17 @magic};

numvarargs = length(varargin);

% now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
optargs(1:numvarargs) = varargin;
% or ...
% [optargs{1:numvarargs}] = varargin{:};

% Place optional args in memorable variable names
[is_entire_Q] = optargs{:};





[J, K] = size(Q);

attr_count = sum(Q, 1);

% If some attributes is required by <=2 item, then should check Theorem 2
if any(attr_count <= 2)
    fprintf('\n Some attribute is required by <=2 items; check Theorem 2 instead.\n')
    return
end


%%%% Condition A: check if the Q-matrix is complete
[is_complete, ind_I_K] = check_complete(Q);


%%%% Condition B: check if each attribute is required by >= 3 items
is_3items = prod(sum(Q, 1) >= 3);


%%%% Condition C: check if the submatrix has K distinct columns
Q_submat = Q(setdiff(1:J, ind_I_K), :);
% count how many distinct columns Q^\star has
num_distinct = size(unique(Q_submat', 'rows'), 1);
is_distinct = (num_distinct == K);

id = is_complete * is_3items * is_distinct;

% if the Q is a submatrix as that in Theorem 2, then do not print the
% message
if is_entire_Q 
    if id==1
        fprintf('\n Under DINA, this Q-matrix IS strictly identifiable by Theorem 1.\n')
    else
        fprintf('\n Under DINA, this Q-matrix IS NOT strictly identifiable by Theorem 1.\n')
    end
end

% mytime = cputime - time0;

end