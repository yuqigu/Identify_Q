% K=3
Q = [eye(3); 1 1 0; 1 0 1; 0 1 1; 1 1 1];

[stric_id, ~, ~, ~] = check_Theorem1(Q)

[generic_id, local_generic_id] = check_Theorem2(Q)

[cond_hold, item_set_1, item_set_2, item_set_3] = check_Theorem4(Q)

% K=8
Q = [eye(8);...
     0 0 1 1 1 0 1 1;...
     0 1 0 1 0 1 1 1;...
     1 0 0 0 1 1 1 1;...
     1 1 1 1 1 1 0 1];
 
 
 check_Theorem1(Q)
 
 check_Theorem2(Q)