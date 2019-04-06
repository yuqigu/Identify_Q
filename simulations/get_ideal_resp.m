function [ir] = get_ideal_resp(Q, Z)

% This function obtains the J * 2^K ideal response matrix given Q-matrix
% and attribute patterns as row vectors of the matrix Z

% @param Q:   Q-matrix of size J * K
% @param Z:   latent attribute profiles of individuals; size N * K

% @return ir: binary ideal response matrix of size N * J



[J, K] = size(Q);
ir = prod(bsxfun(@power, reshape(Z, [1 size(Z,1) K]), reshape(Q, [J 1 K])), 3)';

ir = ir';


end