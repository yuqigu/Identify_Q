function[I] = get_I(D,M)
    
index0 = 2.^((M-1):-1:0)'; % length M
indices = index0;
for d = 2:D
    newindices = bsxfun(@plus, indices', ...
        index0(d:M)); % outer sum, size (M-d+1) * M
    indices = [indices ; newindices(:)];
    [sortedtemp, ordertemp, temp] = unique(indices, 'first');
    indices = indices(sort(ordertemp));
end
I = binary(indices,M);

end
