function X = mat(x)
    % Converts non-symmetric block vector to block diagonal matrix
    
    global n ns n2is nblocks;

    X = zeros(n);
    for j = 1:nblocks
        X(block(j)) = reshape(x(n2is(j)+1:n2is(j+1)), ns(j+1), ns(j+1));
    end
    
end
