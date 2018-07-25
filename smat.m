function X = smat(x)
    % converts symmetric block vector to block diagonal matrix
    
    global n ns nis ntis nblocks;
    
    X = zeros(n);
    start = 1;
    
    % copy over top-right section
    for j = 1:nblocks
        start = 1;
        for i = 1:ns(j+1)
            X(nis(j)+1:nis(j)+i, nis(j)+i) = x(ntis(j)+start:ntis(j)+start+i-1);
            start = start+i;
        end
    end
    
    % scale and fill in bottom-left
    for i = 1:n
        for j = i+1:n
            X(i,j) = X(i,j)/sqrt(2);
            X(j,i) = X(i,j);
        end
    end
end
