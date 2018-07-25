function X = smat_1block(x)
    % converts a symmetric vector to a matrix

    n = floor(sqrt(2*length(x)));
    X = zeros(n);
    start = 1;
    
    % copy over top-right section
    for i = 1:n
        X(1:i, i) = x(start:start+i-1);
        start = start+i;
    end
    
    % scale and fill in bottom-left
    for i = 1:n
        for j = i+1:n
            X(i,j) = X(i,j)/sqrt(2);
            X(j,i) = X(i,j);
        end
    end
end
