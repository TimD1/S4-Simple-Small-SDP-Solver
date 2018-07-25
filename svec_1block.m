function [ x ] = svec_1block( X )
	% converts symmetric matrix to symmetric vector
    
    % allocate storage
    x = zeros(size(X,1)*(size(X,1)+1)/2, 1);
    
    % copy diagonal elements
    for i = 1:size(X,1)
        x(i*(i+1)/2) = X(i,i);
    end
    
    % scale and copy remaining elements
    for i = 1:size(X,1)-1
        x(1+i*(i+1)/2 : (i+1)*(i+2)/2 - 1) = X(1:i,i+1)*sqrt(2);
    end
end
