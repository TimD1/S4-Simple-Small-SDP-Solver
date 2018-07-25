function [ x ] = svec( X )
    % converts symmetric block diagonal matrix to symmetric vector
    
    global ns nis nt nts ntis nblocks;
    
    % allocate storage
    x = zeros(nt, 1);
    
    % copy diagonal elements
    for j = 1:nblocks
        for i = 1:ns(j+1)
            x(ntis(j)+i*(i+1)/2) = ...
                X(nis(j)+i, nis(j)+i);
        end
    end
    
    % copy remaining elements
    for j = 1:nblocks
        start = 1;
        for i = 1:ns(j+1)
            x(ntis(j)+start : ntis(j)+start+i-2 ) = ...
                X(nis(j)+1 : nis(j)+i-1, nis(j)+i)*sqrt(2);
            start = start+i;
        end
    end
end
