function M = skmult(G, K, A)
    % Block-Wise Symmetric Kronecker Product times vector

    % calculate useful indices 
    global nblocks;
    
    M = zeros(size(A));

    % block-wise symmetric kronecker product
    for j = 1:nblocks
        for k = 1:size(A,2)
            M(vblock(j), k) = (1/2)*...
				svec_1block(K(block(j)) *...
				smat_1block(A(vblock(j), k)) *...
				G(block(j))' +...
                G(block(j)) *...
				smat_1block(A(vblock(j), k)) *...
				K(block(j))');
        end
    end
end
