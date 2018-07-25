function X = H(U, P)
    % Block-Wise Symmetrization Operator

    global nis nblocks;

    % apply H block-wise
    X = zeros(nis(end));
    for j = 1:nblocks
        X(block(j)) = (1/2) * (P(block(j))*U(block(j))/P(block(j)) + ...
            (P(block(j))')\(U(block(j))')*(P(block(j))'));
    end
end
