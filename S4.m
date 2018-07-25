function [] = S4(filename)
	%
	% S4: Simple Small SDP Solver
    % Author: Tim Dunn, 2018
	%
    % Usage: S4 <SDP filename>

	% get data from file
	tic;
	if nargin == 0
		filename = './problems/prob1.dat-s';
	end
	file = fopen(filename);
	data = fscanf(file, '%f');

	global n ns nis n2 n2s n2is nt nts ntis nblocks;

	% m = number of constraints
	m = data(1);

	% nblocks = number of blocks
	nblocks = data(2);

	% 1D indices
	% n's = 1D sizes (n) of each nxn block (0 prefix)
	ns = data(3:3+nblocks-1);
	ns = [ 0; abs(ns)];
	
    % n = sum( n's )
	n = sum(ns);
	
    % n_i's = accumulative 1D indices of each block
	nis = zeros(size(ns));
	for i = 1:length(ns)
		nis(i) = sum(ns(1:i));
	end

	% 2D indices
	% n^2's = 2D sizes (n*n) of each nxn block (0 prefix)
	n2s = ns.^2;
	
    % n2 = sum( n^2's )
	n2 = sum(n2s);
	
    % n^2_i's = accumulative 2D indices of each block
	n2is = zeros(size(ns));
	for i = 1:length(ns)
		n2is(i) = sum(n2s(1:i));
	end

	% svec indices
	% nt's = triangle numbers (n*(n+1)/2) for size of svec(nxn block) (0 prefix)
	nts = ns.*(ns+1)/2;
	
    % nt_i's = accumulative triangle-number indices of svec(nxn block) 
	ntis = zeros(size(nts));
	for i = 1:length(nts)
		ntis(i) = sum(nts(1:i));
    end
    
    % nt = sum(nt's) = sum( n_i*(n_i+1)/2 )
	nt = sum(nts);
    

	% read data from file and allocate space
	b = -data(3+nblocks:3+nblocks+m-1);
	data = data(3+nblocks+m:end);
	c = zeros(n2, 1);
	A = zeros(m, n2);

	% read in constraint matrices
	for idx = 1:5:size(data,1)-1
		matno = data(idx);
		blkno = data(idx+1);
		i = data(idx+2);
		j = data(idx+3);
		entry = data(idx+4);
		id1 = (j-1)*ns(blkno+1) + i + n2is(blkno); 
		id2 = (i-1)*ns(blkno+1) + j + n2is(blkno);
		if matno == 0
			c(id1) = -entry;
			c(id2) = -entry;
		else
			A(matno, id1) = -entry;
			A(matno, id2) = -entry;
		end
	end

	% convert A and c to svec versions
	c = svec(mat(c));
	A_temp = A;
	A = zeros(m, nt);
	for i = 1:m
		A(i,:) = svec(mat(A_temp(i,:)));
	end

	% set ksi
	ksi = zeros(nblocks,1);
	temp = zeros(m, 1);
	for j = 1:nblocks
		for k = 1:m
			temp(k) = (1 + abs(b(k))) / (1 + norm(A(k, vblock(j))));
		end
		ksi(j) = ns(j+1) * max(temp);
	end

	% set eta
	eta = zeros(nblocks, 1);
	temp = zeros(m,1);
	for j = 1:nblocks
		for k = 1:m
			temp(k) = norm(A(k, vblock(j)));
		end
		eta(j) = (1/sqrt(ns(j+1))) * (1 + max( max(temp), norm(c(vblock(j)))));
	end

	% set initial iterates using ksi and eta
	X = zeros(n);
	for j = 1:nblocks
		X(block(j)) = eye(ns(j+1))*ksi(j);
	end
	x = svec(X);
	y = zeros(m,1);
	Z = zeros(n);
	for j = 1:nblocks
		Z(block(j)) = eye(ns(j+1))*eta(j);
	end
	z = svec(Z);

	% calculate normA for each block, save sqrt of each
	normsA = zeros(nblocks,1);
	for j = 1:nblocks
		normsA(j) = norm(A(:, vblock(j)), 'fro');
	end
	normsA = max(1, sqrt(normsA));
		
	% calculate normb
	normb = max(1, norm(b));

	% calculate normc for each block, save max
	normc = 1;
	for j = 1:nblocks
		normc = max(normc, norm(c(vblock(j))));
	end

	% scale constraints
	for j = 1:nblocks
		A(:, vblock(j)) = A(:, vblock(j)) / normsA(j);
		c(vblock(j)) = c(vblock(j)) / (normc*normsA(j));
	end
	b = b / normb;

	% scale initial iterates x and z
	for j = 1:nblocks
		x(vblock(j)) = x(vblock(j)) * normsA(j);
		z(vblock(j)) = z(vblock(j)) / (normc*normsA(j));
	end

	% set initial values
	rp = b - A*x;
	Rd = c - z -A'*y;
	relgap = dot(x,z) / (1 + max(abs(dot(c,x)), abs(dot(b,y))));
	pinfeas = norm(rp) / (1 + norm(b));
	dinfeas = norm(Rd) / (1 + norm(c));
	phi = max(pinfeas, dinfeas);
	soln_relgap = 1e-6;
	soln_phi = 1e-6;
	g = 0.9;
	iter = 0;

	% print header
	fprintf('Problem Filename: %s\n', filename);
	fprintf('iter\talpha\tbeta\tsigma\t\tmu\t\te\tpinfeas\t   dinfeas    relgap\t primal\t     dual\t phi\t    gamma rp\t\tRd\n');

	% update intermediate solution until convergence
	while(relgap > soln_relgap || phi > soln_phi)    
		
		% get Cholesky decompositions of X and Z (upper triangular!)
		X = smat(x);
		Q = chol(X);
		Z = smat(z);
		P = chol(Z);
		
		% calculate residuals, relgap, feasability
		rp = b - A*x;
		Rd = c - z -A'*y;
		mu = dot(x,z) / n;
		relgap = dot(x,z) / (1 + max(abs(dot(c,x)), abs(dot(b,y))));
		phi = max(norm(rp)/(1+norm(b)), norm(Rd)/(1+norm(c)));
		Rc = -svec( H(X*Z, P) );
		
		% solve system (predictor)   
		M = A*skmult(X, Z\eye(n), A');
		h = rp + A*(skmult(X, Z\eye(n), Rd) - skmult(P\eye(n), P\eye(n), Rc));
		dy = M\h;
		dz = Rd - A'*dy;
		dx = -x - skmult(X, Z\eye(n), dz);
		
		% calculate alpha (primal step-length)
		alphas = zeros(1,nblocks);
		for j = 1:nblocks
			XDX = (Q(block(j))') \ smat_1block(dx(vblock(j))) / Q(block(j));
			eigs_a = sort(eig(XDX));
			lambda_a = eigs_a(1);
			if lambda_a < 0
				alphas(j) = -1/lambda_a;
			else
				alphas(j) = Inf;
			end
		end
		alpha = min(1, g*min(alphas));
		
		% calculate beta (dual step-length)
		betas = zeros(1,nblocks);
		for j = 1:nblocks
			ZDZ = (P(block(j))') \ smat_1block(dz(vblock(j))) / P(block(j));
			eigs_b = sort(eig(ZDZ));
			lambda_b = eigs_b(1);
			if lambda_b < 0
				betas(j) = -1/lambda_b;
			else
				betas(j) = Inf;
			end
		end
		beta = min(1, g*min(betas));
		
		% calculate exponent (for sigma)
		if mu > 1e-6
			if min(alpha, beta) < 1 / sqrt(3)
				e = 1;
			else
				e = max(1, 3*min(alpha, beta)^2);
			end
		else
			e = 1;
		end
		
		% calculate sigma (centering parameter)
		if dot(x+alpha*dx, z+beta*dz) < 0
			sigma = 0.8;
		else
			frac = dot(x+alpha*dx, z+beta*dz)/dot(x,z);
			sigma = min(1, frac^e);
		end
		
		% solve modified system (corrector)
		Rc = svec(sigma*mu*eye(n)-H(X*Z, P)-H(smat(dx)*smat(dz), P));
		h = rp  + A*(skmult(X, Z\eye(n), Rd) - skmult(P\eye(n), P\eye(n), Rc));
		dy = M\h;
		dz = Rd - A'*dy;
		dx = skmult(P\eye(n), P\eye(n), Rc) - skmult(X, Z\eye(n), dz);
		g = 0.9 + 0.09*min(alpha, beta);
		
		% calculate alpha (primal step-length)
		for j = 1:nblocks
			XDX = (Q(block(j))') \ smat_1block(dx(vblock(j))) / Q(block(j));
			eigs_a = sort(eig(XDX));
			lambda_a = eigs_a(1);
			if lambda_a < 0
				alphas(j) = -1/lambda_a;
			else
				alphas(j) = Inf;
			end
		end
		alpha = min(1, g*min(alphas));
		
		% calculate beta (dual step-length)
		for j = 1:nblocks
			ZDZ = (P(block(j))') \ smat_1block(dz(vblock(j))) / P(block(j));
			eigs_b = sort(eig(ZDZ));
			lambda_b = eigs_b(1);
			if lambda_b < 0
				betas(j) = -1/lambda_b;
			else
				betas(j) = Inf;
			end
		end
		beta = min(1, g*min(betas));
		
		% update solution
		x = x + alpha*dx;
		y = y + beta*dy;
		z = z + beta*dz;
		g = 0.9 + 0.09*min(alpha, beta);
		
		% print results of iteration
		iter = iter + 1;
		primal = dot(c,x);
		dual = dot(b,y);
		pinfeas = norm(rp) / (1 + norm(b));
		dinfeas = norm(Rd) / (1 + norm(c));  
		fprintf('%d\t%.3f\t%.3f\t%.3e\t%.3e   %.3f\t%.2e   %.2e   %.2e   %.2e   %.2e   %.2e   %.2f  %2e  %2e\n', ...
			iter, alpha, beta, sigma, mu, e, pinfeas, dinfeas, relgap, primal, dual, phi, g, norm(rp), norm(Rd));
	end
	disp(' ');
	disp('Solver has converged!');

	% un-scale constraints
	for j = 1:nblocks
		A(:, vblock(j)) = A(:, vblock(j)) * normsA(j);
		c(vblock(j)) = c(vblock(j)) * (normc*normsA(j));
	end
	b = b * normb;

	% un-scale solutions
	for j = 1:nblocks
		x(vblock(j)) = x(vblock(j)) * normb / normsA(j);
		z(vblock(j)) = z(vblock(j)) * (normc*normsA(j));
	end
	y = y * normc;
	X = smat(x);
	Z = smat(z);

	% display results
	disp(['Primal: ', num2str(dot(c,x))]);
	disp(['Dual:   ', num2str(dot(b,y))]);
	time = toc;
	disp(['Time:   ', num2str(time), ' seconds']);
end
