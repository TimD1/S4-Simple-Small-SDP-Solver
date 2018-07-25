function [block_indices] = block(blkno)
	% Calculates indices of a larger matrix which belong to a particular block
	% A(block(i)) is now equivalent to A(nis(i)+1:nis(i+1), nis(i)+1:nis(i+1))

    global n nis;

    block_indices = sub2ind([n, n],...
		ndgrid(nis(blkno)+1:nis(blkno+1), nis(blkno)+1:nis(blkno+1)),...
		ndgrid(nis(blkno)+1:nis(blkno+1), nis(blkno)+1:nis(blkno+1))');
            
return 
