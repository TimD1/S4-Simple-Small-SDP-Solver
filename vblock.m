function [vblock_indices] = vblock(blkno)
	% Calculates indices of a larger vec which belong to a particular block
	% A(vblock(i))  is equivalent to A(ntis(i)+1:ntis(i+1))

    global ntis;

    vblock_indices = ntis(blkno)+1:ntis(blkno+1);
            
return 
