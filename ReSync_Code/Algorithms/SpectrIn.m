%% Spectral Initialization for Rotation Synchronization
%%------------------------------------------------
%% Input Parameters: 
%% Ind: edge_num by 2 "edge indices matrix". Each row is the index of an edge (i,j). that is sorted as (1,2), (1,3), (1,4),... (2,3), (2,4),.... 
%% edge_num is the number of edges.
%% RijMat: 3 by 3 by edge_num tensor that stores the given relative rotations corresponding to Ind


%% Output:
%% R_est: Estimated rotations (3x3xn)



function R_est = SpectrIn(Ind,RijMat)

 % building the graph   
    Ind_i = Ind(:,1);
    Ind_j = Ind(:,2);
    n=max(Ind,[],'all');
    m=size(Ind_i,1);
    AdjMat = sparse(Ind_i,Ind_j,1,n,n); % Adjacency matrix
    AdjMat = full(AdjMat + AdjMat');
    d=3;
    mat_size = ones(1,n)*d;
    cum_ind = [0,cumsum(mat_size)];
    Rij_blk = zeros(n*d);
    for k = 1:m
       i = Ind_i(k); j=Ind_j(k);
       Rij_blk((cum_ind(i)+1):cum_ind(i+1), (cum_ind(j)+1):cum_ind(j+1))= RijMat(:,:,k);    
    end

    Rij_blk = Rij_blk+Rij_blk';

    % Spectral 
    [V,~] = eigs(Rij_blk,d, 'la');
    % comstruct U by reversing the sign of the last column of V so that 
    % the determinants of U and V differ by a sign.
    U = V * diag([ones(1,d-1), -1]);

    R_est = zeros(d,d,n);
    normR = 0;
    for i=1:n
       Ri = V((cum_ind(i)+1):cum_ind(i+1), :); 
       [Ur,Lr,Vr] = svd(Ri);
       S0 = diag([ones(1,d-1),det(Ur*Vr')]);
       normR = normR + norm(S0 - Lr,"fro")^2; % the distance dist(V, SO(d)^n)
       R_est(:,:,i) = Ur*S0*Vr';  % the projection on SO(d)
    end

    Q_est = zeros(d,d,n);
    normQ = 0;
    for i=1:n
       Ri = U((cum_ind(i)+1):cum_ind(i+1), :); 
       [Ur,Lr,Vr] = svd(Ri);
       S0 = diag([ones(1,d-1),det(Ur*Vr')]);
       normQ = normQ + norm(S0 - Lr,"fro")^2;  % the distance dist(U, SO(d)^n)
       Q_est(:,:,i) = Ur*S0*Vr';  % the projection on SO(d)
    end
    
    % choose eithor U or V based on the order between dist(V, SO(d)^n) and dist(U, SO(d)^n)
    if normR > normQ 
        R_est = Q_est;
    end

   fprintf('Initialization completed! \n \n'); 

end