function x = CSpeC(L1,cortex1,Betta,Cs,P_M,k,Epsil,pow_flag,elec_index)
% CSpeC: Convex Spectral Clustering

% inputs:

% L1: leadfield matrix
% cortex1: discretized cortex vertices and faces
% Betta: source contribution measure
% Cs: source covariance matrix wo silence
% P_M: estimated scalp power 
% k: estimated number of silence
% Epsil: upper bound on the power matching error
% pow_flag: 0-> no power constraint
% elec_index: indices of \phi selected electrodes in \mathcal{S}_{elec}^{(r)}

% output:

% x: solution vector of the CSpeC optimization

% 
% 
% Examples: 
% for high-res brain model:
%   x = CSpeC(L1,cortex1,Betta_new,Cs,P_M,k,Epsil,pow_flag,elec_index);
% for low-res brain model:
%   x = CSpeC(L1,cortex1,Betta_new,[],[],k,[],pow_flag,[]);
% 
% See also: README.txt and [1] for more info.

% [1] A. Chamanzar, M. Behrmann, and P. Grover, 
%   "Neural silences can be localized using noninvasive scalp EEG", 
%   To be submitted to Nature BME, 2020.

% Author: Alireza 	Date: 2020/05/02 23:55:48 	Revision: 0.1 
% Copyright: This project is licensed - see the LICENSE.md file for details

%% Relaxed RatioCut partitioning based on a knn graph:

p = size(L1,2);
lambda_contig = logspace(0,2,5) * (ones(p,1)'*Betta);


Err = zeros(size(lambda_contig,2),1);
Err_Cont = zeros(size(lambda_contig,2),1);

Loc_graph = cortex1.vertices;% Contructing the knn graph based on the 3D locations 
% of sources in the discretized brain model:
A = zeros(p,p);
KK = k;
% Betta = Betta/(ones(p,1)'*Betta);
Sim = Loc_graph;% similarity function of the knn graph

for i=1:p
    for j=i:p
        A(i,j) = norm(Sim(i,:)-Sim(j,:),2);
    end
end
A = A + permute(A,[2,1]);

W = zeros(p,p);% Weighted adjacency matrix 
for i=1:p
    [~,inds_H] = sort(A(i,:),'ascend');
    [~,inds_V] = sort(A(:,i),'ascend');
    inds_H = inds_H(1:KK+1);%KK+1 because of self edges
    inds_V = inds_V(1:KK+1);
    inds = unique([inds_H,inds_V']);
    W(i,inds) = A(i,inds);
    W(inds,i) = A(inds,i);
end
W_temp = reshape(W,[],1);
W_temp(W_temp==0) = [];
Sigm = var(W_temp);

W = exp(-(W.^2)/(2*Sigm));
W(W==1) = 0;

D = diag(sum(W,2));%Graph Degree matrix
L = D - W;%Graph Laplacian matrix
X = [];

if (pow_flag && k>1)%for the high-res with power constraints
    for ll = 1:size(lambda_contig,2)
        cvx_begin
        
        variables x(p)
        minimize((Betta'*(ones(p,1)-x))+lambda_contig(ll)*((ones(p,1)-x)'*L*(ones(p,1)-x)))
        
        subject to
        x <= ones(p,1)
        x >= zeros(p,1)
        norm(x,1) <= p-k
        
        for r = 1:length(elec_index)
                a_r = L1(elec_index(r),:);
                A_tilda_sq_r = diag(a_r);
                sum_square_abs(((ones(p,1))'*A_tilda_sq_r*Cs*A_tilda_sq_r*(x))-(P_M(elec_index(r),elec_index(r)))) <= Epsil(r)
        end
        cvx_end
        
        X = cat(2,X,x);
        Err(ll) = (Betta'*(ones(p,1)-x));
        Err_Cont(ll) = ((ones(p,1)-x)'*L*(ones(p,1)-x));
        

    end
elseif (pow_flag && k==1)%for the high-res with power constraints
    cvx_begin
    
    variables x(p)
    minimize((Betta'*(ones(p,1)-x)))
    
    subject to
    x <= ones(p,1)
    x >= zeros(p,1)
    norm(x,1) <= p-k
    
    for r = 1:length(elec_index)
            a_r = L1(elec_index(r),:);
            A_tilda_sq_r = diag(a_r);
            sum_square_abs(((ones(p,1))'*A_tilda_sq_r*Cs*A_tilda_sq_r*(x))-(P_M(elec_index(r),elec_index(r)))) <= Epsil(r)
    end
    cvx_end
    return
elseif (~pow_flag && k>1)%for the low-res wo power constraints
    for ll = 1:size(lambda_contig,2)
        cvx_begin
        
        variables x(p)
        minimize((Betta'*(ones(p,1)-x))+lambda_contig(ll)*((ones(p,1)-x)'*L*(ones(p,1)-x)))
        
        subject to
        x <= ones(p,1)
        x >= zeros(p,1)
        norm(x,1) <= p-k
        
        cvx_end
        
        X = cat(2,X,x);
        Err(ll) = (Betta'*(ones(p,1)-x));
        Err_Cont(ll) = ((ones(p,1)-x)'*L*(ones(p,1)-x));
    end
elseif (~pow_flag && k==1)%for the low-res wo power constraints
    cvx_begin
    
    variables x(p)
    minimize((Betta'*(ones(p,1)-x)))
    
    subject to
    x <= ones(p,1)
    x >= zeros(p,1)
    norm(x,1) <= p-k
    
    cvx_end
    return
end
    

%choosing the best lambda based on the total normalized error of source contribution and contiguity: 
max_Err = max(Err);
max_Err_Cont = max(Err_Cont);

dx = (Err./max_Err).^2 + (Err_Cont./max_Err_Cont).^2;
[~,idx_dx] = min(dx);
indx_Lambda = idx_dx;

x = X(:,indx_Lambda);


end

