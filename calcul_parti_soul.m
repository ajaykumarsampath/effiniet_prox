function [ cur_opt] = calcul_parti_soul( sys,Tree,V,current_opts)

% This function calcualte the particular soultion 
% for the equation Eu+E_d d=0 given by 
% u=-pinv(E)Ed d;
% INPUT------    sys  :  system dynamics 
%               Tree  :  Tree structure 
%     current_demand  :  Current demand estimated
%            inv_hat  :  Previous_vhat;  
% Output---- cur_opt  :  
% current_opts=par_sol_opt;
Nd=size(Tree.stage,1);
nu=size(sys.L1,1);
nv=size(sys.L,2);

vhat=zeros(nu,Nd);
zeta=zeros(nu,Nd);
beta=zeros(nv,Nd);
w=zeros(sys.nx,Nd);

current_demand=current_opts.demand;
prev_vhat=current_opts.prev_vhat;

alpha_bar=zeros(nv,sys.Np);

Wv1=V.Wu*sys.L;

for k=1:sys.Np
    nodes_stage=find(Tree.stage==k-1);
    alpha_bar(:,k)=(V.alpha(k,:)*sys.L)';
    for j=1:length(nodes_stage)
        w(:,nodes_stage(j))=sys.Gd*(Tree.value(nodes_stage(j),:)'+current_demand(k,:)');
        vhat(:,nodes_stage(j))=sys.L1*(current_demand(k,:)'...
           +Tree.value(nodes_stage(j),:)');
    end
end

for k=1:sys.Np
    nodes_stage=find(Tree.stage==k-1);
    if(k==1)
        zeta(:,nodes_stage)=(vhat(:,nodes_stage)-prev_vhat);
        zeta(:,nodes_stage)=-(vhat(:,2)-vhat(:,nodes_stage))...
            +zeta(:,nodes_stage);
        beta(:,nodes_stage)=2*(zeta(:,nodes_stage)'*Wv1)';
    else
        for j=1:length(nodes_stage)
            zeta(:,nodes_stage(j))=Tree.prob(nodes_stage(j))*(vhat(:,nodes_stage(j))...
                -vhat(:,Tree.ancestor(nodes_stage(j))));
            if(k<sys.Np)
                nchild=Tree.children{nodes_stage(j)};
                for l=1:length(nchild)
                    zeta(:,nodes_stage(j))=Tree.prob(nchild(l))*(vhat(:,nchild(l))-...
                        vhat(:,nodes_stage(j)))+zeta(:,nodes_stage(j));
                end
            end
            beta(:,nodes_stage(j))=2*(zeta(:,nodes_stage(j))'*Wv1)';
        end
    end
end

cur_opt=current_opts;
cur_opt.beta=beta;
cur_opt.vhat=vhat;
cur_opt.alpha_bar=alpha_bar;
cur_opt.w=w;

end

