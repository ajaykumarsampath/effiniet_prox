function [ Ptree_null] = factor_apg_effiniet_null( sys_null,V,Tree)

% This function calculate the matrices factor step of the effiniet 
% cost function. 
%%
Ptree_null=struct('P',cell(1,1),'K',cell(1,1),'Phi',cell(1,1),'omega',cell(1,1),...
    'd',cell(1,1),'f',cell(1,1));

Wv=sys_null.L'*V.Wu*sys_null.L;
Gbar=sys_null.G;
Bbar=sys_null.B;
nv=size(sys_null.L,2); %reduced variable 

Ns=size(Tree.leaves,1); % no of scenarios
%Nd=size(Tree.stage,1)-Ns; % no of non-leave nodes

for k=1:Ns
    Ptree_null.P{Tree.leaves(k)+1,1}=zeros(nv);
end 


for k=sys_null.Np-1:-1:0
    if(k<1)
        Pbar=2*Wv+Ptree_null.P{2};
        Ptree_null.omega{1}=-0.5*(Pbar\eye(nv));
        Ptree_null.Phi{1,1}=Ptree_null.omega{1,1}*Gbar';
        Ptree_null.Theta{1,1}=Ptree_null.omega{1,1}*Bbar';
        Ptree_null.K{1,1}=-2*Ptree_null.omega{1}*Wv;
        Ptree_null.P{1,1}=-Wv*Ptree_null.K{1,1};
        Ptree_null.d{1,1}=Ptree_null.K{1,1}'*Gbar';
        Ptree_null.f{1,1}=Ptree_null.K{1,1}'*Bbar';
    else
        nodes_stage=find(Tree.stage==k-1);
        
        for j=1:length(nodes_stage)
            nchild=Tree.children{nodes_stage(j)};
            prob_sum=Tree.prob(nodes_stage(j))+sum(Tree.prob(nchild));
            if(k==sys_null.Np-1)
                Pbar=Tree.prob(nchild)*Wv;
            else
                Pbar=prob_sum*Wv;
                for l=1:length(nchild)
                    Pbar=Pbar+Ptree_null.P{nchild(l)+1};
                end
            end
            Ptree_null.omega{nodes_stage(j)+1,1}=-0.5*(Pbar\eye(nv));
            Ptree_null.Phi{nodes_stage(j)+1,1}=Ptree_null.omega{nodes_stage(j)+1,1}*Gbar';
            Ptree_null.Theta{nodes_stage(j)+1,1}=Ptree_null.omega{nodes_stage(j)+1,1}*Bbar';
            Ptree_null.K{nodes_stage(j)+1,1}=-2*Tree.prob(nodes_stage(j))*...
                Ptree_null.omega{nodes_stage(j)+1}*Wv;
            Ptree_null.P{nodes_stage(j)+1,1}=-Tree.prob(nodes_stage(j))*Wv*...
                Ptree_null.K{nodes_stage(j)+1,1};
            Ptree_null.d{nodes_stage(j)+1,1}=Ptree_null.K{nodes_stage(j)+1,1}'*Gbar';
            Ptree_null.f{nodes_stage(j)+1,1}=Ptree_null.K{nodes_stage(j)+1,1}'*Bbar';
        end
    end
    
end

Ptree_null.Bbar=Bbar;
Ptree_null.Gbar=Gbar;
end



