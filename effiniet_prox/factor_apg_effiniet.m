function [ Ptree] = factor_apg_effiniet( sys,V,Tree)

% This function calculate the matrices factor step of the effiniet
% cost function.
%%
Ptree=struct('P',cell(1,1),'K',cell(1,1),'Phi',cell(1,1),'omega',cell(1,1),...
    'd',cell(1,1),'f',cell(1,1));

Wv=sys.L'*V.Wu*sys.L;
nv=size(sys.L,2); %reduced variable

Ns=size(Tree.leaves,1); % no of scenarios
%Nd=size(Tree.stage,1)-Ns; % no of non-leave nodes

if(sys.cell)
    %Gbar=sys.G*sys.L;
    Bbar=sys.B*sys.L;
    
    for k=1:Ns
        Ptree.P{Tree.leaves(k)+1,1}=zeros(nv);
        Ptree.Gbar{Tree.leaves(k)+1,1}=sys.G{Tree.leaves(k)+1,1}*sys.L;
    end
    
    for k=sys.Np-1:-1:0
        if(k<1)
            Ptree.Gbar{1,1}=sys.G{1,1}*sys.L;
            Pbar=2*Wv+Ptree.P{2};
            Ptree.omega{1}=-0.5*(Pbar\eye(nv));
            Ptree.Phi{1,1}=Ptree.omega{1,1}*Ptree.Gbar{1,1}';
            Ptree.Theta{1,1}=Ptree.omega{1,1}*Bbar';
            Ptree.K{1,1}=-2*Ptree.omega{1}*Wv;
            Ptree.P{1,1}=-Wv*Ptree.K{1,1};
            Ptree.d{1,1}=Ptree.K{1,1}'*Ptree.Gbar{1,1}';
            Ptree.f{1,1}=Ptree.K{1,1}'*Bbar';
        else
            nodes_stage=find(Tree.stage==k-1);
            
            for j=1:length(nodes_stage)
                nchild=Tree.children{nodes_stage(j)};
                prob_sum=Tree.prob(nodes_stage(j))+sum(Tree.prob(nchild));
                if(k==sys.Np-1)
                    Pbar=Tree.prob(nchild)*Wv;
                else
                    Pbar=prob_sum*Wv;
                    for l=1:length(nchild)
                        Pbar=Pbar+Ptree.P{nchild(l)+1};
                    end
                end
                Ptree.Gbar{nodes_stage(j)+1,1}=sys.G{nodes_stage(j)+1,1}*sys.L;
                Ptree.omega{nodes_stage(j)+1,1}=-0.5*(Pbar\eye(nv));
                Ptree.Phi{nodes_stage(j)+1,1}=Ptree.omega{nodes_stage(j)+1,1}*...
                    Ptree.Gbar{nodes_stage(j)+1,1}';
                Ptree.Theta{nodes_stage(j)+1,1}=Ptree.omega{nodes_stage(j)+1,1}*Bbar';
                Ptree.K{nodes_stage(j)+1,1}=-2*Tree.prob(nodes_stage(j))*...
                    Ptree.omega{nodes_stage(j)+1}*Wv;
                Ptree.P{nodes_stage(j)+1,1}=-Tree.prob(nodes_stage(j))*Wv*...
                    Ptree.K{nodes_stage(j)+1,1};
                Ptree.d{nodes_stage(j)+1,1}=Ptree.K{nodes_stage(j)+1,1}'*...
                    Ptree.Gbar{nodes_stage(j)+1,1}';
                Ptree.f{nodes_stage(j)+1,1}=Ptree.K{nodes_stage(j)+1,1}'*Bbar';
            end
        end
        
    end
    Ptree.Bbar=Bbar;
    %Ptree.Gbar=Gbar;
else
    Gbar=sys.G*sys.L;
    Bbar=sys.B*sys.L;
    
    for k=1:Ns
        Ptree.P{Tree.leaves(k)+1,1}=zeros(nv);
    end
    
    for k=sys.Np-1:-1:0
        if(k<1)
            Pbar=2*Wv+Ptree.P{2};
            Ptree.omega{1}=-0.5*(Pbar\eye(nv));
            Ptree.Phi{1,1}=Ptree.omega{1,1}*Gbar';
            Ptree.Theta{1,1}=Ptree.omega{1,1}*Bbar';
            Ptree.K{1,1}=-2*Ptree.omega{1}*Wv;
            Ptree.P{1,1}=-Wv*Ptree.K{1,1};
            Ptree.d{1,1}=Ptree.K{1,1}'*Gbar';
            Ptree.f{1,1}=Ptree.K{1,1}'*Bbar';
        else
            nodes_stage=find(Tree.stage==k-1);
            
            for j=1:length(nodes_stage)
                nchild=Tree.children{nodes_stage(j)};
                prob_sum=Tree.prob(nodes_stage(j))+sum(Tree.prob(nchild));
                if(k==sys.Np-1)
                    Pbar=Tree.prob(nchild)*Wv;
                else
                    Pbar=prob_sum*Wv;
                    for l=1:length(nchild)
                        Pbar=Pbar+Ptree.P{nchild(l)+1};
                    end
                end
                Ptree.omega{nodes_stage(j)+1,1}=-0.5*(Pbar\eye(nv));
                Ptree.Phi{nodes_stage(j)+1,1}=Ptree.omega{nodes_stage(j)+1,1}*Gbar';
                Ptree.Theta{nodes_stage(j)+1,1}=Ptree.omega{nodes_stage(j)+1,1}*Bbar';
                Ptree.K{nodes_stage(j)+1,1}=-2*Tree.prob(nodes_stage(j))*...
                    Ptree.omega{nodes_stage(j)+1}*Wv;
                Ptree.P{nodes_stage(j)+1,1}=-Tree.prob(nodes_stage(j))*Wv*...
                    Ptree.K{nodes_stage(j)+1,1};
                Ptree.d{nodes_stage(j)+1,1}=Ptree.K{nodes_stage(j)+1,1}'*Gbar';
                Ptree.f{nodes_stage(j)+1,1}=Ptree.K{nodes_stage(j)+1,1}'*Bbar';
            end
        end
        
    end
    Ptree.Bbar=Bbar;
    Ptree.Gbar=Gbar;
end



end

