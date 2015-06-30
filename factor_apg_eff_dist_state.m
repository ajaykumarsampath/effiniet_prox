function [ Ptree] = factor_apg_eff_dist_state( sys,V,Tree)

% This function calculate the matrices factor step of the effiniet
% cost function.
%%
Ptree=struct('P',cell(1,1),'K',cell(1,1),'Phi',cell(1,1),'omega',cell(1,1),...
    'd',cell(1,1),'f',cell(1,1));

Wv=sys.L'*V.Wu*sys.L;
nv=size(sys.L,2); %reduced variable

if(sys.cell)
    %Gbar=sys.G*sys.L;
    Bbar=sys.B*sys.L;
    
    for k=sys.Np:-1:1
        nodes_stage=find(Tree.stage==k-1);
        for j=1:length(nodes_stage)
            if(k==sys.Np)
                Wv_k=Tree.prob(nodes_stage(j))*Wv;
            else
                Wv_k=Tree.prob(nodes_stage(j))*Wv+Ptree.P{nodes_stage(j)+1};
                nchild=Tree.children{nodes_stage(j)};
                for l=1:length(nchild)
                    Wv_k=Wv_k+Tree.prob(nchild(l))*Wv;
                end
            end
            Ptree.Gbar{nodes_stage(j),1}=sys.G{nodes_stage(j),1}*sys.L;
            
            Ptree.omega{nodes_stage(j),1}=-0.5*(Wv_k\eye(nv));
            Ptree.Phi{nodes_stage(j),1}=Ptree.omega{nodes_stage(j),1}*...
                (sys.F{nodes_stage(j)}*Bbar)';
            Ptree.Psi{nodes_stage(j),1}=Ptree.omega{nodes_stage(j),1}*...
                Ptree.Gbar{nodes_stage(j)}';
            Ptree.Theta{nodes_stage(j),1}=Ptree.omega{nodes_stage(j),1}*Bbar';
            
            Ptree.K{nodes_stage(j),1}=-2*Tree.prob(nodes_stage(j))*...
                (Ptree.omega{nodes_stage(j),1}*Wv);
            
            Ptree.d{nodes_stage(j),1}=Ptree.K{nodes_stage(j),1}'*...
                (sys.F{nodes_stage(j)}*Bbar)';
            Ptree.f{nodes_stage(j),1}=Ptree.K{nodes_stage(j),1}'*...
                Ptree.Gbar{nodes_stage(j)}';
            Ptree.g{nodes_stage(j),1}=Ptree.K{nodes_stage(j),1}'*Bbar';   
        end
        
        if(k==1)
            Ptree.P{1}=-Wv*Ptree.K{1};
        else
            Pnodes_stage=find(Tree.stage==k-2);
            for j=1:length(Pnodes_stage)
                nchild=Tree.children{Pnodes_stage(j)};
                for l=1:length(nchild)
                    if(l==1)
                        Ptree.P{Pnodes_stage(j)+1,1}=Tree.prob(nchild(l))*Ptree.K{nchild(l),1};
                    else
                        Ptree.P{Pnodes_stage(j)+1,1}=Ptree.P{Pnodes_stage(j)+1,1}+...
                            Tree.prob(nchild(l))*Ptree.K{nchild(l),1};
                    end
                end 
                Ptree.P{Pnodes_stage(j)+1,1}=-Wv*Ptree.P{Pnodes_stage(j)+1,1};
            end
        end 
    end
    Ptree.Bbar=sys.B*sys.L;
else
    Gbar=sys.G*sys.L;
    Bbar=sys.B*sys.L;
    Fbar=sys.F*Bbar;
     
   for k=sys.Np:-1:1
        nodes_stage=find(Tree.stage==k-1);
        for j=1:length(nodes_stage)
            if(k==sys.Np)
                Wv_k=Tree.prob(nodes_stage(j))*Wv;
            else
                Wv_k=Tree.prob(nodes_stage(j))*Wv+Ptree.P{nodes_stage(j)+1};
                nchild=Tree.children{nodes_stage};
                for l=1:length(nchild)
                    Wv_k=Wv_k+Tree.prob(nchild(l))*Wv_k;
                end
            end
            Ptree.omega{nodes_stage(j),1}=-0.5*(Wv_k\eye(nv));
            Ptree.Phi{nodes_stage(j),1}=Ptree.omega{nodes_stage(j),1}*(Fbar+Gbar)';
            Ptree.Theta{nodes_stage(j),1}=Ptree.omega{nodes_stage(j),1}*Bbar';
            Ptree.K{nodes_stage(j),1}=-2*Tree.prob(nodes_stage(j))*...
                (Ptree.omega{nodes_stage(j),1}*Wv);
            Ptree.d{nodes_stage(j),1}=Ptree.K{nodes_stage(j),1}'*(Fbar+Gbar)';
            Ptree.f{nodes_stage(j),1}=Ptree.K{nodes_stage(j),1}'*Bbar';
        end
        if(k==1)
            Ptree.P{1}=-Wv*Ptree.K{1};
        else
            Pnodes_stage=find(Tree.stage==k-2);
            for j=1:length(Pnodes_stage)
                nchild=Tree.children{Pnodes_stage(j)};
                for l=1:length(nchild)
                    if(l==1)
                        Ptree.P{Pnodes_stage(j)+1,1}=Tree.prob(nchild(l))*Ptree.K{nchild(l),1};
                    else
                        Ptree.P{Pnodes_stage(j)+1,1}=Ptree.P{Pnodes_stage(j)+1,1}+...
                            Tree.prob(nchild(l))*Ptree.K{nchild(l),1};
                    end
                end 
                Ptree.P{Pnodes_stage(j)+1,1}=-Wv*Ptree.P{Pnodes_stage(j)+1,1};
            end
        end 
    end
    Ptree.Bbar=Bbar;
    Ptree.Gbar=Gbar;
end



end



