function [ Z,details] = solve_apg_effinet_modified( sys_box,Tree,Ptree_box,Y,x,opts)

%% This function performs the solve step for the DWNs problem 
% 
% INPUT-----   sys     :  system dynamics 
%             Tree     :  Tree structure 
%                V     :  Cost function
%            Ptree     :  Factor step matrices 
%                y     :  Dual variable 
%        init_cond     :  x_0=p
%                         v_{-1}=q;
%                         \hat{v}_{-1}=\hat{q};
%
% OUTPUT-----   Z   : X states 
%                   : U control/Input
%         details   :

% opts=opts_apg.state;
Nd=size(Tree.stage,1);% total number of nodes including leaves 
Ns=size(Tree.leaves,1);% total number of leaves
nv=size(sys_box.L,2);


Z.X=zeros(sys_box.nx,Nd+1);
Z.U=zeros(sys_box.nu,Nd);
v=zeros(nv,Nd);

q=zeros(sys_box.nx,Nd+1);
r=zeros(nv,Nd+1);
sigma=zeros(nv,Nd+1);
%sum_q=zeros(sys.nx,1);
%sum_r=zeros(sys.nu,1);
y=Y.y;
yt=Y.yt;
nyt=size(yt,1);


for i=1:Ns
    q(:,Tree.leaves(i)+1)=sys_box.A'*sys_box.F{Nd+1+i}(1:nyt,:)'*yt(:,i);
end 

for k=sys_box.Np:-1:1
    nodes_stage=find(Tree.stage==k-1);
    if(k==1)
        sigma(:,1)=r(:,2)+opts.beta(:,1);
        v(:,1)=Ptree_box.Phi{1}*y(:,1)+Ptree_box.Theta{1}*q(:,2)+...
            Ptree_box.omega{1}*sigma(:,1);
        
        r(:,1)=Ptree_box.K{1}*sigma(:,1)+Ptree_box.d{1}*y(:,1)+Ptree_box.f{1}*...
            q(:,2);
        q(:,1)=sys_box.F{1}'*y(:,1)+sys_box.A'*q(:,2);
    else
        Pnodes_stage=find(Tree.stage==k-2);
        for j=1:length(nodes_stage)
            sigma(:,nodes_stage(j))=r(:,nodes_stage(j)+1)+opts.beta(:,nodes_stage(j));
            v(:,nodes_stage(j))=Ptree_box.Phi{nodes_stage(j)}*y(:,nodes_stage(j))+...
                Ptree_box.Theta{nodes_stage(j)}*q(:,nodes_stage(j)+1)+...
                Ptree_box.omega{nodes_stage(j)}*sigma(:,nodes_stage(j));
        end
        
        for j=1:length(Pnodes_stage)
            nchild=Tree.children{Pnodes_stage(j)};
            for l=1:length(nchild)
                if(l==1)
                    r(:,Pnodes_stage(j)+1)=Ptree_box.K{nchild(l)}*sigma(:,nchild(l))+...
                        Ptree_box.d{nchild(l)}*y(:,nchild(l))+Ptree_box.f{nchild(l)}*q(:,nchild(l)+1);
                    q(:,Pnodes_stage(j)+1)=sys_box.F{nchild(l)}'*y(:,nchild(l))...
                        +sys_box.A'*q(:,nchild(l)+1);
                else
                    r(:,Pnodes_stage(j)+1)=Ptree_box.K{nchild(l)}*sigma(:,nchild(l))+...
                        Ptree_box.d{nchild(l)}*y(:,nchild(l))+Ptree_box.f{nchild(l)}*...
                        q(:,nchild(l)+1)+r(:,Pnodes_stage(j)+1);
                    q(:,Pnodes_stage(j)+1)=sys_box.F{nchild(l)}'*y(:,nchild(l))...
                        +sys_box.A'*q(:,nchild(l)+1)+q(:,Pnodes_stage(j)+1);
                end
            end
        end
    end 
end

Z.X(:,1)=x;

for kk=1:Nd-Ns+1
    if(kk==1)
        v(:,kk)=Ptree_box.K{kk,1}*opts.v+v(:,kk);
        Z.U(:,kk)=sys_box.L*v(:,kk)+opts.vhat(:,kk);
        Z.X(:,kk+1)=sys_box.A*Z.X(:,1)+sys_box.B*Z.U(:,kk)+opts.w(:,kk);
    else
        nchild=Tree.children{kk-1};
        for l=1:length(nchild)
            v(:,nchild(l))=Ptree_box.K{nchild(l),1}*v(:,Tree.ancestor(nchild(l),1))+v(:,nchild(l));
            Z.U(:,nchild(l))=sys_box.L*v(:,nchild(l))+opts.vhat(:,nchild(l));
            Z.X(:,nchild(l)+1)=sys_box.A*Z.X(:,kk)+sys_box.B*Z.U(:,nchild(l))+...
                opts.w(:,nchild(l));
        end
    end
end

details.q=q;
end



