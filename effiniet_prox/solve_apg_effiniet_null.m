function [ Z] = solve_apg_effiniet_null( sys_null,Tree,Ptree,Y,x,opts)

% This function performs the solve step for the DWNs problem 
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

% opts=current_state_opt;
Nd=size(Tree.stage,1);% total number of nodes including leaves 
Ns=size(Tree.leaves,1);% total number of leaves
nv=size(sys_null.L,2);


Z.X=zeros(sys_null.nx,Nd+1);
Z.U=zeros(size(sys_null.L,1),Nd+1-Ns);
v=zeros(nv,Nd-Ns+1);

q=zeros(sys_null.nx,Nd+1);
r=zeros(nv,Nd+1);
sigma=zeros(nv,Nd+1);
%sum_q=zeros(sys.nx,1);
%sum_r=zeros(sys.nu,1);
y=Y.y;
yt=Y.yt;
nx=sys_null.nx;
alpha_bar=opts.alpha_bar;


for i=1:Ns
    q(:,Tree.leaves(i)+1)=sys_null.F(1:2*nx,:)'*yt(:,i);
    r(:,Tree.leaves(i)+1)=Ptree.Gbar(1:2*nx,:)'*yt(:,i);
end 

for k=sys_null.Np-1:-1:0
    if(k==0)
        sigma(:,1)=r(:,2)+alpha_bar(:,1)+opts.beta(:,1);
        r(:,1)=Ptree(:,1).K{1,1}'*sigma(:,1)+Ptree.d{1,1}*y(:,1)+...
            Ptree.f{1,1}*q(:,2);
        v(:,1)=Ptree.Phi{1,1}*y(:,1)+Ptree.Theta{1,1}*q(:,2)+...
            Ptree.omega{1,1}*sigma(:,1);
        q(:,1)=sys_null.A'*q(:,2);
    else
        nodes_stage=find(Tree.stage==k-1);
        for j=1:length(nodes_stage)
            nchild=Tree.children{nodes_stage(j)};
            sum_q=q(:,nchild(1)+1);
            sum_r=r(:,nchild(1)+1);
            if(length(nchild)>1)
                for l=2:length(nchild)
                    sum_q=sum_q+q(:,nchild(l)+1);
                    sum_r=sum_r+r(:,nchild(l)+1);
                end
            end
            sigma(:,nodes_stage(j)+1)=sum_r+Tree.prob(nodes_stage(j))*alpha_bar(:,k+1)...
                +opts.beta(:,nodes_stage(j)+1);
            r(:,nodes_stage(j)+1)=Ptree.K{nodes_stage(j)+1,1}'*sigma(:,nodes_stage(j)+1)+...
                Ptree.d{nodes_stage(j)+1,:}*y(:,nodes_stage(j)+1)...
                +Ptree.f{nodes_stage(j)+1,:}*sum_q;
            v(:,nodes_stage(j)+1)=Ptree.Phi{nodes_stage(j)+1,1}*...
                y(:,nodes_stage(j)+1)+Ptree.Theta{nodes_stage(j)+1,1}*sum_q+...
                Ptree.omega{nodes_stage(j)+1,1}*sigma(:,nodes_stage(j)+1);
            q(:,nodes_stage(j)+1)=sys_null.F'*y(:,nodes_stage(j)+1)+sys_null.A'*sum_q;
        end
    end
    
end

Z.X(:,1)=x;
for kk=1:Nd-Ns+1
    if(kk==1)
        v(:,kk)=Ptree.K{kk,1}*opts.v+v(:,kk);
        Z.U(:,kk)=sys_null.L*v(:,kk)+opts.vhat(:,kk);
        Z.X(:,kk+1)=sys_null.A*Z.X(:,1)+Ptree.Bbar*v(:,kk)+opts.w(:,kk);
    else
        v(:,kk)=Ptree.K{kk,1}*v(:,Tree.ancestor(kk-1,1)+1)+v(:,kk);
        nchild=Tree.children{kk-1};
        Z.U(:,kk)=sys_null.L*v(:,kk)+opts.vhat(:,kk);
        for l=1:length(nchild)
            Z.X(:,nchild(l)+1)=sys_null.A*Z.X(:,kk)+Ptree.Bbar*v(:,kk)+...
                opts.w(:,kk);
        end
    end
end

Z.v=v;
end

