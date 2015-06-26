function [ Z,details] = solve_apg_effinet( sys,Tree,Ptree,Y,x,opts)

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
nv=size(sys.L,2);


Z.X=zeros(sys.nx,Nd+1);
Z.U=zeros(sys.nu,Nd+1-Ns);
v=zeros(nv,Nd-Ns+1);

q=zeros(sys.nx,Nd+1);
r=zeros(nv,Nd+1);
sigma=zeros(nv,Nd+1);
%sum_q=zeros(sys.nx,1);
%sum_r=zeros(sys.nu,1);
y=Y.y;
yt=Y.yt;
nyt=size(yt,1);

nx=sys.nx;
alpha_bar=opts.alpha_bar;


for i=1:Ns
    q(:,Tree.leaves(i)+1)=sys.F{Nd-Ns+1+i}(1:nyt,:)'*yt(:,i);
    r(:,Tree.leaves(i)+1)=Ptree.Gbar{Nd-Ns+1+i}(1:nyt,:)'*yt(:,i);
end 

for k=sys.Np-1:-1:0
    if(k==0)
        sigma(:,1)=r(:,2)+alpha_bar(:,1)+opts.beta(:,1);
        r(:,1)=Ptree(:,1).K{1,1}'*sigma(:,1)+Ptree.d{1,1}*y(:,1)+...
            Ptree.f{1,1}*q(:,2);
        v(:,1)=Ptree.Phi{1,1}*y(:,1)+Ptree.Theta{1,1}*q(:,2)+...
            Ptree.omega{1,1}*sigma(:,1);
        q(:,1)=sys.A'*q(:,2);
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
            q(:,nodes_stage(j)+1)=sys.F{nodes_stage(j)+1}'*y(:,nodes_stage(j)+1)+sys.A'*sum_q;
        end
    end
    
end

%{
for kk=Nd-Ns:-1:1
    stage=Tree.stage(kk);
    nchild=Tree.children{kk};
    sum_q=q(:,nchild(1));
    sum_r=r(:,nchild(1));
    if(length(nchild)>1)
        for l=2:lenght(nchild)
            sum_q=sum_q+q(:,nchild(l));
            sum_r=sum_r+r(:,nchild(l));
        end
    end
    sigma(:,kk)=sum_r+Tree.prob(kk)*alpha_bar(:,stage+1)+opts.beta(:,kk);
    r(:,kk)=sigma(:,kk)+Ptree.Gbar'*y((kk-1)*ny+1:kk*ny,1)+Ptree.Bbar'*sum_q;
    v(:,kk)=Ptree.Phi{kk,1}*y((kk-1)*ny+1:kk*ny,1)+Ptree.Theta{kk,1}*sum_q+...
        Ptree.omega{kk,1}*sigma(:,kk);
    if(kk==1)
        q(:,kk)=sys.A'*sum_q;
    else
        q(:,kk)=sys.F'*y((kk-1)*ny+1:kk*ny,1)+sys.A'*sum_q;
    end
end
%}
Z.X(:,1)=x;
for kk=1:Nd-Ns+1
    if(kk==1)
        v(:,kk)=Ptree.K{kk,1}*opts.v+v(:,kk);
        Z.U(:,kk)=sys.L*v(:,kk)+opts.vhat(:,kk);
        Z.X(:,kk+1)=sys.A*Z.X(:,1)+sys.B*Z.U(:,kk)+opts.w(:,kk);
    else
        v(:,kk)=Ptree.K{kk,1}*v(:,Tree.ancestor(kk-1,1)+1)+v(:,kk);
        nchild=Tree.children{kk-1};
        Z.U(:,kk)=sys.L*v(:,kk)+opts.vhat(:,kk);
        for l=1:length(nchild)
            Z.X(:,nchild(l)+1)=sys.A*Z.X(:,kk)+sys.B*Z.U(:,kk)+...
                opts.w(:,kk);
        end
    end
end

details.v=v;
end

