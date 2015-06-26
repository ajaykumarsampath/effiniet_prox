function [ Dual_hessian,details] = dual_hessian_calculate_box( sys,Tree,V,apg_opts)
% This function calculate the dual hessian for the effiniet problem
%
% INPUT ----    sys      :
%               Tree     :
%               V        :
%               apg_opts :
%
% OUTPUT -----  Dual_hessian :
%               details      :

%%
nx=sys.nx;
nu=sys.nu;
ne=size(apg_opts.E,1);
Nd=size(Tree.stage,1);
Ns=size(Tree.leaves,1);
nz=nx+nu;
nz1=2*nx+nu;
%Np=sys.Np;
Bbar=zeros(Nd*(nx+ne),nu+Ns*nx+nz*(Nd-Ns));
%Bbar=zeros(Nd*nx,nu+Ns*nx+nz*(Nd-Ns));
H=zeros(nz*Nd);
ny=size(sys.F{1},1);
A=zeros(nz+(Nd-Ns)*ny,nz*Nd);
epsilon=0e-4*eye(nz);

%{
for k=1:Np
    if(k==1)
        Bbar(1:nx,1:nz)=[sys.B -eye(nx)];
        Bbar(Nd*nx+(k-1)*ne+1:Nd*nx+k*ne,1:nu)=apg_opts.E;
    else
        nodes_stage=find(Tree.stage==k-2);
        
        for j=1:length(nodes_stage)
            nchild=Tree.children{nodes_stage(j)};
            for i=1:length(nchild)
                Bbar((nchild(i)-1)*nx+1:nchild(i)*nx,...
                    (nodes_stage(j)-1)*nz+nu+1:nodes_stage(j)*nz+nu)=[sys.A sys.B];
                Bbar((nchild(i)-1)*nx+1:nchild(i)*nx,...
                   (nchild(i)-1)*nz+nu+1:(nchild(i)-1)*nz+nu+nx)=-eye(nx);
                Bbar(Nd*nx+(nchild(i)-1)*ne+1:Nd*nx+nchild(i)*ne,...
                    nodes_stage(j)*nz+1:nodes_stage(j)*nz+nu)=apg_opts.E;
            end
        end
    end
end
%}

for i=1:Nd
    if(i==1)
        Bbar(1:nx,1:nz)=[sys.B -eye(nx)];
        Bbar(Nd*nx+(i-1)*ne+1:Nd*nx+i*ne,1:nu)=apg_opts.E;
    else
        Bbar((i-1)*nx+1:i*nx,nu+(i-2)*nz+1:nu+(i-2)*nz+nz1)=...
            [sys.A sys.B -eye(nx)];
        Bbar(Nd*nx+(i-1)*ne+1:Nd*nx+i*ne,(i-1)*nz+1:(i-1)*nz+nu)=apg_opts.E;
    end 
end
for i=1:Nd
    if(i==1)
        H(1:nu,1:nu)=2*V.Wu+epsilon(1:nu,1:nu);
        H(1:nu,nz+1:nz+nu)=-V.Wu;
        H(nu+1:nz,nu+1:nz)=epsilon(nu+1:nz,nu+1:nz);
    elseif(i==Nd)
        H((Nd-1)*nz+1:(Nd-1)*nz+nu,(Nd-2)*nz+1:(Nd-2)*nz+nu)=-Tree.prob(Nd-1,1)*V.Wu;
        H((Nd-1)*nz+1:(Nd-1)*nz+nu,(Nd-1)*nz+1:(Nd-1)*nz+nu)=Tree.prob(Nd-1,1)*...
            V.Wu+epsilon(1:nu,1:nu);
        H((Nd-1)*nz+nu+1:Nd*nz,(Nd-1)*nz+nu+1:Nd*nz)=epsilon(nu+1:nz,nu+1:nz);
    else
        H((i-1)*nz+1:(i-1)*nz+nu,(i-2)*nz+1:(i-2)*nz+nu)=-Tree.prob(i-1,1)*V.Wu;
        H((i-1)*nz+1:(i-1)*nz+nu,(i-1)*nz+1:(i-1)*nz+nu)=...
            (Tree.prob(i-1,1)+Tree.prob(i))*V.Wu+epsilon(1:nu,1:nu);
        H((i-1)*nz+1:(i-1)*nz+nu,i*nz+1:i*nz+nu)=-Tree.prob(i,1)*V.Wu;
        H((i-1)*nz+nu+1:i*nz,(i-1)*nz+nu+1:i*nz)=epsilon(nu+1:nz,nu+1:nz);
    end
end 
%}
Kbar=[H Bbar';Bbar zeros(size(Bbar,1))];

L=null(Bbar);
Wbar=L'*H*L;


if(sys.cell)
    for i=1:Nd+1
        if(i==1)
            A(1:nu,1:nu)=sys.G{i}(nx+1:nz,:);
        elseif(i==Nd+1)
            A(nu+(i-2)*ny+1:nz+(i-2)*ny,nu+(i-2)*nz+1:(i-1)*nz)=...
                sys.F{i}(1:nx,:);
        else
            A(nu+(i-2)*ny+1:nu+(i-1)*ny,nu+(i-2)*nz+1:nu+(i-1)*nz)=...
                [sys.F{i} sys.G{i}];
        end
    end
else
    for i=1:Nd+1
        if(i==1)
            A(1:nu,1:nu)=sys.G(nx+1:nz,:);
        elseif(i==Nd+1)
            A(nu+(i-2)*ny+1:nz+(i-2)*ny,nu+(i-2)*nz+1:(i-1)*nz)=...
                sys.F(1:nx,:);
        else
            A(nu+(i-2)*ny+1:nu+(i-1)*ny,nu+(i-2)*nz+1:nu+(i-1)*nz)=...
                [sys.F sys.G];
        end
    end
end

Abar=A*L;
K1=Wbar\(eye(size(Wbar,1)));
details.cond_primal=cond(K1);
details.min_eig=min(eig(K1));


Dual_hessian=Abar*K1*Abar';
details.condition_number=cond(Dual_hessian);
details.norm=norm(Dual_hessian);

end



