function [ Dual_hessian,details] = dual_hessian_calculate( sys,Tree,V,apg_opts)
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
Bbar=zeros(Nd*(nx+ne),nu+Ns*nx+nz*(Nd-Ns));
%Bbar=zeros(Nd*nx,nu+Ns*nx+nz*(Nd-Ns));
H=zeros(nz*Nd);
ny=size(sys.F{1},1);
A=zeros(2*nz+(Nd-Ns)*ny,nz*Nd);
epsilon=0e-4*eye(nz);
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
        H((Nd-1)*nz+1:(Nd-1)*nz+nu,(Nd-2)*nz+1:(Nd-2)*nz+nu)=-V.Wu;
        H((Nd-1)*nz+1:(Nd-1)*nz+nu,(Nd-1)*nz+1:(Nd-1)*nz+nu)=V.Wu+epsilon(1:nu,1:nu);
        H((Nd-1)*nz+nu+1:Nd*nz,(Nd-1)*nz+nu+1:Nd*nz)=epsilon(nu+1:nz,nu+1:nz);
    else
        H((i-1)*nz+1:(i-1)*nz+nu,(i-2)*nz+1:(i-2)*nz+nu)=-V.Wu;
        H((i-1)*nz+1:(i-1)*nz+nu,(i-1)*nz+1:(i-1)*nz+nu)=2*V.Wu+epsilon(1:nu,1:nu);
        H((i-1)*nz+1:(i-1)*nz+nu,i*nz+1:i*nz+nu)=-V.Wu;
        H((i-1)*nz+nu+1:i*nz,(i-1)*nz+nu+1:i*nz)=epsilon(nu+1:nz,nu+1:nz);
    end
end 

Kbar=[H Bbar';Bbar zeros(size(Bbar,1))];

L=null(Bbar);
Wbar=L'*H*L;


if(sys.cell)
    for i=1:Nd+1
        if(i==1)
            A(1:2*nu,1:nu)=sys.G{i}(2*nx+1:2*nz,:);
        elseif(i==Nd+1)
            A(2*nu+(i-2)*ny+1:2*nz+(i-2)*ny,nu+(i-2)*nz+1:(i-1)*nz)=...
                sys.F{i}(1:2*nx,:);
        else
            A(2*nu+(i-2)*ny+1:2*nu+(i-1)*ny,nu+(i-2)*nz+1:nu+(i-1)*nz)=...
                [sys.F{i} sys.G{i}];
        end
    end
else
    for i=1:Nd+1
        if(i==1)
            A(1:2*nu,1:nu)=sys.G(2*nx+1:2*nz,:);
        elseif(i==Nd+1)
            A(2*nu+(i-2)*ny+1:2*nz+(i-2)*ny,nu+(i-2)*nz+1:(i-1)*nz)=...
                sys.F(1:2*nx,:);
        else
            A(2*nu+(i-2)*ny+1:2*nu+(i-1)*ny,nu+(i-2)*nz+1:nu+(i-1)*nz)=...
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
%{
Abar_2=eye(nz*Nd);
Dual_hessian_1=Abar_2*L*K1*(Abar_2*L)';
details.condition_number_1=cond(Dual_hessian_1);
details.norm_1=norm(Dual_hessian_1);

Hbar=Kbar\(eye(size(Kbar,1)));
K11=Hbar(1:size(H,1),1:size(H,1));

details.cond_primal_2=cond(K11);
details.min_eig_2=min(eig(K11));
Dual_hessian_2=A*K11*A';
details.condition_number_2=cond(Dual_hessian_2);
details.norm_2=norm(Dual_hessian_2);
%}


end


