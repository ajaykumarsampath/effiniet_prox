function [ Dual_hessian,details] = dual_hessian_calculate_u( sys,Tree,V,apg_opts)
% This function calculate the dual hessian for the effiniet problem
%
% INPUT---- sys
%%
nx=sys.nx;
nu=sys.nu;
ne=size(apg_opts.E,1);
Nd=size(Tree.stage,1);
Ns=size(Tree.leaves,1);
nz=nx+nu;
%nz1=2*nx+nu;
Bbar=zeros(Nd*(nx+ne),Nd*nu);
%Bbar=zeros(Nd*nx,Nd*nu);
H=zeros(nu*Nd);
ny=size(sys.F{1},1);
A=zeros(2*nz+(Nd-Ns)*ny,nz*Nd);
epsilon=0e-4*eye(nz);

for i=1:Nd
    Bbar((i-1)*nx+1:Nd*nx,(i-1)*nu+1:i*nu)=kron(ones(Nd-i+1,1),sys.B);
end 

Bbar(Nd*nx+1:end,:)=kron(eye(Nd),apg_opts.E);

for i=1:Nd
    if(i==1)
        H(1:nu,1:nu)=2*V.Wu+epsilon(1:nu,1:nu);
        H(1:nu,nu+1:nu+nu)=-V.Wu;
    elseif(i==Nd)
        H((Nd-1)*nu+1:(Nd-1)*nu+nu,(Nd-2)*nu+1:(Nd-2)*nu+nu)=-V.Wu;
        H((Nd-1)*nu+1:(Nd-1)*nu+nu,(Nd-1)*nu+1:(Nd-1)*nu+nu)=V.Wu+epsilon(1:nu,1:nu);
    else
        H((i-1)*nu+1:(i-1)*nu+nu,(i-2)*nu+1:(i-2)*nu+nu)=-V.Wu;
        H((i-1)*nu+1:(i-1)*nu+nu,(i-1)*nu+1:(i-1)*nu+nu)=2*V.Wu+epsilon(1:nu,1:nu);
        H((i-1)*nu+1:(i-1)*nu+nu,i*nu+1:i*nu+nu)=-V.Wu;
    end
end 

%PP=H+Bbar'*Bbar;
L=null(Bbar);
Wbar=L'*H*L;
K1=Wbar\(eye(size(Wbar,1)));
details.cond_primal=cond(K1);
details.min_eig=min(eig(K1));

%{
Kbar=[H Bbar';Bbar zeros(size(Bbar,1))];
Hbar=Kbar\(eye(size(Kbar,1)));
K1=Hbar(1:size(H,1),1:size(H,1));
cond(K1)
min(eig(K1))
L=null(Bbar);
Wu=zeros(nu*Nd,Nd*nu);
for i=1:Nd
    if(i==1)
        Wu(1:nu,1:nu)=2*V.Wu;
        Wu(1:nu,nu+1:2*nu)=-V.Wu;
    elseif(i==Nd)
        Wu((Nd-1)*nu+1:Nd*nu,(Nd-1)*nu+1:Nd*nu)=V.Wu;
        Wu((Nd-1)*nu+1:Nd*nu,(Nd-2)*nu+1:(Nd-1)*nu)=-V.Wu;
    else
        Wu((i-1)*nu+1:i*nu,(i-2)*nu+1:(i-1)*nu)=-V.Wu;
        Wu((i-1)*nu+1:i*nu,(i-1)*nu+1:i*nu)=2*V.Wu;
        Wu((i-1)*nu+1:i*nu,i*nu+1:(i+1)*nu)=-V.Wu;
    end 
end 
dual_hessian=L'*(Wu\L);
%Kbar=[Wu Bbar';Bbar zeros(size(Bbar,1))];
%cond(Kbar)
%\eye(size(Kbar,1))
%hessian=inv(Kbar);
%cond(hessian)
%}

end

