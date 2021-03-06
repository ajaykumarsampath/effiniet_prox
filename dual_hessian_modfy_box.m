function [ Dual_hessian,details] = dual_hessian_modfy_box( sys,Tree,V,apg_opts)
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
Np=sys.Np;
nz=nx+nu;
nz1=2*nx+nu;
%Np=sys.Np;
Bbar=zeros(Np*(nx+ne),Np*nz);
H=zeros(nz*Np);
ny=size(sys.F{1},1);
A=zeros(Np*nz,Np*nz);
epsilon=0e-4*eye(nz);


for i=1:Np
    if(i==1)
        Bbar(1:nx,1:nz)=[sys.B -eye(nx)];
        Bbar(Np*nx+(i-1)*ne+1:Np*nx+i*ne,1:nu)=apg_opts.E;
    else
        Bbar((i-1)*nx+1:i*nx,nu+(i-2)*nz+1:nu+(i-2)*nz+nz1)=...
            [sys.A sys.B -eye(nx)];
        Bbar(Np*nx+(i-1)*ne+1:Np*nx+i*ne,(i-1)*nz+1:(i-1)*nz+nu)=apg_opts.E;
    end 
end

for i=1:Np
    if(i==Np)
        H((i-1)*nz+1:(i-1)*nz+nu,(i-1)*nz+1:(i-1)*nz+nu)=Tree.prob(i)*V.Wu;
        H((i-1)*nz+1:(i-1)*nz+nu,(i-2)*nz+1:(i-2)*nz+nu)=-Tree.prob(i)*V.Wu;
        H((i-1)*nz+nu+1:i*nz,(i-1)*nz+nu+1:i*nz)=epsilon(nu+1:nz,nu+1:nz);
    else
        if(i>1)
            H((i-1)*nz+1:(i-1)*nz+nu,(i-2)*nz+1:(i-2)*nz+nu)=-Tree.prob(i)*V.Wu;
        end
        H((i-1)*nz+1:(i-1)*nz+nu,(i-1)*nz+1:(i-1)*nz+nu)=(Tree.prob(i)+Tree.prob(i+1))*V.Wu;
        H((i-1)*nz+1:(i-1)*nz+nu,i*nz+1:i*nz+nu)=-Tree.prob(i)*V.Wu;
        H((i-1)*nz+nu+1:i*nz,(i-1)*nz+nu+1:i*nz)=epsilon(nu+1:nz,nu+1:nz);
    end
end

%}
%Kbar=[H Bbar';Bbar zeros(size(Bbar,1))];

L=null(Bbar);
Wbar=L'*H*L;



if(sys.cell)
    for i=1:Np
        A((i-1)*nz+1:i*nz,(i-1)*nz+1:i*nz)=blkdiag(sys.G{i},...
            sys.F{i});
    end
else
    for i=1:Np
        A((i-1)*nz+1:i*nz,(i-1)*nz+1:i*nz)=blkdiag(sys.G,...
            sys.F);
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




