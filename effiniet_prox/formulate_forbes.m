function [ prob ] = formulate_forbes( sys,V,Tree,forbes_opts )
% This function formulates  problem for forbes 

%%
nx=sys.nx;
nu=sys.nu;
ne=size(forbes_opts.E,1);
Nd=size(Tree.stage,1);
Ns=size(Tree.leaves,1);
nz=nx+nu;
nz1=2*nx+nu;
ny=size(sys.F{1},1);

Bbar=zeros(Nd*(nx+ne),nu+Ns*nx+nz*(Nd-Ns));
bbar=zeros(Nd*(nx+ne),1);

H=zeros(nz*Nd);
%A=zeros(2*nz+(Nd-Ns)*ny,nz*Nd);
q=zeros(nz*Nd,1);

lb=zeros(Nd*nz,1);
ub=zeros(Nd*nz,1);

epsilon=0e-4*eye(nz);

% function f
% <Hz,z>+<q,z> such that Bbar z =bbar;

% Bbar z =bbar ;
for i=1:Nd
    if(i==1)
        Bbar(1:nx,1:nz)=[sys.B -eye(nx)];
        Bbar(Nd*nx+(i-1)*ne+1:Nd*nx+i*ne,1:nu)=forbes_opts.E;
        bbar(1:nx,1)=-sys.A*forbes_opts.x-forbes_opts.w(:,1);
        bbar(Nd*nx+(i-1)*ne+1:Nd*nx+i*ne,1)=-forbes_opts.Ed*forbes_opts.demand(1,:)';
    else
        Bbar((i-1)*nx+1:i*nx,nu+(i-2)*nz+1:nu+(i-2)*nz+nz1)=...
            [sys.A sys.B -eye(nx)];
        Bbar(Nd*nx+(i-1)*ne+1:Nd*nx+i*ne,(i-1)*nz+1:(i-1)*nz+nu)=forbes_opts.E;
        bbar((i-1)*nx+1:i*nx,1)=-forbes_opts.w(:,i);
        bbar(Nd*nx+(i-1)*ne+1:Nd*nx+i*ne,1)=-forbes_opts.Ed*forbes_opts.demand(i,:)';
    end 
end

% <Hz,z>+<q,z> ;
for i=1:Nd
    if(i==1)
        H(1:nu,1:nu)=2*V.Wu+epsilon(1:nu,1:nu);
        H(1:nu,nz+1:nz+nu)=-V.Wu;
        H(nu+1:nz,nu+1:nz)=epsilon(nu+1:nz,nu+1:nz);
        q(1:nu,1)=V.alpha(1,:)'-2*V.Wu'*forbes_opts.uprev;
    elseif(i==Nd)
        H((Nd-1)*nz+1:(Nd-1)*nz+nu,(Nd-2)*nz+1:(Nd-2)*nz+nu)=-V.Wu;
        H((Nd-1)*nz+1:(Nd-1)*nz+nu,(Nd-1)*nz+1:(Nd-1)*nz+nu)=V.Wu+epsilon(1:nu,1:nu);
        H((Nd-1)*nz+nu+1:Nd*nz,(Nd-1)*nz+nu+1:Nd*nz)=epsilon(nu+1:nz,nu+1:nz);
        q((Nd-1)*nz+1:(Nd-1)*nz+nu,1)=V.alpha(Nd,:)';
    else
        H((i-1)*nz+1:(i-1)*nz+nu,(i-2)*nz+1:(i-2)*nz+nu)=-V.Wu;
        H((i-1)*nz+1:(i-1)*nz+nu,(i-1)*nz+1:(i-1)*nz+nu)=2*V.Wu+epsilon(1:nu,1:nu);
        H((i-1)*nz+1:(i-1)*nz+nu,i*nz+1:i*nz+nu)=-V.Wu;
        H((i-1)*nz+nu+1:i*nz,(i-1)*nz+nu+1:i*nz)=epsilon(nu+1:nz,nu+1:nz);
        q((i-1)*nz+1:(i-1)*nz+nu,1)=V.alpha(i,:)';
    end
end 

H=2*H;

for i=1:Nd+1
    if(i==1)
        lb(1:nu,1)=sys.umin(1:nu,1);
        ub(1:nu,1)=sys.umax(1:nu,1);
    elseif(i==Nd+1)
        lb((Nd-1)*nz+nu+1:Nd*nz,1)=sys.xmin(1:nx,1);
        ub((Nd-1)*nz+nu+1:Nd*nz,1)=sys.xmax(1:nx,1);
    else
        lb((i-2)*nz+nu+1:(i-1)*nz+nu,1)=[sys.xmin(1:nx,1);sys.umin(1:nu,1)];
        ub((i-2)*nz+nu+1:(i-1)*nz+nu,1)=[sys.xmax(1:nx,1);sys.umax(1:nu,1)];
    end
end
prob.f1=quad_over_affine(H,q,Bbar,bbar);

prob.g=indBox(lb,ub);

prob.A1 = 1;
prob.B = -1;
prob.b = zeros(nz*(Nd),1);





%{
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
%}

end

