function [ sys,V] = system_prox_formation_modified( S,P,Tree,ops_sys)
% This function transform the effiniet problem into
% dual proximal gradient formulation.
%
%   INPUT-----    S  : Effinet system discription
%                 P  : Effinet cost description
%              Tree  : Scenario tree discription
%
%   OUTPUT------ sys : Output system for dual proximal form

% remove the input corresponding to u4;
sys.nu=S.nu-1;
sys.nx=S.nx;
sys.Np=P.Hp;
Nd=size(Tree.stage,1);
Ns=size(Tree.leaves,1);
%sys.F=[S.A;zeros(S.nu,S.nx)];
%sys.G=[S.B;eye(S.nu)];
nx=sys.nx;
nu=sys.nu;
sys.xmin=P.xs;
sys.xmax=S.xmax;

sys.xs=P.beta*sys.xmax;

sys.gamma_min=ops_sys.gamma_min;
sys.gamma_xs=ops_sys.gamma_xs;
sys.gamma_max=ops_sys.gamma_max;

% Null spaces and particular soultion format
sys.L=null(S.E(:,[1:3 5:6]));
sys.L1=-pinv(S.E(:,[1:3 5:6]))*(S.Ed);

sys.A=S.A;
sys.B=S.B(:,[1:3 5:6])/3600;
sys.Gd=S.Gd/3600;
sys.umin=3600*S.umin([1:3 5:6],1);
sys.umax=3600*S.umax([1:3 5:6],1);

% cost function
V.Wu=P.Wu([1:3 5:6],[1:3 5:6]);

% Normalise the constraints
if(ops_sys.cell)
    sys.F=cell(Nd+1,1);
    sys.G=cell(Nd+1-Ns,1);
    for j=1:Nd+1
        sys.F{j}=[eye(nx);-eye(nx);zeros(2*nu,nx)];
        sys.G{j}=[zeros(2*nx,nu);eye(nu);-eye(nu)];
        if(ops_sys.normalise)
            for i=1:S.nx
                if(sys.xmax(i)>0)
                    sys.F{j}(i,i)=sys.F{j}(i,i)/sys.xmax(i);
                    if(j==Nd+1)
                        sys.xmax(i)=1;
                    end
                end
                
                if(sys.xmin(i)>0)
                    sys.F{j}(S.nx+i,i)=sys.F{j}(S.nx+i,i)/sys.xmin(i);
                    if(j==Nd+1)
                        sys.xmin(i)=1;
                    end
                end
            end
            for i=1:sys.nu
                if(sys.umax(i)>0)
                    sys.G{j}(2*S.nx+i,i)=sys.G{j}(2*S.nx+i,i)/sys.umax(i);
                    if(j==Nd+1)
                        sys.umax(i)=1;
                    end
                end
                
                if(sys.umin(i)>0)
                    sys.G{j}(2*S.nx+S.nu+i,i)=sys.G{j}(2*S.nx+S.nu+i,i)/sys.umin(i);
                    if(j==Nd+1)
                        sys.umin(i)=1;
                    end
                end
            end
        end
    end
else
    sys.F=[eye(S.nx);-eye(S.nx);zeros(2*S.nu,S.nx)];
    sys.G=[zeros(2*S.nx,S.nu);eye(S.nu);-eye(S.nu)];
    if(ops_sys.normalise)
        for i=1:S.nx
            if(sys.xmax(i)>0)
                sys.F(i,i)=sys.F(i,i)/sys.xmax(i);
                sys.xmax(i)=1;
            end
            
            if(sys.xmin(i)>0)
                sys.F(S.nx+i,i)=sys.F(S.nx+i,i)/sys.xmin(i);
                sys.xmin(i)=1;
            end
            
        end
        
        for i=1:S.nu
            if(sys.umax(i)>0)
                sys.G(2*S.nx+i,i)=sys.G(2*S.nx+i,i)/sys.umax(i);
                sys.umax(i)=1;
            end
            
            if(sys.umin(i)>0)
                sys.G(2*S.nx+S.nu+i,i)=sys.G(2*S.nx+S.nu+i,i)/sys.umin(i);
                sys.umin(i)=1;
            end
            
        end
    end
end
sys.cell=ops_sys.cell;
sys.normalized=ops_sys.normalise;

sys.umin=kron(ones(Nd-Ns+1,1),sys.umin);
sys.umax=kron(ones(Nd-Ns+1,1),sys.umax);

sys.xmin=kron(ones(Nd+1,1),sys.xmin);
sys.xmax=kron(ones(Nd+1,1),sys.xmax);
sys.xs=kron(ones(Nd+1,1),sys.xs);
%}

%}
end


