function [ sys_new ] = precondition_calculate( sys,DualHessian,Tree)

% Precondition_calculate calcuate the diagonal preconditioning
% of the system.

%DualHessian=DH_normalized;

nx=size(sys.A,1);
nu=size(sys.B,2);

Nd=size(Tree.stage,1);
Ns=size(Tree.leaves,1);

sys_new=sys;

DH_diag=diag(diag(DualHessian).^(-0.5));
%DH_diag=diag(DualHessian).^(-0.5);
%sys_new.xmin=kron(ones(Nd+1,1),sys.xmin);
%sys_new.xmax=kron(ones(Nd+1,1),sys.xmax);

%sys_new.umin=kron(ones(Nd+1-Ns,1),sys.umin);
%sys_new.umax=kron(ones(Nd+1-Ns,1),sys.umax);

if(sys.cell)
    ny=size(sys.F{1},1);
    for j=1:Nd+1
        if(j==1)
            sys_new.G{j,1}=sys.G{j,1};
            sys_new.G{j,1}(2*sys.nx+1:end,:)=DH_diag(1:2*nu,1:2*nu)*sys.G{j,1}(2*sys.nx+1:end,:);
            sys_new.umax(1:nu,1)=DH_diag(1:nu,1:nu)*sys.umax(1:nu,1);
            sys_new.umin(1:nu,1)=DH_diag(nu+1:2*nu,nu+1:2*nu)*sys.umin(1:nu,1);
        elseif(j==Nd+1)
            sys_new.F{j,1}=sys.F{j,1};
            sys_new.F{j,1}(1:2*sys.nx,:)=DH_diag(ny*(Nd-Ns)+2*nu+1:end,ny*(Nd-Ns)+2*nu+1:end)...
                *sys.F{j,1}(1:2*sys.nx,:);
            sys_new.xmax((Nd+Ns-1)*nx+1:end,1)=DH_diag((Nd-Ns)*ny+2*nu+1:...
                (Nd-Ns)*ny+2*nu+nx,(Nd-Ns)*ny+2*nu+1:...
                (Nd-Ns)*ny+2*nu+nx)*sys.xmax((Nd+Ns-1)*nx+1:end,1);
            sys_new.xmin((Nd+Ns-1)*nx+1:end,1)=DH_diag((Nd-Ns)*ny+2*nu+nx+1:...
                (Nd-Ns)*ny+2*(nu+nx),(Nd-Ns)*ny+2*nu+nx+1:...
                (Nd-Ns)*ny+2*(nu+nx))*sys.xmin((Nd+Ns-1)*nx+1:end,1);
        else
            %sys_new.F{j,1}=sys.F{j,1};
            sys_new.F{j,1}(1:2*sys.nx,:)=DH_diag((j-2)*ny+2*nu+1:(j-1)*ny,(j-2)*ny+2*nu+1:(j-1)*ny)...
                *sys.F{j,1}(1:2*sys.nx,:);
            sys_new.G{j,1}(2*sys.nx+1:end,:)=DH_diag((j-1)*ny+1:(j-1)*ny+2*nu,(j-1)*ny+1:(j-1)*ny+2*nu)...
                *sys.G{j,1}(2*sys.nx+1:end,:);
            
            sys_new.xmax((j-2)*nx+1:(j-1)*nx,1)=DH_diag((j-2)*ny+2*nu+1:(j-2)*ny+2*nu+nx,...
                (j-2)*ny+2*nu+1:(j-2)*ny+2*nu+nx)*sys.xmax((j-2)*nx+1:(j-1)*nx,1);
            sys_new.xmin((j-2)*nx+1:(j-1)*nx,1)=DH_diag((j-2)*ny+2*nu+nx+1:(j-1)*ny,...
                (j-2)*ny+2*nu+nx+1:(j-1)*ny)*sys.xmin((j-2)*nx+1:(j-1)*nx,1);
            
            sys_new.umax((j-1)*nu+1:j*nu,1)=DH_diag((j-1)*ny+1:(j-1)*ny+nu,...
                (j-1)*ny+1:(j-1)*ny+nu)*sys.umax((j-1)*nu+1:j*nu,1);
            sys_new.umin((j-1)*nu+1:j*nu,1)=DH_diag((j-1)*ny+nu+1:(j-1)*ny+2*nu,...
                (j-1)*ny+nu+1:(j-1)*ny+2*nu)*sys.umin((j-1)*nu+1:j*nu,1);
            
        end
    end
else
    ny=size(sys.F,1);
    sys_new.F=cell(Nd+1,1);
    sys_new.G=cell(Nd+1-Ns,1);
    for j=1:Nd+1
        if(j==1)
            sys_new.G{j,1}=DH_diag(1:2*nu,1:2*nu)*sys.G(2*sys.nx+1:end,:);
            sys_new.umax(1:nu,1)=DH_diag(1:nu,1:nu)*sys.umax(1:nu,1);
            sys_new.umin(1:nu,1)=DH_diag(nu+1:2*nu,nu+1:2*nu)*sys.umin(1:nu,1);
        elseif(j==Nd+1)
            sys_new.F{j,1}=DH_diag(ny*(Nd-Ns)+2*nu+1:end,ny*(Nd-Ns)+2*nu+1:end)...
                *sys.F(1:2*sys.nx,:);
            sys_new.xmax((Nd+Ns-1)*nx+1:end,1)=DH_diag((Nd-Ns)*ny+2*nu+1:...
                (Nd-Ns)*ny+2*nu+nx,(Nd-Ns)*ny+2*nu+1:...
                (Nd-Ns)*ny+2*nu+nx)*sys.xmax((Nd+Ns-1)*nx+1:end,1);
            sys_new.xmin((Nd+Ns-1)*nx+1:end,1)=DH_diag((Nd-Ns)*ny+2*nu+nx+1:...
                (Nd-Ns)*ny+2*(nu+nx),(Nd-Ns)*ny+2*nu+nx+1:...
                (Nd-Ns)*ny+2*(nu+nx))*sys.xmin((Nd+Ns-1)*nx+1:end,1);
        else
            sys_new.F{j,1}=DH_diag((j-1)*ny+2*nu+1:j*ny,(j-1)*ny+2*nu+1:j*ny)...
                *sys.F(1:2*sys.nx,:);
            sys_new.G{j,1}=DH_diag(j*ny+1:j*ny+2*nu,j*ny+1:j*ny+2*nu)...
                *sys.G(2*sys.nx+1:end,:);
            
            sys_new.xmax((j-1)*nx+1:j*nx,1)=DH_diag((j-1)*ny+2*nu+1:(j-1)*ny+2*nu+nx,...
                (j-1)*ny+2*nu+1:(j-1)*ny+2*nu+nx)*sys.xmax((j-1)*nx+1:j*nx,1);
            sys_new.xmin((j-1)*nx+1:end,1)=DH_diag((j-1)*ny+2*nu+nx+1:j*ny,...
                (j-1)*ny+2*nu+nx+1:j*ny)*sys.xmin((j-1)*nx+1:j*nx,1);
            
            sys_new.umax(j*nu+1:(j+1)*nu,1)=DH_diag(j*ny+1:j*ny+nu,...
                j*ny+1:j*ny+nu)*sys.umax(j*nu+1:(j+1)*nu,1);
            sys_new.umin(j*nu+1:(j+1)*nu,1)=DH_diag(j*ny+nu+1:j*ny+2*nu,...
                j*ny+nu+1:j*ny+2*nu)*sys.umin(j*nu+1:(j+1)*nu,1);
            
        end
    end
end

end

