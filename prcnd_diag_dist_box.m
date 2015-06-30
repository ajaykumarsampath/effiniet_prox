function [ sys_new ] = prcnd_diag_dist_box( sys,DualHessian,opts)

% Precondition_calculate calcuate the diagonal preconditioning
% of the system. 
% INPUT -----      sys    : 
%          DualHessian    :
%                 Tree    :
%
% OUTPUT -----   sys_new  :
%%
% DualHessian=DH_normalized_not;

if(nargin<3)
    opts.dual_hessian=1;
end 

nx=size(sys.A,1);
nu=size(sys.B,2);
nz=2*nx+nu;
Np=sys.Np;

if(opts.dual_hessian)
    DH_diag=diag(diag(DualHessian).^(-0.5));
else
    DH_diag=DualHessian;
end

%DH_diag=diag(diag(DualHessian).^(-0.5));

sys_new=sys;

if(sys.cell)
    for j=1:Np
        sys_new.F{j,1}=sys.F{j,1};
        sys_new.G{j,1}=sys.G{j,1};
        sys_new.G{j,1}=DH_diag((j-1)*nz+1:(j-1)*nz+nu,(j-1)*nz+1:(j-1)*nz+nu)*...
            sys_new.G{j,1};
        sys_new.F{j,1}=DH_diag((j-1)*nz+nu+1:j*nz,(j-1)*nz+nu+1:j*nz)*...
            sys_new.F{j,1};
        
        sys_new.umax((j-1)*nu+1:j*nu,1)=DH_diag((j-1)*nz+1:(j-1)*nz+nu,...
            (j-1)*nz+1:(j-1)*nz+nu)*sys_new.umax((j-1)*nu+1:j*nu,1);
        sys_new.umin((j-1)*nu+1:j*nu,1)=DH_diag((j-1)*nz+1:(j-1)*nz+nu,...
            (j-1)*nz+1:(j-1)*nz+nu)*sys_new.umin((j-1)*nu+1:j*nu,1);
        
        sys_new.xmax((j-1)*nx+1:j*nx,1)=DH_diag((j-1)*nz+nu+1:(j-1)*nz+nu+nx,...
            (j-1)*nz+nu+1:(j-1)*nz+nu+nx)*sys_new.xmax((j-1)*nx+1:j*nx);
        sys_new.xmin((j-1)*nx+1:j*nx,1)=DH_diag((j-1)*nz+nu+1:(j-1)*nz+nu+nx,...
            (j-1)*nz+nu+1:(j-1)*nz+nu+nx)*sys_new.xmin((j-1)*nx+1:j*nx);
        sys_new.xs((j-1)*nx+1:j*nx,1)=DH_diag((j-1)*nz+nu+nx+1:j*nz,...
            (j-1)*nz+nu+nx+1:j*nz)*sys_new.xmin((j-1)*nx+1:j*nx);
    end
else
    sys_new.F=cell(Np,1);
    sys_new.G=cell(Np,1);
    for j=1:Np
        sys_new.F{j,1}=sys.F;
        sys_new.G{j,1}=sys.G;
        sys_new.G{j,1}(nx+1:nz,:)=DH_diag((j-1)*nz+1:(j-1)*nz+nu,(j-1)*nz+1:(j-1)*nz+nu)*...
            sys_new.G{j,1};
        sys_new.F{j,1}(1:nx,:)=DH_diag((j-1)*nz+nu+1:j*nz,(j-1)*nz+nu+1:j*nz)*...
            sys_new.F{j,1};
        
        sys_new.umax((j-1)*nu+1:j*nu,1)=DH_diag((j-1)*nz+1:(j-1)*nz+nu,...
            (j-1)*nz+1:(j-1)*nz+nu)*sys_new.umax((j-1)*nu+1:j*nu,1);
        sys_new.umin((j-1)*nu+1:j*nu,1)=DH_diag((j-1)*nz+1:(j-1)*nz+nu,...
            (j-1)*nz+1:(j-1)*nz+nu)*sys_new.umin((j-1)*nu+1:j*nu,1);
        
        sys_new.xmax((j-1)*nx+1:j*nx,1)=DH_diag((j-1)*nz+nu+1:(j-1)*nz+nu+nx,...
            (j-1)*nz+nu+1:(j-1)*nz+nu+nx)*sys_new.xmax((j-1)*nx+1:j*nx);
        sys_new.xmin((j-1)*nx+1:j*nx,1)=DH_diag((j-1)*nz+nu+1:(j-1)*nz+nu+nx,...
            (j-1)*nz+nu+1:(j-1)*nz+nu+nx)*sys_new.xmin((j-1)*nx+1:j*nx);
        sys_new.xs((j-1)*nx+1:j*nx,1)=DH_diag((j-1)*nz+nu+nx+1:j*nz,...
            (j-1)*nz+nu+nx+1:j*nz)*sys_new.xmin((j-1)*nx+1:j*nx);
    end
end

end






