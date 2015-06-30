function [ sys_new,details] = prcnd_sys_dist_box(sys_box,Tree,V,apg_opts)

% This function preconditions the system with scenario trees.
%
% INPUT-----    sys_box      :
%               Tree         :
%               V            :
%               apg_opts     :
%
% OUTPUT-----   sys_new      :
%
%%
nx=sys_box.nx;
nu=sys_box.nu;
nz=2*nx+nu;
Nd=size(Tree.stage,1);
Np=sys_box.Np;

sys_temp.A=sys_box.A;
sys_temp.B=sys_box.B;
sys_temp.nx=sys_box.nx;
sys_temp.nu=sys_box.nu;
sys_temp.Np=sys_box.Np;

sys_temp.F=cell(Np,1);
sys_temp.G=cell(Np,1);
sys_temp.cell=1;

Tree_new.stage=(0:Np-1)';
Tree_new.leaves=Np;
Tree_new.prob=ones(Np,1);

for i=1:Np
    sys_temp.F{i}=sys_box.F{1};
    sys_temp.G{i}=sys_box.G{1};
end
sys_temp.umin=kron(ones(Np,1),sys_box.umin(1:nu,1));
sys_temp.umax=kron(ones(Np,1),sys_box.umax(1:nu,1));

sys_temp.xmin=kron(ones(Np,1),sys_box.xmin(1:nx,1));
sys_temp.xmax=kron(ones(Np,1),sys_box.xmax(1:nx,1));
sys_temp.xs=kron(ones(Np,1),sys_box.xs(1:nx,1));

[DH_nml,dts_nml]=dual_hessian_dist_box(sys_temp,Tree_new,V,apg_opts);
details.norm_act=dts_nml.norm;

if(apg_opts.exact_prcnd)
    [sys_prcnd,details.pre_cnd]=precond_modify_sdp_box(sys_temp,DH_nml);
    
    [DH_nml_prcnd,dts_nml_prcnd]=dual_hessian_modfy_box(sys_prcnd,Tree_new,V,apg_opts);
    
    details.norm=dts_nml_prcnd.norm;
    
    prob=sqrt(Tree.prob);
    
    sys_new=sys_box;
    
    for j=1:Nd
        
        sys_new.G{j,1}(nx+1:end,:)=prob(j)*details.pre_cnd.Eg1*sys_box.G{j,1}(nx+1:end,:);
        sys_new.umax((j-1)*nu+1:j*nu,1)=prob(j)*details.pre_cnd.Eg1*sys_box.umax((j-1)*nu+1:j*nu,1);
        sys_new.umin((j-1)*nu+1:j*nu,1)=prob(j)*details.pre_cnd.Eg1*sys_box.umin((j-1)*nu+1:j*nu,1);
        
        sys_new.F{j,1}(1:nx,:)=prob(j)*details.pre_cnd.Ef1*sys_box.F{j,1}(1:nx,:);
        sys_new.xmax((j-1)*nx+1:j*nx,1)=prob(j)*details.pre_cnd.Ef1*sys_box.xmax((j-1)*nx+1:j*nx,1);
        sys_new.xmin((j-1)*nx+1:j*nx,1)=prob(j)*details.pre_cnd.Ef1*sys_box.xmin((j-1)*nx+1:j*nx,1);
        sys_new.xs((j-1)*nx+1:j*nx,1)=prob(j)*details.pre_cnd.Ef1*sys_box.xs((j-1)*nx+1:j*nx,1);
    end
else
    sys_prcnd=prcnd_diag_dist_box(sys_temp,DH_nml);
    
    [DH_nml_prcnd,dts_nml_prcnd]=dual_hessian_dist_box(sys_prcnd,Tree_new,V,apg_opts);
    details.norm=dts_nml_prcnd.norm;
    
    DH_diag=diag(diag(DH_nml).^(-0.5));
    prob=sqrt(Tree.prob);
    %prob=ones(Nd,1);
    sys_new=sys_box;
    
    for j=1:Nd
        
        sys_new.G{j,1}=sys_box.G{j,1};
        sys_new.F{j,1}=sys_box.F{j,1};
        
        k=Tree.stage(j)+1;
        sys_new.G{j,1}=prob(j)*DH_diag((k-1)*nz+1:(k-1)*nz+nu,(k-1)*nz+1:(k-1)*nz+nu)*...
            sys_box.G{j,1};
        sys_new.umax((j-1)*nu+1:j*nu,1)=prob(j)*DH_diag((k-1)*nz+1:(k-1)*nz+nu,(k-1)*nz+1:(k-1)*nz+nu)*...
            sys_box.umax((j-1)*nu+1:j*nu,1);
        sys_new.umin((j-1)*nu+1:j*nu,1)=prob(j)*DH_diag((k-1)*nz+1:(k-1)*nz+nu,(k-1)*nz+1:(k-1)*nz+nu)*...
            sys_box.umin((j-1)*nu+1:j*nu,1);
        
        sys_new.F{j,1}=prob(j)*DH_diag((k-1)*nz+nu+1:k*nz,(k-1)*nz+nu+1:k*nz)*...
            sys_box.F{j,1};
        sys_new.xmax((j-1)*nx+1:j*nx,1)=prob(j)*DH_diag((k-1)*nz+nu+1:(k-1)*nz+nu+nx,...
            (k-1)*nz+nu+1:(k-1)*nz+nu+nx)*sys_box.xmax((j-1)*nx+1:j*nx,1);
        sys_new.xmin((j-1)*nx+1:j*nx,1)=prob(j)*DH_diag((k-1)*nz+nu+1:(k-1)*nz+nu+nx,...
            (k-1)*nz+nu+1:(k-1)*nz+nu+nx)*sys_box.xmin((j-1)*nx+1:j*nx,1);
        sys_new.xs((j-1)*nx+1:j*nx,1)=prob(j)*DH_diag((k-1)*nz+nu+nx+1:k*nz,(k-1)*nz+nu+nx+1:k*nz)*...
            sys_box.xs((j-1)*nx+1:j*nx,1);
    end
end

%%
%sys_temp_prcnd=precondition_calculate_box(sys_temp,DH_nml,Tree_new);



end


