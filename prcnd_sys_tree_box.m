function [ sys_new,details] = prcnd_sys_tree_box(sys_box,Tree,V,apg_opts)

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
%Nd=size(Tree.stage,1);
%Ns=size(Tree.leaves,1);
Np=sys_box.Np;
ny=size(sys_box.F{1},1);

sys_temp.A=sys_box.A;
sys_temp.B=sys_box.B;
sys_temp.nx=sys_box.nx;
sys_temp.nu=sys_box.nu;

sys_temp.F=cell(Np+1,1);
sys_temp.G=cell(Np+1,1);
sys_temp.cell=1;

Tree_new.stage=(0:Np-1)';
Tree_new.leaves=Np;
Tree_new.prob=ones(Np,1);

for i=1:Np+1
   sys_temp.F{i}=sys_box.F{1};
   sys_temp.G{i}=sys_box.G{1};
end
sys_temp.umin=kron(ones(Np,1),sys_box.umin(1:nu,1));
sys_temp.umax=kron(ones(Np,1),sys_box.umax(1:nu,1));

sys_temp.xmin=kron(ones(Np+1,1),sys_box.xmin(1:nx,1));
sys_temp.xmax=kron(ones(Np+1,1),sys_box.xmax(1:nx,1));


[DH_nml,dts_nml]=dual_hessian_calculate_box(sys_temp,Tree_new,V,apg_opts);
details.norm_act=dts_nml.norm;

sys_prcnd=precondition_calculate_box(sys_temp,DH_nml,Tree_new);
[DH_nml_prcnd,dts_nml_prcnd]=dual_hessian_calculate_box(sys_prcnd,Tree_new,V,apg_opts);
details.norm=dts_nml_prcnd.norm;


DH_diag=diag(diag(DH_nml).^(-0.5));

prob=sqrt(Tree.prob);

sys_new=sys_box;

for k=1:Np+1
    
    if(k==1)
        
        sys_new.G{k,1}=sys_box.G{k,1};
        sys_new.G{k,1}(nx+1:end,:)=DH_diag(1:nu,1:nu)*sys_box.G{k,1}(nx+1:end,:);
        sys_new.umax(1:nu,1)=DH_diag(1:nu,1:nu)*sys_box.umax(1:nu,1);
        sys_new.umin(1:nu,1)=DH_diag(1:nu,1:nu)*sys_box.umin(1:nu,1);
    par_sol_opt.demand=3600*DemandData(kk:kk+P.Hu,:);
prev_vhat=3600*sys_box.L1*DemandData(kk-1,:)';

par_sol_opt.prev_vhat=prev_vhat;
    elseif(k==Np+1)
        
        nodes_stage=Tree.leaves;
        
        for j=1:length(nodes_stage)
            sys_new.F{nodes_stage(j)+1,1}(1:nx,:)=prob(nodes_stage(j))*...
                DH_diag(ny*(Np-1)+nu+1:end,ny*(Np-1)+nu+1:end)...
                *sys_box.F{nodes_stage(j)+1,1}(1:nx,:);
            sys_new.xmax(nodes_stage(j)*nx+1:(nodes_stage(j)+1)*nx)=...
                prob(nodes_stage(j))*DH_diag((Np-1)*ny+nu+1:(Np-1)*ny+nu+nx,(Np-1)*ny+nu+1:...
                (Np-1)*ny+nu+nx)*sys_box.xmax(nodes_stage(j)*nx+1:(nodes_stage(j)+1)*nx);
            sys_new.xmin(nodes_stage(j)*nx+1:(nodes_stage(j)+1)*nx)=...
                prob(nodes_stage(j))*DH_diag((Np-1)*ny+nu+1:(Np-1)*ny+nu+nx,(Np-1)*ny+nu+1:...
                (Np-1)*ny+nu+nx)*sys_box.xmin(nodes_stage(j)*nx+1:(nodes_stage(j)+1)*nx);
        end
    else
        nodes_stage=find(Tree.stage==k-2);
        
        for j=1:length(nodes_stage)
            sys_new.F{nodes_stage(j)+1,1}(1:nx,:)=prob(nodes_stage(j))*...
                DH_diag(ny*(k-2)+nu+1:ny*(k-1),ny*(k-2)+nu+1:ny*(k-1))...
                *sys_box.F{nodes_stage(j)+1,1}(1:nx,:);
            sys_new.G{nodes_stage(j)+1,1}(nx+1:end,:)=prob(nodes_stage(j))*...
                DH_diag(ny*(k-1)+1:ny*(k-1)+nu,ny*(k-1)+1:ny*(k-1)+nu)...
                *sys_box.G{nodes_stage(j)+1,1}(nx+1:end,:);
            
            sys_new.xmax(nodes_stage(j)*nx+1:(nodes_stage(j)+1)*nx,1)=prob(nodes_stage(j))*...
                DH_diag((k-2)*ny+nu+1:(k-1)*ny,(k-2)*ny+nu+1:(k-1)*ny)...
                *sys_box.xmax(nodes_stage(j)*nx+1:(nodes_stage(j)+1)*nx,1);
            sys_new.xmin(nodes_stage(j)*nx+1:(nodes_stage(j)+1)*nx,1)=prob(nodes_stage(j))*...
                DH_diag((k-2)*ny+nu+1:(k-1)*ny,(k-2)*ny+nu+1:(k-1)*ny)...
                *sys_box.xmin(nodes_stage(j)*nx+1:(nodes_stage(j)+1)*nx,1);
            
            sys_new.umax(nodes_stage(j)*nu+1:(nodes_stage(j)+1)*nu,1)=prob(nodes_stage(j))*...
                DH_diag((k-1)*ny+1:(k-1)*ny+nu,(k-1)*ny+1:(k-1)*ny+nu)...
                *sys_box.umax(nodes_stage(j)*nu+1:(nodes_stage(j)+1)*nu,1);
            sys_new.umin(nodes_stage(j)*nu+1:(nodes_stage(j)+1)*nu,1)=prob(nodes_stage(j))*...
                DH_diag((k-1)*ny+1:(k-1)*ny+nu,(k-1)*ny+1:(k-1)*ny+nu)...
                *sys_box.umin(nodes_stage(j)*nu+1:(nodes_stage(j)+1)*nu,1);
        end
    end
end 
%%
%sys_temp_prcnd=precondition_calculate_box(sys_temp,DH_nml,Tree_new);



end

