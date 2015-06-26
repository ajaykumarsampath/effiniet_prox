close all;
clear all
clc
%%
kk=5300;
load('dwn');
P.Hp=24;
P.Hu=23;
P.beta=0.1;
Dd = DemandData(:,1); % My data
% Define data
tree_ops.Tmax = 800; % time steps for which data are available
tree_ops.nScen = 250; % no. of scenarios
tree_ops.N = 23; % prediction

load ('svm_error.mat');
SVM_ERROR=svm_demand_error(1:tree_ops.nScen,1:tree_ops.N);
SVM_ERROR=reshape(SVM_ERROR,tree_ops.nScen,1,tree_ops.N);
W = SVM_ERROR(1:tree_ops.nScen,:,1:tree_ops.N);
tree_ops.nw = size(W,2); % no. of component of the stochastic variable

% Compare all demands...
std_demand=zeros(S.nd,1);
range_demand=zeros(S.nd,1);
Demand_ratio=zeros(S.nd,1);

for i=1:S.nd
    temp=(DemandData(:,i)./Dd(:,1));
    std_demand(i,1)=std(temp);
    range_demand(i,1)=max(temp)-min(temp);
    Demand_ratio(i,1)=DemandData(1,i)./Dd(1,1);
end
%Demand_ratio=3600*Demand_ratio;
%% scenario tree generation
tree_ops.ni = ones(tree_ops.N,1); % max number of desired scenarios for each time stage
branching_factor = [3 2 1 1];
tree_ops.ni(1:length(branching_factor)) = branching_factor;
tree_ops.Wscaling = 1; % 1/0: do/don't scale disturbance W
tree_ops.Wfiltering = 0; % 1/0 do/don't filter disturbance W by system dynamics and u=Kx

[old_tree,scen_reduct_details]=Tree_formation(W,tree_ops);

old_tree{1}.value=kron(Demand_ratio',old_tree{1}.value);
old_tree=old_tree{1};
[Tree,new_tree_details]=Tranform_tree(old_tree);

Nscenarios=length(Tree.leaves);
Nstage=length(Tree.stage)+1;
%% 
% transform the effiniet into proximal algorithm 
ops_sys.gamma_min=inf;
ops_sys.gamma_xs=inf;
ops_sys.gamma_max=inf;
ops_sys.normalise=0;
ops_sys.cell=0;
sys_actual=system_prox_formation(S,P,Tree,ops_sys);
ops_sys.cell=1;
ops_sys_poly=0;
sys_NN=system_prox_formation(S,P,Tree,ops_sys);

[sys_box,V]=system_prox_formation_box(S,P,Tree,ops_sys);
%V.Wu=V.Wu/3600;
apg_opts.E=S.E;
apg_opts.Ed=S.Ed;
%%
Ns=size(Tree.leaves,1);
Np=sys_box.Np;

sys_temp.A=sys_box.A;
sys_temp.B=sys_box.B;
sys_temp.nx=sys_box.nx;
sys_temp.nu=sys_box.nu;

sys_temp.F=cell(Np+1,1);
sys_temp.G=cell(Np+1,1);
sys_temp.cell=1;

Tree_new.stage=(0:Np-1)';
Tree_new.leaves=Np;

for i=1:Np+1
   sys_temp.F{i}=sys_box.F{1};
   sys_temp.G{i}=sys_box.G{1};
end
sys_temp.umin=kron(ones(Np,1),sys_box.umin(1:sys_box.nu,1));
sys_temp.umax=kron(ones(Np,1),sys_box.umax(1:sys_box.nu,1));

sys_temp.xmin=kron(ones(Np+1,1),sys_box.xmin(1:sys_box.nx,1));
sys_temp.xmax=kron(ones(Np+1,1),sys_box.xmax(1:sys_box.nx,1));

SI=scenario_index(Tree);
DH_diag=zeros((sys_box.nx+sys_box.nu)*Np,Ns+1);

for i=1:Ns+1
    if(i==Ns+1)
        Tree_new.prob=ones(Np,1);
    else
        Tree_new.prob=Tree.prob(SI{i});
    end
    [DH_nml,dts_nml]=dual_hessian_calculate_box(sys_temp,Tree_new,V,apg_opts);
    details.norm_act=dts_nml.norm;
    
    sys_prcnd=precondition_calculate_box(sys_temp,DH_nml,Tree_new);
    [DH_nml_prcnd,dts_nml_prcnd]=dual_hessian_calculate_box(sys_prcnd,Tree_new,V,apg_opts);
    details.norm=dts_nml_prcnd.norm;
    DH_diag(:,i)=diag(DH_nml).^(-0.5);
end 
