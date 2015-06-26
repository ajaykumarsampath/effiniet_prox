%%
% This function eliminate the equality constraints and transform the water
% network into effinet model.

clear all
clear classes
clc
%%
kk=3320;
load('dwn');
P.Hp=24;
P.Hu=23;
Dd = DemandData(:,1); % My data
% Define data
tree_ops.Tmax = 800; % time steps for which data are available
tree_ops.nScen = 100; % no. of scenarios
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
%Demand_ratio=0*Demand_ratio;
%% scenario tree generation
tree_ops.ni = ones(tree_ops.N,1); % max number of desired scenarios for each time stage
branching_factor = [4 2 1 1];
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
sys_actual=system_prox_modified_box(S,P,Tree,ops_sys);
ops_sys.cell=1;
ops_sys_poly=0;

[sys_box,V]=system_prox_modified_box(S,P,Tree,ops_sys);
%V.Wu=V.Wu/3600;
apg_opts.E=S.E;
apg_opts.Ed=S.Ed;
%% Particular solution calculation
par_sol_opt.demand=3600*DemandData(kk:kk+P.Hu,:);
prev_vhat=3600*sys_box.L1*DemandData(kk-1,:)';

par_sol_opt.prev_vhat=prev_vhat;
%V.alpha=3600*(kron(ones(P.Hp,1),P.alpha1')+P.alpha2(kk:kk+P.Hu,:));
V.alpha=(kron(ones(P.Hp,1),P.alpha1')+P.alpha2(kk:kk+P.Hu,:));
current_state_opt=calParSol_modified(sys_box,Tree,V,par_sol_opt);
current_state_opt.v=3600*[0.0656 0.00 0.0849 0.0934]';

%%
Ptree_box=factor_apg_effinet_modified(sys_box,V,Tree);

x=0.1*(S.xmax-S.xmin)+S.xmin;
Nd=size(Tree.stage,1);
Ns=size(Tree.leaves,1);

Y.y=10*rand(size(sys_box.F{1},1),Nd);
Y.yt=10*rand(sys_box.nx,Ns);
Y.y(:,1)=zeros(9,1);
%%
opts_apg.state=current_state_opt;

opts_apg.x=x;
opts_apg.E=S.E;
opts_apg.Ed=S.Ed;
opts_apg.steps=10;
%opts_apg.lambda=1/dts_prmultiple lines comment in latexcnd.norm_act;

[Z1,details_apg]=solve_apg_effinet_modified(sys_box,Tree,Ptree_box,Y,x,opts_apg.state);
Z.X_APG=Z1.X;
Z.U_APG=Z1.U;
%grad_opts=current_state_opt;
%%
current_state_opt.x=x;
grad_opts.u=sys_box.L*current_state_opt.v+current_state_opt.prev_vhat;
grad_opts.x=current_state_opt.x;
grad_opts.y=Y.y;
grad_opts.yt=Y.yt;
grad_opts.E=S.E;
grad_opts.Ed=S.Ed;
grad_opts.demand=current_state_opt.demand;
dual_grad_yalmip=dual_gradient_modified_yalmip(sys_box,Tree,V,grad_opts);
[results,error]=dual_grad_yalmip{{grad_opts.x,grad_opts.u,Y.y,Y.yt}};
Z.X_yalmip=results{1,1};
Z.U_yalmip=results{1,2};
%%
grad_null_opts.x=current_state_opt.x;
grad_null_opts.demand=current_state_opt.demand;

grad_null_opts.vhat=current_state_opt.vhat;
grad_null_opts.prev_vhat=current_state_opt.prev_vhat;
grad_null_opts.alpha_bar=current_state_opt.alpha_bar;
grad_null_opts.w=current_state_opt.w;
grad_null_opts.v=current_state_opt.v;

dual_grad_yalmip=dual_grad_modif_yalm_null(sys_box,Tree,Ptree_box,V,grad_null_opts);
[results_1,error]=dual_grad_yalmip{{grad_opts.x,grad_null_opts.v,Y.y,Y.yt}};
Z.X_yalmip_null=results_1{1,1};
Z.U_yalmip_null=results_1{1,3};
%%
max(max(abs(Z.X_yalmip_null-Z.X_yalmip)))
max(max(abs(Z.U_yalmip_null-Z.U_yalmip)))

max(max(abs(Z.U_APG-Z.U_yalmip)))