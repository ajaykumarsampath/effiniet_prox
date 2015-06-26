close all;
clear all
clc
%%
kk=4000;
load('dwn');
P.Hp=24;
P.Hu=23;
P.beta=0.1;
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
%ele=0*ele;
%% scenario tree generation
tree_ops.ni = ones(tree_ops.N,1); % max number of desired scenarios for each time stage
branching_factor = [1 1 1 1];
tree_ops.ni(1:length(branching_factor)) = branching_factor;
tree_ops.Wscaling = 1; % 1/0: do/don't scale disturbance W
tree_ops.Wfiltering = 0; % 1/0 do/don't filter disturbance W by system dynamics and u=Kx

[old_tree,scen_reduct_details]=Tree_formation(W,tree_ops);

old_tree{1}.value=kron(Demand_ratio',old_tree{1}.value);
old_tree=old_tree{1};
[Tree,new_tree_details]=Tranform_tree(old_tree);
%% 
% transform the effiniet into proximal algorithm 
ops_sys.gamma_min=inf;
ops_sys.gamma_xs=inf;
ops_sys.gamma_max=inf;
ops_sys.normalise=0;
ops_sys.cell=0;
sys_actual=system_prox_formation(S,P,Tree,ops_sys);
sys_actual_mod=system_prox_formation_modified(S,P,Tree,ops_sys);
ops_sys.cell=1;
ops_sys_poly=0;

sys_NN=system_prox_formation(S,P,Tree,ops_sys);

[sys_box,V]=system_prox_formation_box(S,P,Tree,ops_sys);

apg_opts.E=S.E;
apg_opts.Ed=S.Ed;

% Preconditioner 

[DH_nml,dts_nml]=dual_hessian_calculate_box(sys_box,Tree,V,apg_opts);

%[DH_nml_l,dts_nml_l]=dual_hessian_calculate(sys_NN,Tree,V,apg_opts);
%% Preconditioning 
opts_prcnd.dual_hessian=1;
sys_prcnd=precondition_calculate_box(sys_box,DH_nml,Tree);
[DH_nml_prcnd,dts_nml_prcnd]=dual_hessian_calculate_box(sys_prcnd,Tree,V,apg_opts);
l=eig(DH_nml_prcnd);
[sys_prcnd_act,dts_prcnd]=prcnd_sys_tree_box(sys_box,Tree,V,apg_opts);
%{
tic
[Pre,t]=precond_sdp(DH_nml);
toc
opts_prcnd.dual_hessian=0;
sys_prcnd_act=precondition_calculate_box(sys_box,Pre,Tree,opts_prcnd);
[DH_nml_prcnd_act,dts_nml_prcnd_act]=dual_hessian_calculate_box(sys_prcnd_act,Tree,V,apg_opts);
l_act=eig(DH_nml_prcnd_act);
%}
%% Particular solution calculation
par_sol_opt.demand=3600*DemandData(kk:kk+P.Hu,:);
prev_vhat=3600*sys_box.L1*DemandData(kk-1,:)';

par_sol_opt.prev_vhat=prev_vhat;
%V.alpha=3600*(kron(ones(P.Hp,1),P.alpha1')+P.alpha2(kk:kk+P.Hu,:));
V.alpha=(kron(ones(P.Hp,1),P.alpha1')+P.alpha2(kk:kk+P.Hu,:));
current_state_opt=calParSol(sys_box,Tree,V,par_sol_opt);
current_state_opt.v=3600*[0.0656 0.00 0.0849 0.0934]';

%% Factor step

Ptree_box=factor_apg_effinet(sys_box,V,Tree);

Ptree_prcnd=factor_apg_effinet(sys_prcnd,V,Tree);

Ptree_prcnd_act=factor_apg_effinet(sys_prcnd_act,V,Tree);

%% Solve steps 
opts_apg.state=current_state_opt;

opts_apg.x=0.1*(S.xmax-P.xs)+P.xs;
opts_apg.E=S.E;
opts_apg.Ed=S.Ed;

opts_apg.steps=1000;
opts_apg.lambda=1/dts_nml.norm;
%opts_apg.lambda=2.1492e-05;
[Z_box,details_apg_box]=APG_effinet_box(sys_box,Ptree_box,Tree,V,opts_apg);

opts_apg.steps=1000;
opts_apg.lambda=1/dts_nml_prcnd.norm;
[Z_prcnd,details_apg_prcnd]=APG_effinet_box(sys_prcnd,Ptree_prcnd,Tree,V,opts_apg);

opts_apg.steps=1000;
opts_apg.lambda=1/dts_nml_prcnd.norm;
[Z_prcnd_act,details_apg_prcnd_act]=APG_effinet_box(sys_prcnd_act,...
    Ptree_prcnd_act,Tree,V,opts_apg);

%% Gurobi Algorithm 
apg_opts.u=sys_actual.L*opts_apg.state.v+opts_apg.state.prev_vhat;
apg_opts.x=opts_apg.x;
apg_opts.E=S.E;
apg_opts.Ed=S.Ed;
apg_opts.demand=opts_apg.state.demand;

effinet_apg=effinet_yalmip(sys_actual,Tree,V,apg_opts);
tic
[result,error]=effinet_apg{{apg_opts.x,apg_opts.u}};
toc
Z.eff_X=result{1,1};
Z.eff_U=result{1,2};
%%
figure
for i=1:sys_actual.nx
    subplot(3,1,i)
    plot(Z_box.X(i,:));
    hold all;
    plot(Z_prcnd.X(i,:))
    plot(Z.eff_X(i,:));
    plot(Z_prcnd_act.X(i,:))
end
legend('APG','APG-diag','gurobi','APG-EXACT')

figure
for i=1:sys_actual.nu
    subplot(3,2,i)
    plot(Z_box.U(i,:));
    hold all;
    plot(Z_prcnd.U(i,:))
     plot(Z.eff_U(i,:));
    plot(Z_prcnd_act.U(i,:));
end
legend('APG','APG-diag','gurobi','APG-EXACT')

figure
plot(details_apg_box.jobj)
hold all;
plot(details_apg_prcnd.jobj)
plot(details_apg_prcnd_act.jobj)
legend('APG','APG-diag','APG-EXACT')
