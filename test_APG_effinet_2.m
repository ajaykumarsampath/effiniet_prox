%%
% This function eliminate the equality constraints and transform the water
% network into effinet model.

clear all
clear classes
clc
%%
kk=3300;
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
ele=zeros(S.nd,1);

for i=1:S.nd
    temp=(DemandData(:,i)./Dd(:,1));
    std_demand(i,1)=std(temp);
    range_demand(i,1)=max(temp)-min(temp);
    ele(i,1)=DemandData(1,i)./Dd(1,1);
end
ele=0*ele;
%% scenario tree generation
tree_ops.ni = ones(tree_ops.N,1); % max number of desired scenarios for each time stage
branching_factor = [1 1 1 1];
tree_ops.ni(1:length(branching_factor)) = branching_factor;
tree_ops.Wscaling = 1; % 1/0: do/don't scale disturbance W
tree_ops.Wfiltering = 0; % 1/0 do/don't filter disturbance W by system dynamics and u=Kx

[old_tree,scen_reduct_details]=Tree_formation(W,tree_ops);

old_tree{1}.value=kron(ele',old_tree{1}.value);
old_tree=old_tree{1};
[Tree,new_tree_details]=Tranform_tree(old_tree);
%% APG algorithm 
% transform the effiniet into proximal algorithm 

%
ops_sys.gamma_min=inf;
ops_sys.gamma_xs=inf;
ops_sys.gamma_max=inf;
ops_sys.normalise=0;
ops_sys.cell=0;
sys_actual=system_prox_formation(S,P,Tree,ops_sys);
sys_actual_mod=system_prox_formation_modified(S,P,Tree,ops_sys);
ops_sys.cell=1;
ops_sys_poly=0;
%{
ops_sys_dst.gamma_min=1e6;
ops_sys_dst.gamma_xs=1e6;
ops_sys_dst.gamma_max=1e6;
%}
sys_NN=system_prox_formation(S,P,Tree,ops_sys);

ops_sys.normalise=1;
%sys=sys_old;
sys=system_prox_formation(S,P,Tree,ops_sys);
[sys_mod,V_mod]=system_prox_formation_modified(S,P,Tree,ops_sys);

[sys_box,V]=system_prox_formation_box(S,P,Tree,ops_sys);

ops_sys.cell=0;
ops_sys.normalise=0;

apg_opts.E=S.E;
apg_opts.Ed=S.Ed;

apg_opts_mod.E=S.E(:,[1:3 5:6]);
apg_opts_mod.Ed=S.Ed;

[sys_null,V_null]=system_prox_formation_null(S,P,Tree,ops_sys);
%%
[DH_normalized_not,dts_not_normalized]=dual_hessian_calculate(sys_NN,Tree,V,apg_opts);
[DH_normalized,dts_normalized]=dual_hessian_calculate(sys,Tree,V,apg_opts);
[DH_mod,dts_mod]=dual_hessian_calculate(sys_mod,Tree,V_mod,apg_opts_mod);
% sys_new=sys;
sys_prcnd=precondition_calculate(sys,DH_normalized,Tree);
[DH_normalized_new,dts_normalized_new]=dual_hessian_calculate(sys_prcnd,Tree,V,apg_opts);

%%


% Factor step calculation
Ptree_NN=factor_apg_effiniet(sys_NN,V,Tree);
Ptree=factor_apg_effiniet(sys,V,Tree);
Ptree_prcnd=factor_apg_effiniet(sys_prcnd,V,Tree);

Ptree_mod=factor_apg_effiniet(sys_mod,V_mod,Tree);
% Factor step calculation
Ptree_null=factor_apg_effiniet_null(sys_null,V,Tree);

%% 
%current_demand=DemandData(kk:kk+P.Hu,:);
par_sol_opt.demand=3600*DemandData(kk:kk+P.Hu,:);
prev_vhat=3600*sys_NN.L1*DemandData(kk-1,:)';

par_sol_opt.prev_vhat=prev_vhat;
V.alpha=(kron(ones(P.Hp,1),P.alpha1')+P.alpha2(kk:kk+P.Hu,:));
current_state_opt=calcul_parti_soul(sys_NN,Tree,V,par_sol_opt);
%current_state_opt.v=0.1*rand(size(sys_dst_bf_prcnd.L,2),1);
current_state_opt.v=3600*[0.0656 0.00 0.0849 0.0934]';
opts_apg.state=current_state_opt;

opts_apg.x=0.1*(S.xmax-P.xs)+P.xs;
opts_apg.primal_inf=0.01;
opts_apg.dual_inf=0.01;
opts_apg.E=S.E;
opts_apg.Ed=S.Ed;


%%

V_mod.alpha=(kron(ones(P.Hp,1),P.alpha1([1:3 5:6])')+P.alpha2(kk:kk+P.Hu,[1:3 5:6]));
prev_vhat=3600*sys_mod.L1*DemandData(kk-1,:)';

par_sol_opt.prev_vhat=prev_vhat;
current_state_opt_mod=calcul_parti_soul(sys_mod,Tree,V_mod,par_sol_opt);
current_state_opt_mod.v=3600*[0.0656 0.0849 0.0934]';
opts_apg_mod.state=current_state_opt_mod;

opts_apg_mod.x=0.1*(S.xmax-P.xs)+P.xs;
opts_apg_mod.primal_inf=0.01;
opts_apg_mod.dual_inf=0.01;
opts_apg_mod.E=S.E;
opts_apg_mod.Ed=S.Ed;
%%
opts_apg.steps=4000;
opts_apg.lambda=18e-4;

[Z,details_apg]=APG_effiniet_2(sys,Ptree,Tree,V,opts_apg);

opts_apg.lambda=18e-6;
[Z_prcnd,details_apg_prcnd]=APG_effiniet_2(sys_prcnd,Ptree_prcnd,Tree,V,opts_apg);

opts_apg.steps=4000;
opts_apg.lambda=10e-6;
%[Z,details_apg]=APG_effiniet_2(sys_dst_bf_prcnd,Ptree,Tree,V,opts_apg);
[Z_NN,details_apg_NN]=APG_effiniet_2(sys_NN,Ptree_NN,Tree,V,opts_apg);
%%
opts_apg.steps=8000;
opts_apg.lambda=18e-6;
[Z_null,details_apg_null]=APG_effiniet_null(sys_null,Ptree_null,Tree,V,opts_apg);
%%
opts_apg_mod.steps=8000;
opts_apg_mod.lambda=18e-4*5;
[Z_mod,details_apg_mod]=APG_effiniet_2(sys_mod,Ptree_mod,Tree,V_mod,opts_apg_mod);
%%
%sys.xmax=10*sys.xmax;
apg_opts.u=sys.L*opts_apg.state.v+opts_apg.state.prev_vhat;
apg_opts.x=opts_apg.x;
apg_opts.E=S.E;
apg_opts.Ed=S.Ed;
apg_opts.demand=opts_apg.state.demand;

effiniet_apg=effiniet_yalmip(sys_actual,Tree,V,apg_opts);
tic
[result,error]=effiniet_apg{{apg_opts.x,apg_opts.u}};
toc
Z.eff_X=result{1,1};
Z.eff_U=result{1,2};

%%
apg_opts_mod.u=sys_mod.L*opts_apg_mod.state.v+opts_apg_mod.state.prev_vhat;
apg_opts_mod.x=opts_apg.x;
apg_opts_mod.demand=opts_apg.state.demand;
effiniet_apg_mod=effiniet_yalmip(sys_actual_mod,Tree,V_mod,apg_opts_mod);
tic
[result_mod,error]=effiniet_apg_mod{{apg_opts_mod.x,apg_opts_mod.u}};
toc
Z_mod.eff_X=result_mod{1,1};
Z_mod.eff_U=result_mod{1,2};
%%
figure
plot(details_apg.jobj);
hold all;
plot(details_apg_null.jobj);
figure
plot(Z.X(1,:));
hold all; 
plot(Z_null.X(1,:))
plot(Z.eff_X(1,:));
%%
figure
subplot(311)
plot(Z_mod.X(1,:));
hold all;
plot(Z_mod.eff_X(1,:));

subplot(312)
plot(Z_mod.X(2,:));
hold all;
plot(Z_mod.eff_X(2,:));

subplot(313)
plot(Z_mod.X(3,:));
hold all;
plot(Z_mod.eff_X(3,:));
