%%
% This function eliminate the equality constraints and transform the water
% network into effinet model.

clear all
clear classes
clc
%%
kk=3000;
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
%ele=0*ele;
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
ops_sys_dst.gamma_min=inf;
ops_sys_dst.gamma_xs=inf;
ops_sys_dst.gamma_max=inf;

[sys_dst_bf_prcnd,V]=system_prox_formation(S,P,ops_sys_dst);

sys=sys_dst_bf_prcnd;

% Factor step calculation
Ptree=factor_apg_effiniet(sys_dst_bf_prcnd,V,Tree);
%% 
%current_demand=DemandData(kk:kk+P.Hu,:);
prev_vhat=sys_dst_bf_prcnd.L1*DemandData(kk-1,:)';
par_sol_opt.demand=DemandData(kk:kk+P.Hu,:);
par_sol_opt.prev_vhat=prev_vhat;
V.alpha=(kron(ones(P.Hp,1),P.alpha1')+P.alpha2(kk:kk+P.Hu,:));
current_state_opt=calcul_parti_soul(sys_dst_bf_prcnd,Tree,V,par_sol_opt);
current_state_opt.v=zeros(size(sys_dst_bf_prcnd.L,2),1);

opts_apg.state=current_state_opt;
opts_apg.lambda=0.1;

opts_apg.x=0.5*(sys_dst_bf_prcnd.xmax-sys_dst_bf_prcnd.xs)+sys_dst_bf_prcnd.xs;
opts_apg.primal_inf=0.01;
opts_apg.dual_inf=0.01;
opts_apg.steps=50;
%%
[Z,details_apg]=APG_effiniet(sys_dst_bf_prcnd,Ptree,Tree,V,opts_apg);
%%
apg_opts.u=sys.L*opts_apg.state.v+opts_apg.state.prev_vhat;
apg_opts.x=opts_apg.x;
apg_opts.E=S.E;
apg_opts.Ed=S.Ed;
apg_opts.demand=opts_apg.state.demand;

effiniet_apg=effiniet_yalmip(sys,Tree,V,apg_opts);
result=effiniet_apg{{apg_opts.x,apg_opts.u}};