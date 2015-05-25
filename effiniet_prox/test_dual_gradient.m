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
branching_factor = [1 1 1];
tree_ops.ni(1:length(branching_factor)) = branching_factor;
tree_ops.Wscaling = 1; % 1/0: do/don't scale disturbance W
tree_ops.Wfiltering = 1; % 1/0 do/don't filter disturbance W by system dynamics and u=Kx
%tree_ops.fig=1;
[old_tree,scen_reduct_details]=Tree_formation(W,tree_ops);
%ele=0*ele;
old_tree{1}.value=kron(ele',old_tree{1}.value);
old_tree=old_tree{1};
[Tree,new_tree_details]=Tranform_tree(old_tree);
%{
% testing the correctness of the tree formed 
SI=scenario_index(old_tree);
SI_new=scenario_index(Tree);
for i=1:size(SI,1)
    value(i,1)=old_tree.prob(SI{i}(end,1));
    value(i,2)=Tree.prob(SI_new{(new_tree_details.scenario_order(i))}(end,1));
end
%}
%% APG algorithm 
% transform the effiniet into proximal algorithm 

ops_sys_dst.gamma_min=1000;
ops_sys_dst.gamma_xs=1000;
ops_sys_dst.gamma_max=1000;
ops_sys_dst.normalise=0;
ops_sys_dst.cell=1;

[sys_dst_bf_prcnd,V]=system_prox_formation(S,P,Tree,ops_sys_dst);
prev_vhat=sys_dst_bf_prcnd.L1*3600*DemandData(kk-1,:)';

sys=sys_dst_bf_prcnd;
%current_demand=DemandData(kk:kk+P.Hu,:);

par_sol_opt.demand=3600*DemandData(kk:kk+P.Hu,:);
par_sol_opt.prev_vhat=prev_vhat;

V.alpha=3600*(kron(ones(P.Hp,1),P.alpha1')+P.alpha2(kk:kk+P.Hu,:));
%V.alpha=(kron(ones(P.Hp,1),P.alpha1')+P.alpha2(kk:kk+P.Hu,:));
Ptree=factor_apg_effiniet(sys_dst_bf_prcnd,V,Tree);

current_state_opt=calcul_parti_soul(sys_dst_bf_prcnd,Tree,V,par_sol_opt);


%%
x=0.1*(S.xmax-S.xmin)+S.xmin;
%current_state_opt.v=zeros(size(sys_dst_bf_prcnd.L,2),1);
current_state_opt.v=3600*[0.0656 0.00 0.0849 0.0934]';
Y.y=1000*rand(size(sys_dst_bf_prcnd.F{1},1),(size(Tree.children,1)+1));
Y.yt=1000*rand(2*sys_dst_bf_prcnd.nx,size(Tree.leaves,1),1);
%Z=solve_apg_effiniet(sys_dst_bf_prcnd,Tree,V,Ptree,current_state_opt);
[Z,details_apg]=solve_apg_effiniet(sys_dst_bf_prcnd,Tree,Ptree,Y,x,current_state_opt);
%grad_opts=current_state_opt;
%%
current_state_opt.x=x;
grad_opts.u=sys.L*current_state_opt.v+current_state_opt.prev_vhat;
grad_opts.x=current_state_opt.x;
grad_opts.y=Y.y;
grad_opts.yt=Y.yt;
grad_opts.E=S.E;
grad_opts.Ed=S.Ed;
grad_opts.demand=current_state_opt.demand;
dual_grad_yalmip=dual_gradient_yalmip(sys,Tree,V,grad_opts);
[results_1,error]=dual_grad_yalmip{{grad_opts.x,grad_opts.u,Y.y,Y.yt}};
Z.X_yalmip=results_1{1,1};
Z.U_yalmip=results_1{1,2};
%%
grad_opts_null.x=current_state_opt.x;
grad_opts_null.y=Y.y;
grad_opts_null.yt=Y.yt;
grad_opts_null.E=S.E;
grad_opts_null.Ed=S.Ed;
grad_opts_null.demand=current_state_opt.demand;

grad_opts_null.vhat=current_state_opt.vhat;
grad_opts_null.prev_vhat=current_state_opt.prev_vhat;
grad_opts_null.alpha_bar=current_state_opt.alpha_bar;
grad_opts_null.w=current_state_opt.w;
grad_opts_null.v=current_state_opt.v;

dual_grad_yalmip_null=dual_gradient_yalmip_null(sys,Tree,Ptree,V,grad_opts_null);
[results,error]=dual_grad_yalmip_null{{grad_opts_null.x,grad_opts_null.v,Y.y,Y.yt}};
Z.Xyalmip_null=results{1,1};
Z.Uyalmip_null=results{1,3};
pp_null{1}=Z.X_yalmip-Z.Xyalmip_null;
pp_null_u{1}=Z.U_yalmip-Z.Uyalmip_null;
%max(max(abs(results{1,1}-results_1{1,1})))
%max(max(results{1,3}-results_1{1,2}))

max(max(abs(Z.X-Z.X_yalmip)))

pp_null{2}=Z.X-Z.Xyalmip_null;
pp_null_u{2}=Z.U-Z.Uyalmip_null;
pp_null{3}=Z.X-Z.X_yalmip;
pp_null_u{3}=Z.U-Z.U_yalmip;