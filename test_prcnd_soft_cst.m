close all;
clear all
clc
%%
kk=3875;
load('dwn');
P.Hp=24;
P.Hu=23;
P.beta=0.1;
Dd = DemandData(:,1); % My data
% Define data
tree_ops.Tmax = 800; % time steps for which data are available
tree_ops.nScen = 200; % no. of scenarios
tree_ops.N = 23; % prediction

load ('svm_error.mat');
SVM_ERROR=svm_demand_error(1:tree_ops.nScen,1:tree_ops.N);
SVM_ERROR=reshape(SVM_ERROR,tree_ops.nScen,1,tree_ops.N);
W = SVM_ERROR(1:tree_ops.nScen,:,1:tree_ops.N);
%W=3600*W;
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
%Demand_ratio=36*Demand_ratio;
%% scenario tree generation
tree_ops.ni = ones(tree_ops.N,1); % max number of desired scenarios for each time stage
branching_factor = [5 2 2 1];
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
ops_sys.gamma_xs=inf;
ops_sys.gamma_xbox=inf;
ops_sys.normalise=0;
ops_sys.cell=0;
sys_actual=system_prox_dist_state_box(S,P,Tree,ops_sys);
ops_sys.cell=1;
ops_sys_poly=0;

sys_box_dist=system_prox_dist_state_box(S,P,Tree,ops_sys);

[sys_box,V]=system_prox_modif_state_box(S,P,Tree,ops_sys);

%V.Wu=V.Wu/3600;
apg_opts.E=S.E;
apg_opts.Ed=S.Ed;
apg_opts.exact_prcnd=0;
%% Particular solution calculation
par_sol_opt.demand=3600*DemandData(kk:kk+P.Hu,:);
prev_vhat=3600*sys_box_dist.L1*DemandData(kk-1,:)';

par_sol_opt.prev_vhat=prev_vhat;
%V.alpha=3600*(kron(ones(P.Hp,1),P.alpha1')+P.alpha2(kk:kk+P.Hu,:));
V.alpha=(kron(ones(P.Hp,1),P.alpha1')+P.alpha2(kk:kk+P.Hu,:));
current_state_opt=calParSol_modified(sys_box_dist,Tree,V,par_sol_opt);
current_state_opt.v=3600*[0.0656 0.00 0.0849 0.0934]';
%% Preconditioning 
opts_prcnd.dual_hessian=1;
opts_prcnd.prcnd=0;

[sys_prcnd_dist,dts_prcnd_dist]=prcnd_sys_dist_box(sys_box_dist,Tree,V,apg_opts);
[sys_prcnd,dts_prcnd]=prcnd_sys_modfy_box(sys_box,Tree,V,apg_opts);


%% Factor step

Ptree_prcnd_dist=factor_apg_eff_dist_state(sys_prcnd_dist,V,Tree);
Ptree_prcnd=factor_apg_eff_modif_state(sys_prcnd,V,Tree);
%% Solve steps  
% solve step options 
opts_apg.state=current_state_opt;

opts_apg.x=0.5*(S.xmax-P.xs)+P.xs;
opts_apg.E=S.E;
opts_apg.Ed=S.Ed;
opts_apg.constraints='soft';
tic
opts_apg.steps=200;
opts_apg.lambda=1.3*1/dts_prcnd_dist.norm;
[Z_prcnd_dist,details_apg_prcnd_dist]=APG_effinet_dist_box(sys_prcnd_dist,Ptree_prcnd_dist,Tree,V,opts_apg);
toc

tic
opts_apg.constraints='soft';
opts_apg.steps=200;
opts_apg.lambda=1.3*1/dts_prcnd.norm;
[Z_prcnd,details_apg_prcnd]=APG_effinet_box_modified(sys_prcnd,Ptree_prcnd,Tree,V,opts_apg);
toc
%% Gurobi Algorithm 
apg_opts.u=sys_actual.L*opts_apg.state.v+opts_apg.state.prev_vhat;
apg_opts.x=opts_apg.x;
apg_opts.E=S.E;
apg_opts.Ed=S.Ed;
apg_opts.demand=opts_apg.state.demand;

effinet_apg=effinet_modif_state_yalmip(sys_actual,Tree,V,apg_opts);
tic
[result,error]=effinet_apg{{apg_opts.x,apg_opts.u}};
toc
Z.X=result{1,1};
Z.U=result{1,2};
%% check the equality constraints 

eq_const_box=control_demand_equation(Z_prcnd,Tree,apg_opts);
eq_const_dist=control_demand_equation(Z_prcnd_dist,Tree,apg_opts);
eq_const_yalmip=control_demand_equation(Z,Tree,apg_opts);

%SI=scenario_index(Tree);
figure
plot(max(abs(eq_const_dist)))
hold all;
plot(max(abs(eq_const_box)))
plot(max(abs(eq_const_yalmip)))
legend('APG-dist','APG-box','gurobi')

figure
for i=1:sys_actual.nx
    subplot(3,1,i)
    plot(Z_prcnd_dist.X(i,:));
    hold all;
    plot(Z_prcnd.X(i,:))
    plot(Z.X(i,:));
    plot(S.xmax(i)*ones(Nstage,1));
    plot(S.xmin(i)*ones(Nstage,1));
    legend('APG-dist','APG-box','gurobi')
end

figure
for i=1:sys_actual.nu
    subplot(3,2,i)
    
    plot(Z_prcnd_dist.U(i,:));
    hold all;
    plot(Z_prcnd.U(i,:))
    plot(Z.U(i,:));
    plot(3600*S.umin(i)*ones(Nstage,1));
    legend('APG-dist','APG-box','gurobi')
    %plot(3600*S.umax(i)*ones(Nstage,1));
end
%%
figure
plot(details_apg_prcnd_dist.jobj);
hold all;
plot(details_apg_prcnd.jobj)
legend('APG-dist','APG-box')