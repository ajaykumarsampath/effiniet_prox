%%
% This function eliminate the equality constraints and transform the water
% network into effinet model.

close all;
clear all
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

ops_sys.normalise=1;
%sys=sys_old;out.x1(())
sys=system_prox_formation(S,P,Tree,ops_sys);
[sys_mod,V_mod]=system_prox_formation_modified(S,P,Tree,ops_sys);

[sys_box,V]=system_prox_formation_box(S,P,Tree,ops_sys);

ops_sys.cell=0;
ops_sys.normalise=0;

forbes_opts.E=S.E;
forbes_opts.Ed=S.Ed;

apg_opts.E=S.E;
apg_opts.Ed=S.Ed;

apg_opts_mod.E=S.E(:,[1:3 5:6]);
apg_opts_mod.Ed=S.Ed;

[sys_null,V_null]=system_prox_formation_null(S,P,Tree,ops_sys);
%V.Wu=10*V.Wu;
%V.Wu(5,5)=1000;
%V.Wu(4,4)=1000;
%% 
par_sol_opt.demand=3600*DemandData(kk:kk+P.Hu,:);
prev_vhat=3600*sys_NN.L1*DemandData(kk-1,:)';

par_sol_opt.prev_vhat=prev_vhat;
%V.alpha=3600*(kron(ones(P.Hp,1),P.alpha1')+P.alpha2(kk:kk+P.Hu,:));
V.alpha=(kron(ones(P.Hp,1),P.alpha1')+P.alpha2(kk:kk+P.Hu,:));
current_state_opt=calcul_parti_soul(sys_NN,Tree,V,par_sol_opt);
current_state_opt.v=3600*[0.0656 0.00 0.0849 0.0934]';

opts.state=current_state_opt;

forbes_opts.x=0.1*(S.xmax-P.xs)+P.xs;
forbes_opts.primal_inf=0.01;
forbes_opts.dual_inf=0.01;
forbes_opts.E=S.E;
forbes_opts.Ed=S.Ed;
forbes_opts.w=current_state_opt.w;
forbes_opts.demand=current_state_opt.demand;
forbes_opts.uprev=sys.L*current_state_opt.v+prev_vhat;
%%
clear opt;
prob=formulate_forbes(sys_NN,V,Tree,forbes_opts);

prob.x0=forbes_opts.x;
opt_lbfgs.display = 1;
opt_lbfgs.method = 'lbfgs';
opt_lbfgs.variant='fast';
opt_lbfgs.memory = 10;
opt_lbfgs.tolOpt=1e-5;
%
tic;out_lbfg = forbes(prob,opt_lbfgs);toc

Z.X=zeros(sys.nx,size(Tree.stage,1));
Z.U=zeros(sys.nu,size(Tree.stage,1));
nz=sys.nx+sys.nu;
Z.X(:,1)=prob.x0;
%Z.X1(:,1)=prob.x0;
for i=1:size(Tree.stage,1)
    %Z.U(:,i)=out.z((i-1)*nz+1:(i-1)*nz+sys.nu,1);
    %Z.X(:,i+1)=out.z((i-1)*nz+sys.nu+1:i*nz,1);
    Z.U(:,i)=out_lbfg.x1((i-1)*nz+1:(i-1)*nz+sys.nu,1);
    Z.X(:,i+1)=out_lbfg.x1((i-1)*nz+sys.nu+1:i*nz,1);
end


opt.variant='fast';
opt.method='fbs';
opt.tolOpt=1e-5;
tic;out = forbes(prob,opt);toc

Z.X1=zeros(sys.nx,size(Tree.stage,1));
Z.U1=zeros(sys.nu,size(Tree.stage,1));
nz=sys.nx+sys.nu;
%Z.X(:,1)=prob.x0;
Z.X1(:,1)=prob.x0;
for i=1:size(Tree.stage,1)
    %Z.U(:,i)=out.z((i-1)*nz+1:(i-1)*nz+sys.nu,1);
    %Z.X(:,i+1)=out.z((i-1)*nz+sys.nu+1:i*nz,1);
    Z.U1(:,i)=out.x1((i-1)*nz+1:(i-1)*nz+sys.nu,1);
    Z.X1(:,i+1)=out.x1((i-1)*nz+sys.nu+1:i*nz,1);
end
%%
%sys=sys_NN;
Ptree_NN=factor_apg_effiniet(sys_NN,V,Tree);
opts_apg.state=current_state_opt;

opts_apg.x=forbes_opts.x;
opts_apg.primal_inf=0.01;
opts_apg.dual_inf=0.01;
opts_apg.E=S.E;
opts_apg.Ed=S.Ed;

opts_apg.steps=out.iterations;
opts_apg.lambda=8e-5;

[Z_NN,details_apg_NN]=APG_effiniet_2(sys_NN,Ptree_NN,Tree,V,opts_apg);

%%
opts_apg.state=current_state_opt;
opts_apg.x=forbes_opts.x;
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
%{
apg_opts_mod.u=sys_mod.L*opts_apg_mod.state.v+opts_apg_mod.state.prev_vhat;
apg_opts_mod.x=opts_apg.x;
apg_opts_mod.demand=opts_apg.state.demand;
effiniet_apg_mod=effiniet_yalmip(sys_actual_mod,Tree,V_mod,apg_opts_mod);
tic
[result_mod,error]=effiniet_apg_mod{{apg_opts_mod.x,apg_opts_mod.u}};
toc
Z_mod.eff_X=result_mod{1,1};
Z_mod.eff_U=result_mod{1,2};
%}
%%
figure
for i=1:sys.nx
subplot(3,1,i)
plot(Z.X(1,:));
hold all; 
plot(Z.X1(1,:))
plot(Z.eff_X(1,:));
plot(Z_NN.X(1,:))
%plot(Z_N.X(1,:))
end
%legend('fobes-fast','fbs-fast','gurobi','APG-NN','APG')
legend('fobes-fast','fbs-fast','gurobi','APG-NN')
figure
for i=1:sys.nu
subplot(3,2,i)
plot(Z.U(i,:));
hold all; 
plot(Z.U1(i,:))
plot(Z.eff_U(i,:));
plot(Z_NN.U(i,:));
%plot(Z_N.U(i,:));
end 
legend('fobes-fast','fbs-fast','gurobi','APG-NN')