close all;
clear all
clc
%%
kk=3870;
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
branching_factor = [2 1 1 1];
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
sys_actual=system_prox_modif_state_box(S,P,Tree,ops_sys);
ops_sys.cell=1;
ops_sys_poly=0;
sys_NN=system_prox_modif_state_box(S,P,Tree,ops_sys);

[sys_box,V]=system_prox_modif_state_box(S,P,Tree,ops_sys);
%V.Wu=V.Wu/3600;
apg_opts.E=S.E;
apg_opts.Ed=S.Ed;
apg_opts.exact_prcnd=0;
%% Particular solution calculation
par_sol_opt.demand=3600*DemandData(kk:kk+P.Hu,:);
prev_vhat=3600*sys_box.L1*DemandData(kk-1,:)';

par_sol_opt.prev_vhat=prev_vhat;
%V.alpha=3600*(kron(ones(P.Hp,1),P.alpha1')+P.alpha2(kk:kk+P.Hu,:));
V.alpha=(kron(ones(P.Hp,1),P.alpha1')+P.alpha2(kk:kk+P.Hu,:));
current_state_opt=calParSol_modified(sys_box,Tree,V,par_sol_opt);
current_state_opt.v=3600*[0.0656 0.00 0.0849 0.0934]';
%% Preconditioning 
opts_prcnd.dual_hessian=1;
opts_prcnd.prcnd=0;
if(Nscenarios==1)
    [DH_nml,dts_nml]=dual_hessian_calculate_box(sys_box,Tree,V,apg_opts);
    sys_prcnd_act=precondition_calculate_box(sys_box,DH_nml,Tree);
    [DH_nml_prcnd,dts_nml_prcnd]=dual_hessian_calculate_box(sys_prcnd_act,Tree,V,apg_opts);
    l=eig(DH_nml_prcnd);
end
[sys_prcnd,dts_prcnd]=prcnd_sys_modfy_box(sys_box,Tree,V,apg_opts);

%apg_opts.exact_prcnd=1;
%[sys_ext_prcnd,dts_exct_prcnd]=prcnd_sys_modfy_box(sys_box,Tree,V,apg_opts);

%% Factor step

Ptree_box=factor_apg_eff_modif_state(sys_box,V,Tree);
%Ptree=factor_apg_effinet(sys_box,V,Tree);
Ptree_prcnd=factor_apg_eff_modif_state(sys_prcnd,V,Tree);
if(apg_opts.exact_prcnd)
    Ptree_exct_prcnd=factor_apg_eff_modif_state(sys_ext_prcnd,V,Tree);
end

if(Nscenarios==1)
    Ptree_prcnd_act=factor_apg_eff_modif_state(sys_prcnd_act,V,Tree);
end

%% Solve steps 
opts_apg.state=current_state_opt;

opts_apg.x=0.7*(S.xmax-P.xs)+P.xs;
opts_apg.E=S.E;
opts_apg.Ed=S.Ed;

if(opts_prcnd.prcnd)
    opts_apg.steps=1000;
    opts_apg.lambda=0.08*1/dts_prcnd.norm_act;
    %opts_apg.lambda=2.1492e-05;
    [Z_box,details_apg_box]=APG_effinet_box_modified(sys_box,Ptree_box,Tree,V,opts_apg);
end

opts_apg.steps=200;
opts_apg.lambda=1.2*1/dts_prcnd.norm;
[Z_prcnd,details_apg_prcnd]=APG_effinet_box_modified(sys_prcnd,Ptree_prcnd,Tree,V,opts_apg);

if(apg_opts.exact_prcnd)
    opts_apg.lambda=1/dts_exct_prcnd.norm;
    [Z_exct_prcnd,details_apg_ext_prcnd]=APG_effinet_box_modified(sys_prcnd,Ptree_prcnd,Tree,V,opts_apg);
end

if(Nscenarios==1)
    opts_apg.steps=1000;
    opts_apg.lambda=1/dts_prcnd.norm;
    [Z_prcnd_act,details_apg_prcnd_act]=APG_effinet_box_modified(sys_prcnd_act,...
        Ptree_prcnd_act,Tree,V,opts_apg);
end
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
if(opts_prcnd.prcnd)
    eq_const=control_demand_equation(Z_box,Tree,apg_opts);
end
eq_const_pre=control_demand_equation(Z_prcnd,Tree,apg_opts);
eq_const_yalmip=control_demand_equation(Z,Tree,apg_opts);

%SI=scenario_index(Tree);
figure
if(opts_prcnd.prcnd)
    plot(max(abs(eq_const)))
end
hold all;
plot(max(abs(eq_const_pre)))
plot(max(abs(eq_const_yalmip)))
if(opts_prcnd.prcnd)
    legend('APG','APG-diag','gurobi')
else
    legend('APG-diag','gurobi')
end

figure
for i=1:sys_actual.nx
    subplot(3,1,i)
    if(opts_prcnd.prcnd)
        plot(Z_box.X(i,:));
    end     
    hold all;
    plot(Z_prcnd.X(i,:))
    if(apg_opts.exact_prcnd)
        plot(Z_exct_prcnd.X(i,:))
    end
    plot(Z.X(i,:));
    plot(S.xmax(i)*ones(Nstage,1));
    plot(S.xmin(i)*ones(Nstage,1));
    if(Nscenarios==1)
        plot(Z_prcnd_act.X(i,:))
    end
end
if(Nscenarios==1)
    legend('APG','APG-diag','gurobi','APG-EXACT')
else
    if(opts_prcnd.prcnd)
        legend('APG','APG-diag','gurobi')
    else
        if(apg_opts.exact_prcnd)
            legend('APG-diag','APG-exact','gurobi')
        else
            legend('APG-diag','gurobi')
        end
    end
end
figure
for i=1:sys_actual.nu
    subplot(3,2,i)
    %}
    if(opts_prcnd.prcnd)
        plot(Z_box.U(i,:));
    end
    hold all;
    plot(Z_prcnd.U(i,:))
    if(apg_opts.exact_prcnd)
        plot(Z_exct_prcnd.U(i,:))
    end
    plot(Z.U(i,:));
    %}
    plot(3600*S.umin(i)*ones(Nstage,1));
    %plot(3600*S.umax(i)*ones(Nstage,1));
    if(Nscenarios==1)
        plot(Z_prcnd_act.U(i,:));
    end
end
if(Nscenarios==1)
    legend('APG','APG-diag','gurobi','APG-EXACT')
else
    if(opts_prcnd.prcnd)
        legend('APG','APG-diag','gurobi')
    else
        if(apg_opts.exact_prcnd)
            legend('APG-diag','APG-exact','gurobi')
        else
            legend('APG-diag','gurobi')
        end
    end
end
%%
figure
if(opts_prcnd.prcnd)
    plot(details_apg_box.jobj)
end
hold all;
plot(details_apg_prcnd.jobj)
if(apg_opts.exact_prcnd)
    plot(details_apg_ext_prcnd.jobj);
end
if(Nscenarios==1)
    plot(details_apg_prcnd_act.jobj)
    legend('APG','APG-diag','APG-EXACT')
else
    if(apg_opts.exact_prcnd)
        legend('APG-diag','APG-exact','gurobi')
    else
        legend('APG-diag','gurobi')
    end
end