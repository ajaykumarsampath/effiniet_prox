close all;
clear all
clc
%%
kk=4875;
load('dwn');
%load('dwn_big');
P.Hp=24;
P.Hu=23;
%S.xmin=0.2*S.xmax;
%P.xs=0.7*S.xmax;
%P.beta=0.1;
Dd = DemandData(:,1); % My data
% Define data
tree_ops.Tmax = 800; % time steps for which data are available
tree_ops.nScen = 300; % no. of scenarios
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
branching_factor =  [3 2 1];
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
ops_sys.gamma_xs=1e4;
ops_sys.gamma_xbox=1e9;
ops_sys.normalise=0;
ops_sys.cell=0;
sys_actual=system_prox_dist_state_box(S,P,Tree,ops_sys);
ops_sys.cell=1;
ops_sys_poly=0;

[sys_box_dist,V]=system_prox_dist_state_box(S,P,Tree,ops_sys);

yamip_opts.E=S.E;
yamip_opts.Ed=S.Ed;
yamip_opts.exact_prcnd=0;
%% Particular solution calculation
par_sol_opt.demand=3600*DemandData(kk:kk+P.Hu,:);
prev_vhat=3600*sys_box_dist.L1*DemandData(kk-1,:)';

par_sol_opt.prev_vhat=prev_vhat;
%V.alpha=3600*(kron(ones(P.Hp,1),P.alpha1')+P.alpha2(kk:kk+P.Hu,:));
V.alpha=(kron(ones(P.Hp,1),P.alpha1')+P.alpha2(kk:kk+P.Hu,:));
V.Qe=100*eye(sys_actual.nx);
V.Qs=100*eye(sys_actual.nx);
current_state_opt=calParSol_modified(sys_box_dist,Tree,V,par_sol_opt);
%current_state_opt.v=3600*[0.0656 0.00 0.0849 0.0934]';
current_state_opt.v=3600*rand(size(sys_actual.L,2),1);
%% Preconditioning opts_apg.constraints='soft';
opts_prcnd.dual_hessian=1;
opts_prcnd.prcnd=0;

[sys_prcnd_dist,dts_prcnd_dist]=prcnd_sys_dist_box(sys_box_dist,Tree,V,yamip_opts);

% Factor step

Ptree_prcnd_dist=factor_apg_eff_dist_state(sys_prcnd_dist,V,Tree);
Ptree_prcnd=factor_apg_eff_modif_state(sys_prcnd_dist,V,Tree);

%% Solve steps  
% solve step options 
opts_apg.state=current_state_opt;

%opts_apg.x=0.1*S.xmin;
opts_apg.x=0.5*(S.xmax-P.xs)+P.xs;
opts_apg.E=S.E;
opts_apg.Ed=S.Ed;
opts_apg.constraints='soft';
opts_apg.distance='yes';

opts_apg.steps=300;
opts_apg.lambda=1.3*1/dts_prcnd_dist.norm;

tic
[Z_prcnd_dist,details_apg_prcnd_dist]=APG_effinet_dist_box(sys_prcnd_dist,Ptree_prcnd_dist,Tree,V,opts_apg);
toc
%{
tic
[Z_prcnd,details_apg_prcnd]=APG_effinet_dist_box(sys_prcnd_dist,Ptree_prcnd,Tree,V,opts_apg);
toc
%}
%% Gurobi Algorithm 
yamip_opts.u=sys_actual.L*opts_apg.state.v+opts_apg.state.prev_vhat;
yamip_opts.x=opts_apg.x;
yamip_opts.E=S.E;
yamip_opts.Ed=S.Ed;
yamip_opts.demand=opts_apg.state.demand;

%effinet_apg=effinet_modif_state_yalmip(sys_actual,Tree,V,yamip_opts);
effinet_apg=effinet_modifi_soft_yalmip(sys_actual,Tree,V,yamip_opts);
tic
[result,error]=effinet_apg{{yamip_opts.x,yamip_opts.u}};
toc
Z.X=result{1,1};
Z.U=result{1,2};

%{
effinet_apg_dist=effinet_modifi_distance_yalmip(sys_actual,Tree,V,yamip_opts);
tic
[result,error]=effinet_apg{{yamip_opts.x,yamip_opts.u}};
toc
Z_yalmip_dist.X=result{1,1};
Z_yalmip_dist.U=result{1,2};
%}
%% check the equality constraints 

eq_const_dist=control_demand_equation(Z_prcnd_dist,Tree,yamip_opts);
eq_const_yalmip=control_demand_equation(Z,Tree,yamip_opts);

%SI=scenario_index(Tree);
figure
plot(max(abs(eq_const_dist)))
hold all;
plot(max(abs(eq_const_yalmip)))
legend('APG-dist','gurobi')

figure
for i=1:sys_actual.nx
    subplot(sys_actual.nx,1,i)
    plot(Z_prcnd_dist.X(i,:));
    hold all;
    plot(Z.X(i,:));
    %plot(Z_prcnd.X(i,:));
    %plot(Z_yalmip_dist.X(i,:));
    plot(sys_actual.xmax(i)*ones(Nstage,1));
    plot(sys_actual.xs(i)*ones(Nstage,1));
    plot(sys_actual.xmin(i)*ones(Nstage,1));
    %legend('APG-dist','gurobi','gurobi-distance','max-vol','safety','min-vol')
    legend('APG-dist','gurobi','max-vol','safety','min-vol')
    xlabel('nodes across the tree')
    ylabel('volume in the tanks')
end

figure
for i=1:sys_actual.nu
    subplot(3,2,i)
    plot(Z_prcnd_dist.U(i,:));
    hold all;
    plot(Z.U(i,:));
    %plot(Z_prcnd.U(i,:));
    %plot(Z_yalmip_dist.U(i,:));
    plot(3600*S.umin(i)*ones(Nstage,1));
    %plot(3600*S.umax(i)*ones(Nstage,1));
    %legend('APG-dist','gurobi','gurobi-distance','umin')
    legend('APG-dist','gurobi','umin')
    xlabel('nodes across the tree')
    ylabel('control action in the tanks')
end
%%
%figure
%plot(details_apg_prcnd_dist.jobj);
%legend('APG-dist')

jobj=[0 0];
Z_prcnd_dist.Xs=kron(ones(1,Nstage),sys_actual.xs(1:sys_actual.nx))-Z_prcnd_dist.X;
Z_prcnd_dist.Xs=max(zeros(sys_actual.nx,Nstage),Z_prcnd_dist.Xs);

%APG
for i=1:Nstage-Nscenarios
    if(i==1)
        jobj(1)=jobj(1)+V.alpha(1,:)*Z_prcnd_dist.U(:,1)+(Z_prcnd_dist.U(:,1)-yamip_opts.u)'...
            *V.Wu*(Z_prcnd_dist.U(:,1)-yamip_opts.u)+...
             Z_prcnd_dist.Xs(:,2)'*V.Qs*Z_prcnd_dist.Xs(:,2);
    else
        stage=Tree.stage(i-1)+1;
        nchild=Tree.children{i-1};
        for l=1:length(nchild)
            jobj(1)=jobj(1)+Tree.prob(nchild(l),1)*(V.alpha(stage+1,:)*Z_prcnd_dist.U(:,nchild(l))+...
                (Z_prcnd_dist.U(:,nchild(l))-Z_prcnd_dist.U(:,Tree.ancestor(nchild(l))))'*V.Wu*...
                (Z_prcnd_dist.U(:,nchild(l))-Z_prcnd_dist.U(:,Tree.ancestor(nchild(l)))))+...
                 Z_prcnd_dist.Xs(:,nchild(l)+1)'*V.Qs*Z_prcnd_dist.Xs(:,nchild(l)+1);
        end
    end
end

Z.Xs=kron(ones(1,Nstage),sys_actual.xs(1:sys_actual.nx))-Z.X;
Z.Xs=max(zeros(sys_actual.nx,Nstage),Z.Xs);

%yalmip
for i=1:Nstage-Nscenarios
    if(i==1)
        jobj(2)=jobj(2)+V.alpha(1,:)*Z.U(:,1)+(Z.U(:,1)-yamip_opts.u)'...
            *V.Wu*(Z.U(:,1)-yamip_opts.u)+Z.Xs(:,2)'*V.Qs*Z.Xs(:,2);
    else
        stage=Tree.stage(i-1)+1;
        nchild=Tree.children{i-1};
        for l=1:length(nchild)
            jobj(2)=jobj(2)+Tree.prob(nchild(l),1)*(V.alpha(stage+1,:)*Z.U(:,nchild(l))+...
                (Z.U(:,nchild(l))-Z.U(:,Tree.ancestor(nchild(l))))'*V.Wu*...
                (Z.U(:,nchild(l))-Z.U(:,Tree.ancestor(nchild(l)))))+...
                 Z.Xs(:,nchild(l)+1)'*V.Qs*Z.Xs(:,nchild(l)+1);
        end
    end
end

display('difference of operation cost between gurobi and yalmip', num2str(jobj(1)-jobj(2)))
