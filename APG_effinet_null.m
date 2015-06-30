function [ Z_null,details] = APG_effinet_null(sys_null,Ptree_null,Tree,V,opts)

% This function is implements the APG algorithm to solve the sc-mpc of
% effinet problem with soft constraints on state and hard constraits
% on input.

% INPUT----     sys   :  containts the system dynamics and proximal
%                        algorithm options.
%             Ptree   :  contains the factor step matrices
%                 V   :  cost function
%               opt   :  options for termination
%%  opts=opts_apg;
Ns=length(Tree.leaves);% total scenarios in the tree
Nd=length(Tree.stage);%toal nodes in the tree
% Initalizing the dual varibables
Y.y0=zeros(size(sys_null.F,1),Nd-Ns+1);
Y.yt0=zeros(2*sys_null.nx,Ns);
Y.y1=zeros(size(sys_null.F,1),Nd-Ns+1);
Y.yt1=zeros(2*sys_null.nx,Ns);
W.y=zeros(size(sys_null.F,1),Nd-Ns+1);
W.yt=zeros(2*sys_null.nx,Ns);

epsilon_prm=1;

prm_fes.soft=zeros(sys_null.nx,Nd);
prm_fes.hard=zeros(sys_null.nu,Nd+1-Ns);

theta=[1 1]';

details.term_crit=zeros(1,4);
dual_grad=prm_fes;

opts_prox.lambda=opts.lambda;
opts_prox.xmin=kron(ones(Nd,1),sys_null.xmin);
opts_prox.xs=kron(ones(Nd,1),sys_null.xs);
opts_prox.xmax=kron(ones(Nd,1),sys_null.xmax);

opts_prox.umin=kron(ones(1,Nd-Ns+1),sys_null.umin)-opts.state.vhat;
opts_prox.umax=kron(ones(1,Nd-Ns+1),sys_null.umax)-opts.state.vhat;


%opts_prox.umin=kron(ones(1,Nd-Ns+1),sys_null.umin);
%opts_prox.umax=kron(ones(1,Nd-Ns+1),sys_null.umax);
opts_prox.nx=sys_null.nx;
opts_prox.nu=size(sys_null.umax,1);
opts_prox.gamma_c=sys_null.gamma_max;
opts_prox.gamma_s=sys_null.gamma_xs;

uprev=sys_null.L*opts.state.v+opts.state.prev_vhat;
%{
%g1=[sys.xmax;-sys.xmin;sys.umax;-sys.umin];
grad_opts.u=uprev;
grad_opts.x=opts.x;
grad_opts.E=opts.E;
grad_opts.Ed=opts.Ed;
grad_opts.demand=opts.state.demand;
%V.Wu=100*eye(sys.nu);
dual_grad_yalmip=dual_gradient_yalmip(sys,Tree,V,grad_opts);


%}
%%
%tic
%opts.steps=200;
%opts.lambda=18e-6;
j=1;
jobj=zeros(1,opts.steps-1);
while(j<opts.steps)
    % Step 1: accelerated step
    W.y=Y.y1+theta(2)*(1/theta(1)-1)*(Y.y1-Y.y0);
    
    W.yt=Y.yt1+theta(2)*(1/theta(1)-1)*(Y.yt1-Y.yt0);
    
    
    %step 2: argmin of the lagrangian using dynamic programming
    %}
    Z_null=solve_apg_effinet_null(sys_null,Tree,Ptree_null,W,opts.x,opts.state);
    
    for i=1:Nd-Ns+1
        if(i==1)
            jobj(j)=jobj(j)+V.alpha(1,:)*Z_null.U(:,1)+(Z_null.U(:,1)-uprev)'*V.Wu*(Z_null.U(:,1)-uprev)...
                ;%+W.y(:,1)'*sys_null.G*Z_null.v(:,i);
        else
            stage=Tree.stage(i-1)+1;
            jobj(j)=jobj(j)+Tree.prob(i-1,1)*(V.alpha(stage,:)*Z_null.U(:,i)+(Z_null.U(:,i)-...
                Z_null.U(:,Tree.ancestor(i-1)+1))'*V.Wu*...
                (Z_null.U(:,i)-Z_null.U(:,Tree.ancestor(i-1)+1)))...
                ;%+W.y(:,i)'*(sys_null.F*Z_null.X(:,i)+sys_null.G*Z_null.v(:,i));
        end
    end
    
    for i=1:Ns
        %jobj(j)=jobj(j)+W.yt(:,i)'*...
         %   (sys_null.F(1:2*sys_null.nx,:)*Z_null.X(:,Tree.leaves(i)+1));
    end
    %{
    [results,error]=dual_grad_yalmip{{grad_opts.x,grad_opts.u,W.y,W.yt}};
    Z.X=results{1,1};
    Z.U=results{1,2};
    jobj(j)=results{1,3};
    %}
    %step 3: proximal with respect to g
    t=proximal_effinet_null(Z_null,W,sys_null,opts_prox);
    
    %step 4: update the dual vector
    Y.y0=Y.y1;
    Y.yt0=Y.yt1;
    
    %dual_grad(1:sys.nu,:)=Z.U;
    %dual_grad(sys.nu+1:end,1:Nd-1)=Z.X(:,2:Nd);
    %prm_fes.soft(1:sys.nx,1:Nd+Ns-1)=(Z.X(:,2:Nd+Ns)-t.x(:,1:Nd+1-Ns));
    %prm_fes.hard=(Z.U-t.u);
    
    %
    Y.y1(1:2*sys_null.nx,2:Nd)=W.y(1:2*sys_null.nx,2:Nd)+opts.lambda*...
        (sys_null.F(1:2*sys_null.nx,:)*Z_null.X(:,2:Nd)-t.x(:,1:Nd-Ns));
    Y.y1(2*sys_null.nx+1:end,:)=W.y(2*sys_null.nx+1:end,:)+opts.lambda*(sys_null.G(2*sys_null.nx+1:end,:)*Z_null.v-t.u);
    Y.yt1=W.yt+opts.lambda*(sys_null.F(1:2*sys_null.nx,:)*Z_null.X(:,Nd+1:Nd+Ns)-t.x(:,Nd-Ns+1:Nd+Ns-1));
    %}
    %{
    Y.y1(1:2*sys.nx,2:Nd)=max(0,W.y(1:2*sys.nx,2:Nd)+opts.lambda*...
        (sys.F(1:2*sys.nx,:)*Z.X(:,2:Nd)-kron(ones(1,Nd-1),g1(1:2*sys.nx,1))));
    Y.y1(2*sys.nx+1:end,:)=max(0,W.y(2*sys.nx+1:end,:)+opts.lambda*...
        (sys.G(2*sys.nx+1:end,:)*Z.U-kron(ones(1,Nd),g1(2*sys.nx+1:end,1))));
    Y.yt1=max(0,W.yt+opts.lambda*...
        (sys.F(1:2*sys.nx,:)*Z.X(:,Nd+1:Nd+Ns)-kron(ones(1,Ns),g1(1:2*sys.nx,1))));
    %}
    iter=j;
    details.prm_cst(iter)=0;%primal cost;
    details.dual_cst(iter)=0;% dual cost;
    
    theta(1)=theta(2);
    theta(2)=(sqrt(theta(1)^4+4*theta(1)^2)-theta(1)^2)/2;
    j=j+1;
    if(mod(j,400)==0)
        %theta=[1 1]';
    end
    %%
end
figure
plot(jobj);
details.gpad_solve=toc;
details.W=W;
details.jobj=jobj;
%details.epsilon_prm_avg= epsilon_prm_avg;
%details.epsilon_prm=epsilon_prm;
end





