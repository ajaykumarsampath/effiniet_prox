function [ Z,details] = APG_effinet_dist_box(sys_box,Ptree_box,Tree,V,opts)

% This function is implements the APG algorithm to solve the sc-mpc of
% effinet problem with soft constraints on state and hard constraits
% on input.

% INPUT----     sys   :  containts the system dynamics and proximal
%                        algorithm options.
%             Ptree   :  contains the factor step matrices
%                 V   :  cost function
%               opt   : 
%                       Termination options: maximum number of steps etc   
%                       Type of constraints: soft or hard constaints
%                       Intial condition and preious control input.
%% opts=opts_apg;
Ns=length(Tree.leaves);% total scenarios in the tree
Nd=length(Tree.stage);%toal nodes in the tree
% Initalizing the dual varibables

ny_x=size(sys_box.F{1},1);
ny_u=size(sys_box.G{1},1);

Y.x0=zeros(ny_x,Nd);
Y.u0=zeros(ny_u,Nd);
Y.x1=zeros(ny_x,Nd);
Y.u1=zeros(ny_u,Nd);
W.x=zeros(ny_x,Nd);
W.u=zeros(ny_u,Nd);

epsilon_prm=1;

theta=[1 1]';

details.term_crit=zeros(1,4);

opts_prox.lambda=opts.lambda;
opts_prox.xmin=reshape(sys_box.xmin,sys_box.nx,Nd);
opts_prox.xs=reshape(sys_box.xs,sys_box.nx,Nd);
opts_prox.xmax=reshape(sys_box.xmax,sys_box.nx,Nd);

opts_prox.umin=reshape(sys_box.umin,sys_box.nu,Nd);
opts_prox.umax=reshape(sys_box.umax,sys_box.nu,Nd);
opts_prox.nx=sys_box.nx;
opts_prox.nu=sys_box.nu;
opts_prox.gamma_xbox=sys_box.gamma_xbox/opts.lambda;
opts_prox.gamma_xs=sys_box.gamma_xs/opts.lambda;
opts_prox.iter=200;
opts_prox.constraints=opts.constraints;

uprev=sys_box.L*opts.state.v+opts.state.prev_vhat;

prm_feas.x=zeros(ny_x,Nd);
prm_feas.u=zeros(ny_u,Nd);

%%
%tic
j=1;
jobj=zeros(1,opts.steps-1);
while(j<opts.steps)
    
    %% Step 1: accelerated step
    W.x=Y.x1+theta(2)*(1/theta(1)-1)*(Y.x1-Y.x0);
    W.u=Y.u1+theta(2)*(1/theta(1)-1)*(Y.u1-Y.u0);
    
    %step 2: argmin of the lagrangian using dynamic programming
    %}
    Z=solve_apg_eff_dist_state(sys_box,Tree,Ptree_box,W,opts.x,opts.state);
    
    %
    for i=1:Nd-Ns+1
        if(i==1)
            jobj(j)=jobj(j)+V.alpha(1,:)*Z.U(:,1)+(Z.U(:,1)-uprev)'*V.Wu*(Z.U(:,1)-uprev)...
                +W.x(:,1)'*sys_box.F{i}*Z.X(:,2)+W.u(:,1)'*sys_box.G{1}*Z.U(:,1);
        else
            stage=Tree.stage(i-1)+1;
            nchild=Tree.children{i-1};
            for l=1:length(nchild)
                jobj(j)=jobj(j)+Tree.prob(nchild(l),1)*(V.alpha(stage+1,:)*Z.U(:,nchild(l))+...
                    (Z.U(:,nchild(l))-Z.U(:,Tree.ancestor(nchild(l))))'*V.Wu*...
                    (Z.U(:,nchild(l))-Z.U(:,Tree.ancestor(nchild(l)))))...
                    +W.x(:,nchild(l))'*sys_box.F{nchild(l)}*Z.X(:,nchild(l)+1)+...
                    W.u(:,nchild(l))'*sys_box.G{nchild(l)}*Z.U(:,nchild(l));
            end
        end
    end
    
    %step 3: proximal with respect to g
    if(strcmp(opts.distance,'yes'))
        t=prox_effinet_dist_box(Z,W,sys_box,Tree,opts_prox);
    else
        t=prox_effinet_modif_state_box(Z,W,sys_box,Tree,opts_prox);
    end
    
    
    %step 4: update the dual vector
    Y.x0=Y.x1;
    Y.u0=Y.u1;
    
    if(sys_box.cell)
        for i=1:Nd
            Y.x1(:,i)=W.x(:,i)+opts.lambda*(sys_box.F{i,1}*Z.X(:,i+1)-t.x(:,i));
            Y.u1(:,i)=W.u(:,i)+opts.lambda*(sys_box.G{i,1}*Z.U(:,i)-t.u(:,i));
            prm_feas.x(:,i)=sys_box.F{i,1}*Z.X(:,i+1)-t.x(:,i);
            prm_feas.u(:,i)=sys_box.G{i,1}*Z.U(:,i)-t.u(:,i);
        end
    else
        Y.x1=W.x+opts.lambda*(sys_box.F*Z.X(:,2:Nd)-t.x);
        Y.u1=W.u+opts.lambda*(sys_box.G(sys_box.nx+1:end,:)*Z.U-t.u);
    end
    
    iter=j;
    details.prm_cst(iter)=0;%primal cost;
    details.dual_cst(iter)=0;% dual cost;
    
    theta(1)=theta(2);
    theta(2)=(sqrt(theta(1)^4+4*theta(1)^2)-theta(1)^2)/2;
    %if(mod(j,300)==1)
        %if(jobj(j)<jobj(j-1))
            %theta=[1 1]';
        %end
    %end
     j=j+1;
end
%figure
%plot(jobj);
details.gpad_solve=toc;
details.W=W;
details.jobj=jobj;
details.prm_feas=prm_feas;
%details.epsilon_prm_avg= epsilon_prm_avg;
%details.epsilon_prm=epsilon_prm;
end







