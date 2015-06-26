function [ Z,details] = APG_effinet_box(sys_box,Ptree_box,Tree,V,opts)

% This function is implements the APG algorithm to solve the sc-mpc of
% effinet problem with soft constraints on state and hard constraits
% on input.

% INPUT----     sys   :  containts the system dynamics and proximal
%                        algorithm options.
%             Ptree   :  contains the factor step matrices
%                 V   :  cost function
%               opt   :  options for termination
%% opts=opts_apg;
Ns=length(Tree.leaves);% total scenarios in the tree
Nd=length(Tree.stage);%toal nodes in the tree
% Initalizing the dual varibables
Y.y0=zeros(size(sys_box.F{1},1),Nd-Ns+1);
Y.yt0=zeros(sys_box.nx,Ns);
Y.y1=zeros(size(sys_box.F{1},1),Nd-Ns+1);
Y.yt1=zeros(sys_box.nx,Ns);
W.y=zeros(size(sys_box.F{1},1),Nd-Ns+1);
W.yt=zeros(sys_box.nx,Ns);

epsilon_prm=1;

prm_fes.soft=zeros(sys_box.nx,Nd);
prm_fes.hard=zeros(sys_box.nu,Nd+1-Ns);

theta=[1 1]';

details.term_crit=zeros(1,4);
dual_grad=prm_fes;

opts_prox.lambda=opts.lambda;
opts_prox.xmin=reshape(sys_box.xmin(sys_box.nx+1:end,1),sys_box.nx,Nd);
opts_prox.xs=reshape(sys_box.xs(sys_box.nx+1:end,1),sys_box.nx,Nd);
opts_prox.xmax=reshape(sys_box.xmax(sys_box.nx+1:end,1),sys_box.nx,Nd);

opts_prox.umin=reshape(sys_box.umin,sys_box.nu,Nd-Ns+1);
opts_prox.umax=reshape(sys_box.umax,sys_box.nu,Nd-Ns+1);
opts_prox.nx=sys_box.nx;
opts_prox.nu=sys_box.nu;
opts_prox.gamma_c=sys_box.gamma_max;
opts_prox.gamma_s=sys_box.gamma_xs;

uprev=sys_box.L*opts.state.v+opts.state.prev_vhat;

prm_feas.x=zeros(sys_box.nx,Nd);
prm_feas.u=zeros(sys_box.nu,Nd-Ns+1);

%g1=[sys.xmax(1:sys.nx);-sys.xmin(1:sys.nx);sys.umax(1:sys.nu);sys.umin(1:sys.nu)];
%{
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
j=1;
jobj=zeros(1,opts.steps-1);
while(j<opts.steps)
    
    %% Step 1: accelerated step
    W.y=Y.y1+theta(2)*(1/theta(1)-1)*(Y.y1-Y.y0);
    
    W.yt=Y.yt1+theta(2)*(1/theta(1)-1)*(Y.yt1-Y.yt0);
    
    %step 2: argmin of the lagrangian using dynamic programming
    %}
    Z=solve_apg_effinet(sys_box,Tree,Ptree_box,W,opts.x,opts.state);
    
    for i=1:Nd-Ns+1
        if(i==1)
            jobj(j)=jobj(j)+V.alpha(1,:)*Z.U(:,1)+(Z.U(:,1)-uprev)'*V.Wu*(Z.U(:,1)-uprev)...
                +W.y(:,1)'*sys_box.G{i}*Z.U(:,i);
        else
            stage=Tree.stage(i-1)+1;
            jobj(j)=jobj(j)+Tree.prob(i-1,1)*(V.alpha(stage+1,:)*Z.U(:,i)+(Z.U(:,i)-...
                Z.U(:,Tree.ancestor(i-1)+1))'*V.Wu*...
                (Z.U(:,i)-Z.U(:,Tree.ancestor(i-1)+1)))...
                +W.y(:,i)'*(sys_box.F{i}*Z.X(:,i)+sys_box.G{i}*Z.U(:,i));
        end
    end
    
    for i=1:Ns
        jobj(j)=jobj(j)+W.yt(:,i)'*...
            (sys_box.F{Tree.leaves(i)}(1:sys_box.nx,:)*Z.X(:,Tree.leaves(i)+1));
    end
    %{
    [results,error]=dual_grad_yalmip{{grad_opts.x,grad_opts.u,W.y,W.yt}};
    Z.X=results{1,1};
    Z.U=results{1,2};
    jobj(j)=results{1,3};
    %}
    %step 3: proximal with respect to g
    t=proximal_effinet_box(Z,W,sys_box,opts_prox);
    
    %step 4: update the dual vector
    Y.y0=Y.y1;
    Y.yt0=Y.yt1;
    
    
    if(sys_box.cell)      
        for i=1:Nd-Ns+2
            if(i==1)
                prm_feas.u(:,1)=(sys_box.G{i,1}(sys_box.nx+1:end,:)*Z.U(:,i)-t.u(:,i));
                Y.y1(sys_box.nx+1:end,i)=W.y(sys_box.nx+1:end,i)+opts.lambda*...
                   (sys_box.G{i,1}(sys_box.nx+1:end,:)*Z.U(:,i)-t.u(:,i));
            elseif(i>Nd-Ns+1)
                for k=1:Ns
                    Y.yt1(:,k)=W.yt(:,k)+opts.lambda*(sys_box.F{i+k-1,1}(1:sys_box.nx,:)*...
                        Z.X(:,Nd-Ns+1+k)-t.x(:,Nd-Ns+k));
                    prm_feas.x(:,i+k-1)=sys_box.F{i+k-1,1}(1:sys_box.nx,:)*...
                        Z.X(:,Nd-Ns+1+k)-t.x(:,Nd-Ns+k);
                end
                %Y.yt1=W.yt+opts.lambda*(sys_box.F{i,1}(1:sys_box.nx,:)*Z.X(:,Nd+1:Nd+Ns)-t.x(:,Nd-Ns+1:Nd+Ns-1));
                %prm_feas.x(:,i-1)=sys_box.F{i,1}(1:sys_box.nx,:)*Z.X(:,Nd+1:Nd+Ns)-t.x(:,Nd-Ns+1:Nd+Ns-1);
            else
                Y.y1(1:sys_box.nx,i)=W.y(1:sys_box.nx,i)+opts.lambda*...
                    (sys_box.F{i,1}(1:sys_box.nx,:)*Z.X(:,i)-t.x(:,i-1));
                Y.y1(sys_box.nx+1:end,i)=W.y(sys_box.nx+1:end,i)...
                    +opts.lambda*(sys_box.G{i,1}(sys_box.nx+1:end,:)*Z.U(:,i)-t.u(:,i));
                prm_feas.x(:,i-1)=sys_box.F{i,1}(1:sys_box.nx,:)*Z.X(:,i)-t.x(:,i-1);
                prm_feas.u(:,i)=sys_box.G{i,1}(sys_box.nx+1:end,:)*Z.U(:,i)-t.u(:,i);
            end
        end
    else
        Y.y1(1:sys_box.nx,2:Nd)=W.y(1:sys_box.nx,2:Nd)+opts.lambda*...
            (sys_box.F(1:sys_box.nx,:)*Z.X(:,2:Nd)-t.x(:,1:Nd-Ns));
        Y.y1(sys_box.nx+1:end,:)=W.y(sys_box.nx+1:end,:)+opts.lambda*...
            (sys_box.G(sys_box.nx+1:end,:)*Z.U-t.u);
        Y.yt1=W.yt+opts.lambda*(sys_box.F(1:sys_box.nx,:)*...
            Z.X(:,Nd+1:Nd+Ns)-t.x(:,Nd-Ns+1:Nd+Ns-1));
    end

    %{
    %
    if(sys.cell)
        for i=1:Nd+1
            if(i==1)
                Y.y1(2*sys.nx+1:end,1)=max(0,W.y(2*sys.nx+1:end,1)+opts.lambda*...
                    (sys.G{1}(2*sys.nx+1:end,:)*Z.U(:,i)-g1(2*sys.nx+1:end,1)));
            elseif(i==Nd+1)
                Y.yt1=max(0,W.yt+opts.lambda*...
                    (sys.F{i}(1:2*sys.nx,:)*Z.X(:,Nd+1:Nd+Ns)-g1(1:2*sys.nx,1)));
            else
                Y.y1(1:2*sys.nx,i)=max(0,W.y(1:2*sys.nx,i)+opts.lambda*...
                    (sys.F{i}(1:2*sys.nx,:)*Z.X(:,i)-g1(1:2*sys.nx,1)));
                Y.y1(2*sys.nx+1:end,i)=max(0,W.y(2*sys.nx+1:end,i)+opts.lambda*...
                    (sys.G{i}(2*sys.nx+1:end,:)*Z.U(:,i)-g1(2*sys.nx+1:end,1)));
            end
        end
    else
        Y.y1(1:2*sys.nx,2:Nd)=max(0,W.y(1:2*sys.nx,2:Nd)+opts.lambda*...
            (sys.F(1:2*sys.nx,:)*Z.X(:,2:Nd)-kron(ones(1,Nd-1),g1(1:2*sys.nx,1))));
        Y.y1(2*sys.nx+1:end,:)=max(0,W.y(2*sys.nx+1:end,:)+opts.lambda*...
            (sys.G(2*sys.nx+1:end,:)*Z.U-kron(ones(1,Nd),g1(2*sys.nx+1:end,1))));
        Y.yt1=max(0,W.yt+opts.lambda*...
            (sys.F(1:2*sys.nx,:)*Z.X(:,Nd+1:Nd+Ns)-kron(ones(1,Ns),g1(1:2*sys.nx,1))));
    end

    %}
    iter=j;
    details.prm_cst(iter)=0;%primal cost;
    details.dual_cst(iter)=0;% dual cost;
    
    theta(1)=theta(2);
    theta(2)=(sqrt(theta(1)^4+4*theta(1)^2)-theta(1)^2)/2;
    %if(mod(j,300)==1)
        %if(jobj(j)<jobj(j-1))
            theta=[1 1]';
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





