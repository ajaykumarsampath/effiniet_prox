function [ Z,details] = APG_effiniet_2(sys,Ptree,Tree,V,opts)

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
Y.y0=zeros(size(sys.F{1},1),Nd-Ns+1);
Y.yt0=zeros(2*sys.nx,Ns);
Y.y1=zeros(size(sys.F{1},1),Nd-Ns+1);
Y.yt1=zeros(2*sys.nx,Ns);
W.y=zeros(size(sys.F{1},1),Nd-Ns+1);
W.yt=zeros(2*sys.nx,Ns);

epsilon_prm=1;

prm_fes.soft=zeros(sys.nx,Nd);
prm_fes.hard=zeros(sys.nu,Nd+1-Ns);

theta=[1 1]';

details.term_crit=zeros(1,4);
dual_grad=prm_fes;

opts_prox.lambda=opts.lambda;
opts_prox.xmin=reshape(-sys.xmin(sys.nx+1:end,1),sys.nx,Nd);
opts_prox.xs=reshape(sys.xs(sys.nx+1:end,1),sys.nx,Nd);
opts_prox.xmax=reshape(sys.xmax(sys.nx+1:end,1),sys.nx,Nd);

opts_prox.umin=reshape(-sys.umin,sys.nu,Nd-Ns+1);
opts_prox.umax=reshape(sys.umax,sys.nu,Nd-Ns+1);
opts_prox.nx=sys.nx;
opts_prox.nu=sys.nu;
opts_prox.gamma_c=sys.gamma_max;
opts_prox.gamma_s=sys.gamma_xs;

uprev=sys.L*opts.state.v+opts.state.prev_vhat;
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
%opts.steps=8000;
%opts.lambda=1e-3;
j=1;
jobj=zeros(1,opts.steps-1);
while(j<opts.steps)
    
    %% Step 1: accelerated step
    W.y=Y.y1+theta(2)*(1/theta(1)-1)*(Y.y1-Y.y0);
    
    W.yt=Y.yt1+theta(2)*(1/theta(1)-1)*(Y.yt1-Y.yt0);
    
    
    %step 2: argmin of the lagrangian using dynamic programming
    %}
    [Z,details_apg]=solve_apg_effiniet(sys,Tree,Ptree,W,opts.x,opts.state);
    
    %
    %results=dual_grad_yalmip{{grad_opts.x,grad_opts.u,W.y,W.yt}};
    %Z.X=results{1,1};
    %Z.U=results{1,2};
    %}
    for i=1:Nd-Ns+1
        if(i==1)
            jobj(j)=jobj(j)+V.alpha(1,:)*Z.U(:,1)+(Z.U(:,1)-uprev)'*V.Wu*(Z.U(:,1)-uprev)...
                +W.y(:,1)'*sys.G{i}*Z.U(:,i);
        else
            stage=Tree.stage(i-1)+1;
            jobj(j)=jobj(j)+Tree.prob(i-1,1)*(V.alpha(stage,:)*Z.U(:,i)+(Z.U(:,i)-...
                Z.U(:,Tree.ancestor(i-1)+1))'*V.Wu*...
                (Z.U(:,i)-Z.U(:,Tree.ancestor(i-1)+1)))...
                +W.y(:,i)'*(sys.F{i}*Z.X(:,i)+sys.G{i}*Z.U(:,i));
        end
    end
    
    for i=1:Ns
        jobj(j)=jobj(j)+W.yt(:,i)'*...
            (sys.F{Tree.leaves(i)}(1:2*sys.nx,:)*Z.X(:,Tree.leaves(i)+1));
    end
    %{
    [results,error]=dual_grad_yalmip{{grad_opts.x,grad_opts.u,W.y,W.yt}};
    Z.X=results{1,1};
    Z.U=results{1,2};
    jobj(j)=results{1,3};
    %}
    %step 3: proximal with respect to g
    t=proximal_effinet_2(Z,W,sys,opts_prox);
    
    %step 4: update the dual vector
    Y.y0=Y.y1;
    Y.yt0=Y.yt1;
    
    %dual_grad(1:sys.nu,:)=Z.U;
    %dual_grad(sys.nu+1:end,1:Nd-1)=Z.X(:,2:Nd);
    %prm_fes.soft(1:sys.nx,1:Nd+Ns-1)=(Z.X(:,2:Nd+Ns)-t.x(:,1:Nd+1-Ns));
    %prm_fes.hard=(Z.U-t.u);
    %}
    if(sys.cell)      
        for i=1:Nd+1
            if(i==1)
               Y.y1(2*sys.nx+1:end,i)=W.y(2*sys.nx+1:end,i)+opts.lambda*...
                   (sys.G{i,1}(2*sys.nx+1:end,:)*Z.U(:,i)-t.u(:,i));
            elseif(i==Nd+1)
                Y.yt1=W.yt+opts.lambda*(sys.F{i,1}(1:2*sys.nx,:)*Z.X(:,Nd+1:Nd+Ns)-t.x(:,Nd-Ns+1:Nd+Ns-1));
            else
                Y.y1(1:2*sys.nx,i)=W.y(1:2*sys.nx,i)+opts.lambda*...
                    (sys.F{i,1}(1:2*sys.nx,:)*Z.X(:,i)-t.x(:,i-1));
                Y.y1(2*sys.nx+1:end,i)=W.y(2*sys.nx+1:end,i)...
                    +opts.lambda*(sys.G{i,1}(2*sys.nx+1:end,:)*Z.U(:,i)-t.u(:,i));
            end
        end
    else
        Y.y1(1:2*sys.nx,2:Nd)=W.y(1:2*sys.nx,2:Nd)+opts.lambda*(sys.F(1:2*sys.nx,:)*Z.X(:,2:Nd)-t.x(:,1:Nd-Ns));
        Y.y1(2*sys.nx+1:end,:)=W.y(2*sys.nx+1:end,:)+opts.lambda*(sys.G(2*sys.nx+1:end,:)*Z.U-t.u);
        Y.yt1=W.yt+opts.lambda*(sys.F(1:2*sys.nx,:)*Z.X(:,Nd+1:Nd+Ns)-t.x(:,Nd-Ns+1:Nd+Ns-1));
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
    %if(j>2)
     %   if(jobj(j)<jobj(j-1))
     %       theta=[1 1]';
     %   end
    %end
     j=j+1;
end
figure
plot(jobj);
details.gpad_solve=toc;
details.W=W;
details.jobj=jobj;
%details.epsilon_prm_avg= epsilon_prm_avg;
%details.epsilon_prm=epsilon_prm;
end




