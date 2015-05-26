function [ Z,details] = APG_effiniet(sys,Ptree,Tree,V,opts)

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
Y.y0=zeros(size(sys.F,1),Nd-Ns+1);
Y.yt0=zeros(2*sys.nx,Ns);
Y.y1=zeros(size(sys.F,1),Nd-Ns+1);
Y.yt1=zeros(2*sys.nx,Ns);
W.y=zeros(size(sys.F,1),Nd-Ns+1);
W.yt=zeros(2*sys.nx,Ns);

epsilon_prm=1;

prm_fes.soft=zeros(sys.nx,Nd);
prm_fes.hard=zeros(sys.nu,Nd+1-Ns);

theta=[1 1]';

details.term_crit=zeros(1,4);
dual_grad=prm_fes;

opts_prox.lambda=opts.lambda;
opts_prox.xmin=kron(ones(Nd,1),sys.xmin);
opts_prox.xs=kron(ones(Nd,1),sys.xs);
opts_prox.xmax=kron(ones(Nd,1),sys.xmax);

opts_prox.umin=kron(ones(1,Nd-Ns+1),sys.umin);
opts_prox.umax=kron(ones(1,Nd-Ns+1),sys.umax);
opts_prox.nx=sys.nx;
opts_prox.nu=sys.nu;
opts_prox.gamma_c=sys.gamma_max;
opts_prox.gamma_s=sys.gamma_xs;

uprev=sys.L*opts.state.v+opts.state.prev_vhat;
%}
grad_opts.u=uprev;
grad_opts.x=opts.x;
grad_opts.E=opts.E;
grad_opts.Ed=opts.Ed;
grad_opts.demand=opts.state.demand;
dual_grad_yalmip=dual_gradient_yalmip(sys,Tree,V,grad_opts);
%}
%%
%tic
j=1;
while(j<opts.steps)
    % Step 1: accelerated step
    W.y=Y.y1+theta(2)*(1/theta(1)-1)*(Y.y1-Y.y0);
    
    W.yt=Y.yt1+theta(2)*(1/theta(1)-1)*(Y.yt1-Y.yt0);
    
    
    %step 2: argmin of the lagrangian using dynamic programming
    %Z=solve_apg_effiniet(sys,Tree,Ptree,W,opts.x,opts.state);
    
    
    results=dual_grad_yalmip{{grad_opts.x,grad_opts.u,W.y,W.yt}};
    Z.X=results{1,1};
    Z.U=results{1,2};
    %step 3: proximal with respect to g
    t=proximal_effinet_2(Z,W,sys,opts_prox);
    
    %step 4: update the dual vector
    Y.y0=Y.y1;
    Y.yt0=Y.yt1;
    
    %dual_grad(1:sys.nu,:)=Z.U;
    %dual_grad(sys.nu+1:end,1:Nd-1)=Z.X(:,2:Nd);
    %prm_fes.soft(1:sys.nx,1:Nd+Ns-1)=(Z.X(:,2:Nd+Ns)-t.x(:,1:Nd+1-Ns));
    %prm_fes.hard=(Z.U-t.u);
    
    Y.y1(1:2*sys.nx,2:Nd)=W.y(1:2*sys.nx,2:Nd)+opts.lambda*(sys.F(1:2*sys.nx,:)*Z.X(:,2:Nd)-t.x(:,1:Nd-Ns));
    Y.y1(2*sys.nx+1:end,:)=W.y(2*sys.nx+1:end,:)+opts.lambda*(sys.G(2*sys.nx+1:end,:)*Z.U-t.u);
    Y.yt1=W.yt+opts.lambda*(sys.F(1:2*sys.nx,:)*Z.X(:,Nd+1:Nd+Ns)-t.x(:,Nd-Ns+1:Nd+Ns-1));
    
    iter=j;
    details.prm_cst(iter)=0;%primal cost;
    details.dual_cst(iter)=0;% dual cost;
    
    theta(1)=theta(2);
    theta(2)=(sqrt(theta(1)^4+4*theta(1)^2)-theta(1)^2)/2;
    j=j+1;
    %%
end
%%
%termination criteria
%{
    if(j==1)
        prm_avg_next=prm_fes;
        epsilon_prm_avg=max(max(max(prm_fes.hard)),max(max(prm_fes.soft)));
    else
        prm_avg_next.hard=(1-theta(2))*prm_avg_next.hard+theta(2)*prm_fes.hard;
        prm_avg_next.soft=(1-theta(2))*prm_avg_next.soft+theta(2)*prm_fes.soft;
        epsilon_prm_avg=max(max(max(prm_avg_next.hard)),max(max(prm_avg_next.soft)));
    end
    if epsilon_prm_avg<=opts.primal_inf %average_primal feasibility less
        
        details.term_crit(1,2)=1;
        details.iterate=j;
        j=10*opts.steps;
    else
        epsilon_prm=max( max(max(prm_fes.soft)), ...
            max(max(prm_fes.hard)) );
        if(epsilon_prm<=opts.primal_inf) % primal feasibility of the iterate
            if (min(min(min(W.y)),min( W.yt))>0)
                sum=0;
                for i=1:Nd-Ns+1
                    if(i==1)
                        sum=sum-W.y(:,i)'*[zeros(sys.nx,1);prm_fes.hard(:,i)];
                    else
                        sum=sum-W.y(:,i)'*[prm_fes.soft(:,i-1);prm_fes.hard(:,i)];
                    end
                    
                end
                for i=1:Ns
                    sum=sum-W.yt(:,i)'*(prm_fes.soft(:,Nd+i-Ns));
                end
                if sum<=opts.dual_gap %condition 29. dual gap
                    details.term_crit(1,2)=1;
                    details.iterate=j;
                    j=10*opts.steps;
                else
                    prm_cst=0;%primal cost;
                    for i=1:Nd-Ns+1
                        if(i==1)
                            prm_cst=prm_cst+Tree.prob(i,1)*((Z.U(:,i)-uprev)'*V.Wu*(Z.U(:,i)-uprev)+...
                                V.alpha(i,:)*Z.U(:,i));
                        else
                            prm_cst=prm_cst+Tree.prob(i,1)*((Z.U(:,i)-Z.U(:,Tree.ancestor(i)))'*...
                                V.R*(Z.U(:,i)-Z.U(:,Tree.ancestor(i)))+V.alpha(i,:)*Z.U(:,i));
                        end
                        
                    end
                    for i=1:Ns
                        prm_cst=prm_cst+Tree.prob(Tree.leaves(i))*(Z.X(:,Tree.leaves(i))'*V.Vf{i,1}*...
                            Z.X(:,Tree.leaves(i)));
                    end
                    if sum<=ops.dual_gap*prm_cst/(1+ops.dual_gap) %condition 30 dual gap
                        details.term_crit(1,3)=1;
                        details.iterate=j;
                        j=10*opts.steps;
                    else
                        %step 4: theta update
                        theta(1)=theta(2);
                        theta(2)=(sqrt(theta(1)^4+4*theta(1)^2)-theta(1)^2)/2;
                        j=j+1;
                    end
                end
            else
                %
                prm_cst=0;%primal cost;
                dual_cst=0;% dual cost;
                for i=1:Nd-Ns
                    prm_cst=prm_cst+Tree.prob(i,1)*(Z.X(:,i)'*V.Q*Z.X(:,i)+Z.U(:,i)'*V.R*Z.U(:,i));
                    dual_cst=dual_cst+Y.y1(i,:)*(prm_fes(:,i)-g_nodes(:,i));
                end
                for i=1:Ns
                    prm_cst=prm_cst+Tree.prob(Tree.leaves(i))*(Z.X(:,Tree.leaves(i))'*V.Vf{i,1}*...
                        Z.X(:,Tree.leaves(i)));
                    dual_cst=dual_cst+Y.yt1{i,:}*(prm_fes_term{i,1}-g_nodes_term{i,1});
                end
                if (-dual_cst<=ops.dual_gap*max(dual_cst,1)) %condtion 27 (dual gap)
                    details.term_crit(1,4)=1;
                    details.iterate=j;
                    j=10*ops.steps;
                else
                %
                    %step 4: theta update
                    theta(1)=theta(2);
                    theta(2)=(sqrt(theta(1)^4+4*theta(1)^2)-theta(1)^2)/2;
                    j=j+1;
                %end
            %end
            else
            %step 4: theta update
            theta(1)=theta(2);
            theta(2)=(sqrt(theta(1)^4+4*theta(1)^2)-theta(1)^2)/2;
            j=j+1;
%}
%end
%end
%details.epsilon_prm_avg(iter)=epsilon_prm_avg;
%details.epsilon_prm(iter)=epsilon_prm;

%end
%%
details.gpad_solve=toc;
details.W=W;
%details.epsilon_prm_avg= epsilon_prm_avg;
%details.epsilon_prm=epsilon_prm;
end


