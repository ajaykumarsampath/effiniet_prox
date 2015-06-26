function [ t ] = prox_effinet_modif_state_box(Z,W,sys,Tree,opts_prox)
% This function calcualte the proximal for the state and input constraints
% soft constraits on state,
% hard constraints on input.
%
% INPUT-----   Z   :  state and control
%              W   :  dual variable
%      opts_prox   :  proximal options.
%


if(sys.cell)
    constraints.x=W.y(1:opts_prox.nx,:)/opts_prox.lambda;
    constraints.u=W.y(opts_prox.nx+1:end,:)/opts_prox.lambda;
    Nd=size(Tree.stage,1);
    
    
    for i=1:Nd
        constraints.x(:,i)=constraints.x(:,i)+sys.F{i,1}(1:sys.nx,:)*Z.X(:,i+1);
        constraints.u(:,i)=constraints.u(:,i)+sys.G{i,1}(sys.nx+1:end,:)*Z.U(:,i);
    end 
    
    opts_prox.xmax=reshape(opts_prox.xmax,sys.nx,size(constraints.x,2));
    opts_prox.xmin=reshape(opts_prox.xmin,sys.nx,size(constraints.x,2));
    opts_prox.xs=reshape(opts_prox.xs,sys.nx,size(constraints.x,2));
    % hard constraints on input u
    t.u=zeros(sys.nu,size(constraints.u,2));
    t.u=min(opts_prox.umax,constraints.u(1:sys.nu,:));
    t.u=max(opts_prox.umin,t.u);
    
    % soft constraints on input x
    xtemp=reshape(constraints.x,sys.nx*Nd,1);
    %t.x=proximal_state_admm(xtemp,opts_prox);
    t.x=proximal_effinet_state(xtemp,opts_prox);
    t.x=reshape(t.x,sys.nx,Nd);
    
    %t1.x=zeros(sys.nx,size(constraints.x,2));
    %t1.x=min(opts_prox.xmax,constraints.x(1:sys.nx,:));
    %t1.x=max(opts_prox.xs,t1.x);
    %t1.x=max(opts_prox.xmin,t1.x);
    %max(max(abs(t.x-t1.x)))
else
    constraints.x=[W.y(1:opts_prox.nx,2:end) W.yt]/opts_prox.lambda+sys.F(1:sys.nx,:)*Z.X(:,2:end);
    constraints.u=W.y(opts_prox.nx+1:end,:)/opts_prox.lambda+sys.G(sys.nx+1:end,:)*Z.U;
    
    opts_prox.xmax=reshape(opts_prox.xmax,sys.nx,size(constraints.x,2));
    opts_prox.xmin=reshape(opts_prox.xmin,sys.nx,size(constraints.x,2));
    % hard constraints on input u
    t.u=zeros(sys.nu,size(constraints.u,2));
    t.u=min(opts_prox.umax,constraints.u(1:sys.nu,:));
    t.u=min(-opts_prox.umin,t.u);
    
    % hard constraints on input u
    t.x=zeros(sys.nx,size(constraints.x,2));
    t.x=min(opts_prox.xmax,constraints.x(1:sys.nx,:));
    t.x=min(-opts_prox.xmin,t.x);
end

end

