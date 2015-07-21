function [ t ] = prox_effinet_dist_box(Z,W,sys,Tree,opts_prox)
% This function calcualte the proximal for the state and input constraints
% soft constraits on state,
% hard constraints on input.
%
% INPUT-----   Z   :  state and control
%              W   :  dual variable
%      opts_prox   :  proximal options.
%%

distance=[0;0];
if(sys.cell)
    constraints.x=W.x/opts_prox.lambda;
    constraints.u=W.u/opts_prox.lambda;
    Nd=size(Tree.stage,1);
    
    
    for i=1:Nd
        constraints.x(:,i)=constraints.x(:,i)+sys.F{i,1}*Z.X(:,i+1);
        constraints.u(:,i)=constraints.u(:,i)+sys.G{i,1}*Z.U(:,i);
    end
    
    
    % hard constraints on input u
    t.u=zeros(sys.nu,size(constraints.u,2));
    t.u=min(opts_prox.umax,constraints.u(1:sys.nu,:));
    t.u=max(opts_prox.umin,t.u);
    
    % soft constraints on input x
    opts_prox.xmax=reshape(opts_prox.xmax,sys.nx*size(constraints.x,2),1);
    opts_prox.xmin=reshape(opts_prox.xmin,sys.nx*size(constraints.x,2),1);
    opts_prox.xs=reshape(opts_prox.xs,sys.nx*size(constraints.x,2),1);
    
    t.x=zeros(2*sys.nx,size(constraints.x,2));
    xtemp=reshape(constraints.x(1:sys.nx,:),sys.nx*Nd,1);
    
    ProjSet1=min(opts_prox.xmax,xtemp);
    ProjSet1=max(opts_prox.xmin,ProjSet1);
    
    distance(1)=norm(xtemp-ProjSet1,2);
    if(distance(1)>opts_prox.gamma_xbox)
        xtemp=xtemp+opts_prox.xbox*(ProjSet1-xtemp)/distance(1);
    else
        xtemp=ProjSet1;
    end
    
    t.x(1:sys.nx,:)=reshape(xtemp,sys.nx,Nd);
    xtemp=reshape(constraints.x(sys.nx+1:2*sys.nx,:),sys.nx*Nd,1);
    
    %ProjSet2=min(opts_prox.xmax,xtemp);
    ProjSet2=max(opts_prox.xs,xtemp);
    
    distance(2)=norm(xtemp-ProjSet2,2);
    if(distance(2)>opts_prox.gamma_xs)
        xtemp=xtemp+opts_prox.gamma_xs*(ProjSet2-xtemp)/distance(2);
    else
        xtemp=ProjSet2;
    end
    
    t.x(sys.nx+1:2*sys.nx,:)=reshape(xtemp,sys.nx,Nd);
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

