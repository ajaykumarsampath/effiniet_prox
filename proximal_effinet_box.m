function [ t ] = proximal_effinet_box(Z,W,sys,opts_prox)
% This function calcualte the proximal for the state and input constraints
% soft constraits on state,
% hard constraints on input.
%
% INPUT-----   Z   :  state and control
%              W   :  dual variable
%      opts_prox   :  proximal options.
%

%dim=(size(W.y,2)+size(W.yt,2)-1)*opts_prox.nx;
if(sys.cell)
    hard_const.x=[W.y(1:opts_prox.nx,2:end) W.yt]/opts_prox.lambda;
    hard_const.u=W.y(opts_prox.nx+1:end,:)/opts_prox.lambda;
    Nd=size(hard_const.x,2);
    Ns=size(W.yt,2);
    for i=1:Nd
        hard_const.x(:,i)=hard_const.x(:,i)+sys.F{i+1,1}(1:sys.nx,:)*Z.X(:,i+1);
        if(i<Nd-Ns+1)
            hard_const.u(:,i)=hard_const.u(:,i)+sys.G{i,1}(sys.nx+1:end,:)*Z.U(:,i);
        end
    end 
    
    opts_prox.xmax=reshape(opts_prox.xmax,sys.nx,size(hard_const.x,2));
    opts_prox.xmin=reshape(opts_prox.xmin,sys.nx,size(hard_const.x,2));
    % hard constraints on input u
    t.u=zeros(sys.nu,size(hard_const.u,2));
    t.u=min(opts_prox.umax,hard_const.u(1:sys.nu,:));
    t.u=max(opts_prox.umin,t.u);
    
    % hard constraints on input u
    t.x=zeros(sys.nx,size(hard_const.x,2));
    t.x=min(opts_prox.xmax,hard_const.x(1:sys.nx,:));
    t.x=max(opts_prox.xmin,t.x);
else
    hard_const.x=[W.y(1:opts_prox.nx,2:end) W.yt]/opts_prox.lambda+sys.F(1:sys.nx,:)*Z.X(:,2:end);
    hard_const.u=W.y(opts_prox.nx+1:end,:)/opts_prox.lambda+sys.G(sys.nx+1:end,:)*Z.U;
    
    opts_prox.xmax=reshape(opts_prox.xmax,sys.nx,size(hard_const.x,2));
    opts_prox.xmin=reshape(opts_prox.xmin,sys.nx,size(hard_const.x,2));
    % hard constraints on input u
    t.u=zeros(sys.nu,size(hard_const.u,2));
    t.u=min(opts_prox.umax,hard_const.u(1:sys.nu,:));
    t.u=min(-opts_prox.umin,t.u);
    
    % hard constraints on input u
    t.x=zeros(sys.nx,size(hard_const.x,2));
    t.x=min(opts_prox.xmax,hard_const.x(1:sys.nx,:));
    t.x=min(-opts_prox.xmin,t.x);
end


%{
% soft constraints on state
p=[1/2 1/2];
yq=reshape(soft_const,dim,1);
prox_details.y=zeros(dim,2);
prox_details.x=zeros(dim,2);
prox_details.gamma_c=opts_prox.gamma_c/(opts_prox.lambda*p(1));
prox_details.gamma_xs=opts_prox.gamma_s/(opts_prox.lambda*p(2));
i=1;
imax=1000;
while(i<imax)
    
    prox_details.x(:,1)=yq-prox_details.y(:,1);
    prox_details.proj_c=min(prox_details.x(:,1),opts_prox.xmax);
    prox_details.proj_c=max(prox_details.proj_c,opts_prox.xmin);
    
    if(norm(prox_details.x(:,1)-prox_details.proj_c,2)>prox_details.gamma_c)
        prox_details.x(:,1)=prox_details.x(:,1)+prox_details.gamma_c*...
            (prox_details.proj_c-prox_details.x(:,1))/(norm(prox_details.x(:,1)-prox_details.proj_c,2));
    else
        prox_details.x(:,1)=prox_details.proj_c;
    end
    
    prox_details.x(:,2)=yq-prox_details.y(:,2);
    prox_details.proj_xs=max(prox_details.x(:,2),opts_prox.xs);
    if(norm(prox_details.x(:,2)-prox_details.proj_xs,2)>prox_details.gamma_xs)
        prox_details.x(:,2)=prox_details.x(:,2)+prox_details.gamma_xs*...
            (prox_details.proj_xs-prox_details.x(:,2))/(norm(prox_details.x(:,2)-prox_details.proj_xs,2));
    else
        prox_details.x(:,2)=prox_details.proj_xs;
    end
    prox_details.y=prox_details.y+prox_details.x-kron([1 1],mean(prox_details.x')');
    if(max(abs(prox_details.x(:,1)-prox_details.x(:,2)))<0.001)
        prox_details.iter=i;
        i=imax+200;
        %i=i+1;
    else
        i=i+1;
    end
end
t.x=reshape(prox_details.x(:,1),opts_prox.nx,dim/opts_prox.nx);
%}
end

