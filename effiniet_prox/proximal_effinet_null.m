function [ t ] = proximal_effinet_null(Z,W,sys,opts_prox)
% This function calcualte the proximal for the state and input constraints
% soft constraits on state,
% hard constraints on input.
%
% INPUT-----   Z   :  state and control
%              W   :  dual variable
%      opts_prox   :  proximal options.
%

%dim=(size(W.y,2)+size(W.yt,2)-1)*opts_prox.nx;
hard_const.x=[W.y(1:2*opts_prox.nx,2:end) W.yt]/opts_prox.lambda+sys.F(1:2*sys.nx,:)*Z.X(:,2:end);
hard_const.u=W.y(2*opts_prox.nx+1:end,:)/opts_prox.lambda+sys.G(2*sys.nx+1:end,:)*Z.v;

opts_prox.xmax=reshape(opts_prox.xmax,sys.nx,size(hard_const.x,2));
opts_prox.xmin=reshape(opts_prox.xmin,sys.nx,size(hard_const.x,2));
% hard constraints on input u
t.u=zeros(2*opts_prox.nu,size(hard_const.u,2));
t.u=min(opts_prox.umax,hard_const.u(1:opts_prox.nu,:));
t.u(opts_prox.nu+1:2*opts_prox.nu,:)=min(-opts_prox.umin,hard_const.u(opts_prox.nu+1:2*opts_prox.nu,:));

% hard constraints on input u
t.x=zeros(2*sys.nx,size(hard_const.x,2));
t.x=min(opts_prox.xmax,hard_const.x(1:sys.nx,:));
t.x(sys.nx+1:2*sys.nx,:)=min(-opts_prox.xmin,hard_const.x(sys.nx+1:2*sys.nx,:));

end



