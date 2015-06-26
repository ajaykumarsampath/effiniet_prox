function [ sys_new,pre_cnd] = precond_modify_sdp_box( sys,DH_nml)
%%
%   Calculate the optimal preconditioner for the dual accelerated gradient
%   method
%   Input  -----   DH_nml     : Dual Hessian.
%                  sys        : Original system.
%
%   Output -----   pre_cnd    : Preconditioned matrixs
%                               Ef1---  constraints on state
%                               Eg1---  constraints on control
%                               Eft1--- constraints on terminal constraints.
%                  sys_new    : preconditioned system.
%                  opts       : const_type---box
%
sdpvar_opts=sdpsettings('solver','sedumi','verbose',0,'cachesolvers',1);

q=rank(DH_nml);
ny=size(DH_nml,1);
[U,S,V]=svd(DH_nml);
sd=sqrt(diag(S));
sd(q+1:end,1)=zeros(ny-q,1);
Q1=U*diag(sd);
Qt=Q1(:,1:q)';


Np=sys.Np;
Fsdp=sdpvar(sys.nx,sys.nx,'diagonal');
Gsdp=sdpvar(sys.nu,sys.nu,'diagonal');


for i=1:Np
    if(i==1)
        P=Gsdp;
        P=blkdiag(P,Fsdp);
    else
        P=blkdiag(P,Gsdp,Fsdp);
    end
end

t=sdpvar(1);
const=(eye(q)<=Qt*P*Qt'<=t*eye(q));
const=const+(1e-3*eye(ny)<=P);

diagnos=solvesdp(const,t,sdpvar_opts);
if diagnos.problem==0
    disp('feasible');
else
    disp('Infeasible');
end
%Pre1=value(P);
%eig=value(t);

Gvalue=value(Gsdp);
Fvalue=value(Fsdp);
%Ftvalue=value(Ftsdp);

pre_cnd.Ef1=sqrt(Fvalue);
pre_cnd.Eg1=sqrt(Gvalue);
sys_new=sys;

for j=1:Np
    sys_new.G{j,1}(sys.nx+1:end,:)=pre_cnd.Eg1*sys.G{j,1}(sys.nx+1:end,:);
    sys_new.F{j,1}(1:sys.nx,:)=pre_cnd.Ef1*sys.F{j,1}(1:sys.nx,:);
end

end

