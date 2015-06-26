function [ Pre,t  ] = precond_sdp( DH_nml )
sdpvar_opts=sdpsettings('solver','sedumi','verbose',0,'cachesolvers',1);

q=rank(DH_nml);
ny=size(DH_nml,1);
[U,S,V]=svd(DH_nml);
sd=sqrt(diag(S));
sd(q+1:end,1)=zeros(ny-q,1);
Q1=U*diag(sd);
Qt=Q1(:,1:q)';
%norm(DH_nml-Qt'*Qt)
%Q=sdpvar(q,ny);
P=sdpvar(ny,ny,'diagonal');
t=sdpvar(1);
const=(eye(q)<=Qt*P*Qt'<=t*eye(q));
const=const+(1e-3*eye(ny)<=P);

%precond_optimiser=optimizer(const,t,sdpvar_opts,{Q},{t,P});
%diagnos=optimize(const,t,sdpvar_opts);
diagnos=solvesdp(const,t,sdpvar_opts);
if diagnos.problem==0
    disp('feasible');
else 
    disp('Infeasible');
end
Pre1=value(P);
t=value(t);

Pre=sqrt(Pre1);

end

