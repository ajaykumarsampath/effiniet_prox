function [ sys,V] = system_prox_dist_state_box( S,P,Tree,ops_sys)
% This function transform the effiniet problem into
% dual proximal gradient formulation.
%
%   INPUT-----    S  : Effinet system discription
%                 P  : Effinet cost description
%              Tree  : Scenario tree discription
%
%   OUTPUT------ sys : Output system for dual proximal form

sys.nu=S.nu;
sys.nx=S.nx;
sys.Np=P.Hp;
Nd=size(Tree.stage,1);
Ns=size(Tree.leaves,1);

sys.xmin=S.xmin;
sys.xmax=S.xmax;
sys.xs=P.xs;
%sys.xs=sys.xmin;

sys.gamma_xbox=ops_sys.gamma_xbox;
sys.gamma_xs=ops_sys.gamma_xs;

% Null spaces and particular soultion format
sys.L=null(S.E);
sys.L1=-pinv(S.E)*(S.Ed);

sys.A=S.A;
sys.B=S.B/3600;
sys.Gd=S.Gd/3600;
sys.umin=3600*S.umin;
sys.umax=3600*S.umax;

% cost function
V.Wu=P.Wu;

% Normalise the constraints
if(ops_sys.cell)
    sys.F=cell(Nd+1,1);
    sys.G=cell(Nd+1-Ns,1);
    for j=1:Nd
        sys.F{j}=[eye(S.nx);eye(S.nx)];
        sys.G{j}=eye(S.nu);
    end
else
    sys.F=[eye(S.nx);eye(S.nx)];
    sys.G=eye(S.nu);
end
sys.cell=ops_sys.cell;
sys.normalized=ops_sys.normalise;

sys.umin=kron(ones(Nd,1),sys.umin);
sys.umax=kron(ones(Nd,1),sys.umax);

sys.xmin=kron(ones(Nd,1),sys.xmin);
sys.xmax=kron(ones(Nd,1),sys.xmax);
sys.xs=kron(ones(Nd,1),sys.xs);
%}

%}
end



