function [ sys_null,V] = system_prox_formation_null( S,P,Tree,ops_sys)
% This function transform the effiniet problem into 
% dual proximal gradient formulation.
%
%   INPUT-----    S  : Effinet system discription 
%                 P  : Effinet cost description
%              Tree  : Scenario tree discription 
% 
%   OUTPUT------ sys : Output system for dual proximal form 

% Null spaces and particular soultion format
sys_null.L=null(S.E);
sys_null.L1=-pinv(S.E)*(S.Ed);


sys_null.nu=size(sys_null.L,2);
sys_null.nx=S.nx;
sys_null.Np=P.Hp;
Nd=size(Tree.stage,1);
Ns=size(Tree.leaves,1);

%sys.F=[S.A;zeros(S.nu,S.nx)];
%sys.G=[S.B;eye(S.nu)];

sys_null.xmin=P.xs;
sys_null.xmax=S.xmax;
sys_null.xs=P.beta*sys_null.xmax;

sys_null.umin=S.umin*3600;
sys_null.umax=S.umax*3600;

sys_null.gamma_min=ops_sys.gamma_min;
sys_null.gamma_xs=ops_sys.gamma_xs;
sys_null.gamma_max=ops_sys.gamma_max;

sys_null.A=S.A;
sys_null.B=S.B*sys_null.L/3600;
sys_null.Gd=S.Gd/3600;

% cost function
V.Wu=P.Wu;

% Normalised the constraints
if(ops_sys.cell)
    sys_null.F=cell(Nd+1,1);
    sys_null.G=cell(Nd+1-Ns,1);
    for j=1:Nd+1
        sys_null.F{j}=[eye(S.nx);-eye(S.nx);zeros(2*S.nu,S.nx)];
        sys_null.G{j}=[zeros(2*S.nx,sys_null.nu);sys_null.L;-sys_null.L];
        if(ops_sys.normalise)
            for i=1:S.nx
                if(sys_null.xmax(i)>0)
                    sys_null.F{j}(i,i)=sys_null.F{j}(i,i)/sys_null.xmax(i);
                    sys_null.xmax(i)=1;
                end
                
                if(sys_null.xmin(i)>0)
                    sys_null.F{j}(S.nx+i,i)=sys_null.F{j}(S.nx+i,i)/sys_null.xmin(i);
                    sys_null.xmin(i)=1;
                end
            end
            for i=1:S.nu
                if(sys_null.umax(i)>0)
                    sys_null.G{j}(2*S.nx+i,i)=sys_null.G{j}(2*S.nx+i,i)/sys_null.umax(i);
                    sys_null.umax(i)=1;
                end
                
                if(sys_null.umin(i)>0)
                    sys_null.G{j}(2*S.nx+S.nu+i,i)=sys_null.G{j}(2*S.nx+S.nu+i,i)/sys_null.umin(i);
                    sys_null.umin(i)=1;
                end
            end
        end
    end
else
sys_null.F=[eye(S.nx);-eye(S.nx);zeros(2*S.nu,S.nx)];
sys_null.G=[zeros(2*S.nx,sys_null.nu);sys_null.L;-sys_null.L];
    if(ops_sys.normalise)
        for i=1:S.nx
            if(sys_null.xmax(i)>0)
                sys_null.F(i,i)=sys_null.F(i,i)/sys_null.xmax(i);
                sys_null.xmax(i)=1;
            end
            
            if(sys_null.xmin(i)>0)
                sys_null.F(S.nx+i,i)=sys_null.F(S.nx+i,i)/sys_null.xmin(i);
                sys_null.xmin(i)=1;
            end
            
        end
        
        for i=1:S.nu
            if(sys_null.umax(i)>0)
                sys_null.G(2*S.nx+i,i)=sys_null.G(2*S.nx+i,i)/sys_null.umax(i);
                sys_null.umax(i)=1;
            end
            
            if(sys_null.umin(i)>0)
                sys_null.G(2*S.nx+S.nu+i,i)=sys_null.G(2*S.nx+S.nu+i,i)/sys_null.umin(i);
                sys_null.umin(i)=1;
            end
            
        end
    end
end

sys_null.cell=ops_sys.cell;
sys_null.normalized=ops_sys.normalise;

end

