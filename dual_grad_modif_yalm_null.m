function [ zoptimiser ] = dual_grad_modif_yalm_null( sys,Tree,Ptree,V,opts)
%  
% Generate the yalmip optimizer for dual gradient calculation 

% opts=grad_opts_null;
default_options=sdpsettings('solver','gurobi','verbose',0,'cachesolvers',1);

Nd=size(Tree.stage,1);
Ns=size(Tree.leaves,1);

nv=size(sys.L,2);
ny=size(sys.F{1},1);

xinit=sdpvar(sys.nx,1);
vinit=sdpvar(nv,1);
X=sdpvar(sys.nx,Nd+1);
v=sdpvar(nv,Nd);
U=sdpvar(sys.nu,Nd);

Y=sdpvar(ny,Nd);
Yt=sdpvar(sys.nx,Ns);

const=(X(:,1)==xinit);
Jobj=0;
Wv=sys.L'*V.Wu*sys.L;

Wv1=V.Wu*sys.L;

vhat=opts.vhat;
prev_vhat=opts.prev_vhat;

for i=1:Nd-Ns+1
    if(i==1)
        Jobj=Jobj+opts.alpha_bar(:,1)'*v(:,1)+(v(:,1)-vinit)'*Wv*(v(:,1)-vinit)...
            +2*(vhat(:,1)-prev_vhat)'*Wv1*(v(:,1)-vinit)+Y(:,1)'*Ptree.Gbar{1}*v(:,i);
        const=const+(U(:,1)==sys.L*v(:,1)+opts.vhat(:,1));
        const=const+(X(:,2)==sys.A*X(:,1)+sys.B*U(:,1)+opts.w(:,1));
    else
        stage=Tree.stage(i-1)+1;
        
        nchild=Tree.children{i-1};
        
        for l=1:length(nchild)
            Jobj=Jobj+Tree.prob(nchild(l),1)*(opts.alpha_bar(:,stage+1)'*v(:,nchild(l))+...
                (v(:,nchild(l))-v(:,Tree.ancestor(nchild(l))))'*Wv*...
                (v(:,nchild(l))-v(:,Tree.ancestor(nchild(l)))))...
                +Y(:,nchild(l))'*(sys.F{nchild(l)}*X(:,i)+Ptree.Gbar{nchild(l)}*v(:,nchild(l)));
            
            const=const+(U(:,nchild(l))==sys.L*v(:,nchild(l))+opts.vhat(:,nchild(l)));
            
            const=const+(X(:,nchild(l)+1)==sys.A*X(:,i)+sys.B*U(:,nchild(l))+...
                sys.Gd*(opts.demand(stage+1,:)+Tree.value(nchild(l),:))');
        end
    end   
end

for i=1:Ns
    Jobj=Jobj+Yt(:,i)'*...
        (sys.F{Nd+1+i}(1:sys.nx,:)*X(:,Tree.leaves(i)+1));
end

zoptimiser=optimizer(const,Jobj,default_options,{xinit,vinit,Y,Yt},{X,v,U,Jobj});

end

