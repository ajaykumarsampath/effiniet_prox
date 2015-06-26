function [ zoptimiser] = dual_gradient_modified_yalmip( sys,Tree,V,opts)
%  
% Generate the yalmip optimizer for dual gradient calculation 

% opts=grad_opts;
default_options=sdpsettings('solver','gurobi','verbose',0,'cachesolvers',1);

Nd=size(Tree.stage,1);
Ns=size(Tree.leaves,1);

%ny=size(sys.F,1);

xinit=sdpvar(sys.nx,1);
uinit=sdpvar(sys.nu,1);
X=sdpvar(sys.nx,Nd+1);
U=sdpvar(sys.nu,Nd);
Y=sdpvar(size(sys.F{1},1),Nd);
Yt=sdpvar(sys.nx,Ns);

const=(X(:,1)==xinit);
Jobj=0;
for i=1:Nd-Ns+1
    if(i==1)
        Jobj=Jobj+V.alpha(1,:)*U(:,1)+(U(:,1)-uinit)'*V.Wu*(U(:,1)-uinit)...
            +Y(:,1)'*sys.G{i}*U(:,i);
        const=const+(opts.E*U(:,1)+opts.Ed*opts.demand(1,:)'==0);
        const=const+(X(:,2)==sys.A*X(:,1)+sys.B*U(:,1)+sys.Gd*opts.demand(1,:)');
    else
        stage=Tree.stage(i-1)+1;
        
        nchild=Tree.children{i-1};
        
        for l=1:length(nchild)
            Jobj=Jobj+Tree.prob(nchild(l),1)*(V.alpha(stage+1,:)*U(:,nchild(l))+...
                (U(:,nchild(l))-U(:,Tree.ancestor(nchild(l))))'*V.Wu*...
                (U(:,nchild(l))-U(:,Tree.ancestor(nchild(l)))))...
                +Y(:,nchild(l))'*(sys.F{nchild(l)}*X(:,i)+sys.G{nchild(l)}*U(:,nchild(l)));
            
            const=const+(X(:,nchild(l)+1)==sys.A*X(:,i)+sys.B*U(:,nchild(l))+...
                sys.Gd*(opts.demand(stage+1,:)+Tree.value(nchild(l),:))');
            
            const=const+(opts.E*U(:,nchild(l))+opts.Ed*(opts.demand(stage+1,:)+...
                Tree.value(nchild(l),:))'==0);
        end
    end   
end

for i=1:Ns
    Jobj=Jobj+Yt(:,i)'*...
        (sys.F{Nd+1+i}(1:sys.nx,:)*X(:,Tree.leaves(i)+1));
end

zoptimiser=optimizer(const,Jobj,default_options,{xinit,uinit,Y,Yt},{X,U,Jobj});
end


