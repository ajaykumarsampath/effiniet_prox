function [ zoptimiser ] = effiniet_yalmip( sys,Tree,V,opts)
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
U=sdpvar(sys.nu,Nd-Ns+1);


const=(X(:,1)==xinit);
Jobj=0;
for i=1:Nd-Ns+1
    if(i==1)
        Jobj=Jobj+V.alpha(1,:)*U(:,1)+(U(:,1)-uinit)'*V.Wu*(U(:,1)-uinit);
        const=const+(opts.E*U(:,1)+opts.Ed*opts.demand(1,:)'==0);
        const=const+(X(:,2)==sys.A*X(:,1)+sys.B*U(:,1)+sys.Gd*opts.demand(1,:)');
        const=const+(sys.umin(1:sys.nu,1)<=U(:,1)<=sys.umax(1:sys.nu,1));
    else
        stage=Tree.stage(i-1)+1;
        Jobj=Jobj+Tree.prob(i-1,1)*(V.alpha(stage+1,:)*U(:,i)+(U(:,i)-U(:,Tree.ancestor(i-1)+1))'*V.Wu*...
            (U(:,i)-U(:,Tree.ancestor(i-1)+1)));
        const=const+(sys.umin((i-1)*sys.nu+1:i*sys.nu,1)<=U(:,i)<=sys.umax((i-1)*sys.nu+1:i*sys.nu,1));
        const=const+(sys.xmin((i-1)*sys.nx+1:i*sys.nx,1)<=X(:,i)<=sys.xmax((i-1)*sys.nx+1:i*sys.nx,1));
        
        nchild=Tree.children{i-1};
        for l=1:length(nchild)
            const=const+(X(:,nchild(l)+1)==sys.A*X(:,i)+sys.B*U(:,i)+sys.Gd*...
                (opts.demand(stage+1,:)+Tree.value(nchild(l),:))');
            const=const+(opts.E*U(:,i)+opts.Ed*(opts.demand(stage+1,:)+...
                Tree.value(nchild(l),:))'==0);
        end
    end   
end

for i=1:Ns
    const=const+(sys.xmin(1:sys.nx,1)<=X(:,Tree.leaves(i)+1)<=sys.xmax(1:sys.nx,1));
end

zoptimiser=optimizer(const,Jobj,default_options,{xinit,uinit},{X,U,Jobj});
end


