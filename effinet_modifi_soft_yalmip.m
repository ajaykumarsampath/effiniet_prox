function [ zoptimiser ] = effinet_modifi_soft_yalmip( sys,Tree,V,opts)
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
Xs=sdpvar(sys.nx,Nd+1);
U=sdpvar(sys.nu,Nd);


const=(X(:,1)==xinit)+(Xs(:,1)==zeros(sys.nx,1));
Jobj=0;
for i=1:Nd-Ns+1
    if(i==1)
        Jobj=Jobj+V.alpha(1,:)*U(:,1)+(U(:,1)-uinit)'*V.Wu*(U(:,1)-uinit)+...
            Xs(:,2)'*V.Qs*Xs(:,2);
        const=const+(opts.E*U(:,1)+opts.Ed*opts.demand(1,:)'==0);
        const=const+(X(:,2)==sys.A*X(:,1)+sys.B*U(:,1)+sys.Gd*opts.demand(1,:)');
        const=const+(sys.umin(1:sys.nu,1)<=U(:,1)<=sys.umax(1:sys.nu,1));
        const=const+(sys.xmin(1:sys.nx,1)<=X(:,2)<=sys.xmax(1:sys.nx,1));
        const=const+(zeros(sys.nx,1)<=Xs(:,2));
        const=const+(sys.xs(1:sys.nx,1)<=Xs(:,2)+X(:,2));
    else
        stage=Tree.stage(i-1)+1;        
        
        nchild=Tree.children{i-1};
        
        for l=1:length(nchild)
            Jobj=Jobj+Tree.prob(nchild(l),1)*(V.alpha(stage+1,:)*U(:,nchild(l))+(U(:,nchild(l))-...
                U(:,Tree.ancestor(nchild(l))))'*V.Wu*(U(:,nchild(l))-U(:,Tree.ancestor(nchild(l)))))...
                +Xs(:,nchild(l))'*V.Qs*Xs(:,nchild(l));
            
            const=const+(X(:,nchild(l)+1)==sys.A*X(:,i)+sys.B*U(:,nchild(l))+sys.Gd*...
                (opts.demand(stage+1,:)+Tree.value(nchild(l),:))');
            
            const=const+(sys.umin((nchild(l)-1)*sys.nu+1:nchild(l)*sys.nu,1)<=...
                U(:,nchild(l))<=sys.umax((nchild(l)-1)*sys.nu+1:nchild(l)*sys.nu,1));
            
            const=const+(sys.xmin((nchild(l)-1)*sys.nx+1:nchild(l)*sys.nx,1)<=...
                X(:,nchild(l)+1)<=sys.xmax((nchild(l)-1)*sys.nx+1:nchild(l)*sys.nx,1));
            
            const=const+(zeros(sys.nx,1)<=Xs(:,nchild(l)+1));
            
            const=const+(sys.xs((nchild(l)-1)*sys.nx+1:nchild(l)*sys.nx,1)...
                <=Xs(:,nchild(l)+1)+X(:,nchild(l)+1));
            
            const=const+(opts.E*U(:,nchild(l))+opts.Ed*(opts.demand(stage+1,:)+...
                Tree.value(nchild(l),:))'==0);
            
        end
    end   
end

zoptimiser=optimizer(const,Jobj,default_options,{xinit,uinit},{X,U,Xs,Jobj});

end



