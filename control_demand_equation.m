function [eq_const] = control_demand_equation(Z,Tree,apg_opts)

% This function calculate the control-demand equation.  
%

Nd=size(Tree.stage,1);
Ns=size(Tree.children,1);
ne=size(apg_opts.E,1);

eq_const=zeros(ne,Nd);
for i=1:Nd-Ns+1
    if(i==1)
        eq_const(:,i)=apg_opts.E*Z.U(:,1)+apg_opts.Ed*apg_opts.demand(1,:)';
    else
        stage=Tree.stage(i-1)+2;
        nchild=Tree.children{i-1};
        for l=1:length(nchild)
            eq_const(:,nchild(l))=apg_opts.E*Z.U(:,nchild(l))+apg_opts.Ed*(apg_opts.demand(stage,:)+...
                Tree.value(nchild(l),:))';
        end
    end   
end
end

