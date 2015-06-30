function [prox,details] = proximal_state_admm(x,opts_prox)

% This function calculate the proximal of the distance function from 2
% sets

dim=size(x,1);
z=zeros(dim,3);
z(:,1)=x;
z(:,2)=x;
z(:,3)=x;
y=zeros(dim,3);
i=0;

gamma=opts_prox.lambda;
WtSet1=opts_prox.gamma_xbox/opts_prox.lambda;
WtSet2=opts_prox.gamma_xs/opts_prox.lambda;
opts_prox.xmax=reshape(opts_prox.xmax,dim,1);
opts_prox.xmin=reshape(opts_prox.xmin,dim,1);
opts_prox.xs=reshape(opts_prox.xs,dim,1);

%zavg=mean(z')';
zavg=zeros(dim,1);
while(i<opts_prox.iter)
    
    z(:,1)=zavg-y(:,1);
    ProjSet1=min(opts_prox.xmax,z(:,1));
    ProjSet1=max(opts_prox.xmin,ProjSet1);
    
   
    
    distance(1)=norm(z(:,1)-ProjSet1,2);
    if(distance(1)>WtSet1)
        z(:,1)=z(:,1)+WtSet1*(ProjSet1-z(:,1))/distance(1);
    else
        z(:,1)=ProjSet1;
    end
    
    z(:,2)=zavg-y(:,2);
    ProjSet2=max(opts_prox.xs,z(:,2));
    distance(2)=norm(z(:,2)-ProjSet2,2); 
    if(distance(2)>WtSet2)
        z(:,2)=z(:,2)+WtSet2*(ProjSet2-z(:,2))/distance(2);
    else
        z(:,2)=ProjSet2;
    end
    
    z(:,3)=zavg-y(:,3);
    z(:,3)=(gamma*x+y(:,3))/(1+gamma);
    
    zavg=mean(z')';
    %t=xt(:,1)-xt(:,2);
    y=y+z-kron([1 1 1],zavg);
    t=z(:,1)-z(:,2);
    if(max(abs(t))<1e-3)
        details.iter=i;
        i=opts_prox.iter+100;
    else
        i=i+1;
    end
end
if(i==opts_prox.iter)
   % max(abs(t))
end
z2=min(opts_prox.xmax,x);
z2=max(opts_prox.xmin,z2);
max(max(abs(z(:,2)-z2)))
details.tolerance=max(abs(t));
prox=z(:,1);
end




