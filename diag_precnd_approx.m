function [ Precond ] = diag_precnd_approx( sys_box,Ptree_box)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

nx=sys_box.nx;
nu=sys_box.nu;
nz=nx+nu;
ny=size(sys_box.F{1},1);

Nd=size(Ptree_box.P,1);
Ns=Nd-size(Ptree_box.K,1);

DH_approx=zeros(nz+(Nd-Ns)*ny,nz+(Nd-Ns)*ny);

for i=1:Nd+1
    if(i==1)
        DH_approx(1:nu,1:nu)=sys_box.L*Ptree_box.omega{1};
        sys_box.L*Ptree_box.Phi{6}
    elseif(i==Nd+1)
    end
end 
end

