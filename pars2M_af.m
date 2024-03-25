function [M_af Mss Ms Mrl Mrp Mrh Mt] = pars2M_af(pars,M_reg)

if ~exist('M_reg','var')
  M_reg = eye(4);
end

Mt=eye(4,4);
Mt(1,4)=pars(1);
Mt(2,4)=pars(2);
Mt(3,4)=pars(3);

Ms=eye(4,4);
Ms(1,1)=pars(4);
Ms(2,2)=pars(5);
Ms(3,3)=pars(6);

Mrl=eye(4,4);
Mrl(2,2) = cos(pars(7));
Mrl(2,3) = sin(pars(7));
Mrl(3,2) = -sin(pars(7));
Mrl(3,3) = cos(pars(7));

Mrp=eye(4,4);
Mrp(1,1) = cos(pars(8));
Mrp(1,3) = sin(pars(8));
Mrp(3,1) = -sin(pars(8));
Mrp(3,3) = cos(pars(8));

Mrh=eye(4,4);
Mrh(1,1) = cos(pars(9));
Mrh(1,2) = sin(pars(9));
Mrh(2,1) = -sin(pars(9));
Mrh(2,2) = cos(pars(9));

Mss = eye(4,4);
Mss(1,2)=pars(10);
Mss(1,3)=pars(11);
Mss(2,3)=pars(12);

M_af = Mss*Ms*Mrl*Mrp*Mrh*Mt*M_reg;

