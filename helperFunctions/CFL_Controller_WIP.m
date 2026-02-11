function CFL = CFL_Controller(CFL0, U, c)

u = U(:,:,2); 
v = U(:,:,3);
CFLMax = 0.7;

lambda = (abs(u) + c)/dx + (abs(v) + c)/dy;





end