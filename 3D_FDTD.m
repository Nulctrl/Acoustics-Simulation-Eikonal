c = 343; %speed of sound
fmax = 1000; %Hz
PPW = 6; %points per wavelength at fmax
duration = 0.1; %seconds
refl_coeff = 0.9; %reflection coefficient

Bx = 10; By = 8; Bz = 6 ; 
x_in = Bx*0.3; y_in = By*0.5; z_in = Bz*0.5; %source

draw = true; %to plot or not
apply_rigid = true; %apply rigid boundaries
apply_loss = true; %apply loss 

if (draw)
   %a mask convenient for plotting
   draw_mask = NaN*in_mask;
   draw_mask(in_mask) = 1;
end

if (apply_loss)
   assert(apply_rigid);
end

%Grid spacing, time step, sample rate
dx = c/fmax/PPW;
dt = sqrt(1/3)*dx/c;
SR = 1/dt;
fprintf('sample rate = %.3f Hz\n',SR) 
fprintf('Δx = %.5f m \n',dx) 

Nx = ceil(Lx/dx)+2;
Ny = ceil(Ly/dx)+2; 
Nz = ceil(Lz/dx)+2; 
Nt = ceil(duration/dt); 

xv = [0:Nx-1]*dx-0.5*dx;
yv = [0:Ny-1]*dx-0.5*dx;
zv = [0:Nz-1]*dx-0.5*dx; 
[X,Y,Z] = ndgrid(xv,yv,zv);

in_mask = false(Nx,Ny,Nz);
in_mask(X(:)>=0 & Y(:)>=0 & Z(:)>=0 & X(:)<Bx & Y(:)<By & Z(:)<Bz) = true;

if (apply_rigid)
   K_map = zeros(Nx,Ny,Nz);
   K_map(2:Nx-1,2:Ny-1,2:Nz-1) = K_map(2:Nx-1,2:Ny-1,2:Nz-1) + in_mask(3:Nx,2:Ny-1,2:Nz-1);
   K_map(2:Nx-1,2:Ny-1,2:Nz-1) = K_map(2:Nx-1,2:Ny-1,2:Nz-1) + in_mask(1:Nx-2,2:Ny-1,2:Nz-1);
   K_map(2:Nx-1,2:Ny-1,2:Nz-1) = K_map(2:Nx-1,2:Ny-1,2:Nz-1) + in_mask(2:Nx-1,3:Ny,2:Nz-1);
   K_map(2:Nx-1,2:Ny-1,2:Nz-1) = K_map(2:Nx-1,2:Ny-1,2:Nz-1) + in_mask(2:Nx-1,1:Ny-2,2:Nz-1);
   K_map(2:Nx-1,2:Ny-1,2:Nz-1) = K_map(2:Nx-1,2:Ny-1,2:Nz-1) + in_mask(2:Nx-1,2:Ny-1,3:Nz);
   K_map(2:Nx-1,2:Ny-1,2:Nz-1) = K_map(2:Nx-1,2:Ny-1,2:Nz-1) + in_mask(2:Nx-1,2:Ny-1,1:Nz-2);
   K_map(~in_mask) = 0;
   ib = find(K_map(:)>0 & K_map(:)<6);
   Kib = K_map(ib);
   clear K_map;
end

u0 = zeros(Nx,Ny,Nz);
u1 = zeros(Nx,Ny,Nz);
u2 = zeros(Nx,Ny,Nz);
u3 = zeros(Nx,Ny,Nz);

%set up an excitation signal
u_in = zeros(Nt,1);
Nh = ceil(5*SR/fmax);
u_in(1:Nh) = 0.5-0.5*cos(2*pi*(0:Nh-1)'./Nh);
u_in(1:Nh) = u_in(1:Nh).*sin(2*pi*(0:Nh-1)'./Nh);

inx = round(x_in/dx+0.5)+1;
iny = round(y_in/dx+0.5)+1;
inz = round(z_in/dx+0.5)+1;
assert(in_mask(inx,iny,inz));

if (apply_loss)
   assert(abs(refl_coeff)<=1.0);
   g = (1-refl_coeff)/(1+refl_coeff);
   lf = 0.5*sqrt(0.5)*g; 
end
bb = 0;

%fdtd update
for nt=0:Nt+1
   u0(2:Nx-1,2:Ny-1,2:Nz-1) = in_mask(2:Nx-1,2:Ny-1,2:Nz-1).*((u1(3:Nx,2:Ny-1,2:Nz-1) + u1(1:Nx-2,2:Ny-1,2:Nz-1) + u1(2:Nx-1,3:Ny,2:Nz-1) + u1(2:Nx-1,1:Ny-2,2:Nz-1) + u1(2:Nx-1,2:Ny-1,3:Nz) + u1(2:Nx-1,2:Ny-1,1:Nz-2))/3 - u2(2:Nx-1,2:Ny-1,2:Nz-1));
   if (apply_rigid)
     u0(ib) = (2-Kib/3).*u1(ib) + (u1(ib+1) + u1(ib-1) + u1(ib+Nx) + u1(ib-Nx) + u1(ib+Nx*Ny) + u1(ib-Nx*Ny))/3 - u2(ib);
       if (apply_loss)
          u0(ib) = (u0(ib) + lf*(6-Kib).*u2(ib))./(1+lf.*(6-Kib));
       end
   end

   u0(inx,iny,inz) = u0(inx,iny,inz) + u_in(nt+1);
   
   if(nt+1==2)
       disp(u0(inx,iny,inz));
   end
   
   %plotting
   if (draw)
         u1g = squeeze(u1(inx,:,:));
      if nt==0
         figure('name','float');
         u_draw = (u1g.*(draw_mask)).';
         hh = imagesc(xv,yv,u_draw,'cdatamapping','scaled');
         set(gca,'ydir','normal');
         axis equal;
         xlabel('x');
         ylabel('y');
         colorbar;
         colormap(bone);
         xlim([min(xv) max(xv)]);
         ylim([min(yv) max(yv)]);
      else
         umax = max(abs(u1g(:)))+eps;
         ustart=u1g(inx,iny);
         u_draw = (u1g.*(draw_mask)).';
         set(hh,'cdata',u_draw);
         clim([-umax umax]);
         drawnow;
      end
   end

   ustart=u1(inx,iny,inz);
   u2 = u1; 
   u1 = u0;
end

% This 3D implementation is a modification based on Brian Hamilton's 2D FDTD tutorial