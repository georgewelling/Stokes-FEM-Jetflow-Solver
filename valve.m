
  load("jetflow/jetflow3.mat")
  quad=10;
  etype = 'triangle';
  nonlinear_tol = 1e-6; % setting the nonlinear tolerance
  iter = 0; % initializing the interation counter
  residual_vector_norm = 1; % initializing the residual vector norm

%% Making the remanining Domain function space
   Omega2.e=fem_get_basis(Omega2.p, quad, etype);
   Omega.e=fem_get_basis(Omega.p, quad, etype);

%% Making velocity function space
  vel = Omega2;  
  vel.name = 'vel';
  vel.e=fem_get_basis(vel.p, quad, etype); 

% Initializing the coefficients for the velocity field to zero (no nodes x number of unknowns per node)
  vel.u = zeros(vel.dm * size(vel.x,1), 1);
  n_vel = size(vel.u,1); 
%% Making pressure function space 
  pres = Omega;
  pres.name = 'pres'; 
  pres.dm = 1; 
  pres.p  = 1; 
  pres.e=fem_get_basis(pres.p, quad, etype); 

% Initializing the coefficients for the pressure field to zero (no nodes x number of unknowns per node)
  pres.u = zeros(size(pres.x,1),1); 
  n_pres = size(pres.u,1);
%% Starting the Newton-Raphson iteration 
for k = 0.05:0.05:1
while (iter == 0) || (residual_vector_norm > nonlinear_tol)
  % Making a FEM object storing my unknown variables (must do this every iteration) 
    vars.vel = vel;
    vars.pres = pres;
    % Making the residual vector and applying boundary conditions (velocity)
    R1 = fem_assemble_block_residual(@stokes_momentum_eq,Omega,vel,vars);
    for i = 1:size(vel.b,1)
    for j = 2:4
        nn = vel.dm * (vel.b(i,j) - 1);
        switch vel.b(i,5)
            case 1  % inlet
                R1(nn+1:nn+2) = -[(-k*200) - vel.u(nn+1); 0 - vel.u(nn+2)];
            case 3  % side outlets
                if vel.x(vel.b(i,j),2) > 0
                    R1(nn+1:nn+2) = -[0 - vel.u(nn+1);  (k*10) - vel.u(nn+2)];
                else
                    R1(nn+1:nn+2) = -[0 - vel.u(nn+1); (-k*10) - vel.u(nn+2)];
                end
            case 4  % walls
                R1(nn+1:nn+2) = -[0 - vel.u(nn+1); 0 - vel.u(nn+2)];
            % case 2 (outlet): skip, natural BC
        end
     end
    end


  % Making the residual vector and applying boundary conditions (on species B) 
    R2 = fem_assemble_block_residual(@stokes_mass_eqn, Omega, pres, vars);
  % Making the global residual vector + computing the norm 
    R = [R1; R2];
    residual_vector_norm = norm(R,2);
    disp(['Current residual (iter=' num2str(iter) '): ' num2str(residual_vector_norm)])
    if(residual_vector_norm < nonlinear_tol)
      continue
    end

  % Creating the discrete matrix operators corresponding to operators a, b, and c
    disp('        constructing the Jacobian blocks ...')
    A  = fem_assemble_block_matrix_perturbation(@stokes_momentum_eq, Omega, vel, vel, vars); 
    B  = fem_assemble_block_matrix_perturbation(@stokes_momentum_eq, Omega, vel, pres, vars); 
    C  = fem_assemble_block_matrix_perturbation(@stokes_mass_eqn, Omega, pres, vel, vars); 
    D = sparse(size(pres.u,1),size(pres.u,1));

  % Editing block matrices for dirichlet conditions 
    for i = 1:size(vel.b,1)
      if(vel.b(i,5) ==1 || vel.b(i,5) == 3 || vel.b(i,5) == 4)
          for j = 2:4
              nn = vel.dm * (vel.b(i,j) - 1);
              A(nn+1:nn+2,:) = 0;
              A(nn+1:nn+2, nn+1:nn+2) = eye(2);
              B(nn+1:nn+2,:) = 0;
          end
      end
    end

  % Composing the Jacobian from our Jacobian blocks
    disp(['        assembly of the Global Jacobian ...'])
    J = [ A B; C D];

  % Apply Newton Raphson 
    disp(['        solving for the NR update ...'])
    U = J \ R;
    vel.u(1:n_vel) = vel.u(1:n_vel) - U(1:n_vel);
    pres.u(1:n_pres) = pres.u(1:n_pres) - U(n_vel+1:end);

  % Update the iteration counter 
    iter = iter + 1;

    disp(['  ']) 
end
end


%% Plots 
figure(1);
quiver(vel.x(:,1), vel.x(:,2), vel.u(1:vel.dm:n_vel), vel.u(2:vel.dm:n_vel), 2); % the '2' scales arrow size
title('Velocity field');
xlabel('x [mm]');
ylabel('y [mm]');
axis equal tight;
grid on;

figure(2);
trisurf(pres.t, pres.x(:,1), pres.x(:,2), pres.u,'Facecolor','interp','LineStyle','none')
view(2);                 
axis equal tight;
colorbar;
colormap(turbo);
title('Pressure field [mmHg]');
xlabel('x [mm]');
ylabel('y [mm]');



%% Pressure drop calculations 
inlet_nodes = unique(pres.b(pres.b(:,4)==1,2:3));
mean_inlet = mean(pres.u(inlet_nodes));

outlet_nodes = unique(pres.b(pres.b(:,4)==2,2:3));
mean_outlet = mean(pres.u(outlet_nodes));

pressure_diff = mean_inlet - mean_outlet; 
fprintf("Results:\n");
fprintf("Pressure difference: %.4f mmHg\n", pressure_diff);

