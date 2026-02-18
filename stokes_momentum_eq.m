function [Re] = stokes_momentum_eq(e, testsp, teste, vars)

  % Problem parameters
  mu = 4e-2;    % viscosity
  rho = 1e-3;   % density
  
  % Getting the local nodal values of velocity and pressure
  vl(:,1) = vars.vel.u (vars.vel.dm * (vars.vel.t(e,:)-1) + 1);
  vl(:,2) = vars.vel.u (vars.vel.dm * (vars.vel.t(e,:)-1) + 2);
  pl(:,1) = vars.pres.u(vars.pres.t(e,:)); 

  % evaulating the weighted sum of velocity and pressure variables at quadrature points
  v = [vars.vele.y(:,:) * vl(:,1) vars.vele.y(:,:) * vl(:,2)];
  p = vars.prese.y(:,:) * pl(:,1);

% quadratic triangle has 6 nodes, therefore 6 basis function 
% we have 12 vals in the first column as we evaluate it at all quad points
 grad_v = zeros(size(vars.vele.y, 1), 2, 2);
 for i = 1:2 
     for j = 1:2
         grad_v(:,i,j) = vars.vele.dy(:,:,j)*vl(:,i);
     end
 end

% v•(gradv)
vgrad_v = zeros(size(v));
for i = 1:2        
    vgrad_v(:,i) = v(:,1) .* grad_v(:,i,1) + v(:,2) .* grad_v(:,i,2);
end
  
% Get number of local basis functions
ne = size(testsp.t, 2);
Re = zeros(testsp.dm * ne, 1);

% Loop over test functions
for i = 1:ne
    % Index for x-component
    vei_x = testsp.dm * (i-1) + 1;
    % Index for y-component
    vei_y = testsp.dm * (i-1) + 2;
    
    % Pressure terms: -p * div(w)
 
    
    % x-component 
    % x-term 1: rho*(v•gradv)_x * w_x
    term1_x = rho * sum(teste.gw .* vgrad_v(:,1).* teste.y(:,i));
    
    % x-term 2: mu(nabla v _x : nabla w_x)
    term2_x = mu * sum(teste.gw .* (grad_v(:,1,1) .* teste.dy(:,i,1) + ...
                                     grad_v(:,1,2) .* teste.dy(:,i,2)));
    
    % x-term 3: -p(nabla • w) for x-component 
     term3_x = -sum(teste.gw .* p .* teste.dy(:,i,1));
    
    
    % y-component
    % y-term 1: rho*(v•gradv)_y * w_y 
    term1_y = rho * sum(teste.gw .* vgrad_v(:,2) .* teste.y(:,i));
    
    % y-term 2: mu(nabla v_y : nabla w_y)
    term2_y = mu * sum(teste.gw .* (grad_v(:,2,1) .* teste.dy(:,i,1) + grad_v(:,2,2) .* teste.dy(:,i,2)));
    
    % y-term 3: -p(nabla · w) for y-component 
    term3_y = -sum(teste.gw .* p .* teste.dy(:,i,2));
    
    vei = (testsp.dm * (i - 1) + 1):(testsp.dm * i);
    % computing k v . w - p div w ...
    Re(vei) = [term1_x + term2_x + term3_x, term1_y+term2_y+term3_y];
end

end
