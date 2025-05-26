%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conforming P1 code for  u_t = div(\nabla u/ \sqrt{\epsi^2 + |nabla u|^2) + f with  Neuman BC %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
format long;
%% Final times
T=1e-2;
itermax=1000;
tol=1e-6;

%% The meshes are available at https://github.com/jdroniou/HHO-Lapl-OM
meshes={'mesh1_1.mat';'mesh1_2.mat';'mesh1_3.mat';'mesh1_4.mat';'mesh1_5.mat';'mesh1_6.mat'}; 
nbmeshes=size(meshes,1);

%% Initiations
% Errors 
MAXL2error = zeros(nbmeshes, 1);
L1W11error = zeros(nbmeshes, 1);
L2W11error = zeros(nbmeshes, 1);
MAXW11error = zeros(nbmeshes, 1);

MAXL2norm = zeros(nbmeshes, 1);
L1W11norm = zeros(nbmeshes, 1);
L2W11norm = zeros(nbmeshes, 1);
MAXW11norm = zeros(nbmeshes, 1);

% Order of convergence 
ocMAXL2error = zeros(nbmeshes - 1, 1);
ocL1W11error = zeros(nbmeshes - 1, 1);
ocL2W11error = zeros(nbmeshes - 1, 1);


%% Test case number 
ucase = 1;

% Min relax
relaxmin = 1e-5;

%% Loop over each mesh in the sequence
for imesh=1:nbmeshes
    %% Load mesh here!
    loadmesh=strcat('load matlab_meshes/', meshes{imesh});
	disp(loadmesh);
    eval(loadmesh);
    disp('mesh loaded');

    %% Compute real centers of mass and mesh size
    cg=gravity_centers(ncell, cell_v, vertex, area);
    h(imesh)=max(abs(diam));

    %% Time steps k = O(h^2)
%    Ndt(imesh) = ceil(T/h(imesh));
    Ndt(imesh) = 5e2;
    dt = T/Ndt(imesh);

    %% Epsilon
    % epsilon = sqrt(h(imesh)); % epsilon is oder of root h
    epsilon = min(dt, h(imesh));%1e-6;

    %% Initial condition
    U_pre = test_cases(0,vertex, ucase)';
 
    %% Assemble Mass and Stifness Matrices
    M = assemble_mass_system(area, ncell, nvert, cell_v);

    %% Print the scheme solution and exact solution at initial condition to view in Paraview
    if imesh == nbmeshes 
    write_solution_vtk(U_pre, strcat('VTKout/p1_c_diffusion_solution0'), ncell, nedge, nvert, cell_v, cell_n, cell_e, vertex);
    write_solution_vtk(test_cases(0,vertex, ucase)', strcat('VTKout/exact_solution0'), ncell, nedge, nvert, cell_v, cell_n, cell_e, vertex);
    end
     
    %% Error norms initiation
    L2error = zeros(Ndt(imesh), 1);
    W11error = zeros(Ndt(imesh), 1); 
    L2norm = zeros(Ndt(imesh),1);
    W11norm = zeros(Ndt(imesh),1);

    %% Time stepping starts here!
    disp('time starts!')
    ITER = 0;
    Res = 0;
    for idt = 1 : Ndt(imesh)
        b = assemble_source(cell_v, ncell, nvert, area, cg, idt * dt, epsilon, ucase); 
        rhs = M * U_pre + dt * b;
        iter = 0;
        res = 1;

        %% Newton
        X_pre = U_pre; %test_cases(idt*dt,vertex, ucase)'; % initial guess of the solution
        relax = 0.1;
        res_pre = 1.;
        while ( iter < itermax && res > tol && relax > relaxmin )
            iter;
            [Jacobian, Fu] = assemble_jacobian_diffusion(cell_v,ncell,nvert,vertex, X_pre, epsilon);
            % S = assemble_diffusion_system(cell_v,ncell,nvert,vertex, X_pre, epsilon);
            Aglob = M + dt .* Jacobian;
            RHS = rhs - ( M * X_pre + dt .* Fu);
            delta_X = Aglob\RHS;

            %% Residues
            res = norm(RHS, inf); % residual of nonlinear eq at previous step: ||F(X_pre) - rhs||
            % fprintf("res %e, relax %e, idt/Ndt %d/%d\n", res, relax, idt, Ndt(imesh));
            
            if (iter == 1 | res < 1.5 * res_pre)
                iter = iter + 1;  
                X = X_pre +  relax * delta_X; % under-relaxation
                X_pre = X; %current Newton iteration solution is updated
                res_pre = res;
            else
                relax = relax/2;
            end
        end % end nonlinear iterations
        if (iter==itermax | relax <= relaxmin)
            res
            iter
            relax
            error('no convergence')
        end
        % Newton loop end: X outputed here is the solution at the current time step idt*dt
        iter;
        ITER = ITER + iter;
        Res = Res + abs(res);
        U_pre = X; % current time solution is retireing

        %% Errors and norm of the exact solution
        % [L2error(idt), W11error(idt)] = compute_norms(X - test_cases(idt * dt, vertex, ucase)', area, vertex, cell_v);
        % [L2norm(idt), W11norm(idt)] = compute_norms(test_cases(idt * dt, vertex, ucase)', area, vertex, cell_v);

        %% Norm of the approximate solution
        [L2norm(idt), W11norm(idt)] = compute_norms(X, area, vertex, cell_v);

        %% Create files for visualising the approximate and exact solution 
        if imesh == nbmeshes
        write_solution_vtk(X,strcat('VTKout/p1_c_diffusion_solution',num2str(idt)),ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex); % Print the scheme solution at current time in Paraview
        write_solution_vtk(test_cases(idt*dt, vertex, ucase)',strcat('VTKout/exact_solution',num2str(idt)),ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex); % Print the exact solution at current time in Paraview
        end

    end 
    
    %% Time norms: L infinity, L1 norm, and L2 norm
    MAXL2error(imesh) = max(L2error);
    L1W11error(imesh) = dt * sum(W11error);
    L2W11error(imesh) = sqrt(dt * sum(W11error.^2)); 

    MAXL2norm(imesh) = max(L2norm);
    L1W11norm(imesh) = dt * sum(W11norm);
    L2W11norm(imesh) = sqrt(dt * sum(W11norm.^2));
    MAXW11norm(imesh) = max(W11norm);

    %% Relative errors
    MAXL2error(imesh) = MAXL2error(imesh)/MAXL2norm(imesh);
    L1W11error(imesh) = L1W11error(imesh)/L1W11norm(imesh);
    L2W11error(imesh) = L2W11error(imesh)/L2W11norm(imesh);

end 

%% Orders of convergence
ocL2error=zeros(nbmeshes-1,1);
ocH10error=zeros(nbmeshes-1,1);
for imesh = 1:nbmeshes-1
    ocMAXL2error(imesh)=log(MAXL2error(imesh)/MAXL2error(imesh+1)) / log(h(imesh)/h(imesh+1));
    ocL1W11error(imesh)=log(L1W11error(imesh)/L1W11error(imesh+1)) / log(h(imesh)/h(imesh+1));
    ocL2W11error(imesh)=log(L2W11error(imesh)/L2W11error(imesh+1)) / log(h(imesh)/h(imesh+1));
end 
% display list of order of convergence
% MAXL2error
% L1W11error;
% L2W11error
% ocMAXL2error
% ocL1W11error;
% ocL2W11error

% Plot of the MAXL2 error of u and L1W11 error Gradu with respect of h in a logarithmic scale
% subplot(2,1,1);
plot(h,MAXL2norm,'*-');
axis([0 0.5 0 1]);
hold on
plot(h,L1W11norm,'*-');
% loglog(h,h); %compute the vector in order to have a slope equal to 1
title('Apriori estimate','Interpreter','latex')
xlabel('mesh size $h$','Interpreter','latex')
% ylabel('$\left\|u\right\|_{L^{\infty}(L^2)}$','Interpreter','latex')
legend('$\left\|u\right\|_{L^{\infty}(L^2)}$', '$\left\| \nabla u \right\|_{L^{1}(L^1)}$', 'interpreter', 'latex','Location','southeast');
grid on

% subplot(2,1,2);
% loglog(h,L1W11norm,'*-');
% axis([0.01 1 0.0001 1]);
% hold on
% loglog(h,h);
% title('Apriori estimates','Interpreter','latex')
% xlabel('$h$','Interpreter','latex')
% ylabel('$\left\| \nabla u \right\|_{L^{\infty}(L^1)}$','Interpreter','latex')
% legend('L^1W^{1,1} norm','Slope equal to 1','Location','southeast')
% grid on

