function out = stratton_chu(E, H, X, Y, Z, n, k, E_0, rho)
    %{
    This functions calculates the farfield radiation pattern for a
    metasurface from the E and H fields over a nearfield plane parallel to
    the metasurface.  It only calculates the value for a single angle,
    given by the argument rho, and therefore needs to be called once for
    each angle that you want to calculate the farfield for.
    
    
    E - Electric field as extracted from your simulation
        Array of num_measurement_points x 3 complex double, with the
        columns being Ex, Ey and Ez.
    H - Magnetic field as extracted from your simulation
        Array of num_measurement_points x 3 complex double, with the
        columns being Hx, Hy and Hz.
    
    X - }  coordinates of each point on the measurement surface.  
    Y -  } Should be arranged as per the output of meshgrid()
    Z - }  
    
    n - The normal to the measurement surface, with positive values
        pointing away from the metasurface.
    k - Scalar wave vector
    E_0 - Background incident field strength
    rho - unit vector pointing to the observation point
    %}
    
    % Impedance of free space
    eta = 119.9169832*pi; % Ohms
    
    rho_dot_x = rho(1)*X+rho(2)*Y+rho(3)*Z;
    phase = exp(1i*k*rho_dot_x);
    
    %{ 
    % Slower version
    nE = zeros(size(E,1), size(E,2));
    rnH = zeros(size(H,1), size(H,2));
    for i=1:size(E,1)
        %nE(i,:) = cross(n, E(i,:));
        nE(i,:) = [(n(2)*E(i,3)-n(3)*E(i,2)),...
            (n(3)*E(i,1)-n(1)*E(i,3)),...
            (n(1)*E(i,2)-n(2)*E(i,1))];
        %rnH(i,:) = cross(rho, cross(n, H(i,:)));
        rnH(i,:) = [(rho(2)*(n(1)*H(i,2)-n(2)*H(i,1)) - rho(3)*(n(3)*H(i,1)-n(1)*H(i,3))),...
            (rho(3)*(n(2)*H(i,3)-n(3)*H(i,2)) - rho(1)*(n(1)*H(i,2)-n(2)*H(i,1))),...
            (rho(1)*(n(3)*H(i,1)-n(1)*H(i,3)) - rho(2)*(n(2)*H(i,3)-n(3)*H(i,2)))];
    end
    %}
    % cross product: n X E
    nE = [(n(2)*E(:,3)-n(3)*E(:,2)),...
        (n(3)*E(:,1)-n(1)*E(:,3)),...
        (n(1)*E(:,2)-n(2)*E(:,1))];
    % cross product: r X (n X H)
    rnH = [(rho(2)*(n(1)*H(:,2)-n(2)*H(:,1)) - rho(3)*(n(3)*H(:,1)-n(1)*H(:,3))),...
        (rho(3)*(n(2)*H(:,3)-n(3)*H(:,2)) - rho(1)*(n(1)*H(:,2)-n(2)*H(:,1))),...
        (rho(1)*(n(3)*H(:,1)-n(1)*H(:,3)) - rho(2)*(n(2)*H(:,3)-n(3)*H(:,2)))];

    % Convert the E and H fields into matrices which match the
    % cartesian coordiantes of the points at which the measurements
    % were taken.
    nE_matrix = zeros([size(X) 3]);
    rnH_matrix = zeros([size(X) 3]);
    for i = 1:3
        nE_matrix(:,:,i) = reshape(nE(:,i), [size(X,2), size(X,1)]);
        rnH_matrix(:,:,i) = reshape(rnH(:,i), [size(X,2), size(X,1)]);
    end
    
    I = simpsonvec((nE_matrix+eta*rnH_matrix).*phase, size(X,1), X(1,1), Y(1,1),...
        X(end,end), Y(end,end));

    
    out = 1i*k/(4*pi*E_0)*cross(rho,I);
    
end

function out = simpsonvec(F, num_points, ax, ay, bx, by)
    %{
    Use Simpson's rule to integrate the function F over 
    num_points from (ax,ay) to (bx,by).
    %}

    cx = ones(num_points,1);
    cy = ones(num_points,1);
    hx = (bx-ax)/(num_points-1);
    hy = (by-ay)/(num_points-1);

    cx(2:2:end-1) = 4;
    cx(3:2:end-2) = 2;
    cy(2:2:end-1) = 4;
    cy(3:2:end-2) = 2;

    S = cy*cx';

    out = zeros(1,3);
    for i=1:3
        out(i) = (hx*hy/9)*sum(sum(S.*F(:,:,i)));
    end
    
end
