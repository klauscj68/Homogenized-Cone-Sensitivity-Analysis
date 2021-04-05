function [Rst0_sig] = Rst_0(pts_eul,...
                            R_b,R_t,H,...
                            nu,epsilon_0,...
                            n_Rst0)
%Compute the initial concentration of Rst_sig for time t = 0
%   pts_eul are the spatial nodes of the eulerian cone mesh.
%   The second two rows of inputs are geometric parameters defining the 
%    cone 
%   n_Rst0 is defined by the data sheet.  If a number, a uniform volumic
%    concentration of Rst0 is given so that its integral over the volume is
%    n_Rst0. If it takes the value NaN, then the user-defined space varying
%    function in section "Nonuniform concentration spatial distribution for
%    Rst0_sig" is used below

%% Initial Parameter
n_pts = size(pts_eul,2);

switch isnan(n_Rst0)
    case true
        %% Nonuniform concentration spatial distribution for Rst0_sig
        % 10 single photon isomerization on 10 equispaced discs
        eta0 = .5*nu*epsilon_0;
        
        % Each column is x-y-z coordinate of intended photon
        n_photons = 10;
        %  Photons deposited from [.3*H ... .75*H] in units of .05*H 
        Phi = [zeros(2,n_photons);
               H*(.3 + .05*(cumsum(ones(1,n_photons)) - 1))];
           
        photon_wt = [9 9.5 10 10.4 10.8 11.2 11.6 12 12.3 12.5];
        
        % If you know the number of distinct z's in pts_eul then that
        % divides into size of pts_eul according to how many pts are in a
        % level. Then at each horizontal level find the three nodes closest
        % to the horizontal projection of Phi.  Join it to the three nodes
        % periodically paired to it above or below.
        
        z_s = unique(pts_eul(3,:));
        n_z = size(z_s,2);
        n_pts = size(pts_eul,2);
        
        assert(mod(n_pts,n_z) == 0,'n_z does not divide n_pts');
        n_cross = n_pts/n_z;
        
        % Find first z_s index to be at or above photon
        for i=n_photons:-1:1
           photon_z(i) = find(Phi(3,i) <= z_s,1);  
        end
        
        % Find closest 3 nodes at the level to photons horizontal
        % projection
        for i=n_photons:-1:1
            ram = pts_eul(1:2,...
                          n_cross*(photon_z(i) -1)+1:n_cross*photon_z(i));
            
            [~,ram] = sort(vecnorm(ram - Phi(1:2,i)));
            ram = ram(1:3) + n_cross*(photon_z(i) -1);
            ram = ram';
            Dirac(:,i) = [ram - n_cross;...
                          ram];
        end
        
        % Create Rst0_sig vector
        Rst0_sig = zeros(n_pts,1);
        for j=1:n_photons
            Rst0_sig(Dirac(:,j)) = photon_wt(j)*1;
        end
        
        Rst0_sig = eta0*Rst0_sig;
        
    case false
         %% Compute a uniform distribution with n_Rst0 photons
        % Cone Volume
        V = pi * H * (-R_t + R_b) ^ 2 / 0.3e1 - ...
            pi * H * R_b * (-R_t + R_b) + pi * H * R_b ^ 2;
        
        % Volume to surface conversion factor
        eta0 = .5*nu*epsilon_0;
        
        % Uniform volumic density
        Rst0_vol = n_Rst0/V;
        
        % Uniform surface density
        Rst0_sig = eta0*Rst0_vol;
        Rst0_sig = repmat(Rst0_sig,n_pts,1);
        
end

end