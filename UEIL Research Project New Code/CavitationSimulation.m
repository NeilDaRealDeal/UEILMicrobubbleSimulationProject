% To run a simulation, first set up the regular K-wave simulation,
% then plug the parameters into this method along with the microbubble
% location and parameters.

% Requirements: 
%  - KWaveArray Class for Initial Source is Used
%  - <FaH>
% Assumptions:
%  - <FaH>
% Note: Don't Forget to Employ Iteration
% Note: Emissions are in the form of time vs. index (of microbubble or of cube)
% Note: After coding this class, code the example (based off of the bowl trasnducer k-wave example)
% Note: The simulation is not deterministic: the microbubble positions are
% generated with random noise
% Note to Self: CODE ANALYTIC SOLUTION AND KWAVE ON GPU AS WELL IN ADDITION
% TO BUBBLESIM ON GPU
% NOTE TO SELF: BS_ModifiedRaleigh and BS_SimulateOscillation UEILFotisNeil
% functions and reset their original scripts, updating their references as
% well and defining more UEILFotisNeil functions as necessary!!!
% Lots of room for optimization + correction!!!
% NOTE TO SELF: Once you finish, check that the simulation and experiments
% work when using a different device


classdef CavitationSimulation < handle
   properties
      %Microbubble positional parameters
      microbubble_count
      microbubble_grid_locations_x %Generated (pixel address)
      microbubble_grid_locations_y %Generated (pixel address)
      microbubble_grid_locations_z %Generated (pixel address)
      microbubble_positions_x %Generated (m)
      microbubble_positions_y %Generated (m)
      microbubble_positions_z %Generated (m)
      microbubbles_in_each_pixel %Generated (array of integers corresponding to the ordered microbubbles for each pixel in order)
      microbubbles_per_pixel %Generated (array of integers coressponding to the number of microbubbles per pixel)
      measurement_cube_pixelated_sidelength %Odd
      measurement_cube_grid_locations_x %Generated (for the center)
      measurement_cube_grid_locations_y %Generated (for the center)
      measurement_cube_grid_locations_z %Generated (for the center)
      mb_fluid_region_matrix %Nx by Ny by Nz binary matrix where 1s indicate the region
      sphere_approximation_pixelated_radius
      microbubble_cube_based_locations %Generated (array of integers corresponding to the ordered cubes for each microbubble in order)
      microbubbles_in_each_cube %Generated (array of structs of arrays of integers corresponding to the ordered microbubbles for each cube in order)
      microbubbles_per_cube %Generated (array of integers corresponding to the number of microbubbles per cube)

      %Individual microbubble parameters
      gas_model
      radius
      thickness
      shear
      viscocity
      liquid_name
      artifical_scale_up

      %Simulation parameters
      kgrid
      medium
      true_kwavearray %Remove later if possible
      true_source_signal %Remove later if possible
      true_source
      true_sensor
      input_args
      source_amp
      current_sim_source %Generated
      sim_sensor %Generated
      bubblesim_driven_dt %Note: Scattered dt sampling rate should be adjusted accordingly to kgrid.dt
                          %Note: We assme this to be smaller than kgrid.dt
      %Results
      microbubble_radii_by_cube %Generated
      microbubble_scattered_pulses_by_cube %Generated
      microbubble_radii %Generated
      microbubble_scattered_pulses %Generated
      sensor_data %Generated
   end
   methods
       function generate_grid_cube_positions(obj) %Note: May leave some edge pixels un-cubed
           X = obj.kgrid.Nx;
           Y = obj.kgrid.Ny;
           Z = obj.kgrid.Nz;
           s = obj.measurement_cube_pixelated_sidelength;
           position_matrix_x = zeros(floor(X/s), floor(Y/s), floor(Z/s));
           position_matrix_y = zeros(floor(X/s), floor(Y/s), floor(Z/s));
           position_matrix_z = zeros(floor(X/s), floor(Y/s), floor(Z/s));
           for i = 1:floor(X/s)
               for j = 1:floor(Y/s)
                   for k = 1:floor(Z/s)
                       position_matrix_x(i, j, k) = s*(i-1) + (s+1)/2;
                       position_matrix_y(i, j, k) = s*(j-1) + (s+1)/2;
                       position_matrix_z(i, j, k) = s*(k-1) + (s+1)/2;
                   end
               end
           end
           obj.measurement_cube_grid_locations_x = zeros(numel(position_matrix_x), 1);
           obj.measurement_cube_grid_locations_y = zeros(numel(position_matrix_y), 1);
           obj.measurement_cube_grid_locations_z = zeros(numel(position_matrix_z), 1);
           for i = 1:numel(position_matrix_x)
               obj.measurement_cube_grid_locations_x(i) = position_matrix_x(i);
               obj.measurement_cube_grid_locations_y(i) = position_matrix_y(i);
               obj.measurement_cube_grid_locations_z(i) = position_matrix_z(i);
           end
       end
       function set_up_sim_sensor_and_cubes(obj)
            obj.generate_grid_cube_positions();
            obj.sim_sensor = zeros(obj.kgrid.Nx, obj.kgrid.Ny, obj.kgrid.Nz);
            for i = 1:numel(obj.measurement_cube_grid_locations_x)
                obj.sim_sensor(obj.measurement_cube_grid_locations_x(i), obj.measurement_cube_grid_locations_y(i), obj.measurement_cube_grid_locations_z) = 1;
            end
       end
       function data = run_sensing_sim(obj) %return data in terms of a matrix with a time and grid cube index axes
           obj.set_up_sim_sensor_and_cubes();
           data = kspaceFirstOrder3D(obj.kgrid, obj.medium, obj.true_source, obj.sim_sensor, obj.input_args{:}, ...
            'DataCast', 'gpuArray-single', ...
            'PlotScale', [-1, 1] * obj.source_amp);
       end
       function generate_microbubble_possandlocs(obj) %grid with small noise, entirely random (DOUBLE CHECK THIS FUNCTION)
           %Note: This should set every generated quantity in the
           %microbubble positional parameters, except for the grid cube
           %positions.
           sum = 0;
           for i = 1:numel(obj.mb_fluid_region_matrix)
               if (obj.mb_fluid_region_matrix(i))
                    sum = sum + 1;
               end
           end
           region_pixels_indices = zeros(1, sum);
           region_pixels_locs = zeros(3, sum);
           
           count = 0;
           for i = 1:numel(obj.mb_fluid_region_matrix)
               if (obj.mb_fluid_region_matrix(i))
                    count = count + 1;
                    region_pixels_indices(count) = i;
                    region_pixels_locs(:, count) = CavitationSimulation.get_coords(i, obj.kgrid.Nx, obj.kgrid.Ny);
               end
           end
           
           obj.microbubble_grid_locations_x = zeros(1, obj.microbubble_count); %
           obj.microbubble_grid_locations_y = zeros(1, obj.microbubble_count); %
           obj.microbubble_grid_locations_z = zeros(1, obj.microbubble_count); %
           obj.microbubble_positions_x = zeros(1, obj.microbubble_count);
           obj.microbubble_positions_y = zeros(1, obj.microbubble_count);
           obj.microbubble_positions_z = zeros(1, obj.microbubble_count);
           obj.microbubbles_per_pixel = zeros(1, obj.kgrid.Nx*obj.kgrid.Ny*obj.kgrid.Nz);
           obj.microbubble_cube_based_locations = zeros(1, obj.microbubble_count);
           obj.microbubbles_per_cube = zeros(1, numel(obj.measurement_cube_grid_locations_x));
    
           for i = 1:obj.kgrid.Nx*obj.kgrid.Ny*obj.kgrid.Nz
               obj.microbubbles_in_each_pixel(i).indices = [];
           end
           for i = 1:numel(obj.measurement_cube_grid_locations_x)
               obj.microbubbles_in_each_cube(i).indices = [];
           end
           
           rand_ints = randi([1 sum], 1, obj.microbubble_count);
            
           for i = 1:obj.microbubble_count %REVIEW THIS
               rand_int = rand_ints(i);
               pixel_loc = region_pixels_locs(:, rand_int);
               pixel_ind = region_pixels_indices(rand_int);
               obj.microbubble_grid_locations_x(i) = pixel_loc(1);
               obj.microbubble_grid_locations_y(i) = pixel_loc(2);
               obj.microbubble_grid_locations_z(i) = pixel_loc(3);
               obj.microbubbles_in_each_pixel(pixel_ind).indices(numel(obj.microbubbles_in_each_pixel(pixel_ind).indices)+1) = i;
               obj.microbubbles_per_pixel(pixel_ind) = obj.microbubbles_per_pixel(pixel_ind) + 1;
               cube_coindex = [floor((pixel_loc(1)-1)/obj.measurement_cube_pixelated_sidelength)+1, floor((pixel_loc(2)-1)/obj.measurement_cube_pixelated_sidelength)+1, floor((pixel_loc(3)-1)/obj.measurement_cube_pixelated_sidelength)+1];
               cube_index = CavitationSimulation.get_ind(cube_coindex(1), cube_coindex(2), cube_coindex(3), floor(obj.kgrid.Nx/obj.measurement_cube_pixelated_sidelength), floor(obj.kgrid.Ny/obj.measurement_cube_pixelated_sidelength));
               obj.microbubble_cube_based_locations(i) = cube_index;
               obj.microbubbles_in_each_cube(cube_index).indices(numel(obj.microbubbles_in_each_cube(cube_index).indices) + 1) = i;
               obj.microbubbles_per_cube(cube_index) = obj.microbubbles_per_cube(cube_index) + 1;
               obj.microbubble_positions_x(i) = obj.kgrid.x(pixel_loc(1), pixel_loc(2), pixel_loc(3)) + (rand()-0.5)*obj.kgrid.dx;
               obj.microbubble_positions_y(i) = obj.kgrid.y(pixel_loc(1), pixel_loc(2), pixel_loc(3)) + (rand()-0.5)*obj.kgrid.dy;
               obj.microbubble_positions_z(i) = obj.kgrid.z(pixel_loc(1), pixel_loc(2), pixel_loc(3)) + (rand()-0.5)*obj.kgrid.dz;
           end
       end
       function data = calculate_cubical_emissions_from_cubical_received(obj, cubical_received)
           %Emissions per microbubble coming from each cube, not scaled by the number of microbubbles in the cube
           %Also doesn't account for inverse square
           kwave_dt = obj.kgrid.dt;
           bubblesim_dt_driv = obj.bubblesim_driven_dt;
           sampling_rate = 1e6/kwave_dt;
           cube_count = size(cubical_received, 1);
           time = 1:(size(cubical_received, 2) * kgrid_dt/bubblesim_dt_driv);
           received_for_bubblesim = zeros(cube_count, length(time));
           for i = 1:cube_count
               received_for_bubblesim(i, :) = interp(cubical_received(i, :), length(time)/size(cubical_received, 2));
           end
           %Bubblesim Simulation
           returnData = zeros(size(cubical_received));
           for i = 1:cube_count
                [particle2,pulse2,linear2,simulation2,graph2] = UEILFotisNeilBubblesimCallBackCustomPulse(obj.gas_model, obj.radius, obj.thickness, obj.shear, obj.viscocity, obj.liquid_name, '30', '2.5', '10', string(sampling_rate), time.', received_for_bubblesim(i, :).');
                returnData(i, :) = simulation2.pr;
           end
           returnData = returnData * obj.artifical_scale_up;
           data = returnData;
       end
       function data = calculate_cubical_radii_from_cubical_received(obj, cubical_received)
           kwave_dt = obj.kgrid.dt;
           bubblesim_dt_driv = obj.bubblesim_driven_dt;
           sampling_rate = 1e6/kwave_dt;
           cube_count = size(cubical_received, 1);
           time = 1:(size(cubical_received, 2) * kgrid_dt/bubblesim_dt_driv);
           received_for_bubblesim = zeros(cube_count, length(time));
           for i = 1:cube_count
               received_for_bubblesim(i, :) = interp(cubical_received(i, :), length(time)/size(cubical_received, 2));
           end
           %Bubblesim Simulation
           returnData = zeros(size(cubical_received));
           for i = 1:cube_count
                [particle2,pulse2,linear2,simulation2,graph2] = UEILFotisNeilBubblesimCallBackCustomPulse(obj.gas_model, obj.radius, obj.thickness, obj.shear, obj.viscocity, obj.liquid_name, '30', '2.5', '10', string(sampling_rate), time.', received_for_bubblesim(i, :).');
                returnData(i, :) = interp(simulation2.a(:, 1), length(time)/size(simulation2.a, 1)).';
           end
           %Data Processing of Output <FaH>
           data = returnData;
           %Note: Account for uneven sampling
       end
       function data = calculate_bubble_emissions_from_cubical_emissions(obj, cubical_emissions) %smoothing, returns row for each microbubble
           bubble_emissions = zeros(obj.microbubble_count, size(cubical_emissions, 2));
           cube_count = size(cubical_emissions, 1);
           for i = 1:cube_count
               for j = 1:numel(obj.microbubbles_in_each_cube(i).indices)
                   bubble_emissions(obj.microbubbles_in_each_cube(i).indices(j)) = cubical_emissions(i, :);
               end
           end    
           data = bubble_emissions;
       end
       function source = calculate_source_from_microbubble_emissions(obj, bubble_emissions) %input 1 row per microbubble
           %Account for large sphere approximation, inverse square law squarerooted, 
           %and sum over pixels with microbubble count scale ups
           pixel_emissions = zeros(obj.kgrid.Nx*obj.kgrid.Ny*obj.kgrid.Nz, size(bubble_emissions, 2));
           for i = 1:obj.kgrid.Nx*obj.kgrid.Ny*obj.kgrid.Nz
               pixel_emissions(i, :) = obj.microbubbles_per_pixel * bubble_emissions(obj.microbubbles_in_each_pixel(i).indices(1), :);
           end
           sourc.p_mask = obj.true_source.p_mask;
           original_p = obj.true_source.p;
           for i = 1:obj.kgrid.Nx
               for j = 1:obj.kgrid.Ny
                   for k = 1:obj.kgrid.Nz
                       if (sourc.p_mask(i, j, k)) 
                           continue
                       end
                       for i0 = (i-obj.sphere_approximation_pixelated_radius):(i+obj.sphere_approximation_pixelated_radius)
                           for j0 = (j-obj.sphere_approximation_pixelated_radius):(j+obj.sphere_approximation_pixelated_radius)
                               for k0 = (k-obj.sphere_approximation_pixelated_radius):(k+obj.sphere_approximation_pixelated_radius)
                                    if (((i-i0)^2 + (j-j0)^2 + (k-k0)^2 <= obj.sphere_approximation_pixelated_radius^2)&&~all(pixel_emissions(cavitationSimulation.getInd(i0, j0, k0, obj.kgrid.Nx, obj.kgrid.Ny), :) == 0)) 
                                        sourc.p_mask(i, j, k) = 1;
                                    end
                               end
                           end
                       end
                   end
               end
           end
           mask_ind = find(sourc.p_mask);
           original_found = find(obj.true_source.p_mask);
           scaling = (str2double(obj.radius)*(1e-6))/(obj.sphere_approximation_pixelated_radius * obj.kgrid.dx); %Assume dx = dy = dz
           for i = 1:length(mask_ind)
               pressures = zeros(1, size(bubble_emissions, 2));
               index = mask_ind(i);
               indices = CavitationSimulation.get_coords(index, obj.kgrid.Nx, obj.kgrid.Ny);
               if (obj.true_source.p_mask(indices(1), indices(2), indices(3)) == 1)
                   pressures = pressures + original_p(find(original_found == index), :);
               end
               i = indices(1);
               j = indices(2);
               k = indices(3);
               for i0 = (i-obj.sphere_approximation_pixelated_radius):(i+obj.sphere_approximation_pixelated_radius)
                    for j0 = (j-obj.sphere_approximation_pixelated_radius):(j+obj.sphere_approximation_pixelated_radius)
                        for k0 = (k-obj.sphere_approximation_pixelated_radius):(k+obj.sphere_approximation_pixelated_radius)
                            if (((i-i0)^2 + (j-j0)^2 + (k-k0)^2 <= obj.sphere_approximation_pixelated_radius^2)&&~all(pixel_emissions(cavitationSimulation.getInd(i0, j0, k0, obj.kgrid.Nx, obj.kgrid.Ny), :) == 0))
                                pressures = pressures + pixel_emissions(cavitationSimulation.getInd(i0, j0, k0, obj.kgrid.Nx, obj.kgrid.Ny), :)*scaling*obj.artificial_scale_up;
                            end
                        end
                    end
               end
               sourc.p(i, :) = pressures;
           end
           obj.current_sim_source = sourc;
           source = sourc;
       end
       function data = iteration(obj, input_data)
           %input data is cubical emissions, 
           % output is cubical emissions after smoothing, 
           % calculating analytical (from microbubble to microbubble (not cube)), 
           % smoothing for the cubes, 
           % SCALING UP PIXEL emissions by the number of microbubbles in each PIXEL, 
           % calculating source with spheres, calculating k-wave, 
           % subtracting without boundary, and recalculating cubical_emissions from received
           data = 0;
       end
       function data = last_iteration(obj, input_data) %only difference is that the last one sets the emissions, calculated the radii, and sets the radii
           data = 0;
       end
       function data = run_final_sim(obj, bubble_emissions) 
           data = 0;
       end
       function data = run(obj) %check for convergence on iterations
            %ALSO OUTPUT RADII IN ADDITION TO PRESSURE AT SENSORS
            data = 0;
       end
       function sidelength = determine_convergent_cube_pixelated_sidelength(obj) %keep re-running run with different sidelengths until convergence is reached
           sidelength = 0;
       end
   end 
   methods(Static)
       function n = get_ind(x, y, z, Nx, Ny)
           n = x + (y-1)*Nx + (z-1)*Nx*Ny;
       end
       function loc = get_coords(n, Nx, Ny)
           x = mod(n-1, Nx)+1;
           y = mod((n-1-X)/Nx, Ny) + 1;
           z = ((((n-1-X)/Nx) - Y)/Ny) + 1;
           loc = [x; y; z];
       end
       function returnval = rescaleCubicalReceived(cubicalRecieved)
           %Make Non-Static OR DELETE
       end
   end
end