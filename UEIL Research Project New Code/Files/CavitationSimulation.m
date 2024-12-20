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


classdef CavitationSimulation
   properties
      %Microbubble positional parameters
      microbubble_count
      microbubble_grid_locations_x %Generated
      microbubble_grid_locations_y %Generated
      microbubble_grid_locations_z %Generated
      microbubble_positions_x %Generated
      microubbble_positions_y %Generated
      microbubble_positions_z %Generated
      measurement_cube_pixelated_sidelength %Odd
      measurement_cube_grid_locations_x %Generated
      measurement_cube_grid_locations_y %Generated
      measurement_cube_grid_locations_z %Generated
      mb_fluid_region_matrix
      gaussian_positional_noise
      sphere_approximation_pixelated_radius
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
      %Results
      microbubble_radii_by_cube %Generated
      microbubble_scattered_pulses_by_cube %Generated
      microbubble_radii %Generated
      microbubble_scattered_pulses %Generated
      sensor_data %Generated
   end
   methods
       function generate_microbubble_possandlocs(obj) %grid with small noise, not entirely random
           
       end
       function generate_grid_cube_positions(obj) 
           
       end
       function data = run_sensing_sim(obj) %return data in terms of a matrix with a time and grid cube index axes
           obj.generate_grid_cube_positions();
           data = 0;
       end
       function data = calculate_cubical_emissions_from_cubical_recieved(obj, cubical_recieved)
           data = 0;
       end
       function data = calculate_cubical_radii_from_cubical_recieved(obj, cubical_recieved)
           data = 0;
       end
       function data = calculate_bubble_emissions_from_cubical_emissions(obj, cubical_emissions) %smoothing
           data = 0;
       end
       function source = calculate_source_from_microbubble_emissions(obj, microbubble_emissions)
           source = 0;
       end
       function data = iteration(obj, input_data) %input data is cubical emissions, output is cubical emissions after smoothing, calculating analytical (from microbubble to cube), calculating source with spheres, calculating k-wave, subtracting without boundary, and recalculating cubical_emissions from recieved
           data = 0;
       end
       function data = last_iteration(obj, input_data) %only difference is that the last one sets the emissions, calculated the radii, and sets the radii
           data = 0;
       end
       function data = run_final_sim(obj, microbubble_emissions)
           data = 0;
       end
       function data = run(obj) %check for convergence on iterations
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
   end
end