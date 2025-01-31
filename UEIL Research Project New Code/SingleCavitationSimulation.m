% To run a simulation, first set up the regular K-wave simulation,
% then plug the parameters into this method along with the microbubble
% location and parameters.

% Assumtions for all models: Same type of microbubbles, translational 
% force of wave is negligible, motion of microbubbles resulting from 
% blood flow is negligible, RP model holds perfectly, etc.

% IMPORTANT: We will run 3d simulations, and although any sensor can
% be used, a kWaveArray transducer source must be utilized.

%Note: We will model the microbubble as a very small octohedron as a sphere
%approximation: potentially improve upon this later
%Note: If Bubblesim Errors, Try Running Startup First
%And Also Try addpaths

% Notes for Next Steps:
% Multiple Bubble Simulation
% Nonspherical Oscillations with Proteus
% Check Validity of Convolution Approach and Sampling rate of 10, and If Wrong, Can Increase kwave.dt and Bubblesim dt 
% Increase simulation length to ensure tapering off

classdef SingleCavitationSimulation
   properties
       microbubble_location_x %x_position of microbubble in grid_units
       microbubble_location_y %y_position of microbubble in grid_units
       microbubble_location_z %z_position of microbubble in grid_units
       microbubble_true_location_x %x_position of microbubble in meters
       microbubble_true_location_y %y_position of microbubble in meters
       microbubble_true_location_z %z_position of microbubble in meters
       gas_model %String, same as UEILFotisNeilBubblesimCallBack
       radius %String, same as UEILFotisNeilBubblesimCallBack
       thickness %String, same as UEILFotisNeilBubblesimCallBack
       shear %String, same as UEILFotisNeilBubblesimCallBack
       viscocity %String, same as UEILFotisNeilBubblesimCallBack
       liquid_name %String, same as UEILFotisNeilBubblesimCallBack
       kgrid %Same as array transducer example
       medium %Same as array transducer example
       true_kwavearray %Same as array transducer example
       true_source_signal %Same as array transducer example
       true_source %Same as array transducer example
       true_sensor %Same as array transducer example
       source_amp %Same as array transducer example

       input_args %Same as array transducer example
       deltax %double meters, used to determine the size of the octahedral bubble (make as small as possible)
       n %int number of points in bubble source, equals 6 for octahedron
       artificial_scale_up;
       ts; %true source, set in code

       outputsave1; %Driven
       outputsave2; %Time
       outputsave3; %Scattered (with new time)
       outputsave4; %Emitted
       outputsave5; %Sim_source_signal
       outputsave6; %sourc.p
       outputsave7; %sim_kwavearray
       outputsave8; %sourc.p_mask
       outputsave9; %sourc
       outputsave10;
   end
   methods(Static)
      function sourc = set_up_true_source(obj)
          sourc.p_mask = obj.true_kwavearray.getArrayBinaryMask(obj.kgrid);
          sourc.p = obj.true_kwavearray.getDistributedSourceSignal(obj.kgrid, obj.true_source_signal);
      end
      function sour = set_up_sim_source(obj, emitted)
          sim_kwavearray = obj.true_kwavearray;
          sim_source_signal = obj.true_source_signal;
          sim_kwavearray.addCustomElement([obj.microbubble_true_location_x + obj.deltax, obj.microbubble_true_location_x - obj.deltax, obj.microbubble_true_location_x, obj.microbubble_true_location_x, obj.microbubble_true_location_x, obj.microbubble_true_location_x; obj.microbubble_true_location_y, obj.microbubble_true_location_y, obj.microbubble_true_location_y + obj.deltax, obj.microbubble_true_location_y - obj.deltax, obj.microbubble_true_location_y, obj.microbubble_true_location_y; obj.microbubble_true_location_z, obj.microbubble_true_location_z, obj.microbubble_true_location_z, obj.microbubble_true_location_z, obj.microbubble_true_location_z + obj.deltax, obj.microbubble_true_location_z - obj.deltax], (4*obj.deltax^3)/3, 3, 'bubble');
          %Reflip later
          %sim_source_signal(size(sim_source_signal, 1) + 1, 1 : length(emitted)) = emitted;
          %sim_source_signal = flip(sim_source_signal, 1);
          %ave('test_sim_source_signal.mat', 'sim_source_signal')
          sourc.p_mask = sim_kwavearray.getArrayBinaryMask(obj.kgrid);
          sourc.p = sim_kwavearray.getDistributedSourceSignal(obj.kgrid, sim_source_signal);
          %save('test_sim_kwavearray.mat', 'sim_kwavearray')
          sour.p_mask = sourc.p_mask;
          sour.p = sourc.p;
          %save('test_sourc.mat', 'sourc')
      end
      function n = get_ind(x, y, z, Nx, Ny)
          n = x + (y-1)*Nx + (z-1)*Nx*Ny;
      end    
      function sourc = set_up_sim_source_smart(obj, emitted)
          sourc.p_mask = obj.true_source.p_mask;
          sourc.p_mask(obj.microbubble_location_x, obj.microbubble_location_y, obj.microbubble_location_z) = 1;
          %base_matrix = zeros(obj.kgrid.Nx, obj.kgrid.Ny, obj.kgrid.Nz);
          %base_matrix(obj.microbubble_location_x, obj.microbubble_location_y, obj.microbubble_location_z) = 1;
          mask_ind = find(sourc.p_mask);
          microbubble_ind = SingleCavitationSimulation.get_ind(obj.microbubble_location_x, obj.microbubble_location_y, obj.microbubble_location_z, obj.kgrid.Nx, obj.kgrid.Ny);
          index = find(mask_ind == microbubble_ind);
          sourc.p(1:(index - 1), :) = obj.true_source.p(1:(index-1), :);
          sourc.p(index, :) = emitted;
          sourc.p((index+1):length(mask_ind), :) = obj.true_source.p(index:(length(mask_ind)-1), :);
          % sourc.p = obj.true_source.p; 
          % sourc.p(size(obj.true_source.p, 1) + 1, :) = sourc.p(550, :);
          % sourc.p(550, :) = emitted;
          % for t = 1:length(emitted)
          %     sourc.p(:, :, :, t) = obj.true_source.p(:, :, :, t) + emitted(t)*base_matrix;
          % end
      end
      function sourc = set_up_sim_source_smart2(obj, emitted)
          sim_kwavearray = obj.true_kwavearray;
          sim_source_signal = obj.true_source_signal;
      end
      function driven = convert_sensed_to_driven(obj, sensed) %output driven (just processing, especially conversion to intervals of 10^-8 from 10^-7)
        %Note: For now, we assume that sensed time resolution is less than
        %the driven.
        kwaverecipdtscale = 7; %Double Check This
        bubblesimrecipdtscale = 8;
        drive = zeros(1, 10^(bubblesimrecipdtscale - kwaverecipdtscale) * length(sensed));
        scale_up = 10^(bubblesimrecipdtscale - kwaverecipdtscale);
        scale_down = 10^(kwaverecipdtscale - bubblesimrecipdtscale);
        for r = 1:length(drive)
            if ceil(scale_down*r) == floor(scale_down*r)
                drive(r) = sensed(floor(scale_down*r));
            elseif floor(scale_down*r) == 0
                drive(r) = (mod(r, scale_up)*sensed(ceil(scale_down*r))/scale_up);
            else
                drive(r) = ((mod(r, scale_up)*sensed(ceil(scale_down*r)) + (scale_up - mod(r, scale_up))*sensed(floor(scale_down*r)))/scale_up);
            end
        end
        driven = drive;
      end
      function scattered = convert_driven_to_scattered(obj, driven) %output scattered - note: driven is in intervals of 10^-8
          time = zeros(1, length(driven));
          for r = 1:length(time)
              time(r) = (r-1)*(10^(-8));
          end
          %If below still doesn't work, shift to a sampling rate of 100!!!
          disp(size(driven));
          %disp("Driven:");
          %disp(driven);
          %save('test_driven.mat', 'driven')
          %disp(size(time));
          %disp("Time");
          %disp(time);
          %save('test_time.mat', 'time')
          [particle2,pulse2,linear2,simulation2,graph2] = UEILFotisNeilBubblesimCallBackCustomPulse(obj.gas_model, obj.radius, obj.thickness, obj.shear, obj.viscocity, obj.liquid_name, '30', '2.5', '10', '10', time.', driven.');
          scattered = simulation2.pr * obj.artificial_scale_up;
          %disp("scattered:")
          %disp(scattered);
          %save('test_scattered.mat', 'scattered')
      end
      function emitted = convert_scattered_to_emitted(obj, scattered) %output emitted (just processing)
          emitted_from_microbubble = scattered.'/obj.n;
          emitted = emitted_from_microbubble*str2num(obj.radius)*(1e-6)*2/obj.kgrid.dx;
      end

      function sensed = run_sensing_sim(obj) %Output bubble sensed pulse (Note: just sense bubble)
         true_source = SingleCavitationSimulation.set_up_true_source(obj);
         %obj.ts = true_source;
         now_sensor.mask = zeros(obj.kgrid.Nx, obj.kgrid.Ny, obj.kgrid.Nz);
         now_sensor.mask(obj.microbubble_location_x, obj.microbubble_location_y, obj.microbubble_location_z) = 1; 
         [data] = kspaceFirstOrder3D(obj.kgrid, obj.medium, obj.ts, now_sensor, obj.input_args{:}, ...
            'DataCast', 'gpuArray-single', ...
            'PlotScale', [-1, 1] * obj.source_amp); %Note: After input_args, inputs are temporary
         sensed = data;
      end

      function sensor_data = run_true_sim(obj, sensed) %Output true sensor data
          bubble_emission = SingleCavitationSimulation.convert_scattered_to_emitted(obj, SingleCavitationSimulation.convert_driven_to_scattered(obj, SingleCavitationSimulation.convert_sensed_to_driven(obj, sensed)));
          %save('test_emission.mat', 'bubble_emission')
          sim_source = SingleCavitationSimulation.set_up_sim_source_smart(obj, bubble_emission);
          disp("sim_source")
          [data] = kspaceFirstOrder3D(obj.kgrid, obj.medium, sim_source, obj.true_sensor, obj.input_args{:}, ...
            'DataCast', 'gpuArray-single', ...
            'PlotScale', [-1, 1] * obj.source_amp); %Note: After input_args, inputs are temporary
          sensor_data = data;
      end
      
      function sensor_data = run(obj) %Output true sensor data
         sensed = SingleCavitationSimulation.run_sensing_sim(obj);
         sensor_data = SingleCavitationSimulation.run_true_sim(obj, sensed);
      end
      % 
      % function output_data = run_testing(obj)
      %     sensed = SingleCavitationSimulation.run_sensing_sim(obj);
      %     sensor_data = SingleCavitationSimulation.run_true_sim(obj, sensed);
      %     output_data = [sensed, sensor_data];
      % end
   end
end