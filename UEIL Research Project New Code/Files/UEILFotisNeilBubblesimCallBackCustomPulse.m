%For now just try custom amplitude but don't try the file stuff... try and
%overrule

function [particle,pulse,linear,simulation,graph]=UEILFotisNeilBubblesimCallBackCustomPulse(gas_model, radius, thickness, shear, viscocity, liquid_name, amplitude, cycles, frequency, sampling_rate, t, p)
action = 'calculate';
% function [particle,pulse,linear,simulation,graph]=BS_BubblesimCallback(action)
%
% Bubblesim: Simulate bubble response
% Callback functions
% Display and/or calculate bubble response
%
% action: 'display':   Display pulse
%         'calculate': Simulate result

% Lars Hoff, NTNU, Dept. of Telecommunications
% Trondheim, Norway

BS_WriteFunctionname
BS_WriteMessage('');
BS_WriteMessage('No warnings', 'warning');

switch(action)
 case 'display',
  [particle, pulse, linear, simulation, graph] = DisplayCustomPulse(gas_model, radius, thickness, shear, viscocity, liquid_name, amplitude, cycles, frequency, sampling_rate, t, p);
 case 'calculate',
  [particle, pulse, linear, simulation, graph] = DisplayCustomPulse(gas_model, radius, thickness, shear, viscocity, liquid_name, amplitude, cycles, frequency, sampling_rate, t , p);
  [particle, pulse, linear, simulation, graph] = ...
               BS_SimulateResult(particle, pulse, linear, simulation, graph);
 otherwise,
  error('CalculateBubblesim: Unknown calculation option')
end
return


%=== DISPLAY PULSES ========================================
function [particle,pulse,linear,simulation,graph] = DisplayCustomPulse(gas_model, radius, thickness, shear, viscocity, liquid_name, amplitude, cycles, frequency, sampling_rate, t, p)
BS_WriteFunctionname
BS_WriteMessage('Plotting input pulses')

particle  = [];  % Particle parameters
pulse     = [];  % Incoming pulse parameters
simulation= [];  % Simulation results
graph     = [];  % Plotting parameters

BS_WriteMessage('Reading input parameters')
[particle,pulse,simulation,graph]= ...
          UEILFotisNeilBubblesimFromBatch(gas_model, radius, thickness, shear, viscocity, liquid_name, amplitude, cycles, frequency, sampling_rate, particle,pulse,simulation,graph,'read');

BS_WriteMessage('Defining input pulse')
[particle, pulse, graph]= BS_DefinePulse(particle,pulse,graph);

pulse(1).t = t;
pulse(1).p = p;
pulse(2).t = t;
pulse(2).p = -1*p;

graph.tmax = max(pulse(1).t)/2;   % Display settings. Can be zoomed in or out
graph.fmax = 4*pulse(1).f0;

graph = BS_DefineGraphs(graph);

%[particle,pulse,simulation,graph]= ...
          %UEILFotisNeilBubblesimFromBatch(gas_model, radius, thickness, shear, viscocity, liquid_name, amplitude, cycles, frequency, sampling_rate, particle,pulse,simulation,graph,'write');

BS_WriteMessage('Calculating linear results ')
[linear]= BS_LinearOscillation (particle, pulse(1).t, pulse(1).p );

BS_WriteMessage('Plotting input pulses')
BS_PlotInitialPulses ( particle, pulse, linear, graph ); %Double Check, Not Essential For Now

set (graph.figure, 'Name', sprintf('Results  %s' , graph.title) );

BS_WriteMessage('Input pulses ready')
drawnow

return


%=== SIMULATE BUBBLE OSCILLATION ===========================
function [particle, pulse, linear, simulation, graph] = ...
   BS_SimulateResult(particle,pulse,linear,simulation,graph);

MessageWindow= gcf;
set (MessageWindow, 'Pointer', 'Watch'); % Display waiting indicator
drawnow;

%=== Simulate response =====================================
if (pulse(1).invert)
  N=2;
  simulation(2)=simulation(1);
else
  N=1;
end
 
for k=1:N
  tic
  BS_WriteMessage(sprintf('Simulating ODE numerically. Pulse %d', k) );
  [t, a, ps]= BS_SimulateOscillation (particle, pulse(k), simulation(k) ); %Check this, Important
  [tr,pr,fs]= BS_ConstantSampleRate( t, ps, pulse(k).fs ); %Check this, Important

  simulation(k).t = t;    % [s]    Time vector from ODE solver,  uneven sampling
  simulation(k).a = a;    % [m]    Radius and velocity, uneven sampling
  simulation(k).p = ps;   % [Pa]   Scattered sound pressure, uneven sampling
  simulation(k).fs= fs;   % [1/s]  Sample rate, after resampling to constant rate
  simulation(k).tr= tr;   % [s]    Time vector, resampled to constant rate
  simulation(k).pr= pr;   % [Pa]   Scattered pressure, resampled to constant rate

  simulation(k).etime = toc;
end

set (MessageWindow, 'Pointer', 'Arrow'); % Remove waiting indicator

%=== Plot results ==========================================
BS_WriteMessage('Plotting results')
BS_PlotSimulation ( particle, pulse, linear, simulation, graph ); %Double check, not essential right now

save ( graph.resultfile , 'particle', 'pulse', 'simulation' );
BS_WriteMessage(sprintf('Finished. Results saved to file %s', graph.resultfile) );

return

