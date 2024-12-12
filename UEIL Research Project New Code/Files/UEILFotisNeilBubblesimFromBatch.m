function [particle,pulse,simulation,graph] = UEILFotisNeilBubblesimFromBatch(gas_model, radius, thickness, shear, viscocity, liquid_name, amplitude, cycles, frequency, sampling_rate, particle,pulse,simulation,graph,action) 
% function [particle,pulse,simulation,graph] =BS_GUIparameters(particle,pulse,simulation,graph,action) 
%
% Bubblesim: Simulate bubble response
% Read or write particle and pulse parameters from/to GUI
%

% Lars Hoff, NTNU, Dept. of Telecommunications
% Trondheim, Norway

global nano micro milli centi kilo Mega

BS_WriteFunctionname

%=== IDENTIFY HANDLES =============================================

%--- Particle
% input.odesolver     = 'Stiff. Variable order';   
% input.ode           = 'R-P with Radiation Damping';   
% input.radius        = '1.5'; 
% input.shellthickness= '0.3';
% input.shellG        = '20'; 
% input.shelleta      = '0.6'; 
% input.liquid        = 'Blood'; 
% 
% %--- Pulse
% input.pulsetype= 'Hanning';
% input.A        = '0.5';
% input.Nc       = '4';
% input.f0       = '2.25';
% input.fs       = '100';
% 
% %--- Simulation
% input.invert         = false;
% input.thermaldamping = 'adiabatic';
% 
% %--- Plotting
% input.plotlinear     = true;
% input.displayprogress= true;             
% input.plot(1)        = true;
% input.plot(2)        = true;
% input.plot(3)        = true;
% input.plot(4)        = false;
% input.plot(5)        = true;
% input.plot(6)        = false;

switch action
 %======== READ PARAMETERS FROM GUI ========================
 case 'read',  
  simulation.solver.name        = 'Stiff. Variable order';
  simulation.solver.command     = 'ode15s';
  simulation.model.name         = 'R-P with Radiation Damping';
  simulation.model.ode          = 'ModifiedRayleigh';
  simulation.thermaldamping.name = gas_model; 
  simulation.thermaldamping.command = gas_model;

  particle.a0= micro*str2num( radius ); 
  particle.ds= nano *str2num( thickness );
  particle.Gs= Mega *str2num( shear ); 
  particle.es=       str2num( viscocity ); 
  particle.liquid.name =  liquid_name;
  particle =         BS_PhysicalConstants(particle);

  pulse(1).envelope.name = 'Hanning';
  pulse(1).envelope.command = {'hanning', 0};
  pulse(1).A       = Mega*str2num( amplitude );
  pulse(1).Nc      =      str2num( cycles );
  pulse(1).f0      = Mega*str2num( frequency );
  pulse(1).fs      = Mega*str2num( sampling_rate );
  pulse(1).invert  = false;

  simulation.displayprogress= true;  
  graph.plotlinear          = true; 
  graph.include(1) = true;
  graph.include(2) = true;
  graph.include(3) = true;
  graph.include(4) = false;
  graph.include(5) = true;
  graph.include(6) = false;
  
  % for k=1:length(input.plot)
  %   graph.include(k) = get( input.plot(k), 'value');
  % end
   
%=== WRITE PARAMETERS TO GUI ===============================
 case 'write'
  % set (input.radius,         'string', particle.a0 *1e6  );
  % set (input.shellthickness, 'string', particle.ds *1e9  );
  % set (input.shellG,         'string', particle.Gs *1e-6 );
  % set (input.shelleta,       'string', particle.es       );  
  % set (input.A,              'string', pulse(1).A/Mega   );
  % set (input.Nc,             'string', pulse(1).Nc       );
  % set (input.f0,             'string', pulse(1).f0/Mega  );
  % set (input.fs,             'string', pulse(1).fs/Mega  );
end

%== DEBUGGING: Display values ========================================
% $$$ input
% $$$ graph
% $$$ particle
% $$$ pulse
% $$$ simulation

return










% particle  = [];  % Particle parameters
% pulse     = [];  % Incoming pulse parameters
% simulation= [];  % Simulation results
% graph     = [];  % Plotting parameters
% 
% BS_WriteMessage('Reading input parameters')
% [particle,pulse,simulation,graph]= ...
%           BS_GUIparameters(particle,pulse,simulation,graph,'read');
% 
% BS_WriteMessage('Defining input pulse')
% [particle, pulse, graph]= BS_DefinePulse(particle,pulse,graph);
% 
% graph.tmax = max(pulse(1).t)/2;   % Display settings. Can be zoomed in or out
% graph.fmax = 4*pulse(1).f0;
% 
% graph = BS_DefineGraphs(graph);
% 
% [particle,pulse,simulation,graph]= ...
%           BS_GUIparameters(particle,pulse,simulation,graph,'write');
% 
% BS_WriteMessage('Calculating linear results ')
% [linear]= BS_LinearOscillation (particle, pulse(1).t, pulse(1).p );
% 
% BS_WriteMessage('Plotting input pulses')
% BS_PlotInitialPulses ( particle, pulse, linear, graph );
% 
% set (graph.figure, 'Name', sprintf('Results  %s' , graph.title) );
% 
% BS_WriteMessage('Input pulses ready')
% drawnow
% 
% return
