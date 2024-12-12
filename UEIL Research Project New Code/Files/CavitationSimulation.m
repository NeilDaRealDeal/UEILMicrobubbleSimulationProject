% To run a simulation, first set up the regular K-wave simulation,
% then plug the parameters into this method along with the microbubble
% location and parameters.

% Assumtions for all models: Same type of microbubbles, translational 
% force of wave is negligible, motion of microbubbles resulting from 
% blood flow is negligible, RP model holds perfectly, etc.
% wave scattering can be considered solely additive

% Note for future: code a "with transducer as a sensor" case

classdef CavitationSimulation
   properties
      Value {mustBeNumeric}
   end
   methods
      function r = roundOff(obj)
         r = round([obj.Value],2);
      end
      function r = multiplyBy(obj,n)
         r = [obj.Value]*n;
      end
   end
end