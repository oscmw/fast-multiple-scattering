% Example script for fast acoustic transmission multiple-scattering simulation
%
%   example_transmission computes the far field ff for the point-source
%   problem (39) in [1]. The exact solution is known for this problem, and
%   so the error is computed and printed in the command window. The far
%   field is computed using fast_multiple_scattering_transmission.m, which
%   implements the algorithm in [1].
%
%   The geometry of the multiple scatterers is set in
%   set_multiple_scattering_configuration.m. The acoustic transmission
%   problem parameters are set in
%   set_multiple_scattering_material_constants.m.
%
%   Changing 'POINTSOURCE' to 'PLANEWAVE' in Line 37 computes the far field
%   induced by a plane wave.
%   
% Algorithm: this code uses the fast Stage 3 algorithm in [1] and the
% TMATROM package from http://www.romapp.org.
%
% See also: check_transmission, example_sound_soft, plane_wave.
%
% References:
%
% [1] A fast algorithm for the two-dimensional Helmholtz transmission
% problem with large multiple scattering configurations, S. C. Hawkins and
% M. Ganesh, J. Acoust. Soc. Am 2024.
%
% Stuart C. Hawkins - 7 May 2024

% Copyright 2024 Stuart C. Hawkins and M. Ganesh.
% 	
% This file is part of fast-multiple-scattering.
% 
% fast-multiple-scattering is free software: you can redistribute it
% and/or modify	it under the terms of the GNU General Public License as
% published by the Free Software Foundation, either version 3 of the License,
% or (at your option) any later version.
% 
% fast-multiple-scattering is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with fast-multiple-scattering. If not, see <http://www.gnu.org/licenses/>.


clear all
close all

% create points at which to compute the far field
tp = linspace(0,2*pi,1000);

% call the function to compute the far field
ff = fast_multiple_scattering_transmission(tp,'POINTSOURCE');

% plot the far field
figure
plot(tp,10*log10(abs(ff).^2),'b-')
xlabel('observation angle')
ylabel('10*log10(|u|^2)')

