% Fast multiple scattering code to accompany the paper [1].
%
% Note: This code requires TMATROM, which is available at 
% http://www.romapp.org. 
%
% Validation for scattering by circles (using plane wave scattering):
%
%   check_sound_soft.m
%   check_transmission.m
%
% Validation for circular and noncircular scatterers (using point source
% problem (39) in Reference [1]).
%
%   example_sound_soft.m
%   example_transmission.m
%
% Functions to solve the multiple scattering problem:
%
%   fast_multiple_scattering_soft.m
%   fast_multiple_scattering_transmission.m
%
% Functions to set the geometry and material properties:
%
%   set_multiple_scattering_configuration.m
%   set_multiple_scattering_material_constants.m
%
% Function for solving single-scatterer transmission problems:
%
%   solverNystromPenetrable.m
%
% Function for computing radiating wavefunction expansions from far field
% data (see Reference [2]):
%
%   farfield2radiatingwavefunctionexpansion.m
% 
%
% References:
%
% [1] A fast algorithm for the two-dimensional Helmholtz transmission
% problem with large multiple scattering configurations, S. C. Hawkins and
% M. Ganesh, J. Acoust. Soc. Am 2024.
%
% [2] Approximation of radiating waves in the near-field: error estimates
% and application to a class of inverse problems, J. Barkhan, M. Ganesh
% and S. C. Hawkins, J. Comput. Appl. Math. vol 401 pp 113769.

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


%
% Stuart C. Hawkins - 7 May 2024

