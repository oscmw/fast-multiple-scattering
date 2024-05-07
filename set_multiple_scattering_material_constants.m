% Set material properties for acoustic transmission examples.
%
%   [refin,rho0,rho1] = set_multiple_scattering_material_constants() sets
%   the refractive index refin and interior density rho1 of the scatterers,
%   and the density rho0 of the exterior medium.
%
%   This function is used inside the fast-multiple-scattering examples to
%   set the geometry.
%
%  Example: water drops in air
%
%   refin = 0.226666666666667;
%   rho0 = 1.25;
%   rho1 = 1000;
%
% See also: check_sound_soft, check_transmission.
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


function [refin,rho0,rho1] = set_multiple_scattering_material_constants()

% refractive index
refin = 0.226666666666667;

% exterior density
rho0 = 1.25;

% density inside scatterer
rho1 = 1000;