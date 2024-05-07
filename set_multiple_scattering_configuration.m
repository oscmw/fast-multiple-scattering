% Set scatterer properties for examples.
%
%   [N,R,kwave,geom] = set_multiple_scattering_configuration() sets the
%   wavenumber kwave and geometry geom for N scatterers equally spaced on a
%   circle of radius R centred at the origin.
%
%   This function is used inside the fast-multiple-scattering examples to
%   set the geometry. Changing this changes the geometry for all examples.
%
% Example: elliptical scatterers with width 2 and height 1
%
%   geom = obstacleEllipse(2,1)
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


function [N,R,kwave,geom] = set_multiple_scattering_configuration()

% number of scatterers
N = 100;

% radius of circle on which scatterers sit
R = 2048;

% shape of the scatterer - circle with radius 1
geom = obstacleCircle(1);

% wavenumber
kwave = 0.091993928362805;
