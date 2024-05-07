% Test script for fast acoustic sound-soft multiple-scattering simulation
%
%   check_transmission computes the far field for multiple scattering of the
%   plane wave plane_wave(0,kwave) where kwave is the wavenumber. The far
%   field is computed using two algorithms: 
%
%      1) fast_multiple_scattering_soft.m which implements the 
%         algorithm in [1];
%
%      2) MieSolver, which implements the Mie series based algorithm in [2].
%
%   When the scatterers are unit circles the far fields should agree.
%
%   The geometry of the multiple scatterers is set in
%   set_multiple_scattering_configuration.m.
%
% Algorithm: this code uses the fast Stage 3 algorithm in [1], the
% TMATROM package from http://www.romapp.org, and the MieSolver package
% from http://www.miesolver.org.
%
% See also: fast_multiple_scattering_soft, plane_wave.
%
% References:
%
% [1] A fast algorithm for the two-dimensional Helmholtz transmission
% problem with large multiple scattering configurations, S. C. Hawkins and
% M. Ganesh, J. Acoust. Soc. Am 2024.
%
% [2] Algorithm 1009: MieSolver - an object-oriented Mie series software
% for wave scattering by cylinders, S. C. Hawkins, ACM Trans. Math. Softw. 
% vol. 46, 19:1--19:28, 2020.
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

%--------------------------------------------------
% fast multiple scattering version
%--------------------------------------------------

% create points at which to compute the far field
tp = linspace(0,2*pi,1000);

% call the function to compute the far field
f = fast_multiple_scattering_soft(tp,'PLANEWAVE');

%--------------------------------------------------
% MieSolver [2] version
%--------------------------------------------------

% N, R, geom, kwave are set by an external function
[N,R,kwave,geom] = set_multiple_scattering_configuration();

% plane wave with direction (1,0)
inc = plane_wave(0,kwave);

% set scatterer positions
pos = R * exp(1i*2*pi*(0:N-1)/N);

% initialise MieSolver object
p = MieSolver(inc);

% add the scatterers
for j=1:N
    s{j} = scatterer(pos(j),1,'SOFT');
    p.addScatterer(s{j});
end

% solve the scattering problem
p.solve

% get the far field
g = p.getFarfield(exp(1i*tp));

%--------------------------------------------------
% finish off
%--------------------------------------------------

% print the error
fprintf('Error vs MieSolver: %0.2e\n',max(abs(f(:)-g(:))));

% plot the far fields
figure
plot(tp,10*log10(abs(f).^2),'b-',tp,10*log10(abs(g).^2),'r-')
xlabel('observation angle')
ylabel('10*log10(|u|^2)')
legend('fast-multiple-scattering','MieSolver')
