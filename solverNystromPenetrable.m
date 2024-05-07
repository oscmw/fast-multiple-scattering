% Nystrom solver for acoustic transmission problem
%
%  This code requires TMATROM, which is available at http://www.romapp.org. 
%
%  s = solverNystromPenetrable(k,inc,b,refin,rho0,rh1) creates a solver
%  object with wavenumber k and incident fields specified in the incident
%  field inc. The scatterer boundary is specified by obstacle object b. The
%  refractive index of the scatterer is refin, the exterior density is
%  rho0, and the density in the scatterer is rho1.
%
%  s = solverNystrom(k,[],...) creates a solver object with wavenumber k
%  and defers setting the incident fields.
%
%  s.setIncidentField(inc) sets the incident field as specified in the cell
%  array inc.
%
%  s.setup(n) sets up the Nystrom discretisation of an integral equation
%  for the sound soft scattering problem using 2n points on the boundary.
%
%  s.solve() solves the Nystrom system for the incident field inc.
%
%  val = s.getFarField(theta) computes the far field of inc{1} at angles
%  specified by theta.
%
%  val = s.getFarField(theta,index) computes the far field of inc{index} at 
%  angles specified by theta. The far field of inc{index(k)} is val(:,k).
%
%  val = s.getCrossSection(self) computes the cross section for the
%  incident field inc{1}.
%
%  val = s.getCrossSection(self,index) computes the cross section for the
%  incident fields inc{index}.
%
%  s.plotFarField(theta) plots the far field for the incident field inc{1}
%  at points with polar angle theta.
%
%  s.plotFarField(theta,index) plots the far field for the incident field 
%  inc{index} at points with polar angle theta.
%
%  s.plotFarField(theta,index,opts) plots the far field for the incident field 
%  inc{index} at points with polar angle theta, with plot linetype
%  specified by opts.
%
%  obj = s.plotFarField(...) returns the handle of the plot.
%
%  s.plotCrossSection(...) plots the cross section in decibels. See
%  plotFarField for the syntax.
%
%  r = s.getRadius() returns the approximate radius r of the obstacle b.
%
%  s.visualize() plots the obstacle b.
%
% Algorithm: Uses the transmission formulation in [1] discretised using the
%   Nystrom scheme described in [2]
%
% References:
%
% [1] A fast algorithm for the two-dimensional Helmholtz transmission
% problem with large multiple scattering configurations, S. C. Hawkins and
% M. Ganesh, J. Acoust. Soc. Am 2024.
%
% [2] Boundary integral equations in time-harmonic acoustic scattering, R.
% Kress, Mathl. Comput. Modelling, vol 15, 229--243, 1991.
%
% See also: solver, solverNystrom, solverNystromRobin, mfsSolver.
%
% Stuart C. Hawkins - 5 May 2024

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


classdef solverNystromPenetrable < solver
    
    properties
        coupling_parameter
        nystrom_parameter
        matrix
        cof
        scatterer
        refin
        rho0
        rho1
    end
    
    methods
        
        %-----------------------------------------
        % constructor
        %-----------------------------------------
        
        function self = solverNystromPenetrable(kwave,incidentField,scatterer,refin,rho0,rho1)
            
            %  call parent constructor
            self = self@solver(kwave,incidentField);
            
            % set coupling parameter
            self.coupling_parameter = self.kwave;
            
            % check that the scatterer is of the correct type
            if ~isa(scatterer,'obstacle')
                error('scatterer must be of class obstacle')
            end
            
            % set the scatterer
            self.scatterer = scatterer;
            
            % set some things to empty... this will highlight if
            % setup etc not run
            self.matrix = [];
            self.cof = [];
            
            % store the physical parameters
            self.refin = refin;
            self.rho0 = rho0;
            self.rho1 = rho1;
            
        end
        
        %===============================================================
        % methods defined in the solver class that must be overloaded 
        %===============================================================
        
        %-----------------------------------------
        % setup
        %-----------------------------------------
        
        function setup(self,nystrom_parameter)
            
            % store Nystrom parameter
            self.nystrom_parameter = nystrom_parameter;
            
            % compute the Nystrom matrix
            self.matrix = self.nystrom_matrix(nystrom_parameter);
            
        end
        
        %-----------------------------------------
        % solve
        %-----------------------------------------
        
        function solve(self)
            
            if isempty(self.matrix)
                
                error('Must call setup() first')
                
            end
            
            % setup the right hand side
            b = self.nystrom_rhs(self.nystrom_parameter);
            
            % solve the system
            self.cof = self.matrix \ b;
            
        end
        
        %-----------------------------------------
        % get far field
        %-----------------------------------------
        
        function val = getFarField(self,points,index)
            
            % set default for index
            if nargin<3
                index = 1;
            end
            
            % check that solve has been run
            if isempty(self.cof)
                
                error('Must run solve() first')
                
            end
            
            % compute the far field
            val = self.nystrom_far_field(points,index);
            
            
        end
        
        %===============================================================
        % additional methods that might be useful to the user
        %===============================================================
        
        %-----------------------------------------
        % visualize
        %-----------------------------------------
        
        % Just a wrapper for visualise.

        function varargout = visualize(varargin)

            [varargout{1:nargout}] = visualise(varargin{:});           
            
        end
        
        %-----------------------------------------
        % visualise
        %-----------------------------------------
        
        function varargout = visualise(self,varargin)

            % call the visualize method of the obstacle
            obj = self.scatterer.visualise(varargin{:});
            
            % return handle if required
            if nargout > 0
                varargout{1} = obj;
            end
            
        end
        
        %-----------------------------------------
        % get radius
        %-----------------------------------------
        
        function val = getRadius(self)
            
            % evaluate the geometry parametrisation at lots of points and
            % take maximum
            t = 2*pi*(0:999)/1000;
            [x,y,qx,qy] = self.scatterer.geom(t);
            val = sqrt(max(qx.^2+qy.^2));
            
        end
        
    end % end methods
    
    methods(Access=private)
        
        %===============================================================
        % methods that are not for use by the user
        %===============================================================

        %-----------------------------------------
        % setup Nystrom matrix
        %-----------------------------------------
        
        function A = nystrom_matrix(self,n)
            
            % set constant
            C = 0.57721566490153286060651209008240243104215933593992;
            
            %-----------------------------------------------------------
            % set up quadrature for nonsingular integrals
            %-----------------------------------------------------------
            
            pp=pi*(0:2*n-1)/n;
            pp=pp(:);
            pw=pi*ones(2*n,1)/n;
            
            %-----------------------------------------------------------
            % set up quadrature for singular integrals
            %-----------------------------------------------------------
            
            % initialize quadrature matrix
            R=zeros(2*n,2*n);
            
            % setup the quadrature points... see Colton and Kress,
            % Inverse EM and Acoustic Scattering Theory, 3rd Ed, p78
            R(1,:)=-(pi/n^2)*cos(n*(pp(1)-pp));
            
            for m=1:n-1
                R(1,:)=R(1,:) - (2*pi/(n*m))*cos(m*(pp(1)-pp.'));
            end
            
            for i=2:2*n
                
                R(i,:) = circshift(R(1,:),[0,i-1]);
                
            end
            
            %-----------------------------------------------------------
            % set up matrix
            %-----------------------------------------------------------
            
            % compute the matrix that gives the tangential derivative of
            % the surface potential cf p236 in [2]
            D = self.setupDerivativeMatrix();
            
            % compute average density
            rho = (self.rho1+self.rho0)/2;
            
            % matrix is of the form [AA,BB;EE + FF*D,DD]... initialize submatrices
            AA = zeros(2*n,2*n);
            BB = zeros(2*n,2*n);
            DD = zeros(2*n,2*n);
            EE = zeros(2*n,2*n);
            FF = zeros(2*n,2*n);
            
            % apply mapping to quadrature points
            [x,y,qx,qy,dqx,dqy,ddqx,ddqy,nqx,nqy,jac]=self.scatterer.geom(pp);
            
            % Compute the matrix... see [2] and also Colton and Kress,
            % Inverse EM and Acoustic Scattering Theory, 3rd Ed, p78
            % Note: hypersingular operator appears only in the form 
            % T0-T1 and so the Cauchy singularity is cancelled and we don't
            % need to include it.
            %
            % Note that the hypersingular operator is split into
            % T = T0 + T1*D following (2.7) in [2].
            
            %- - - - - - - - - - - - - - - - - - - - - - - - -
            % compute the off-diagonal part
            %- - - - - - - - - - - - - - - - - - - - - - - - -
            
            for i=1:2*n
                
                % (d,e) is the p-q in G(p,q)
                d = qx(i)-qx;
                e = qy(i)-qy;
                
                % compute |p-q|
                r=sqrt( d.^2 + e.^2 );
                
                % set dummy value for diagonal
                r(i)=1000;
                
                % compute various dot products that appear in the kernels
                dp = d.*nqx + e.*nqy;
                dpt = d*nqx(i) + e*nqy(i);
                dd = nqx*nqx(i) + nqy*nqy(i);
                dtang = (dqx(i)*d+dqy(i)*e)/jac(i);

                % compute exterior Green's function kernels
                single_kernel0 = 0.25*1i*besselh(0,self.kwave*r).*jac;                
                double_kernel0 = 0.25*1i*self.kwave*besselh(1,self.kwave*r).*dp.*jac./r;                
                transpose_kernel0 = -0.25*1i*self.kwave*besselh(1,self.kwave*r).*dpt.*jac./r;                
                hyp1_kernel0 = 0.25i*self.kwave^2*dd.*besselh(0,self.kwave*r).*jac;                
                hyp2_kernel0 = -0.25i*self.kwave*besselh(1,self.kwave*r).*dtang./r;
                
                % compute interior Green's function kernels                
                single_kernel1 = 0.25*1i*besselh(0,self.refin*self.kwave*r).*jac;
                double_kernel1 = 0.25*1i*self.refin*self.kwave*besselh(1,self.refin*self.kwave*r).*dp.*jac./r;                                
                transpose_kernel1 = -0.25*1i*self.refin*self.kwave*besselh(1,self.refin*self.kwave*r).*dpt.*jac./r;                
                hyp1_kernel1 = 0.25i*(self.refin*self.kwave)^2*dd.*besselh(0,self.refin*self.kwave*r).*jac;                
                hyp2_kernel1 = -0.25i*self.refin*self.kwave*besselh(1,self.refin*self.kwave*r).*dtang./r;

                % compute singular part of exterior Green's function kernels                
                single_sing0 = -(1/(4*pi))*besselj(0,self.kwave*r).*jac;
                double_sing0 = -(1/(4*pi))*self.kwave*besselj(1,self.kwave*r).*dp.*jac./r;
                transpose_sing0 = (1/(4*pi))*self.kwave*besselj(1,self.kwave*r).*dpt.*jac./r;
                hyp1_sing0 = -(1/(4*pi))*self.kwave^2*dd.*besselj(0,self.kwave*r).*jac;
                hyp2_sing0 = (0.25/pi)*self.kwave*besselj(1,self.kwave*r).*dtang./r;
                
                % compute singular part of interior Green's function kernels                
                single_sing1 = -(1/(4*pi))*besselj(0,self.refin*self.kwave*r).*jac;
                double_sing1 = -(1/(4*pi))*self.refin*self.kwave*besselj(1,self.refin*self.kwave*r).*dp.*jac./r;
                transpose_sing1 = (1/(4*pi))*self.refin*self.kwave*besselj(1,self.refin*self.kwave*r).*dpt.*jac./r;
                hyp1_sing1 = -(1/(4*pi))*(self.refin*self.kwave)^2*dd.*besselj(0,self.refin*self.kwave*r).*jac;
                hyp2_sing1 = (0.25/pi)*self.refin*self.kwave*besselj(1,self.refin*self.kwave*r).*dtang./r;

                % compute smooth part of exterior Green's function kernels
                single_smooth0 = single_kernel0 - single_sing0.*log(4* (sin(0.5*(pp(i)-pp))).^2 );
                double_smooth0 = double_kernel0 - double_sing0.*log(4* (sin(0.5*(pp(i)-pp))).^2 );                
                transpose_smooth0 = transpose_kernel0 - transpose_sing0.*log(4* (sin(0.5*(pp(i)-pp))).^2 );
                hyp1_smooth0 = hyp1_kernel0 - hyp1_sing0.*log(4* (sin(0.5*(pp(i)-pp))).^2 );
                hyp2_smooth0 = hyp2_kernel0 - hyp2_sing0.*log(4* (sin(0.5*(pp(i)-pp))).^2 );

                % compute smooth part of interior Green's function kernels
                single_smooth1 = single_kernel1 - single_sing1.*log(4* (sin(0.5*(pp(i)-pp))).^2 );
                double_smooth1 = double_kernel1 - double_sing1.*log(4* (sin(0.5*(pp(i)-pp))).^2 );
                transpose_smooth1 = transpose_kernel1 - transpose_sing1.*log(4* (sin(0.5*(pp(i)-pp))).^2 );
                hyp1_smooth1 = hyp1_kernel1 - hyp1_sing1.*log(4* (sin(0.5*(pp(i)-pp))).^2 );
                hyp2_smooth1 = hyp2_kernel1 - hyp2_sing1.*log(4* (sin(0.5*(pp(i)-pp))).^2 );

                % assemble the Nystrom matrices for the individual
                % operators using (3.5) in [2].
                S0 = (single_smooth0.*pw  + single_sing0.*R(:,i));
                S1 = (single_smooth1.*pw  + single_sing1.*R(:,i));
                K0 = (double_smooth0.*pw  + double_sing0.*R(:,i));
                K1 = (double_smooth1.*pw  + double_sing1.*R(:,i));
                D0 = (transpose_smooth0.*pw  + transpose_sing0.*R(:,i));
                D1 = (transpose_smooth1.*pw  + transpose_sing1.*R(:,i));
                TA0 = (hyp1_smooth0.*pw  + hyp1_sing0.*R(:,i));
                TA1 = (hyp1_smooth1.*pw  + hyp1_sing1.*R(:,i));
                TB0 = (hyp2_smooth0.*pw  + hyp2_sing0.*R(:,i));
                TB1 = (hyp2_smooth1.*pw  + hyp2_sing1.*R(:,i));
                
                % assemble the Nystrom matrix for (13) in [1]
                AA(i,:) = self.rho0 * K0 - self.rho1 * K1;
                BB(i,:) = self.rho0^2 * S0 - self.rho1^2 * S1;
                DD(i,:) = self.rho1 * D1 - self.rho0 * D0;
                EE(i,:) = TA1 - TA0; 
                FF(i,:) = TB1 - TB0; 
                                
            end
            
            %- - - - - - - - - - - - - - - - - - - - - - - - -
            % compute the diagonal part
            %- - - - - - - - - - - - - - - - - - - - - - - - -
            
            % compute exterior Green's function kernels
            single_smooth_diagonal0= 0.5*(0.5i - C/pi - (1/(2*pi))*log( 0.25*self.kwave^2*jac.^2 ) ).*jac.*pw;
            double_smooth_diagonal0=-0.5*(dqx.*ddqy-dqy.*ddqx)./(2*n*(dqx.^2+dqy.^2));
            transpose_smooth_diagonal0=-0.5*(dqx.*ddqy-dqy.*ddqx)./(2*n*(dqx.^2+dqy.^2));
            hyp1_smooth_diagonal0 = 0.5*self.kwave^2*(0.5i - C/pi - (0.5/pi)*log( 0.25*self.kwave^2*jac.^2 ) ).*jac.*pw;
            hyp2_smooth_diagonal0 = -0.5*(0.5/pi)*(dqx.*ddqx+dqy.*ddqy)./(dqx.^2+dqy.^2).*pw./jac;

            % compute interior Green's function kernels
            single_smooth_diagonal1= 0.5*(0.5i - C/pi - (1/(2*pi))*log( 0.25*(self.refin*self.kwave)^2*jac.^2 ) ).*jac.*pw;
            double_smooth_diagonal1=-0.5*(dqx.*ddqy-dqy.*ddqx)./(2*n*(dqx.^2+dqy.^2));
            transpose_smooth_diagonal1=-0.5*(dqx.*ddqy-dqy.*ddqx)./(2*n*(dqx.^2+dqy.^2));
            hyp1_smooth_diagonal1 = 0.5*(self.refin*self.kwave)^2*(0.5i - C/pi - (0.5/pi)*log( 0.25*(self.refin*self.kwave)^2*jac.^2 ) ).*jac.*pw;
            hyp2_smooth_diagonal1 = -0.5*(0.5/pi)*(dqx.*ddqx+dqy.*ddqy)./(dqx.^2+dqy.^2).*pw./jac;

            % compute singular part of exterior Green's function kernels
            single_sing_diagonal0 = -(1/(4*pi))*besselj(0,0).*jac;
            double_sing_diagonal0 = zeros(size(double_smooth_diagonal0));
            transpose_sing_diagonal0 = zeros(size(double_smooth_diagonal0));
            hyp1_sing_diagonal0 = -(0.25/pi)*self.kwave^2*ones(size(jac)).*jac;
            hyp2_sing_diagonal0 = zeros(size(hyp2_smooth_diagonal0));

            % compute singular part of interior Green's function kernels
            single_sing_diagonal1 = -(1/(4*pi))*besselj(0,0).*jac;            
            double_sing_diagonal1 = zeros(size(double_smooth_diagonal1));
            transpose_sing_diagonal1 = zeros(size(double_smooth_diagonal1));
            hyp1_sing_diagonal1 = -(0.25/pi)*(self.refin*self.kwave)^2*ones(size(jac)).*jac;
            hyp2_sing_diagonal1 = zeros(size(hyp2_smooth_diagonal1));
                            
            % loop along diagonal
            for i=1:2*n
                % assemble the Nystrom matrice entries for the individual
                % operators using (3.5) in [2].                
                S0 = (single_smooth_diagonal0(i)  + single_sing_diagonal0(i)*R(i,i));
                S1 = (single_smooth_diagonal1(i)  + single_sing_diagonal1(i)*R(i,i));
                K0 = (double_smooth_diagonal0(i)  + double_sing_diagonal0(i)*R(i,i));
                K1 = (double_smooth_diagonal1(i)  + double_sing_diagonal1(i)*R(i,i));
                D0 = (transpose_smooth_diagonal0(i)  + transpose_sing_diagonal0(i)*R(i,i));
                D1 = (transpose_smooth_diagonal1(i)  + transpose_sing_diagonal1(i)*R(i,i));
                TA0 = (hyp1_smooth_diagonal0(i)  + hyp1_sing_diagonal0(i)*R(i,i));
                TA1 = (hyp1_smooth_diagonal1(i)  + hyp1_sing_diagonal1(i)*R(i,i));
                TB0 = (hyp2_smooth_diagonal0(i)  + hyp2_sing_diagonal0(i)*R(i,i));
                TB1 = (hyp2_smooth_diagonal1(i)  + hyp2_sing_diagonal1(i)*R(i,i));

                % assemble the Nystrom matrix entries for (13) in [1]
                AA(i,i) = rho + self.rho0 * K0 - self.rho1 * K1;
                BB(i,i) = self.rho0^2 * S0 - self.rho1^2 * S1;
                DD(i,i) = rho + self.rho1 * D1 - self.rho0 * D0;
                EE(i,i) = TA1 - TA0; 
                FF(i,i) = TB1 - TB0; 
            end
            
            % finally, assemble the matrix for (13) in [1] from the blocks
            % computed above
            A = [AA,BB;EE+FF*D,DD];
            
        end
        
        %-----------------------------------------------------------
        % set up right hand side
        %-----------------------------------------------------------
        
        function rhs = nystrom_rhs(self,n)
            
            % set up quadrature points
            pp=pi*(0:2*n-1)/n;
            pp=pp(:);
            
            % apply mapping to quadrature points
            [x,y,qx,qy,dqx,dqy,ddqx,ddqy,nqx,nqy,jac]=self.scatterer.geom(pp);
            
            % get qmap in complex format for the incident field evaluation
            qz = qx + 1i*qy;
            
            % setup the right hand side vectors
            for k=1:self.numIncidentField
                
                % evaluate the incident field and its gradient
                p = self.incidentField{k}.evaluate(qz);
                [px,py] = self.incidentField{k}.evaluateGradient(qz);
                
                a(:,k) = -self.rho0*p;
                
                % compute Robin trace on boundary
                b(:,k) = (nqx.*px+nqy.*py);
            end
            

            rhs = [a;b];
            
            
        end
        
        %-----------------------------------------------------------
        % compute far field
        %-----------------------------------------------------------
        
        function val = nystrom_far_field(self,points,index)
            
            % make sure points is a column vector
            points = points(:);
            
            % initialize temporary array... do this because it is
            % convenient to compute the far field at each point for all
            % right hand sides in one go. This gives the transpose of val
            tmp = zeros(length(index),length(points));
            
            % get coordinates (p,q) of observation point on the unit circle
            p=cos(points);
            q=sin(points);

            % set up quadrature for nonsingular integrals           
            pp=pi*(0:2*self.nystrom_parameter-1)/self.nystrom_parameter;
            pp=pp(:);
            pw=pi*ones(2*self.nystrom_parameter,1)/self.nystrom_parameter;

            % apply parametrisation of boundary
            [x,y,qx,qy,dqx,dqy,ddqx,ddqy,nqx,nqy,jac] = self.scatterer.geom(pp);

            % loop through point
            for j = 1:length(points)
                
                % commpute kernel using details on p75 of Colton and Kress,
                % Inverse EM and Acoustic Scattering Theory, 3rd Ed
                
                d = p(j)-qx;
                e = q(j)-qy;
                
                r=sqrt( d.^2 + e.^2 );
                
                ndp = p(j).*nqx + q(j).*nqy;
                dp = p(j).*qx + q(j).*qy;
                
                kernelS = ((1+1i)/sqrt(2))/sqrt(8*pi*self.kwave)*self.rho0*jac.*exp(-1i*self.kwave*dp);
                kernelK = ((1+1i)/sqrt(2))/sqrt(8*pi*self.kwave)*jac.*(-1i*self.kwave*ndp).*exp(-1i*self.kwave*dp);
                
                % loop through the right hand sides
                for k = 1:length(index)
                
                    % compute the far field for each right hand side using
                    % (9) in [1]
                    tmp(k,j) = sum(pw.*kernelS.*self.cof(2*self.nystrom_parameter+1:end,index(k))) ...
                        + sum(pw.*kernelK.*self.cof(1:2*self.nystrom_parameter,index(k)));
                
                end
                
            end
            
            % take transpose... do this because it is
            % convenient to compute the far field at each point for all
            % right hand sides in one go. This gives the transpose of val
            val = tmp.';

        end
        
        %-----------------------------------------------------------
        % matrix D such that D * f gives the tangential derivative of
        % the function at the Nystrom points, where f gives the values of
        % the function at the Nystrom points.
        %-----------------------------------------------------------

        % See p236 in [2]
            
        function matrix = setupDerivativeMatrix(self)
                       
            % set up quadrature points
            pp=pi*(0:2*self.nystrom_parameter-1)/self.nystrom_parameter;
            pp=pp(:);
            
            % set up matrix of quadrature points
            dt = repmat(pp,1,2*self.nystrom_parameter) - repmat(pp.',2*self.nystrom_parameter,1);
            
            % setup matrix of indices
            k = repmat((0:2*self.nystrom_parameter-1).',1,2*self.nystrom_parameter);
            
            % set diagonal of dt to dummy value
            dt = dt - diag(diag(dt)) + diag(ones(2*self.nystrom_parameter,1));
            
            % compute matrix
            matrix = 0.5*(-1).^(k-k.').*cot(0.5*dt);

            % correct the diagonal
            matrix = matrix - diag(diag(matrix));
            
        end
        
        
    end % end methods
    
end