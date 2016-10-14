%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program name: 	CD1 (convection-diffusion I)
% Summary:		Numerical solution of 1-D convection-diffusion equation with Dirichlet
% 			boundary conditions:
%			dv/dx - (1/Pe)d^2v/dx^2 = 0,  v(0) = a,   v(1) = b
% 			Cell centered grid: Central schemes I/II for convection, local grid refinement
%
% Remarks:		The program uses function 'exact_solution.m'
% 			Theory is given in Sect. 2.3 of
%			Computational Fluid Dynamics
%			Lecture Notes by P. Wesseling
%			Department of Applied Mathematical Analysis, ITS, TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		% .......................Coefficients for the differential operator..................%
a = 	0.2;	%  Value left Dirichlet boundary condition
b = 	1.0;	%  Value right Dirichlet boundary condition
pe = 	40;	%  Peclet number

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		% .......................Parameter for the numerical model..........................%
scheme =  1;	% Enter 1,2 or something else for upwind, central scheme (2.20)
		% or central scheme (2.21), rspectively.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	The program can use a uniform mesh, or local grid refinement in the boundary layer
%	where the refinement zone is scaled with the
%
% 	Definition of cell numbering
%
%                              |<---- del----->|
%     x=0       m1 cells       |    m2 cells  x=1
% grid |---o---|---o---|---o---|-o-|-o-|-o-|-o-|
%      1   1   2   2   3   3   4 4 5 5 6 6 7 7 8


		% .......................Parameters for the grid..........................%
m1 = 	  9;	% Number of cells outside refinement zone
m2 =      9;	% Number of cells inside refinement zone
refinement = 0;	% Enter 1 for for local grid refinement or something else for
		% uniform grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if refinement == 1,
	del = m2/pe;		% del = size of refinement zone, scales with Peclet number
else,
	del = m2/(m1+m2);	% No refinement: add the refinement zone to the normal grid:
				% it occupies the fraction m2/(m1+m2) of the original domain.
end


x = zeros(m1+m2+1,1);		   % x contains coordinates of cell boundaries
j = 1:(m1+1);
x(j) = (j-1)*(1-del)/m1;
j = (m1+2):(m1+m2+1);
x(j) = x(m1+1) + (j-1-m1)*del/m2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for the MATLAB newbies:
% for a vector X,DIFF(X) is [X(2)-X(1)  X(3)-X(2) ... X(n)-X(n-1)] (intrinsic matlab function)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dx = diff(x);			   % dx contains cell sizes: dx(j) = x(j+1)-x(j)
n = length(x)-1;
y = (x(1:n) + x(2:n+1))/2;	   % y contains coordinates of cell centers
dy = [dx(1)/2; diff(y); dx(n)/2];  % dy contains distances between cell centers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficients in numerical flux are called beta:
% F_{j+1/2} = beta0(j)v(j) + beta1(j+1)v(j+1), see equation (2.27)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

betadiffusion = 1./(pe*dy);	% Contribution of diffusion term

if scheme == 1			% Upwind scheme

  beta0 = 1 + betadiffusion(2:n+1);
  beta1 =   - betadiffusion(1:n);
  gamma0 = a*(1+(2/pe)/dx(1));			% Boundary corrections for
  gamma1 = -b*(2/pe)/dx(n);			% right-hand side

elseif scheme == 2		% Central scheme I ( averaging)

  beta0 = 1/2 + betadiffusion(2:n+1);
  beta1 = 1/2 - betadiffusion(1:n);
  beta0(n) = betadiffusion(n+1);
  beta1(1) = - betadiffusion(1);
  gamma0 = (1+(2/pe)/dx(1))*a;			% Boundary corrections for
  gamma1 = (1-(2/pe)/dx(n))*b;			% right-hand side

else				% Central scheme II (linear interpolation)

  beta0 = 0.5*dx(2:n)./dy(2:n) +  betadiffusion(2:n);
  beta1 = 0.5*dx(1:n-1)./dy(2:n) -  betadiffusion(2:n);
  beta0 = [beta0; betadiffusion(n+1)];
  beta1 = [ - betadiffusion(1); beta1];
  gamma0 = (1+(2/pe)/dx(1))*a;			% Boundary corrections for
  gamma1 = (1-(2/pe)/dx(n))*b;			% right-hand side

end

f = zeros(size(y));		% Right-hand side
f(1) = f(1) + gamma0; f(n) = f(n) - gamma1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Numerical scheme gives system Av = f.
% 	Diagonals of tridiagonal matrix A are stored in
%	[-beta0 beta0-beta1 beta1]:
% 	A(j,j-1) = -beta0(j-1);
%	A(j,j) = beta0(j) - beta1(j);
%	A(j,j+1) = beta1(j+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = spdiags([-beta0 beta0-beta1 beta1], -1:1, n, n);
numerical_solution = A\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of norm pf global error: difference between exact and numerical solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error = numerical_solution - exact_solution(y,pe,a,b);
n = length(y); norm1 = norm(error,1)/n
norm2 = norm(error,2)/sqrt(n)
norminf = norm(error,inf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the numerical solution (markers) together with the exact solution (line)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf, figure(1), hold on
plot(y,numerical_solution,'*','markersize',13)
fplot('exact_solution',[0 1],[],[],'-',pe,a,b)
h = findobj; set(h(3),'FontSize',18)
if scheme == 1,
	tt = 'Upwind scheme';
elseif scheme == 2,
	tt = 'Central scheme I: averaging';
else,
	tt = 'Central scheme II: linear interpolation';
end
title(['Pe=',num2str(pe),',  ',num2str(n),' cells, ',tt, ' |e|_2 =', num2str(norm2)],'fontsize',20)
%title(['Pe=',num2str(pe),',  Pe_h=',num2str(pe*min(dx)),',  ',tt, ' |e|_2 =', num2str(norm2)],'fontsize',20)
box on
