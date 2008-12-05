function simul(dr)
% function simul(dr)
% computes simulations
%  
% INPUTS
%   dr: structure of decision rules for stochastic simulations
%  
% OUTPUTS
%   ...
% ALGORITHM
%   ...
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 1996-2007 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

global M_ options_ oo_ 
global ys0_ ct_

if size(M_.lead_lag_incidence,2)-nnz(M_.lead_lag_incidence(M_.maximum_endo_lag+1,:)) > 0
  mess = ['DYNARE: error in model specification : variable ' M_.endo_names(find(M_.lead_lag_incidence(M_.maximum_lag+1,:)==0),:)] ;
  mess = [mess ' doesn''t appear as current variable.'] ; 
  error (mess) ;
end

options_ = set_default_option(options_,'simul_algo',0);
options_ = set_default_option(options_,'dynatol',0.00001);
options_ = set_default_option(options_,'maxit',10);
options_ = set_default_option(options_,'slowc',1);
options_ = set_default_option(options_,'timing',0);
options_ = set_default_option(options_,'gstep',1e-2);
options_ = set_default_option(options_,'scalv',1);
if ~isfield(options_,'periods') & ~isempty(options_.periods)
  options_.periods = options_.periods
end
options_ = set_default_option(options_,'periods',0);
if options_.periods == 0
  error('SIMUL: number of periods for the simulation isn''t specified')
end
options_.periods = options_.periods;
ct_=0;

if options_.simul_algo == 0
  if ~ options_.initval_file
      make_ex_;
      make_y_;
  end

  if isempty(options_.scalv) | options_.scalv == 0
    options_.scalv = oo_.steady_state ;
  end

  options_.scalv= 1 ;

  if M_.maximum_endo_lag ==1 & M_.maximum_endo_lead <= 1
    sim1 ;
  else
    simk ;
  end
else
  set_default_option('replic',1);
  set_default_option('simul_seed',1);
  if isfield(dr,'ghxx')
    set_default_option('order',2);
  else
    set_defaut_option('order',1);
  end
  oo_.endo_simul=simult(oo_.steady_state,dr,options_);
end

dyn2vec;

% 6/18/01 MJ added dyn2vec if 40 variables or less
% 01/16/03 MJ use dyn2vec whatever the number of variables
% 02/18/03 MJ added oo_.steady_state for calling simult
% 05/24/03 MJ added options_ and options_.periods









