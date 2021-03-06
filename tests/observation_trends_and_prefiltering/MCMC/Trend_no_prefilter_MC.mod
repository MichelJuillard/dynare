@#include "../Trend_model_no_prefilter_common.inc"

estimation(order=1,datafile='../AR1_trend_data_with_constant',mh_replic=2000,
            mode_compute=4,first_obs=1,smoother,mh_nblocks=1,mh_jscale=0.3,
            filtered_vars, filter_step_ahead = [1,2,4],
            mcmc_jumping_covariance='MCMC_jump_covar',forecast=100,prefilter=0) P_obs Y_obs junk2;
            
load('../AR1_trend_data_with_constant');
@#include "../Trend_load_data_common.inc" 

loaded_par=load('../orig_params');
if max(abs((M_.params-loaded_par.orig_params)./loaded_par.orig_params))>0.03
    error('Parameter estimates do not match')
end
y_forecast_100_periods=loaded_par.orig_params(strmatch('const_y',loaded_par.param_names,'exact'))+(options_.first_obs+options_.nobs-1+options_.forecast)*loaded_par.orig_params(strmatch('g_y',loaded_par.param_names,'exact'));
p_forecast_100_periods=loaded_par.orig_params(strmatch('const_p',loaded_par.param_names,'exact'))+(options_.first_obs+options_.nobs-1+options_.forecast)*loaded_par.orig_params(strmatch('g_p',loaded_par.param_names,'exact'));

@#include "../Trend_diagnostics_MCMC_common.inc"