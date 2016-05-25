#' Use quantile mapping to correct bias
#'
#' The function is a simple wrapper around quantile mapping functions provided in the \code{\link{qmap}} package. Quantile mapping is done on each variable included in the input data.table, the mapping is specific for each month.
#'
#' @param obs_ctrl observed data for the control period
#' @param sim_ctrl simulated data for the control period
#' @param sim_scen scenario data
#' @param ... other arguments passed to fitQmap
#'
#' @return data.table with corrected scenario simulation
#' @export correct
#' @details The correction is calibrated by comparison of \code{obs_ctrl} and \code{sim_ctrl}. Quantile mapping based on empirical quantiles is considered. The correction for wet.day frequency is applied for precipitation.
#' @examples
#' data(basin_PT)
#'
#' # correct future simulation
#' cor_future = correct(obs_ctrl = basin_PT[['obs_ctrl']], sim_ctrl = basin_PT[['sim_ctrl']], sim_scen =  basin_PT[['sim_scen']])
#'
#' # get the corrected control simulation
#' cor_future = correct(obs_ctrl = basin_PT[['obs_ctrl']], sim_ctrl = basin_PT[['sim_ctrl']])
correct = function(obs_ctrl, sim_ctrl, sim_scen = copy(sim_ctrl), ...){

  mx = melt(obs_ctrl, id.vars = 'DTM', value.name = 'obs')
  my = melt(sim_ctrl, id.vars = 'DTM', value.name = 'sim')
  mz = melt(sim_scen, id.vars = 'DTM', value.name = 'pred')
  m = mx[my, on = c('DTM', 'variable')]

  mon = 1
  v = 'PR'
  for (v in m[, unique(variable)]){
    for (mon in 1:12){
      f = m[variable==v & month(DTM) == mon, fitQmap(obs, sim, method = 'QUANT', wet.day = v == 'PR', ...)]
      mz[variable==v & month(DTM) == mon , cor:=doQmap(pred, f)]
    }
  }

  copy(dcast(mz[, .(DTM, variable, cor)], DTM ~ variable, value.var = 'cor'))
}


