#' musica - An R package for multiscale climate model assessment
#'
#' The package provides functions for flexible assessment of climate model bias and changes at multiple time scales. See documentation for \code{\link{decomp}}, \code{\link{compare}} and \code{\link{vcompare}}.
#'
#' @section Package options:
#' Following option(s) are available:
#' \describe{
#'  \item{additive_variables}{At several places the package compares values. The character vector \code{additive_values} specifies for which variables difference should be used for comparison instead of ratio. Defaults to \code{additive_values = "TAS"}. See \code{\link{options}} for setting or examining options.}
#' }
#'
#'
#' @author Martin Hanel \email{hanel@@fzp.czu.cz}
#' @references Hanel, M., Kozin, R. (2016) Bias correction for hydrological modelling, submitted.
#' @docType package
#' @name musica-package
NULL


.onLoad <- function(libname, pkgname) {
  op <- options()
  op.musica <- list(
    'additive_variables' = 'TAS'
  )
  toset <- !(names(op.musica) %in% names(op))
  if(any(toset)) options(op.musica[toset])

  invisible()
}




#' @export
dif = function(x, y, variable){
  if (variable %in% getOption('additive_variables')) return(x-y)
  r = x/y
  r[y==0] = 0
  return(r)
}

#' @export
rev_dif = function(x, y, variable){
  if (variable %in% getOption('additive_variables')) return(x+y)
  r = x*y
  #r[y==0] = 0
  return(r)
}

#' @export
rev_difv = function(x, v){
  if (v %in% getOption('additive_variables')) {return(sum(x))} else {return(prod(x))}

}

#' @export
prob = function(x){
  (rank(x)-.3) / (length(x)+.4)
}

#' @export
vprod = function(x, variable, ...){
  if (variable %in% getOption('additive_variables')) return(sum(x, ...)) else (return(prod(x, ...)))
}

#' @export
period2code = function(periods){
  num = suppressWarnings(as.integer(gsub("[^\\d]+", "", periods, perl = TRUE)))
  if (any(is.na(num))) stop("Invalid scales specification - number of intervals must be provided. See ?decomp for details.")
  let = gsub('day|days', 'D', periods) %>% gsub('month|months', 'M', .) %>% gsub('year|years', 'Y', .) %>% gsub("[\\d]+| ", '', ., perl = TRUE)
  paste0(let, num)
}


code2period = function(code){

  num = suppressWarnings(as.integer(gsub("[^\\d]+", "", code, perl = TRUE)))
  per = gsub("[\\d]+| ", '', code, perl = TRUE) %>% gsub('Y', 'year', .) %>% gsub('M', 'month', .) %>% gsub('D', 'day', .)

  out = paste(num, per)
  names(out) = code
  out[code %in% c('G0', 'G1')] = '0'
  out
}


#' Decomposition of time-series
#'
#' Enhance the input data.table with series of averages over the periods specified in the \code{periods} argument.
#'
#' @param x data.table with columns \code{DTM} (date), \code{variable} and \code{value}. Any number of variables are in principle allowed.
#' @param year_starts A Period object indicating the start of the year, e.g. \code{year_starts = -months(1)} results in climatological seasons
#' @param periods The periods over which the averages will be calculated, see Details
#' @param agg_by Function for specification of the period (season, month) to be additionaly included in output, see Details
#'
#' @details The original time series in daily time step is decomposed into series of averages ove periods specified in \code{periods} argument using keywords "day", "month", "year" preceded by an integer and a space and optionally followed by "s" (the specification is further passed to \code{cut.Date}, see \code{\link{cut.Date}} for details). The periods must be given in order from longest to shortest, the overall mean is always included (and needs not to be included). Shorter periods are always identified within the closest longer periods, i.e. each shorter period is fully included in exactly one longer period. As a result, the averages may be calculated over shorter periods than specified. This is due to varying length of "month" and "year" periods. The actual length used for averaging is included in the output. To make further assessment of the decomposed objects easier, indicator of period within the year (e.g. quarter or month) as specified by \code{agg_by} argument is included in the output.
#'
#'@return data.table with variables:
#'\describe{
#'  \item{variable}{factor indicating the variable}
#'  \item{DTM}{date}
#'  \item{period}{specification of the averaging length with `D` - day(s), `M` - month(s), `Y` - year(s) and `G1` - the overall mean }
#'  \item{value}{value of the variable for given averaging length}
#'  \item{sub_period}{indication of the aggregating scale specified by \code{agg_by} argument}
#'  \item{period_pos}{average date of the interval}
#'  \item{N}{real length of the vectors used for calculating averages}
#'  \item{TS}{averaging length in hours}
#'}
#' @export decomp
#'
#' @examples
#' data(basin_PT)
#' str(basin_PT)
#' basin_PT[['obs_ctrl']]
#' dobs = decomp(basin_PT[['obs_ctrl']], scales = c('1 year', '1 month', '1 day'))
decomp = function(x, year_starts = months(0), periods = c('1 year', '6 months', '3 months', '1 month', '15 days', '1 day'), agg_by = quarter, full_return = FALSE, remove_incomplete = TRUE){

  if (!all(grepl('year|month|day', periods))) stop("periods argument must be specified using 'year', 'month', 'day' keywords only. See ?decomp for details.")
  names(periods) = period2code(periods)
  x = copy(x)

  # shift for arbitrary year start (climatological, hydrological etc.)
  x[, cDTM := as.IDate(DTM %m-% year_starts)]

  if (remove_incomplete){
    x[, N:=.N, by = .(year(cDTM))]
    x = x[N>364]
    x[, N:=NULL]
  }

  for (i in 1:length(periods)){
    bby = if (i>1) {names(periods)[1:(i-1)]} else {NULL}
    x[, (names(periods)[i]) := cut(cDTM, breaks = periods[i]), by = bby]
  }

  wh = names(x)[! names(x) %in% c('DTM', 'cDTM', names(periods))]
  x[, G1:= mean(cDTM)]
  mo = melt(x, measure.vars = wh)

  R = list()
  for (i in 0:length(periods)){

    bby = if (i>0) {names(periods)[i]} else {'G1'}
    ab = agg_by
    #if ((identical(agg_by, quarter) & (tscale(bby) <= tscale('M3'))) | (identical(agg_by, month) & (tscale(bby) <= tscale('M1')))) {agg_by} else {function(dtm){1}}
    R[[i+1]] = if (!full_return) {
      mo[,  .(period = bby, value = mean(value), sub_period = unique(ab(cDTM)), period_pos = unique(prse), N = .N), by = .(variable, prse = eval(parse(text = bby)))]
    } else {
      mo[,  .(dDTM = DTM, period = bby, value = mean(value), sub_period = (ab(cDTM)), period_pos = unique(prse), N = .N), by = .(variable, prse = eval(parse(text = bby)))]
    }

    }

  R = do.call(rbind, R)
  R[, TS:= tscale(period)]
  setnames(R, 'prse', 'DTM')
  R[, DTM := as.IDate(DTM %m+% year_starts)]
  #R[, period_pos := as.IDate(period_pos %m+% year_starts)] # NEVYHODIT ? DTM == period_pos
  R[, period_pos := as.IDate(DTM) + ifelse(period!='G1', tscale(period)/24/2, 0)]
  #R[, excl:= N < 0.9 * (TS/24)]
  #if (year_starts < months(0) ) {return(R[year(DTM)!=min(year(DTM))])}

  #R[, DUPL:=duplicated(value), ]

  copy(R)

}

#' Convert averaging length code to hours
#'
#' Period durations are calculated by the \code{\link{lubridate}} package.
#'
#' @param x Vector of the averaging period codes
#' @param nyears Overall number of years - used for conversion of the overall mean
#'
#' @return numerical vector of durations in hours
#' @export tscale
#'
#' @examples
#' tscale('M1')
#' tscale('G1', nyears = 25)
tscale = function(x, nyears = 30){

  num = suppressWarnings(as.integer(gsub("[^\\d]+", "", x, perl = TRUE)))
  xx = gsub('G', 'year ',  x)
  xx = gsub('Y', 'year ', xx)
  xx = gsub('M', 'month ', xx)
  xx = gsub('D', 'day ', xx)
  xx[grepl('G', x)] = gsub(1, nyears, xx[grepl('G', x)])

  xs = strsplit(xx, ' ')
  xs = data.table(x, do.call(rbind, xs))
  xs[, ts:= unclass(duration(as.double(V2[1]), V1[1]))/3600, by = x]
  return(xs$ts)

}

# name_scales = function(scales){
#   num = suppressWarnings(as.integer(gsub("[^\\d]+", "", scales, perl = TRUE)))
#   if (any(is.na(num))) stop("Invalid scales specification - number of intervals must be provided. See ?decomp for details.")
#   let = gsub('day|days', 'D', scales) %>% gsub('month|months', 'M', .) %>% gsub('year|years', 'Y', .) %>% gsub("[\\d]+| ", '', ., perl = TRUE)
#   paste0(let, num)
# }

#' Compare decomposed variables
#'
#' @param x List of decomposed variables to be compared
#' @param compare_to Decomposed variable used as a reference
#' @param fun Function used for comparison
#' @param wet_int_only (logical) Should only the wet intervals be considered?
#' @param wet_int_thr Numeric value specifying the minimum depth to be consider wet
#' @param exclude_below Some of the intervals might not be of required length, e.g. D10 interval may have less than 10 days available. The \code{exclude_below} argument controls the minimum fraction of the interval that has to be available in order to be considered in the summary statistics.
#'
#' @return data.table summarizing the differences
#' @export compare
#'
#' @examples
#' data(basin_PT)
#' dobs = decomp(basin_PT[['obs_ctrl']])
#'
compare = function(x, compare_to, fun = mean, wet_int_only = TRUE, wet_int_thr = 0.1, exclude_below = 0.9){

  lst = c(x, COMPARE_TO = list(compare_to))
  for (i in 1:length(lst)){
    lst[[i]][, ID:=names(lst)[i]]
  }
  lst = do.call(rbind, lst)
  lst = lst[N >= exclude_below * (TS/24)]
  stat =
    if (wet_int_only){
      lst[!(! variable  %in% getOption('additive_variables') & value<= wet_int_thr) , .(value = fun(value)), by = .(variable, period, TS, sub_period, ID)]
    } else {
      lst[, .(value = fun(value)), by = .(variable, period, TS, sub_period, ID)]
    }
  cstat = dcast.data.table(stat, variable + period + TS + sub_period ~ ID)
  stat = melt(cstat, measure.vars = names(x), variable.name = 'comp')
  stat[, .(DIF = dif(value, COMPARE_TO, variable[1])), by = .(variable, period, TS, sub_period, comp)]

}

#' Assess the relations between two decomposed variables
#'
#' @param x Decomposed object
#' @param fun Function to sumarize dependence (like \code{cor}, \code{cov})
#' @param wet_int_only (logical) Should only the wet intervals be considered?
#' @param wet_int_thr Numeric value specifying the minimum depth to be consider wet
#' @param exclude_below Some of the intervals might not be of required length, e.g. D10 interval may have less than 10 days available. The \code{exclude_below} argument controls the minimum fraction of the interval that has to be available in order to be considered in the summary statistics.
#'
#' @return data.table summarizing the relation
#' @export vcompare
#'
#' @examples
vcompare = function(x, fun = cor, wet_int_only = TRUE, wet_int_thr = 0.1, exclude_below = 0.9){

  lst = x #c(x, COMPARE_TO = list(compare_to))
  for (i in 1:length(lst)){
    lst[[i]][, ID:=names(lst)[i]]
  }
  lst = do.call(rbind, lst)

  vars = lst[, unique(variable)]
  cmb = t(combn(vars, 2))
  map = data.table(TYP = rep(names(x),  each = nrow(cmb)), cmb)

  B = list()
  for (i in 1:nrow(map)){
    yy = dcast(lst[ID == map[i, TYP] & N >= exclude_below * (TS/24), ], ... ~ variable, value.var = 'value')
    cy = yy[, .(ID = map[i, TYP], V1 = map[i, V1], V2 = map[i, V2], TS = tscale(period), value = fun(eval(parse(text = map[i, V1])),  eval(parse(text = map[i, V2]))), VARS = map[i, paste(V1, V2, sep = ' x ')]), by = .(period, sub_period) ]
    B[[i]] = cy
  }

  do.call(rbind, B)

}



#' Convenience function for calculation of quantiles
#'
#' @param p Specification of the quantile
#' @param ... other arguments passed to \code{\link{quantile}}
#'
#' @return function calculating the p-th quantile
#' @export Q
#'
#' @examples
Q = function(p, ...){
  function(x)quantile(x, p, ...)
}

#' @export
run_bil = function(dta, b = NULL, vars = 'RM'){
  if (is.null(b)) {b = bil.new(type = 'd', data = dta[, .(DTM, P = PR, T = TAS)]) } else {
    bil.set.values(b, dta[, .(DTM, P = PR, T = TAS)])
  }
  bil.pet(b)
  d = data.frame(data.frame(bil.run(b))[, vars])
  names(d) = vars
  d
}

#' @export
month2sea = function(dtm, year_starts){
  i = ((1:12-1) + month(as.Date('1970-01-01') + year_starts)) %% 12
  i[i==0] = 12
  id = c('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N' , 'D')[i]
  id = rep(c(paste0(id[1:3], collapse = ''), paste0(id[1:3+3], collapse = ''), paste0(id[1:3+6], collapse = ''), paste0(id[1:3+9], collapse = '')), each = 3)
  i = month(dtm %m-% year_starts)#((month(dtm)-1) + month(as.Date('1970-01-01') + year_starts)) %% 12
  i[i==0] = 12
  id[i]
}

#' @export
sscale2sea = function(sub_scale, year_starts){

  i = ((1:12-1) + month(as.Date('1970-01-01') + year_starts)) %% 12
  i[i==0] = 12
  if (length(unique(sub_scale)) == 12) return( (c('DJF', 'DJF', 'MAM', 'MAM', 'MAM', 'JJA', 'JJA', 'JJA', 'SON', 'SON', 'SON', 'DJF')[i])[sub_scale]  )

  id = c('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N' , 'D')[i]
  id = c(paste0(id[1:3], collapse = ''), paste0(id[1:3+3], collapse = ''), paste0(id[1:3+6], collapse = ''), paste0(id[1:3+9], collapse = ''))
  # i = month(dtm %m-% year_starts)#((month(dtm)-1) + month(as.Date('1970-01-01') + year_starts)) %% 12
  # i[i==0] = 12
  id[sub_scale]
}

#' @export
period2code = function(periods){
  num = suppressWarnings(as.integer(gsub("[^\\d]+", "", periods, perl = TRUE)))
  if (any(is.na(num))) stop("Invalid scales specification - number of intervals must be provided. See ?decomp for details.")
  let = gsub('day|days', 'D', periods) %>% gsub('month|months', 'M', .) %>% gsub('year|years', 'Y', .) %>% gsub("[\\d]+| ", '', ., perl = TRUE)
  paste0(let, num)
}

#' @export
code2period = function(code){

  num = suppressWarnings(as.integer(gsub("[^\\d]+", "", code, perl = TRUE)))
  per = gsub("[\\d]+| ", '', code, perl = TRUE) %>% gsub('Y', 'year', .) %>% gsub('M', 'month', .) %>% gsub('D', 'day', .)

  out = paste(num, per)
  names(out) = code
  out[code %in% c('G0', 'G1')] = '0'
  out
}

