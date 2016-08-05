library(musica)

data("basin_PT")

scen = basin_PT$sim_scen
ctrl = basin_PT$sim_ctrl
obs = basin_PT$obs_ctrl

rev_difv = function(x, v){
  if (v %in% getOption('additive_variables')) {return(sum(x))} else {return(prod(x))}

}

abb = month
#dta = rbind(data.table(RID = 'SCEN', scen), data.table(RID = 'CTRL', ctrl), data.table(RID = 'OBS', obs), data.table(RID = 'MODIF', obs))
#dta = rbind(data.table(RID = 'TO', obs), data.table(RID = 'FROM', ctrl), data.table(RID = 'NEWDATA', scen), data.table(RID = 'MODIF', scen))
dta = rbind(data.table(RID = 'TO', obs), data.table(RID = 'FROM', ctrl), data.table(RID = 'NEWDATA', ctrl), data.table(RID = 'MODIF', ctrl))
mdta = melt(dta, id.vars = c('RID', 'DTM'), value.name = 'ORIG')

mdta[, G1 := mean(ORIG), by = .(RID, variable)]
mdta[, Y1 := dif(mean(ORIG), G1, variable[1]), by = .(RID, variable, year(DTM))]
mdta[, M1 := dif(mean(ORIG), rev_dif(G1, Y1, variable[1]), variable[1]), by = .(RID, variable, year(DTM), month(DTM))]
mdta[, D1 := dif(ORIG, rev_difv(c(G1, Y1, M1), variable[1]), variable[1]), by = .(RID, variable, DTM)]
mdta[, GRP := month(DTM)]
#mdta[, C := rev_difv(c(G1, Y1, M1, D1), variable[1]), by = .(RID, variable, DTM)]

pr = c('G1', 'Y1', 'M1', 'D1')

for (i in 1:length(pr)){

}
