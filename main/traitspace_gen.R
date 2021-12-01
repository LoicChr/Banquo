##############
#
#     New traitspace objects - linear models
#
##############

load("data/data_ready_all.Rdata")
source("lib/traitspace_and_banquo.R")

## Trait model
multi.mod <- lm(cbind(poros, sla, height) ~ DAYSUB, data = na.omit(data), weights = Cover)
P_S_E_tr <- traitspace(multi.mod, env_p, pdf_species, N = 500, avg.out = T)

save(list = c("P_S_E_tr", 'env_p', "obs.comm", 'pdf_species',"multi.mod","t.avg"), file = "data/traitspace_objs_all.Rdata")

