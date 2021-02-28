source("mfbavart_func.R")

# -----------------------------------------------------------------------------
# load dataset
library(alfred) # load dataset from FRED (Fed St. Louis)

variables <- c("CPIAUCSL","UNRATE","GDPC1")
out <- lapply(variables, get_alfred_series,observation_start = "1980-01-01",observation_end = "2020-12-01",realtime_start = "2020-12-01",realtime_end = "2020-12-01")
alfred_to_ts <- function(x, freq){
  ts(x[,3],start=c(1980,1),frequency=freq)
}
mf_list <- mapply(alfred_to_ts, x = out, freq = c(12, 12, 4))
names(mf_list) <- variables
log_diff <- function(x) {
  freq <- frequency(x)
  100 * freq * diff(log(x))}
mf_list[c("CPIAUCSL", "GDPC1")] <- lapply(mf_list[c("CPIAUCSL", "GDPC1")], log_diff)
mf_list <- mapply(window, x = mf_list, start = list(c(1980, 4), c(1980, 4), c(1980, 2)))
data <- mf_list

est_obj <- mfbavart(data,itr="grw",fhorz=2,prior.sig=c(200,0.75))

Yq <- est_obj$Yq
GDPC1_post <- t(apply(Yq,c(2,3),quantile,probs=c(0.16,0.5,0.84),na.rm=T)[,,"GDPC1"])
ts.plot(GDPC1_post)
