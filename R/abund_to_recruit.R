
#' Convert abundance estimates to a daily recruitment profile
#' for one region/year
#'
#' @param abund An abundance record as returned from clean_data
#' @param phi=0.6
#'
#' @return a dataframe with estimated daily abundance profile.  Note that these
#' are not trimmed to nonnegative values.
#' @export
#'
#' @examples
#' library(readxl)
#' abund <- clean_data(read_excel("annualabundances/2015.xls"))
#' upu <- filter(abund,region=="Uncompahgre Upper")
#' abund2recruit(upu)
#' qplot(date,recruit,data=abund2recruit(upu))
#' all_recruit_profile <- abund %>% group_by(region) %>% do(abund2recruit(.))
#' qplot(date,recruit,data=all_recruit_profile,color=region,geom="line")
abund2recruit <- function(abund,phi=0.6){
  nab <- nrow(abund)
  stopifnot(!is.unsorted(abund$date))
  rdates <- seq(from=abund$date[1], to=abund$date[nab],by="day")
  ab_ind <- abund$date - abund$date[1]+1

  # survival matrix
  smat <- outer(seq(length(rdates)),seq(length(rdates)),
                function(i,j)ifelse(i-j<0,0,phi^(i-j)))
  # gradient of objective function
  bmat <- outer(seq(length(rdates)),seq(length(rdates)),
                function(i,j)ifelse(abs(i-j)==1,-1,ifelse(i==j,2,0)))
  # fix the corners
  bmat[1,1] <- 1
  bmat[length(rdates),length(rdates)] <- 1
  # constraint matrix
  a <- smat[ab_ind,]
  kmat <- cbind(rbind(bmat,a),rbind(t(a),matrix(0,nrow = nab,ncol=nab)))
  bvect <- c(rep(0,length(rdates)),abund$est)
  data.frame(date=rdates,
             recruit=solve(kmat,bvect)[1:length(rdates)]
  )
}

#' Calculate total recruitment for one region/year
#'
#' @param abund An abundance record as returned from clean_data
#' @param phi=0.6
#'
#' @return The sum of the nonnegative recruitments
#' @export
#'
#' @examples
#' abund <- clean_data(read_excel("annualabundances/2015.xls"))
#' upu <- filter(abund,region=="Uncompahgre Upper")
#' total_recruit(upu)
#' all_recruit <- abund %>% group_by(region) %>% do(data_frame(total=total_recruit(.)))

total_recruit <- function(abund,phi=0.6){
  rc <-abund2recruit(abund,phi)$recruit
  sum(ifelse(rc<0,0,rc))
}

#' Make a bootstrap sample from an abundance record for one region/year
#'
#' @param abund A dataframe with columns est, lower, upper
#'
#' @return a similarly formatted dataframe with estimates varied by
#' log normal distributions based on the original estimate and upper
#' confidence level
#' @export
#'
#' @examples
#' abund <- clean_data(read_excel("annualabundances/2015.xls"))
#' upu <- filter(abund,region=="Uncompahgre Upper")
#' makebootstrap(upu)
#' apply(replicate(1000, makebootstrap(upu)$est),1,
#'  function(dat)quantile(dat,probs=c(0.025,0.975)))
makebootstrap <- function(abund){
  makelnparams <- function(est,upper){
    list(lmu=log(est),lsd=log(upper/est)/qnorm(0.975))
  }
  params <- makelnparams(abund$est,abund$upper)
  goodparams <- abund$est >0 # if est=0 we return est=0
  bootabund <- rep(0,nrow(abund))
  bootabund[goodparams] <- rlnorm(sum(goodparams),meanlog=params$lmu[goodparams],
                                  sdlog=params$lsd[goodparams])
  abund$est <- bootabund
  abund
}

#' Calculate a 95% confidence interval for a single region/year using
#' bootstrapping.
#'
#' @param nboot number of replications (typically 1000)
#' @param abund an abundance record for one region/year
#' @param phi=0.6
#'
#' @return a dataframe with the estimate and upper and lower confidence
#' limits
#' @export
#'
#' @examples
#' abund <- clean_data(read_excel("annualabundances/2015.xls"))
#' upu <- filter(abund,region=="Uncompahgre Upper")
#' boot_ci(1000,upu)
#' abund %>% group_by(region) %>% do(boot_ci(1000,.))

boot_ci<- function(nboot,abund,phi=0.6){
  boots <- replicate(nboot,total_recruit(makebootstrap(abund),phi=phi))
  ci <- quantile(boots,probs = c(0.025,0.975))
  year <- as.numeric(strftime(abund$date[1], "%Y"))
  region <- abund[[1,"region"]]
  data_frame(year=year,
             region=region,
             est=total_recruit(abund),
             lcl=ci[1],
             ucl=ci[2]
  )
}
