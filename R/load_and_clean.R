
#' Take an excel file recording annual observations and build a data frame
#' appropriate for converting to distance software
#'
#' @param filename
#'
#' @return A data frame with appropriate columns
#' @export
#'
#' @examples
#'
#' library(readxl)
#' read_ufb_datasheet("Distance Input 2015.xlsx")
#' write.csv(di,row.names=FALSE,na="",file="distanceinput2015.csv", eol = "\r\n")
read_ufb_datasheet <- function(filename){
  sheets <- excel_sheets(filename)
  # for consitency use title case for all regions
  # check that the input file has sheets in proper order
  regionlup <- c("Uncompahgre", "Uncompahgre","Uncompahgre","Redcloud","Redcloud","Baldy Chato","Baldy Chato")
  names(regionlup) <- sheets
  sublup <- c("Upper","Middle","Lower","Upper","Lower","North","South")
  names(sublup) <- sheets
  lup_names <- c("Baldy Chato North", "Baldy Chato South", "Redcloud Lower", "Redcloud Upper",
                 "Uncompahgre Lower", "Uncompahgre Middle", "Uncompahgre Upper" )
  linelengths <- c(385, 127, 371, 474, 457, 751, 583)
  names(linelengths) <- lup_names
  area <- c(23199, 8190, 21360, 29490, 28440, 41490, 32379)
  names(area) <- lup_names
  tometers <- function(dist)ifelse(dist==1,0.25,
                                 ifelse(dist==2,0.75,
                                        ifelse(dist==3,1.25,
                                               ifelse(dist==4,1.75,
                                                      ifelse(dist==5,2.5,
                                                             ifelse(dist==6,4.0,NA))))))


  read_sheet <- function(sname){
    dat <- read_excel(filename,sheet=sname,col_types = c("text","text","text","text","text"))[1:4]
    dat <- dat %>% filter(!is.na(DATE),DATE != "DATE") # remove blank lines

    dat <- dat%>%
      mutate(DIST = as.numeric(DIST),
             OBS = trimws(OBS,which="both"), # remove spare white space in observer
             date=as.Date(as.numeric(DATE),origin = as.Date("1899-12-30")),
             region=regionlup[sname],
             sub=sublup[sname],
             perp=tometers(DIST),
             linelength=linelengths[paste(region,sub)],
             area=area[paste(region,sub)])

    dat
  }
  # return all sheets as a data frame
  bind_rows(lapply(sheets,read_sheet))
}

#' Format a data frame produced by read_ufb_datasheet for importing
#' into distance software
#'
#' @param dat
#'
#' @return a data frame with appropriate naming and date conventions for
#' distance software
#' @export
#'
#' @examples
#' library(readxl)
#' ds <- read_ufb_datasheet("Distance Input 2015.xlsx")
#' dsf <- format_for_distance(ds)
#' write.csv(dsf,row.names=FALSE,na="",file="distanceinput2015.csv", eol = "\r\n")
format_for_distance <- function(dat){
  dat %>% transmute("Region*date"=strftime(date,format="%m/%d/%y"),
                    "Observation*location"=LINE,
                    "Observation*observer"=OBS,
                    "Observation*Perp distance"=perp,
                    "Line transect*line length"=linelength,
                    "Region*area"=area,
                    "Region*label"=region,
                    "Region*subregion"=sub)
}

#format_for_distance(read_ufb_datasheet("Distance Input 2015.xlsx"))

#' Find daily butterfly counts for each regioin
#'
#' @param di A dataframe produced by read_ufb_datasheet
#'
#' @return A dataframe with daily observation counts for each region.
#' This is appropriate for generating label columns when joined with
#' daily abundance estimates from distance software.  Note that the
#' two subregions for Baldy Chato are reversed in the output of distance
#' so they need to be manually corrected.
#' @export
#'
#' @examples
#' tt <- read_ufb_datasheet("Distance Input 2015.xlsx")
#' make_counts(tt)

make_counts <- function(di){
  di %>% group_by(region,sub,date)  %>%
    summarize(count=sum(perp>0,na.rm=T))%>% ungroup() %>%
    arrange(desc(region),desc(sub))
}

#' clean_data prepares data from Distance software for analysis
#'
#' @param df A dataframe obtained from read_excel on an annual abundance
#' file.  The format for this file is: date region sub_region est lower upper
#' where date is in ISO: YYYY-MM-DD, region is main region (title case) and
#' sub_region is title case.
#'
#'
#'
#' @return A cleaned version with the dates converted from integers
#' and a single region column.  Note that region and sub_region should
#' be in title case.
#' @export
#'
#' @examples
#' library(readxl)
#' clean_data(read_excel("filename"))
clean_data <- function(df){
  # origin should be 1900-01-01.  I don't know what is wrong.
  df$date <- as.Date(df$date,origin=as.Date("1899-12-30"))
  df$region <- paste(df$region,df$sub_region)
  df$sub_region <- NULL
  df
}


