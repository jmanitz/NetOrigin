
#' Reads a data file as provided by 'Deutsche Bahn' (for internal use).
#'
#' @param file character with path and file name containing the variables for  'stationID', 'date', 'hour', 'minutes', and 'delay'
#' @return \code{data.frame} with variables \code{'node'}, \code{'time'}, \code{'delay'}
#'
#' @family data_handling
#' @importFrom utils read.table
#' @export
read_DB_data <- function(file){
    dat <- read.table(file, header=TRUE, sep=';', na.strings=c('NULL') ,
                      colClasses=c('character','character','numeric','numeric','numeric'))
    names(dat)[c(1,5)] <- c('node','delay')
  # specify date+time
    dat$time <- strftime(paste(dat[,2], ' ', dat[,3],':',dat[,4],sep=''), format='%Y-%m-%d %H:%M')

  # set negative delays to be zero
    dat$delay[which(dat$delay < 0)] <- 0

  # return value
    ret <- dat[,c('node', 'time', 'delay')]
    return(ret)
}

#' convert individual event information to aggregated information per network node
#'
#' @param dat \code{data.frame} with variables \code{'node'}, \code{'time'}, \code{'delay'}, events data with single events with count magnitude
#' @param from character in \code{\link{strftime}} format, e.g. \code{"2014-06-12 16:15"}, data is subsetted accordingly before aggregation
#' @param cumsum logical indicating whether data is aggregated by cumulative sum, default is \code{TRUE}
#' @return \code{data.frame} of dimension \code{(TxK)}, where \code{T} is the number of observation times and \code{K} the number of network nodes. Thus, each row represents a snapshot of the spreading process at a specific observation time with the event magnitude observed at the network nodes. Rownames are observation times, colnames are node names.
#'
#' @family data_handling
#' @importFrom stats aggregate reshape
#' @export
aggr_data <- function(dat, from = NULL, cumsum = TRUE){#'%Y-%m-%d %H:%M'

    # remove trains before starting time
    if(!is.null(from)){
        day <- min(as.Date(dat$time)) # day of the events
        from2 <- strftime(paste(day, from), format='%Y-%m-%d %H:%M')
        dat <- dat[!dat$time < from2,]
    }

    # aggregate total delays in a station
    datTot <- aggregate(dat$delay, by=list(dat$node,dat$time), FUN=sum, na.rm=TRUE)
    resTot <- reshape(datTot, direction='wide', timevar='Group.1', idvar='Group.2')
    resTot[is.na(resTot)] <- 0

    # aggregate number of delayed trains in a station
    dat$count <- 1
    datNo <- aggregate(dat$count, by=list(dat$node,dat$time), FUN=sum, na.rm=TRUE)
    resNo <- reshape(datNo, direction='wide', timevar='Group.1', idvar='Group.2')
    resNo[is.na(resNo)] <- 0

    # cumulative aggregation and normalization by the number of trains
    if(cumsum){
        resTotC <- apply(resTot[,-1],MARGIN=2,FUN=cumsum)
        resNoC  <- apply(resNo[,-1] ,MARGIN=2,FUN=cumsum)
        
        # normalize cumulative total 
        relDel <- resTotC/resNoC
    }else{
        # normalize total delay 
        relDel <- resTot[,-1]/resNo[,-1]
    }

    # return value
    relDel[is.nan(as.matrix(relDel))] <- 0
    rownames(relDel) <- resNo[,1] # timevar
    return(relDel)
}

