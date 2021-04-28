#' RedLosses Multiblock data
#'
#' Data from RedLosses project.
#'
#' @docType data
#'
#' @usage data(RedLosses)
#'
#' @format Dataframe gathering Metagenomic, Volatilome and Sensory blocks from RedLosses project.
#'
#' @keywords datasets
#'
#' @references Luong et al. (2021) International Journal of Food Microbiology
#'
#' @examples
#' data(RedLosses)
#' Group  # identifying the three blocks
#' delta  # matrix associated to the path-diagram with directed links
#' J <- rep(1:length(Group),times=Group)
#' metag <- data[,J==1]  # metagenomic block
#' volat <- data[,J==2]  # volatilome block
#' senso <- data[,J==3]  # sensory block
"RedLosses"
