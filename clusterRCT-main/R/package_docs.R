



#'@title Fake Cluster Randomized Trial data
#'
#'@description This dataset is a three-level dataset with students in
#'  schools in districts, with schools randomized to treatment and
#'  control.  It is used to illustrate and test the clusterRCT pacakge.
#'
#'@format A data frame with 1500 rows and 7 variables: \describe{
#'  \item{\code{V.k}}{double District-level covariate (think standardized district-wide
#'   average test score relative to state).}
#'  \item{\code{X.jk}}{double School-level (think percent on Free/Reduced Price Lunch).}
#'  \item{\code{C.ijk}}{double Student-level covariate (think baseline measured SES).}
#'  \item{\code{S.id}}{integer School ID.}
#'  \item{\code{D.id}}{integer District ID.}
#'  \item{\code{Yobs}}{double Observed outcome (think test score).)}
#'  \item{\code{T.x}}{integer Treatment assignment (1 treated, 0 control).} }
#'@details These data were generated via the PUMP package.
"fakeCRT"


#'@title Fake Broken Cluster Randomized Trial data
#'
#'@description This dataset is a version of `fakeCRT` with missing data and some
#' blocks that are all tx or all co.
#'
#'@format A data frame with 1500 rows and 7 variables: \describe{
#'  \item{\code{V.k}}{double District-level covariate (think standardized district-wide
#'  average test score relative to state).}
#'  \item{\code{X.jk}}{double School-level (think percent on Free/Reduced Price Lunch).}
#'  \item{\code{C.ijk}}{double Student-level covariate (think baseline measured SES).}
#'  \item{\code{S.id}}{integer School ID.}
#'  \item{\code{D.id}}{integer District ID.}
#'  \item{\code{Yobs}}{double Observed outcome (think test score).)}
#'  \item{\code{T.x}}{integer Treatment assignment (1 treated, 0 control).} }
#'@details These data were generated via the PUMP package.
"fakeBrokeCRT"




#'@title Fake Cluster Randomized Trial data (set 2)
#'
#'@description This dataset is a three-level dataset with students in
#'  schools in districts, with schools randomized to treatment and
#'  control.  It is used to illustrate and test the clusterRCT pacakge.
#'
#' @seealso \code{\link{fakeCRT}}
#'
#'@format A data frame with 1500 rows and 7 variables: \describe{
#'  \item{\code{S.id}}{integer School ID.}
#'  \item{\code{D.id}}{integer District ID.}
#'  \item{\code{Y1}}{double Potential outcome under treatment.}
#'  \item{\code{Y0}}{double Potential outcome under no treatment.}
#'  \item{\code{Yobs}}{double Observed outcome (think test score).)}
#'  \item{\code{T.x}}{integer Treatment assignment (1 treated, 0 control).} }
#'@details These data were generated via the PUMP package.
"fakeCRT2"

