

#' Compute the p-values for the composite null hypothesis.
#'
#'@param px a numeric vectors of p-values for the first individual hypothesis (e.g., the association between the exposure and mediator).
#'@param py a numeric vectors of p-values for the second individual hypothesis (e.g., the association between the mediator and outcome).
#'
#'@return The p-values for the composite null.
#'
#'@author Yaowu Liu
#'
#'@details Suppose that there are N items (e.g., mediators). Each item has two independent p-values.
#'Our method only calculates the composite null p-values for items whose maximum of the two individual p-values is less than 0.5. Otherwise, NA will be returned.
#'When the composite null p-values are large (e.g., > 0.1), they might not be accurate. Do not take these large p-values seriously (e.g., to draw a QQ plot).
#'
#'@references Liu, Y. (2025). A simple and powerful method for large-scale composite null hypothesis testing with applications in mediation analysis.  \emph{Biometrics}, 81(1), ujaf011.
#'
#'@examples N = 100000  # number of items
#'@examples set.seed(123)
#'@examples
#'@examples ### Generate the individual p-values from z scores. The first 100 items are under the composite alternative
#'@examples x = rnorm(N); y = rnorm(N)
#'@examples x[1:100] = x[1:100] + runif(100,min=4,max = 4)
#'@examples y[1:100] = y[1:100] + runif(100,min=4,max = 4)
#'@examples x[101:200]  = x[101:200] + runif(100,min=4,max = 4)
#'@examples y[201:300]  = y[201:300] + runif(100,min=4,max = 4)
#'@examples px = 2*pnorm(-abs(x)); py = 2*pnorm(-abs(y))
#'@examples
#'@examples ### Compute the p-values for the composite null.
#'@examples pvals_compo = DEI_B(px,py)
#'@examples
#'@examples ### Apply the BH procedure to control FDR. Note there are many NA's in pvals_compo
#'@examples pvals_compo[is.na(pvals_compo)] = 1
#'@examples which(p.adjust(pvals_compo, method = "BH") < 0.05)
#'@export

DEI_B <- function(px,py){
    N = length(px)
    if (length(py)!=N){
        stop("The length of the two sets of p-values must be the same!")
    }

    is.small.x = px < py
    pmin = px; pmin[!is.small.x] = py[!is.small.x]
    pmax = px; pmax[is.small.x] = py[is.small.x]

    pmin.sort = sort(pmin)

    pmax.sort = sort(pmax,index.return = T)
    pmax.sort.index = pmax.sort[["ix"]]
    pmax.sort = pmax.sort[["x"]]

    Pvals.est = rep(NA,N)

    loc.min = 1
    i = 1
    epsilon = pmax.sort[i]
    while (i <= N && epsilon<0.4) {
        if (pmin.sort[loc.min] <= epsilon){
            loc.min = loc.min + 1
        }else{
            Nb = loc.min - i   ### region C is < epsilon instead of < = epsilon
            pval.est = epsilon/(1-2*epsilon)*((Nb/(N))-epsilon)
            if (i>1 && pval.est<Pvals.est[i-1]){
                Pvals.est[i] = Pvals.est[i-1]
            }else{
                Pvals.est[i] = pval.est
            }

            i = i + 1
            epsilon = pmax.sort[i]
        }
    }

    pvals.out = rep(NA,N)
    pvals.out[pmax.sort.index] = Pvals.est

    return(pvals.out)
}
