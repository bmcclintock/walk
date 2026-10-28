#' @title Define specifics for model structure and fitting control
#' @name arg_funs
#' @aliases ctmc_model
#' @aliases ctmc_control
#' @param form Model formula for CTMC model components, residency or movement
#' @param link A link funtion for the residency or movement model. One of `"soft_plus"` or `"log"`,
#' or `"logit"` (residency model only; requires `U > L >= 0`).
#' @param a Scale parameter for the `"soft_plus"` link function. Ignored for `link="log"`. The `a` parameter determines the approximation to 
#' a hard plus function, i.e., as `a` becomes large the soft plus function converges to `g^{-1}(x) = max(0,x)`. For this specification, `a` must be greater than or equal to 1.
#' @param L Lower bound for general logit link
#' @param U Upper bound for general logit link
#' @param q_r A named list for the residency model. Must contain elements: `form`, `link`, and `a`. 
#' The easiest way to construct this is the `ctmc_model` function.
#' @param q_m A named list for the movement portion of the model. See `q_r`.
#' @param p Logical. Should a zero-inflation parameter be used. Necessary when there might be outliers when the location likelihood does not
#' overlap likely locations. 
#' @param delta Vector of probabilities for the initial location. Can be one of `"uniform"` for a uniform distribution over all cells, `"stationary"` to use the 
#' stationary distribution implied by the parameters, or a numeric vector the same length as the number of cells. 
#' @param form The form of the rate matrix entries. Can be one of: `"mult"` for residency times movement, i.e., q_\{ij\} = q_r(i) * pi_m(i,j),
#' and `"sde"` for a Langevin diffusion approximation (similar to an additive model q_\{ij\} = q_r(i) + pi_m(i,j),
#' which is not implemented). 
#' @param norm Should the movement portion be adjusted to sum to 1. This is the parameterization suggested by Hewitt et al. (2023). 
#' @param clip Positive value for the maximum movement rate for all cells. All cells with modeled total movement rates above this value will be truncated to `clip`.
#' @references Hewitt, J., Gelfand, A. E., & Schick, R. S. (2023). Time-discretization approximation enriches continuous-time discrete-space models for animal movement. The Annals of Applied Statistics, 17:740-760.
#' @export
ctmc_model <- function(form=~1, link="soft_plus", a=1, L=0, U=0){
  if(!inherits(form, "formula")) stop("'form' must be a one-sided formula, e.g. ~ 1 or ~ d_c1.")
  if(!is.character(link) || length(link) != 1 || !(link %in% c("soft_plus","log","logit")))
    stop("'link' must be one of \"soft_plus\", \"log\" or \"logit\".")
  if(!is.numeric(a) || length(a) != 1 || !is.finite(a) || a < 1) stop("'a' must be a single number >= 1.")
  if(link=="logit"){
    if(U <= 0 || U <= L) stop("U>0 and U>L must be specified when using the logit link.")
    if(L < 0) stop("L for the logit link must be >= 0.")
  }
  out <- list(
    form = form, link = link, a = a, L=L, U=U
  )
  class(out) <- c("ctmc_mod","list")
  return(out)
}

#' @name arg_funs
#' @importFrom stats terms
#' @export
ctmc_control <- function(q_r = ctmc_model(), q_m=ctmc_model(), 
                         p=FALSE, delta="uniform", form="mult", norm=TRUE, clip=0){
  if(!inherits(q_r, "ctmc_mod") || !inherits(q_m, "ctmc_mod")) stop("'q_r' and 'q_m' must be created with ctmc_model().")
  if(!is.character(form) || length(form) != 1 || !(form %in% c("mult","sde")))
    stop("'form' must be \"mult\" or \"sde\" (\"add\" is not implemented).")
  if(!(is.numeric(delta) || (is.character(delta) && length(delta) == 1 && delta %in% c("uniform","stationary"))))
    stop("'delta' must be \"uniform\", \"stationary\", or a numeric vector with one value per cell.")
  if(!is.logical(p) || length(p) != 1 || is.na(p)) stop("'p' must be TRUE or FALSE.")
  if(q_m$link=="logit") stop("'logit' link not implimented for movement portion of the CTMC model. Please select either 'log' or 'soft_plus'")
  if(clip<0) stop("'clip' must be > 0!")
  if(length(attr(terms(q_r$form),"term.labels"))>0 & form=="sde") stop("The SDE cannot handle q_r forms other than ~1 at this time.")
  if(form=="sde") norm <- FALSE
  out <- list(
    q_r = q_r, q_m=q_m, 
    p=p, delta=delta, form=form, norm=norm, clip=clip
  )
  class(out) <- c("ctmc_cont","list")
  return(out)
}
