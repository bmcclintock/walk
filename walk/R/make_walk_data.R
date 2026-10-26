#' @title Produce design data for use in fitting MMPP movement models
#' @param proc_data A sparse matrix with rows corresponding to time steps and columns corresponding to cells. 
#' The entries are probabilities that the animal is located in the corresponding cell and time.
#' @param cell_data A `SpatRaster` object from the `terra` package.
#' @param grad A character vector of names of `cell_data` layers for which gradient covariates will be constructed.
#' @param get_subset_idx Logical. Should a rate matrix subsetting approximation be used to increase computational speed.  
#' @param rast_mask Raster mask for inaccessible cells when `cell_data` is of type `SpatRaster` from the `terra` package. This is ignored
#' if \code{cell_data} is an `POLYGON` data frame from the `sf` package.
#' @param directions Neighborhood structure, one of `"rook"` or `"queen"`.
#' @param mask_zero Logical. If `TRUE` (default) cells where `rast_mask` is `0` are dropped as
#' well as cells where it is `NA`. A mask built with
#' `terra::rasterize(land, r, field = 0, background = 1)` contains no `NA`s, so with
#' `mask_zero = FALSE` it would drop nothing.
#' @param crop_to_data Optional buffer distance (map units). If supplied, `cell_data` (and
#' `rast_mask`) are cropped to the extent of the cells with non-zero observation likelihood plus
#' this buffer before the rate-matrix data are built.
#' @param debug Debugging level: 1-3 mainly for package developers.
#' @param ... Ignored arguments.
#' @import dplyr 
#' @importFrom Matrix rowSums
#' @importFrom igraph graph_from_data_frame
#' @export
make_walk_data <- function(proc_data, cell_data, grad=NULL, rast_mask=NULL, 
                           get_subset_idx=TRUE, directions="rook", mask_zero=TRUE,
                           crop_to_data=NULL, debug=0,...){
  
  cells <- obs <- timestamp <- quad <- NULL
  
  if(debug==1) browser()
  cell_name <- NULL
  quad_data <- NULL
  
  if(debug==2) browser()
  
  if(inherits(cell_data, "SpatRaster")){
    if(ncol(proc_data$L) != terra::ncell(cell_data)) stop("'proc_data' was not built on 'cell_data' (different number of cells).")
    if(!is.null(rast_mask) && mask_zero){
      rast_mask <- terra::ifel(rast_mask == 0, NA, rast_mask)
    }
    if(!is.null(crop_to_data)){
      used <- which(Matrix::colSums(proc_data$L) > 0)
      xy_used <- terra::xyFromCell(cell_data, used)
      e <- terra::ext(min(xy_used[,1]) - crop_to_data, max(xy_used[,1]) + crop_to_data,
                      min(xy_used[,2]) - crop_to_data, max(xy_used[,2]) + crop_to_data)
      e <- terra::intersect(e, terra::ext(cell_data))
      new_cells <- terra::cells(cell_data, e)
      old_ce <- terra::xyFromCell(cell_data, new_cells)
      cell_data <- terra::crop(cell_data, e, snap="out")
      if(!is.null(rast_mask)) rast_mask <- terra::crop(rast_mask, e, snap="out")
      # remap columns of L to the cropped raster
      newid <- terra::cellFromXY(cell_data, old_ce)
      keep <- !is.na(newid)
      M <- Matrix::sparseMatrix(i=new_cells[keep], j=newid[keep], x=1,
                                dims=c(ncol(proc_data$L), terra::ncell(cell_data)))
      lost <- sum(proc_data$L[, -new_cells[keep], drop=FALSE])
      if(lost > 0) warning("Cropping removed some observation likelihood mass (", signif(lost,3), ").")
      proc_data$L <- as(proc_data$L %*% M, "dgCMatrix")
    }
    q_list <- make_q_data_rast(cell_data, grad=grad, rast_mask=rast_mask, directions=directions)
    attr(q_list$q_m, "directions") <- directions
    attr(q_list$q_m, "cov_class") <- "SpatRaster"
  } else if(inherits(cell_data, "sf")){
    stop("In the current version cell_data must be a 'SpatRaster' from the terra package. In future versions we hope to incorporate sf or stars data sets.")
    # q_list <- make_q_data_sf(cell_data, cell_name)
    # attr(q_list$q_m, "directions") <- NULL
    # attr(q_list$q_m, "cov_class") <- "sf"
  }
  out <- list()
  out$L <- proc_data$L[,q_list$q_r$cell]
  if(min(Matrix::rowSums(out$L))==0) warning("Some observations have likelihood values of 0 in all ras cells after 'rast_mask' cells are removed!")
  out$times <- proc_data$times
  if(is.null(out$times$seg)) out$times$seg <- 1L
  if(is.null(out$times$id)) out$times$id <- "1"
  if(is.null(attr(out$times, "time_unit"))) warning("'proc_data$times' has no 'time_unit' attribute; 'dt' units are unknown.")
  out$cell_data <- terra::wrap(cell_data[[1]])
  out$q_r <- q_list$q_r
  out$q_m <- q_list$q_m
  class(out) <- c(class(out), "walk_ddl")
  # Check for irreducubility
  g <- graph_from_data_frame(d = out$q_m[,c("from_cellx","cellx")], directed = TRUE)
  scc <- components(g, mode = "strong")
  if(scc$no>1){
    warning("There appears to be ", scc$no, " disconnected groups of cells in the masked habitat.",
            " This can cause fit failure becuase the implied CTMC is not irreducible. ",
            "In addition, the limiting distribution does not exist.")
  }
  if(debug==3) browser()
  return(out)
}




