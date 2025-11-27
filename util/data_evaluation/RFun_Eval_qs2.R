## Simple library of functions to assess and visualize the results of the coupled simulations
## Modified to work with .qs2 files (qs2 package format)

## Time-stamp: "Last modified 2025-11-25"

require(qs2)  ## for reading .qs2 files
require(stringr)

# Note: RedModRphree, Rmufits, and Rcpp functions for DHT/PHT reading are kept
# but you'll need those packages only if you use ReadAllDHT/ReadAllPHT functions

curdir <- dirname(sys.frame(1)$ofile) ##path.expand(".")

print(paste("RFun_Eval_qs2.R is in ", curdir))

## ============================================================================
## NEW: Functions for reading .qs2 simulation outputs
## ============================================================================

## function which reads all simulation results in a given directory (.qs2 format)
ReadRTSims_qs2 <- function(dir) {
    pattern <- "^iter_.*\\.qs2$"
    files_full <- list.files(dir, pattern = pattern, full.names = TRUE)
    files_name <- list.files(dir, pattern = pattern, full.names = FALSE)
    
    if (length(files_full) == 0) {
        warning(paste("No .qs2 files found in", dir, "with pattern", pattern))
        return(NULL)
    }
    
    res <- lapply(files_full, qs2::qs_read)
    names(res) <- gsub(".qs2", "", files_name, perl = TRUE)
    
    return(res[str_sort(names(res), numeric = TRUE)])
}
## Read a single .qs2 file
ReadQS2 <- function(file) {
    if (!file.exists(file)) {
        stop(paste("File not found:", file))
    }
    qs2::qs_read(file)
}

## Extract chemistry field data from .qs2 iteration file
## Assumes structure similar to old .rds format with $C (chemistry) and $T (transport)
ExtractChemistry <- function(qs2_data) {
    if ("C" %in% names(qs2_data)) {
        return(qs2_data$C)
    } else if (is.data.frame(qs2_data)) {
        return(qs2_data)
    } else {
        warning("Could not find chemistry data in expected format")
        return(qs2_data)
    }
}

## Extract transport field data from .qs2 iteration file
ExtractTransport <- function(qs2_data) {
    if ("T" %in% names(qs2_data)) {
        return(qs2_data$T)
    } else {
        warning("Could not find transport data in expected format")
        return(NULL)
    }
}

## ============================================================================
## ORIGINAL: DHT/PHT reading functions (kept for surrogate analysis)
## ============================================================================

# Only load these if needed (requires Rcpp compilation)
if (requireNamespace("Rcpp", quietly = TRUE) && file.exists(paste0(curdir, "/interpret_keys.cpp"))) {
    library(Rcpp)
    sourceCpp(file = paste0(curdir, "/interpret_keys.cpp"))
    
    # Wrapper around previous sourced Rcpp function
    ConvertDHTKey <- function(value) {
      rcpp_key_convert(value)
    }
    
    ConvertToUInt64 <- function(double_data) {
      rcpp_uint64_convert(double_data)
    }
} else {
    if (!requireNamespace("Rcpp", quietly = TRUE)) {
        message("Note: Rcpp not available. DHT/PHT reading functions will not work.")
    }
    # Create dummy functions so the rest of the script doesn't break
    ConvertDHTKey <- function(value) {
      stop("Rcpp not available. Cannot convert DHT keys.")
    }
    ConvertToUInt64 <- function(double_data) {
      stop("Rcpp not available. Cannot convert to UInt64.")
    }
}

## function which reads all successive DHT stored in a given directory
ReadAllDHT <- function(dir, new_scheme = TRUE) {
    files_full <- list.files(dir, pattern="iter.*\\.dht$", full.names=TRUE)
    files_name <- list.files(dir, pattern="iter.*\\.dht$", full.names=FALSE)
    
    if (length(files_full) == 0) {
        warning(paste("No .dht files found in", dir))
        return(NULL)
    }
    
    res <- lapply(files_full, ReadDHT, new_scheme = new_scheme)
    names(res) <- gsub("\\.dht$","",files_name)
    return(res)
}

## function which reads one .dht file and gives a matrix
ReadDHT <- function(file, new_scheme = TRUE) {
    conn <- file(file, "rb")  ## open for reading in binary mode
    if (!isSeekable(conn))
        stop("Connection not seekable")

    ## we first reposition ourselves to the end of the file...
    tmp <- seek(conn, where=0, origin = "end")
    ## ... and then back to the origin so to store the length in bytes
    flen <- seek(conn, where=0, origin = "start")

    ## we read the first 2 integers (4 bytes each) containing dimensions in bytes
    dims <- readBin(conn, what="integer", n=2)

    ## compute dimensions of the data
    tots <- sum(dims)
    ncol <- tots/8
    nrow <- (flen - 8)/tots ## 8 here is 2*sizeof("int")
    buff <- readBin(conn, what="double", n=ncol*nrow)
    ## close connection
    close(conn)
    res <- matrix(buff, nrow=nrow, ncol=ncol, byrow=TRUE)

    if (new_scheme) {
      nkeys <- dims[1] / 8
      keys <- res[, 1:nkeys]

      conv <- apply(keys, 2, ConvertDHTKey)
      res[, 1:nkeys] <- conv
    }

    return(res)
}

## function which reads all successive PHT stored in a given directory
ReadAllPHT <- function(dir, with_info = FALSE) {
    files_full <- list.files(dir, pattern="iter.*\\.pht$", full.names=TRUE)
    files_name <- list.files(dir, pattern="iter.*\\.pht$", full.names=FALSE)
    
    if (length(files_full) == 0) {
        warning(paste("No .pht files found in", dir))
        return(NULL)
    }
    
    res <- lapply(files_full, ReadPHT, with_info = with_info)
    names(res) <- gsub("\\.pht$","",files_name)
    return(res)
}

## function which reads one .pht file and gives a matrix
ReadPHT <- function(file, with_info = FALSE) {
    conn <- file(file, "rb")  ## open for reading in binary mode
    if (!isSeekable(conn))
        stop("Connection not seekable")

    ## we first reposition ourselves to the end of the file...
    tmp <- seek(conn, where=0, origin = "end")
    ## ... and then back to the origin so to store the length in bytes
    flen <- seek(conn, where=0, origin = "start")

    ## we read the first 2 integers (4 bytes each) containing dimensions in bytes
    dims <- readBin(conn, what="integer", n=2)

    ## compute dimensions of the data
    tots <- sum(dims)
    ncol <- tots/8
    nrow <- (flen - 8)/tots ## 8 here is 2*sizeof("int")
    buff <- readBin(conn, what="double", n=ncol*nrow)
    ## close connection
    close(conn)
    res <- matrix(buff, nrow=nrow, ncol=ncol, byrow=TRUE)

    nkeys <- dims[1] / 8
    keys <- res[, 1:nkeys]

    timesteps <- res[, nkeys + 1]
    conv <- apply(keys, 2, ConvertDHTKey)

    ndata <- dims[2] / 8
    fill_rate <- ConvertToUInt64(res[, nkeys + 2])

    buff <- c(conv, timesteps, fill_rate)

    if (with_info) {
      ndata <- dims[2]/8
      visit_count <- ConvertToUInt64(res[, nkeys + ndata])
      buff <- c(buff, visit_count)
    }

    res <- matrix(buff, nrow = nrow, byrow = FALSE)

    return(res)
}

## ============================================================================
## PLOTTING and ANALYSIS functions (work with both .rds and .qs2 data)
## ============================================================================

## Scatter plots of each variable in the iteration 
PlotScatter <- function(sam1, sam2, which=NULL, labs=c("NO DHT", "DHT"), pch=".", cols=3, ...) {
    if ((!is.data.frame(sam1)) & ("T" %in% names(sam1)))
        sam1 <- sam1$C
    if ((!is.data.frame(sam2)) & ("T" %in% names(sam2)))
        sam2 <- sam2$C
    if (is.numeric(which))
        inds <- which
    else if (is.character(which))
        inds <- match(which, colnames(sam1))
    else if (is.null(which))
        inds <- seq_along(colnames(sam1))

    rows <- ceiling(length(inds) / cols)
    par(mfrow=c(rows, cols))
    a <- lapply(inds, function(x) {
        plot(sam1[,x], sam2[,x], main=colnames(sam1)[x],  xlab=labs[1], ylab=labs[2], pch=pch, col="red", ...)
        abline(0,1, col="grey", lwd=1.5)
    })
    invisible()
}

##### Some metrics for relative comparison

## Root Mean Square Error
RMSE <- function(pred, obs)
    sqrt(mean((pred - obs)^2, na.rm = TRUE))

## Using range as norm
RranRMSE <- function(pred, obs)
    sqrt(mean((pred - obs)^2, na.rm = TRUE))/abs(max(pred, na.rm = TRUE) - min(pred, na.rm = TRUE))

## Using max val as norm
RmaxRMSE <- function(pred, obs)
    sqrt(mean((pred - obs)^2, na.rm = TRUE))/abs(max(pred, na.rm = TRUE))

## Using sd as norm
RsdRMSE <- function(pred, obs)
    sqrt(mean((pred - obs)^2, na.rm = TRUE))/sd(pred, na.rm = TRUE)

## Using mean as norm
RmeanRMSE <- function(pred, obs)
    sqrt(mean((pred - obs)^2, na.rm = TRUE))/mean(pred, na.rm = TRUE)

## Using mean as norm
RAEmax <- function(pred, obs)
    mean(abs(pred - obs), na.rm = TRUE)/max(pred, na.rm = TRUE)

## Max absolute error
MAE <- function(pred, obs)
    max(abs(pred - obs), na.rm = TRUE)

## Mean Absolute Percentage Error
MAPE <- function(pred, obs)
    mean(abs((obs - pred) / obs) * 100, na.rm = TRUE)

## workhorse function for ComputeErrors and its use with mapply
AppliedFun <- function(a, b, .fun) {
    # Extract chemistry data if needed
    if (!is.data.frame(a) && "C" %in% names(a)) a <- a$C
    if (!is.data.frame(b) && "C" %in% names(b)) b <- b$C
    
    mapply(.fun, as.list(a), as.list(b))
}

## Compute the diffs between two simulation, iter by iter,
## with a given metric (passed in form of function name to this function)
ComputeErrors <- function(sim1, sim2, FUN=RMSE) {
    if (length(sim1)!= length(sim2)) {
        cat("The simulations do not have the same length, subsetting to the shortest\n")
        a <- min(length(sim1), length(sim2))
        sim1 <- sim1[1:a]
        sim2 <- sim2[1:a]
    }
    if (!is.function(match.fun(FUN))) {
        stop("Invalid function\n")
    }

    t(mapply(AppliedFun, sim1, sim2, MoreArgs=list(.fun=FUN)))
}

## Function to display the error progress between 2 simulations
ErrorProgress <- function(mat, ignore, colors, metric, ...) {
    if (is.null(mat)) {
        stop("Cannot plot: matrix is NULL")
    }
    
    # Convert to matrix if it's a vector or data frame
    if (is.vector(mat)) {
        stop("Cannot plot: input is a vector (need at least 2 columns). Check that your data has multiple columns.")
    }
    if (is.data.frame(mat)) {
        mat <- as.matrix(mat)
    }
    
    if (nrow(mat) == 0 || ncol(mat) == 0) {
        stop("Cannot plot: matrix is empty")
    }
    
    if (missing(colors))
        colors <- sample(rainbow(ncol(mat)))

    if (missing(metric))
        metric <- "Metric"

    ## if the optional argument "ignore" (a character vector) is
    ## passed, we remove the matching column names
    if (!missing(ignore)) {
        to_remove <- match(ignore, colnames(mat))
        to_remove <- to_remove[!is.na(to_remove)]  # Remove NAs
        if (length(to_remove) > 0) {
            mat <- mat[, -to_remove, drop = FALSE]
            colors <- colors[-to_remove]
        }
    }
    
    yc <- mat[nrow(mat),]
    par(mar=c(5,4,2,8))
    matplot(mat, type="l", lty=1, lwd=2, col=colors, xlab="iteration", ylab=metric, ...)
    mtext(colnames(mat), side = 4, line = 0.5, outer = FALSE, at = yc, adj = 0, col = colors, las=2, cex=0.7)
}

## Function which exports all simulations to ParaView's .vtu 
## Requires package RcppVTK
ExportToParaview <- function(vtu, nameout, results) {
    if (!requireNamespace("RcppVTK", quietly = TRUE)) {
        stop("Package RcppVTK is required for this function")
    }
    require(RcppVTK)
    
    n <- length(results)
    vars <- colnames(results[[1]])
    ## strip eventually present ".vtu" from nameout
    nameout <- sub(".vtu", "", nameout, fixed=TRUE)
    namesteps <- paste0(nameout,  ".", sprintf("%04d",seq(1,n)), ".vtu")
    for (step in seq_along(results)) {
        file.copy(from=vtu, to=namesteps[step], overwrite = TRUE)
        cat(paste("Saving step ", step, " in file ", namesteps[step], "\n"))
        ret <- ExportMatrixToVTU (fin=vtu, fout=namesteps[step], names=colnames(results[[step]]), mat=results[[step]])
    }
    invisible(ret)
}


## Version of Rmufits::PlotCartCellData with the ability to fix the
## "breaks" for color coding of 2D simulations
Plot2DCellData <- function (data, grid, nx, ny, contour = TRUE,
                            nlevels = 12, breaks, palette = "heat.colors",
                            rev.palette = TRUE, scale = TRUE, plot.axes=TRUE, ...) {
    if (!missing(grid)) {
        xc <- unique(sort(grid$cell$XCOORD))
        yc <- unique(sort(grid$cell$YCOORD))
        nx <- length(xc)
        ny <- length(yc)
        if (!length(data) == nx * ny) 
            stop("Wrong nx, ny or grid")
    } else {
        xc <- seq(1, nx)
        yc <- seq(1, ny)
    }
    z <- matrix(round(data, 6), ncol = nx, nrow = ny, byrow = TRUE)
    pp <- t(z[rev(seq(1, nrow(z))), ])

    if (missing(breaks)) {
        breaks <- pretty(data, n = nlevels)
    }
    
    breakslen <- length(breaks)
    colors <- do.call(palette, list(n = breakslen - 1))
    if (rev.palette) 
        colors <- rev(colors)
    if (scale) {
        par(mfrow = c(1, 2))
        nf <- layout(matrix(c(1, 2), 1, 2, byrow = TRUE), widths = c(4, 
            1))
    }
    par(las = 1, mar = c(5, 5, 3, 1))
    image(xc, yc, pp, xlab = "X [m]", ylab = "Y[m]", las = 1, asp = 1,
          breaks = breaks, col = colors, axes = FALSE, ann=plot.axes,
          ...)

    if (plot.axes) {
        axis(1)
        axis(2)
    }
    if (contour) 
        contour(unique(sort(xc)), unique(sort(yc)), pp, breaks = breaks, 
            add = TRUE)
    if (scale) {
        par(las = 1, mar = c(5, 1, 5, 5))
        if (requireNamespace("Rmufits", quietly = TRUE)) {
            Rmufits::PlotImageScale(data, breaks = breaks, add.axis = FALSE, 
                axis.pos = 4, col = colors)
        }
        axis(4, at = breaks)
    }
    invisible(pp)
}

PlotAsMP4 <- function(data, nx, ny, to_plot, out_dir, name,
                      contour = FALSE, scale = FALSE, framerate = 30) {
  sort_data <- data[str_sort(names(data), numeric = TRUE)]
  plot_data <- lapply(sort_data, function(x) {
      if (!is.data.frame(x) && "C" %in% names(x)) {
          return(x$C[[to_plot]])
      } else {
          return(x[[to_plot]])
      }
  })
  pad_size <- ceiling(log10(length(plot_data)))

  dir.create(out_dir, showWarnings = FALSE)
  output_files <- paste0(out_dir, "/", name, "_%0", pad_size, "d.png")
  output_mp4 <- paste0(out_dir, "/", name, ".mp4")

  png(output_files,
    width = 297, height = 210, units = "mm",
    res = 100
  )

  for (i in 1:length(plot_data)) {
    if (requireNamespace("Rmufits", quietly = TRUE)) {
        Rmufits::PlotCartCellData(plot_data[[i]], nx = nx, ny = ny, contour = contour, scale = scale)
    } else {
        Plot2DCellData(plot_data[[i]], nx = nx, ny = ny, contour = contour, scale = scale)
    }
  }
  dev.off()

  ffmpeg_command <- paste(
    "ffmpeg -y -framerate", framerate, "-i", output_files,
    "-c:v libx264 -crf 22", output_mp4
  )

  system(ffmpeg_command)
  message(paste("Created video:", output_mp4))
}

cat("\n=== RFun_Eval_qs2.R loaded successfully ===\n")
cat("New functions for .qs2 files:\n")
cat("  - ReadRTSims_qs2(dir)  : Read all iteration .qs2 files\n")
cat("  - ReadQS2(file)        : Read single .qs2 file\n")
cat("All other functions work as before!\n\n")
