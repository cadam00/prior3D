get_biodiv_raster <- function()
  terra::rast(system.file("external/biodiv_raster.tif", package="prior3D"))
get_depth_raster <- function()
  terra::rast(system.file("external/depth_raster.tif", package="prior3D"))

## Function to read all .asc and .tif files from folder
get_rast <- function(path){
  terra::rast(list.files(path = path,
                         pattern='.asc$|.tif$',
                         all.files=TRUE, full.names=TRUE))
}

## Function to split 2D feature distributions into 3D ones
split_rast <- function(biodiv_raster,
                       depth_raster,
                       breaks,
                       biodiv_df,
                       val_depth_range=TRUE,
                       sep_biodiv_df=","){

  ## Ensure that biodiv_raster is rast file with depths or path of its file
  if (is.character(biodiv_raster)){
    biodiv_raster <- get_rast(biodiv_raster) * 1
  } else if (!is(biodiv_raster, "SpatRaster")) {
    stop("biodiv_raster can be only either character path or SpatRaster class")
  }

  ## Ensure that depth_raster is rast file with depths or path of its file
  if (is.character(depth_raster)){
    depth_raster <- terra::rast(depth_raster)  * 1
  } else if (!is(depth_raster, "SpatRaster")) {
    stop("depth_raster can be only either character or SpatRaster class")
  }
  breaks <- sort(breaks, decreasing = TRUE)
  out_of_analysis <- depth_raster > breaks[1] |
                     depth_raster < breaks[length(breaks)]
  depth_raster[out_of_analysis] <- NA

  if (is.character(biodiv_df)){
    if (file.exists(biodiv_df))
      biodiv_df <- switch(tools::file_ext(biodiv_df),
                          csv = read.csv(biodiv_df, header = TRUE,
                                         sep=sep_biodiv_df),
                          xls = readxl::read_xls(biodiv_df),
                          xlsx = readxl::read_xlsx(biodiv_df)
      )
  } else if (!is(biodiv_df, "data.frame")) {
    stop("biodiv_df can be only either character path or data.frame class")
  }
  possible_names <- c("species_name", "pelagic", "min_z", "max_z", "group")
  possible_names <- possible_names[possible_names %in% names(biodiv_df)]
  if ( !any(c("species_name", "pelagic") %in%  possible_names) ){
    stop("Please provide species_name and pelagic in biodiv_df")
  } else if ( !("species_name" %in%  possible_names) ){
    stop("Please provide species_name in biodiv_df")
  } else if ( !("pelagic" %in%  possible_names) ){
    stop("Please provide pelagic in biodiv_df")
  }
  biodiv_df <- biodiv_df[,possible_names]

  names(biodiv_raster) <- tolower(names(biodiv_raster))
  biodiv_df$species_name <- tolower(biodiv_df$species_name)
  ## Remove species for which we have the 2D distribution
  ## but not their depth behavior
  biodiv_raster <-biodiv_raster[[names(biodiv_raster) %in%
                                   biodiv_df$species_name]]

  ## Crop depth_raster and biodiv_raster at the same extent
  ## We keep only the smaller ones

  if ( any(terra::res(depth_raster) != terra::res(biodiv_raster)) ||
       (terra::crs(depth_raster) != terra::crs(biodiv_raster))){
    warning(
      paste0("Different resolution or coordinates among biodiv_raster and",
                   " depth_raster:\nproject biodiv_raster on depth_raster")
      )
    depth_raster <- terra::project(depth_raster,biodiv_raster, method="average")
  } else {
    ext_depth_raster <- terra::ext(depth_raster)
    ext_biodiv_raster <- terra::ext(biodiv_raster)
    if (any(ext_depth_raster != ext_biodiv_raster)){
      warning(
        "Different extent among biodiv_raster and depth_raster: crop them"
        )
      min_x <- max( c(ext_depth_raster[1], ext_biodiv_raster[1]) )
      max_x <- min( c(ext_depth_raster[2], ext_biodiv_raster[2]) )
      min_y <- max( c(ext_depth_raster[3], ext_biodiv_raster[3]) )
      max_y <- min( c(ext_depth_raster[4], ext_biodiv_raster[4]) )
      exts <- terra::ext(min_x, max_x, min_y, max_y)
      if (exts != ext_depth_raster)
        depth_raster <- crop(depth_raster, exts)
      if (exts != ext_biodiv_raster)
        biodiv_raster <- crop(biodiv_raster, exts)
    }
  }

  ## Create a data.frame with only common names
  ## between 2D distribution data and biodiv_df
  biodiv_df <- merge(data.frame(species_name=names(biodiv_raster)), biodiv_df)
  biodiv_df <- biodiv_df[order(match(biodiv_df$species_name,
                                     names(biodiv_raster))),]
  ## Get vector of number of layers
  n_depths <- length(breaks)-1
  one_n_depths <- 1:n_depths

  ## Classify minimum and maximum depth in depth_raster by user-specified breaks
  depth_raster <- n_depths -
                  terra::classify(depth_raster, breaks, include.lowest=TRUE)

  ## Cast pelagic to logical
  biodiv_df$pelagic <- as.logical(biodiv_df$pelagic)

  ## Get which layers refer to pelagic and which to benthic species
  is_pelagic <- which(biodiv_df$pelagic)
  is_benthic <- which(!biodiv_df$pelagic)

  ## Create a list of rasters where each one correspondes to each depth level
  ## (in the same order)
  ## 1) mask each species class by corresponding depth raster
  ## 2) combine them in the order of initial 2D distribution
  current_rasters <- lapply(one_n_depths, ## For every depth
                            FUN=function(depth_i, is_pelagic, is_benthic,
                                         biodiv_raster, depth_raster,
                                         one_n_depths){
                              if (length(is_pelagic) == 0){
                                rast_i <- terra::mask(
                                  biodiv_raster[[is_benthic]],
                                  depth_raster,
                                  maskvalues=c(NA,one_n_depths[-depth_i]))[[
                                    order(is_benthic)]]
                              } else if (length(is_benthic) == 0){
                                rast_i <- terra::mask(
                                  biodiv_raster[[is_pelagic]],
                                  depth_raster,
                                  maskvalues=c(NA,one_n_depths[
                                    depth_i>one_n_depths]))[[order(is_pelagic)]]
                              } else {
                                benthics <- terra::mask(
                                  biodiv_raster[[is_benthic]],
                                  depth_raster,
                                  maskvalues=c(NA,one_n_depths[-depth_i]))
                                pelagics <- terra::mask(
                                  biodiv_raster[[is_pelagic]],
                                  depth_raster,
                                  maskvalues=c(NA,one_n_depths[
                                    depth_i>one_n_depths]))
                                rast_i <- c(pelagics,benthics)[[
                                  order(c(is_pelagic,is_benthic))]]
                              }
                              return(rast_i)
                            },
                            is_pelagic        = is_pelagic,
                            is_benthic        = is_benthic,
                            biodiv_raster     = biodiv_raster,
                            depth_raster      = depth_raster,
                            one_n_depths      = one_n_depths
  )

  if (!val_depth_range && (("min_z" %in% possible_names) || ("max_z" %in%
                                                            possible_names))){
    if ( !("min_z" %in% possible_names)){
      biodiv_df$min_z <- -Inf
    }
    if ( !("max_z" %in% possible_names)){
      biodiv_df$max_z <- Inf
    }
    ## Classify minimum and maximum depth in biodiv_df by user-specified breaks
    biodiv_df$min_z[is.na(biodiv_df$min_z)] <- -Inf
    biodiv_df$max_z[is.na(biodiv_df$max_z)] <- Inf
    if (any(biodiv_df$min_z > biodiv_df$max_z))
      stop("It should be min_z <= max_z.")
    biodiv_df$min_z <- (n_depths+1) - cut(biodiv_df$min_z, breaks = breaks,
                                          labels = FALSE, include.lowest = TRUE,
                                          ordered_result = TRUE)
    biodiv_df$max_z <- (n_depths+1) - cut(biodiv_df$max_z, breaks = breaks,
                                          labels = FALSE, include.lowest = TRUE,
                                          ordered_result = TRUE)
    nan_rast <- rast(
      nrow=nrow(current_rasters[[1]]),
      ncol=ncol(current_rasters[[1]]),
      resolution=res(current_rasters[[1]]),
      extent=ext(current_rasters[[1]]),
      crs = crs(current_rasters[[1]])
    )

    n_sp <- nrow(biodiv_df)
    one_n_sp <- 1:n_sp

    for (depth_i in one_n_depths){

      remove_these <- which((biodiv_df$min_z<depth_i)|(biodiv_df$max_z>depth_i))

      n_remove <- length(remove_these)

      if (n_remove > 0){

        if (n_sp != n_remove){
          current_rasters[[depth_i]] <-
            current_rasters[[depth_i]][[-remove_these]]
          add(current_rasters[[depth_i]]) <- rep(nan_rast, n_remove)
          current_rasters[[depth_i]] <-
            current_rasters[[depth_i]][[order(c(one_n_sp[-remove_these],
                                                one_n_sp[remove_these]))]]

        } else {
          current_rasters[[depth_i]] <- rep(nan_rast, n_remove)
        }
        names(current_rasters[[depth_i]]) <- names(biodiv_raster)
      }

    }

  }

  names(current_rasters) <- rev(levels(cut(0, breaks = breaks,
                                           include.lowest = TRUE)))

  return(current_rasters)
}

## Function to set problems
## It creates a prioritizr problem
.set_problem <- function(pu_rast, features, penalty, edge_factor,
                         budget, gap, weight_data, threads, portfolio,
                         portfolio_opts, locked_in_raster, locked_out_raster,
                         verbose){

  init_problem <-
    problem(pu_rast, features) |>
    add_binary_decisions() |>
    add_default_solver(gap = gap, verbose = verbose, threads = threads)

  if (any(penalty != 0)){
    init_problem <- init_problem |> add_boundary_penalties(penalty,edge_factor)
  }

  init_problem <- init_problem |> add_max_utility_objective(budget)

  if ( !is.null(weight_data) ){
    init_problem <- init_problem |> add_feature_weights(weight_data)
  }

  if (portfolio == "extra"){
    init_problem <- init_problem |> add_extra_portfolio()
  } else if (portfolio == "top"){
    if (length(portfolio_opts) == 0){
      init_problem <- init_problem |> add_top_portfolio()
    } else {
      init_problem <-
        do.call(add_top_portfolio, c(init_problem, as.list(portfolio_opts)))
    }
  } else if (portfolio == "gap"){
    if (length(portfolio_opts) == 0){
      init_problem <- init_problem |> add_gap_portfolio()
    } else {
      init_problem <-
        do.call(add_gap_portfolio, c(init_problem, as.list(portfolio_opts)))
    }
  } else if (portfolio == "cuts"){
    if (length(portfolio_opts) == 0){
      init_problem <- init_problem |> add_cuts_portfolio()
    } else {
      init_problem <-
        do.call(add_cuts_portfolio, c(init_problem, as.list(portfolio_opts)))
    }
  } else if (portfolio == "shuffle"){
    if (length(portfolio_opts) == 0){
      init_problem <- init_problem |> add_shuffle_portfolio()
    } else {
      init_problem <-
        do.call(add_shuffle_portfolio, c(init_problem, as.list(portfolio_opts)))
    }
  }

  if (!is.null(locked_in_raster)){
    if (global(locked_in_raster, "sum", na.rm=TRUE) == 0){
      locked_in_raster <- NULL
    }
  }

  if (!is.null(locked_in_raster)){
    init_problem <- init_problem |> add_locked_in_constraints(locked_in_raster)
  }

  if (!is.null(locked_out_raster)){
    if (global(locked_out_raster, "sum", na.rm=TRUE) == 0){
      locked_out_raster <- NULL
    }
  }

  if (!is.null(locked_out_raster)){
    init_problem <- init_problem |>
      add_locked_out_constraints(locked_out_raster)
  }

  return(init_problem)
}

## Main function for 3D prioritization
## Solves 3D problem for a single budget
.single_3D <- function(split_features,
                       depth_raster,
                       breaks,
                       biodiv_df = NULL,
                       priority_weights = NULL,
                       budget_percent=0.3,
                       budget_weights="equal",
                       penalty = 0,
                       edge_factor = 0.5,
                       gap=0.1,
                       threads = 1L,
                       portfolio="gap",
                       portfolio_opts=list(number_solutions = 10,
                                           pool_gap = 0.1),
                       locked_in_raster=NULL,
                       locked_out_raster=NULL,
                       verbose = FALSE){

  weight_data <- NULL
  if (!is.null(biodiv_df)){
    group_pos <- which(names(biodiv_df) == "group")
    if (!is.null(priority_weights) && (length(group_pos) == 1)){
      biodiv_df <- biodiv_df[biodiv_df$species_name %in%
                               names(split_features[[1]]),]
      weight_data <- biodiv_df[, c("species_name", "group")]
      if (!all(c("group", "weight") %in% names(priority_weights))) {
        stop('priority_weights should have column names c("group", "weight")')
      }
      priority_weights <- priority_weights[,c("group", "weight")]
      priority_weights <- priority_weights[priority_weights[,1] %in%
                                             weight_data[,2],]
    } else if (length(group_pos) > 1){
      stop('weight_data should have only one group column')
    } else if ((is.null(priority_weights) && (length(group_pos) == 1)) ||
               (!is.null(priority_weights) && (length(group_pos) == 0))){
      stop(
        paste0('Define both priority_weights data frame and "group" column of ",
               "biodiv_df of neither of them')
        )
    }
  }

  split_features <- rev(split_features)
  if ( !is.null(priority_weights) && !is.null(weight_data)){
    names(priority_weights) <- c("group", "weight")
    names(weight_data) <- c("names", "group")
    weight_data <- merge(data.frame("names"=names(split_features[[1]])),
                         weight_data, all.x=TRUE)
    weight_data$id <- 1:dim(weight_data)[1]
    weight_data <-
      merge(weight_data, priority_weights, by="group", all.x=TRUE, all.y=FALSE)
    weight_data <- as.matrix(weight_data[order(weight_data[,2]),4])
  }

  length_split_features <- length(split_features)
  rev_breaks <- rev(breaks)
  pu_rast <-
    subst(depth_raster >= rev_breaks[1] & depth_raster <= rev_breaks[2],
          from=FALSE, to=NA)

  if (!is.null(locked_in_raster)){
    locked_in_raster_i <- mask(locked_in_raster, pu_rast[[1]])
    pu_rast[locked_in_raster_i] <- 0
  } else {
    locked_in_raster_i <- NULL
  }

  if (!is.null(locked_out_raster)){
    locked_out_raster_i <- mask(locked_out_raster, pu_rast[[1]])
  } else {
    locked_out_raster_i <- NULL
  }


  for (i in 2:length_split_features){
    add(pu_rast) <- !terra::allNA(c(pu_rast[[1:i-1]],
                                    subst(depth_raster > rev_breaks[i] &
                                            depth_raster <= rev_breaks[i+1],
                                          from=FALSE, to=NA)))
    pu_rast[[i]] <- terra::subst(pu_rast[[i]], from=FALSE, to=NA)

    if (!is.null(locked_in_raster)){
      add(locked_in_raster_i) <- mask(locked_in_raster, pu_rast[[i]])
      pu_rast[[i]][locked_in_raster_i[[i]]]<- 0
    } else {
      locked_in_raster_i <- NULL
    }

    if (!is.null(locked_out_raster)){
      add(locked_out_raster_i) <- mask(locked_out_raster, pu_rast[[i]])
      pu_rast[[i]][locked_out_raster_i[[i]]]<- 0
    } else {
      locked_out_raster_i <- NULL
    }

  }

  total_budget <-
    budget_percent * as.numeric(global(pu_rast[[length_split_features]],
                                       "sum", na.rm=TRUE))

  if (length(budget_weights) == 1 && is.character(budget_weights)){
    if (budget_weights == "equal"){
      budget_weights <- rep(1, length_split_features)
    } else if (budget_weights == "area"){
      budget_weights <-
        as.numeric(unlist(terra::global(pu_rast, "sum", na.rm=TRUE)))
    } else if (budget_weights == "richness"){
      budget_weights <- rev(sapply(split_features, function(x){
        to_keep <- unlist(global(x, "max", na.rm=TRUE))
        sum((to_keep > 0) & !is.na(to_keep))}))
    }
  }
  if (sum(budget_weights) > 1)
    budget_weights <- budget_weights / sum(budget_weights)
  budget_weights <- rev(budget_weights)

  budget <- budget_weights * as.integer(total_budget)

  budget_new <- budget
  budget_new <- floor(budget)
  indices <- tail(order(budget-budget_new),
                  round(sum(budget)) - sum(budget_new))
  budget_new[indices] <- budget_new[indices] + 1
  budget <- budget_new

  ## Features of first depth
  # Find and remove features with no values at this depth
  # Note that we keep initial feature maps as they are
  to_keep <- unlist(global(split_features[[1]], "max", na.rm=TRUE))
  to_keep <- (to_keep > 0) & !is.na(to_keep)

  if (any(to_keep)){

    # Surplus if the budget given is greater than
    # depth level number of planning units
    surplus <- budget[1] - as.numeric(global(pu_rast[[1]], "sum", na.rm=TRUE))
    if (surplus > 0){
      budget[1] <- budget[1] - surplus
      budget[2] <- budget[2] + surplus
    }

    # Surplus if budget given is greater than
    # depth level number of existing features
    surplus <- budget[1] -
      as.numeric(global(!allNA(split_features[[1]]) &
                          any(split_features[[1]] > 0, na.rm=TRUE),
                        "sum", na.rm=TRUE))
    if (surplus > 0){
      budget[1] <- budget[1] - surplus
      budget[2] <- budget[2] + surplus
    }

    init_problem <- .set_problem(pu_rast[[1]], split_features[[1]][[to_keep]],
                                 penalty, edge_factor,
                                 budget[1], gap, weight_data[to_keep], threads,
                                 portfolio, portfolio_opts,
                                 locked_in_raster_i[[1]],
                                 locked_out_raster_i[[1]], verbose)
    new_cost <- solve(init_problem, force = TRUE)

    if (length(new_cost) > 1){
      # Select the best solution
      new_cost <- new_cost[[1]]
    }
    # Get it ready for the next depth
    # Which cells will be removed is based from the next depth
    new_cost <- terra::subst(new_cost, from=NA, to=0)
  } else {
    budget[2] <- budget[2] + budget[1]
    new_cost <- terra::subst(pu_rast[[1]] * 0, from=NA, to=0)
  }


  ## For the next ones depths
  for (i in 2:length_split_features){

    to_keep <- unlist(global(split_features[[i]], "max", na.rm=TRUE))
    to_keep <- (to_keep > 0) & !is.na(to_keep)
    if (any(to_keep)){


      ## Move solutions of lower levels to the upper ones
      mew_pokemon <- (all(is.na(split_features[[i]][[to_keep]]) |
                            (split_features[[i]][[to_keep]] == 0)) &
                        !is.na(pu_rast[[i]])) * 1e-4
      names(mew_pokemon) <- "mew_pokemon"

      ## 1 - priority map inside the area of current depth
      pu_rast[[i]] <- mask(1 - new_cost, pu_rast[[i]])

      if (!is.null(locked_in_raster_i[[i]])){
        pu_rast[[i]][locked_in_raster_i[[i]]] <- 0
      }

      if (!is.null(locked_out_raster_i[[i]])){
        pu_rast[[i]][locked_out_raster_i[[i]]] <- 0
      }

      # Surplus if budget given is greater than
      # depth level number of planning units
      surplus <- budget[i] - as.numeric(global(pu_rast[[i]], "sum", na.rm=TRUE))
      if ((surplus > 0) && (i != length_split_features)){
        budget[i] <- budget[i] - surplus
        budget[i+1] <- budget[i+1] + surplus
      }

      # Surplus if budget given is greater than
      # depth level number of existing features
      if ( i != length_split_features ){
        surplus <- budget[i] -
          as.numeric(
            global(!allNA(split_features[[i]]) & any(split_features[[i]] > 0,
                            na.rm=TRUE), "sum", na.rm=TRUE))
        if (surplus > 0){
          budget[i] <- budget[i] - surplus
          budget[i+1] <- budget[i+1] + surplus
        }
      }
      mew_bool <- as.numeric(global(mew_pokemon, "max", na.rm=TRUE)) == 0
      if (ifelse(length(mew_bool) > 0, mew_bool, FALSE)) {
        use_features <- split_features[[i]][[to_keep]]
        weight_features <- weight_data[to_keep]
      } else {
        use_features <- c(mew_pokemon,split_features[[i]][[to_keep]])
        if (!is.null(weight_data)){
          weight_features <- c(1, weight_data[to_keep])
        } else {
          weight_features <- weight_data
        }
      }
      init_problem <- .set_problem(pu_rast[[i]], use_features, penalty,
                                   edge_factor, budget[i], gap, weight_features,
                                   threads, portfolio, portfolio_opts,
                                   locked_in_raster_i[[i]],
                                   locked_out_raster_i[[i]], verbose)
      new_cost <- solve(init_problem, force = TRUE)
      if (length(new_cost) > 1){
        # Select the best solution
        new_cost <- new_cost[[1]]
      }
      if (i != length_split_features)
        new_cost <- terra::subst(new_cost, from=NA, to=0)
    } else {
      if (i != length_split_features)
        budget[i+1] <- budget[i+1] + budget[i]
    }
  }

  return(new_cost)

}

## Solves a casual 2D prioritization problem
.prioritize_2D <- function(biodiv_raster,
                           pu_rast,
                           biodiv_df,
                           priority_weights=NULL,
                           budget_percent=0.3,
                           penalty = 0,
                           edge_factor = 0.5,
                           gap=0.1,
                           threads = 1L,
                           portfolio="gap",
                           portfolio_opts=list(number_solutions = 10,
                                               pool_gap = 0.1),
                           locked_in_raster=NULL,
                           locked_out_raster=NULL,
                           verbose = FALSE){

  to_keep <- unlist(global(biodiv_raster, "max", na.rm=TRUE))
  to_keep <- (to_keep > 0) & !is.na(to_keep)
  total_budget <-
    budget_percent * as.numeric(global(pu_rast, "sum", na.rm=TRUE))
  weight_data <- NULL
  group_pos <- which(names(biodiv_df) == "group")
  if (!is.null(priority_weights) && (length(group_pos) == 1)){
    biodiv_df <- biodiv_df[biodiv_df$species_name %in% names(biodiv_raster),]
    weight_data <- biodiv_df[, c("species_name", "group")]
    if (!all(c("group", "weight") %in% names(priority_weights))) {
      stop('priority_weights should have column names c("group", "weight")')
    }
    priority_weights <- priority_weights[,c("group", "weight")]
    priority_weights <-
      priority_weights[priority_weights[,1] %in% weight_data[,2],]
  } else if (length(group_pos) > 1){
    stop('weight_data should have only one group column')
  } else if ((is.null(priority_weights) && (length(group_pos) == 1)) ||
             (!is.null(priority_weights) && (length(group_pos) == 0))){
    stop(paste0('Define both priority_weights data frame and "group" column of",
                " biodiv_df of neither of them'))
  }

  if ( !is.null(priority_weights) && (length(group_pos) > 0)){
    names(priority_weights) <- c("group", "weight")
    names(weight_data) <- c("names", "group")
    weight_data <-
      merge(data.frame("names"=names(biodiv_raster)), weight_data, all.x=TRUE)
    weight_data$id <- 1:dim(weight_data)[1]
    weight_data <-
      merge(weight_data, priority_weights, by="group", all.x=TRUE, all.y=FALSE)
    weight_data <- as.matrix(weight_data[order(weight_data[,2]),4])
  }
  if (!is.null(locked_in_raster)){
    locked_in_raster_i <- mask(locked_in_raster, pu_rast)
    pu_rast[locked_in_raster_i] <- 0
  } else {
    locked_in_raster_i <- NULL
  }

  if (!is.null(locked_out_raster)){
    locked_out_raster_i <- mask(locked_out_raster, pu_rast)
  } else {
    locked_out_raster_i <- NULL
  }

  init_problem <- .set_problem(pu_rast, biodiv_raster[[to_keep]], penalty,
                               edge_factor, total_budget, gap,
                               weight_data=weight_data[to_keep],
                               threads, portfolio, portfolio_opts,
                               locked_in_raster_i, locked_out_raster_i, verbose)
  new_cost <- solve(init_problem, force = TRUE)

  if (length(new_cost) > 1){
    # Select the best solution
    new_cost <- new_cost[[1]]
  }

  return(new_cost)
}

## Function to evaluate prioritization solution over 3D feature distributions
evaluate_3D <- function(solution, split_features){

  n_depths <- length(split_features)


  absolute_held <- total_amount <-
    matrix(0, nrow=nlyr(split_features[[1]]), ncol=n_depths)

  for (i in 1:n_depths){
    absolute_held[,i] <-
      as.numeric(unlist(global(
        split_features[[i]] * solution,"sum",na.rm=TRUE)))
    total_amount[,i]  <-
      as.numeric(unlist(global(split_features[[i]],"sum",na.rm=TRUE)))
  }
  relative_held_raw <- absolute_held / total_amount
  overall_available <- absolute_held / base::rowSums(total_amount, na.rm=TRUE)
  depth_overall_available <- base::colMeans(overall_available, na.rm=TRUE)
  relative_held <- base::colMeans(relative_held_raw, na.rm=TRUE)
  overall_held <- base::rowSums(absolute_held, na.rm=TRUE) /
                  base::rowSums(total_amount, na.rm=TRUE)
  names(overall_held) <- names(split_features[[1]])
  return(list(relative_held_raw=relative_held_raw, relative_held=relative_held,
              overall_held=overall_held, overall_available=overall_available,
              depth_overall_available=depth_overall_available,
              absolute_held=absolute_held, total_amount = total_amount))
}

## Function to compute Jaccard coefficient between two
## binary SpatRaster objects(=intersection/union)
terra_jaccard <- function(x, y){
  sum_solutions <- sum(x, y, na.rm=TRUE)
  as.numeric(terra::global(sum_solutions == 2, "sum", na.rm=TRUE) /
               terra::global(sum_solutions > 0, "sum", na.rm=TRUE))
}

## Function to sum list of SpatRaster objects
sumrast <- function(x, normalize=TRUE){
  x <- sum(terra::rast(x), na.rm=TRUE)
  if (normalize){
    max_x <- as.numeric(terra::global(x, "max", na.rm=TRUE))
    min_x <- as.numeric(terra::global(x, "min", na.rm=TRUE))
    x <- (x-min_x)/(max_x-min_x)
  }
  return(x)
}

## Function to plot sum list of SpatRaster objects
plot_sumrast <- function(x, normalize=TRUE, add_lines=TRUE, ...){
  terra::plot(sumrast(x, normalize=normalize), ...)
  if (add_lines){
    maps::map(add=TRUE)
  }
}

## Function to plot the output of prioritize_3D
plot_3D <- function(x, to_plot="all", add_lines=TRUE){
  if (!is.null(attr(x, "from.function"))){
    if (attr(x, "from.function") == "prioritize_3D"){
      if (to_plot=="all"){
        if (length(x$budget_percent) > 1){
          cols <- rev(viridis(length(x$depth_levels_names)))
          layout(mat = matrix(c(1,2,3),nrow = 3,ncol = 1, byrow = TRUE),
                 heights = c(0.4, 0.4, 0.2))
          plot(classify(sumrast(x$solution3D), c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
                        include.lowest=TRUE, brackets=TRUE),
               col = rev(terrain.colors(255)),
               main="3D\nSum plot", all_levels=TRUE)
          if (add_lines) maps::map(add=TRUE)
          matplot(x$budget_percent, x$relative_helds3D, type="l", col=cols,
                  lwd=3, lty=1, axes=FALSE, xlim = rev(range(x$budget_percent)),
                  main="Relative held per depth", xlab="%PU used",
                  ylab="Relative held")
          axis(1, at = seq(0, 1, 0.1), labels = rep("",11), las=1)
          axis(1, at = seq(1, 0, -0.2), labels = seq(100, 0, -20), las=1)
          axis(2, at = seq(0, 1, 0.1), labels = rep("",11), las=1)
          axis(2, at = seq(0, 1, 0.2), labels = seq(0, 100, 20), las=1)
          plot(1, type = "n", axes=FALSE, xlab="", ylab="")
          legend(x = "top",inset = 0,
                 legend = rev(x$depth_levels_names),
                 col=cols, lwd=3, lty=1, cex=1.1, box.col = NA, horiz = TRUE)
          layout(mat = matrix(1))
        } else {
          cols <- rev(viridis(length(x$depth_levels_names)))
          layout(mat = matrix(c(1,2),nrow = 2,ncol = 1, byrow = TRUE))
          plot(x$solution3D[[1]], main="3D\nSum plot",
               col = rev(terrain.colors(255)))
          if (add_lines) maps::map(add=TRUE)
          barplot(x$relative_helds3D,
                  main = "Relative held", ylim=c(0,1), yaxt='n')
          axis(2, at = seq(0, 1, 0.2), labels = seq(0, 100, 20), las=1)
          layout(mat = matrix(1))
        }

      } else if (to_plot == "maps"){
        if (length(x$budget_percent) > 1){
          layout(mat = matrix(1))
          plot(classify(sumrast(x$solution3D), c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
                        include.lowest=TRUE, brackets=TRUE),
               main="3D\nSum plot", all_levels=TRUE,
               col = rev(terrain.colors(255)))
          if (add_lines) maps::map(add=TRUE)
        } else {
          layout(mat = matrix(1))
          plot(x$solution3D[[1]], main="3D\nSum plot",
                        col = rev(terrain.colors(255)))
          if (add_lines) maps::map(add=TRUE)
        }

      } else if (to_plot == "relative_held"){
        if (length(x$budget_percent) > 1){
          cols <- rev(viridis(length(x$depth_levels_names)))
          layout(mat = matrix(c(1,2),nrow = 2,ncol = 1,byrow = TRUE),
                 heights = c(0.7, 0.3))
          matplot(x$budget_percent, x$relative_helds3D, type="l", col=cols,
                  lwd=3, lty=1, axes=FALSE, xlim = rev(range(x$budget_percent)),
                  main="3D\nRelative held per depth", xlab="%PU used",
                  ylab="Relative held")
          axis(1, at = seq(0, 1, 0.1), labels = rep("",11), las=1)
          axis(1, at = seq(1, 0, -0.2), labels = seq(100, 0, -20), las=1)
          axis(2, at = seq(0, 1, 0.1), labels = rep("",11), las=1)
          axis(2, at = seq(0, 1, 0.2), labels = seq(0, 100, 20), las=1)
          plot(1, type = "n", axes=FALSE, xlab="", ylab="")
          legend(x = "top",inset = 0,
                 legend = rev(x$depth_levels_names),
                 col=cols, lwd=3, lty=1, cex=1.1, box.col = NA, horiz = TRUE)
          layout(mat = matrix(1))
        } else {
          layout(mat = matrix(1))
          barplot(x$relative_helds3D, main = "3D\nRelative held", ylim=c(0,1),
                  yaxt='n')
          axis(2, at = seq(0, 1, 0.2), labels = seq(0, 100, 20), las=1)
        }
      } else {
        stop("to_plot must be one of 'all', 'maps', 'relative_held'")
      }
    }
  }
}

## Function to plot the output of Compare_2D_3D
plot_Compare_2D_3D <- function(x, to_plot="all", add_lines=TRUE){
  if (!is.null(attr(x, "from.function"))){
    if (attr(x, "from.function") == "Compare_2D_3D"){
      if (to_plot=="all"){
        if (length(x$budget_percent) > 1){
          cols <- rev(viridis(length(x$depth_levels_names)))
          layout(mat = matrix(c(1,2,3, 4, 5, 5),nrow = 3,ncol = 2,byrow = TRUE),
                 heights = c(0.4, 0.4, 0.2))
          plot(classify(sumrast(x$solution2D), c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
                        include.lowest=TRUE, brackets=TRUE),
               col = rev(terrain.colors(255)),
               main="2D\nSum plot", all_levels=TRUE)
          if (add_lines) maps::map(add=TRUE)
          plot(classify(sumrast(x$solution3D), c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
                        include.lowest=TRUE, brackets=TRUE),
               col = rev(terrain.colors(255)),
               main="3D\nSum plot", all_levels=TRUE)
          if (add_lines) maps::map(add=TRUE)
          ylim <- c(min(x$relative_helds2D, x$relative_helds3D),
                    max(x$relative_helds2D, x$relative_helds3D))
          matplot(x$budget_percent, x$relative_helds2D, type="l",
                  col=cols, lwd=3, lty=1, axes=FALSE,
                  xlim = rev(range(x$budget_percent)), ylim = ylim,
                  main="Relative held per depth", xlab="%PU used",
                  ylab="Relative held")
          axis(1, at = seq(0, 1, 0.1), labels = rep("",11), las=1)
          axis(1, at = seq(1, 0, -0.2), labels = seq(100, 0, -20), las=1)
          axis(2, at = seq(0, 1, 0.1), labels = rep("",11), las=1)
          axis(2, at = seq(0, 1, 0.2), labels = seq(0, 100, 20), las=1)
          matplot(x$budget_percent, x$relative_helds3D, type="l", col=cols,
                  lwd=3, lty=1, axes=FALSE, xlim = rev(range(x$budget_percent)),
                  ylim = ylim, main="Relative held per depth", xlab="%PU used",
                  ylab="Relative held")
          axis(1, at = seq(0, 1, 0.1), labels = rep("",11), las=1)
          axis(1, at = seq(1, 0, -0.2), labels = seq(100, 0, -20), las=1)
          axis(2, at = seq(0, 1, 0.1), labels = rep("",11), las=1)
          axis(2, at = seq(0, 1, 0.2), labels = seq(0, 100, 20), las=1)
          plot(1, type = "n", axes=FALSE, xlab="", ylab="")
          legend(x = "top",inset = 0,
                 legend = rev(x$depth_levels_names),
                 col=cols, lwd=3, lty=1, cex=1.1, box.col = NA, horiz = TRUE)
          layout(mat = matrix(1))
        } else {
          cols <- rev(viridis(length(x$depth_levels_names)))
          layout(mat = matrix(c(1,2, 3, 4),nrow = 2,ncol = 2,byrow = TRUE),
                 heights = c(0.4, 0.4, 0.2))
          plot(x$solution2D[[1]], main="2D\n\nSum plot",
               col = rev(terrain.colors(255)))
          if (add_lines) maps::map(add=TRUE)
          plot(x$solution3D[[1]], main="3D\n\nSum plot",
               col = rev(terrain.colors(255)))
          if (add_lines) maps::map(add=TRUE)
          barplot(x$relative_helds2D, main = "2D\nRelative held",
                  ylim=c(0,1), yaxt='n')
          axis(2, at = seq(0, 1, 0.2), labels = seq(0, 100, 20), las=1)
          barplot(x$relative_helds3D, main = "3D\nRelative held",
                  ylim=c(0,1), yaxt='n')
          axis(2, at = seq(0, 1, 0.2), labels = seq(0, 100, 20), las=1)
          layout(mat = matrix(1))
        }

      } else if (to_plot == "maps"){
        if (length(x$budget_percent) > 1){
          layout(mat = matrix(c(1,2),nrow = 1,ncol = 2,byrow = TRUE))
          plot(classify(sumrast(x$solution2D), c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
                        include.lowest=TRUE, brackets=TRUE),
               col = rev(terrain.colors(255)),
               main="2D\nSum plot", all_levels=TRUE)
          if (add_lines) maps::map(add=TRUE)
          plot(classify(sumrast(x$solution3D), c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
                        include.lowest=TRUE, brackets=TRUE),
               col = rev(terrain.colors(255)),
               main="3D\nSum plot", all_levels=TRUE)
          if (add_lines) maps::map(add=TRUE)
          layout(mat = matrix(1))
        } else {
          layout(mat = matrix(c(1,2),nrow = 1,ncol = 2,byrow = TRUE))
          plot(x$solution2D[[1]], main="2D\nSum plot",
               col = rev(terrain.colors(255)))
          if (add_lines) maps::map(add=TRUE)
          plot(x$solution3D[[1]], main="3D\nSum plot",
               col = rev(terrain.colors(255)))
          if (add_lines) maps::map(add=TRUE)
          layout(mat = matrix(1))
        }

      } else if (to_plot == "relative_held"){
        if (length(x$budget_percent) > 1){
          cols <- rev(viridis(length(x$depth_levels_names)))
          layout(mat = matrix(c(1,2, 3, 3),nrow = 2,ncol = 2,byrow = TRUE),
                 heights = c(0.7, 0.3))
          ylim <- c(min(x$relative_helds2D, x$relative_helds3D),
                    max(x$relative_helds2D, x$relative_helds3D))
          matplot(x$budget_percent, x$relative_helds2D, type="l", col=cols,
                  lwd=3, lty=1, axes=FALSE, xlim = rev(range(x$budget_percent)),
                  ylim = ylim, main="2D\nRelative held per depth",
                  xlab="%PU used", ylab="Relative held")
          axis(1, at = seq(0, 1, 0.1), labels = rep("",11), las=1)
          axis(1, at = seq(1, 0, -0.2), labels = seq(100, 0, -20), las=1)
          axis(2, at = seq(0, 1, 0.1), labels = rep("",11), las=1)
          axis(2, at = seq(0, 1, 0.2), labels = seq(0, 100, 20), las=1)
          matplot(x$budget_percent, x$relative_helds3D, type="l", col=cols,
                  lwd=3, lty=1, axes=FALSE, xlim = rev(range(x$budget_percent)),
                  ylim = ylim,main="3D\nRelative held per depth",
                  xlab="%PU used", ylab="Relative held")
          axis(1, at = seq(0, 1, 0.1), labels = rep("",11), las=1)
          axis(1, at = seq(1, 0, -0.2), labels = seq(100, 0, -20), las=1)
          axis(2, at = seq(0, 1, 0.1), labels = rep("",11), las=1)
          axis(2, at = seq(0, 1, 0.2), labels = seq(0, 100, 20), las=1)
          plot(1, type = "n", axes=FALSE, xlab="", ylab="")
          legend(x = "top",inset = 0,
                 legend = rev(x$depth_levels_names),
                 col=cols, lwd=3, lty=1, cex=1.1, box.col = NA, horiz = TRUE)
          layout(mat = matrix(1))
        } else {
          layout(mat = matrix(c(1,2),nrow = 1,ncol = 2,byrow = TRUE))
          barplot(x$relative_helds2D, main = "2D\nRelative held",
                  ylim=c(0,1), yaxt='n')
          axis(2, at = seq(0, 1, 0.2), labels = seq(0, 100, 20), las=1)
          barplot(x$relative_helds3D, main = "3D\nRelative held",
                  ylim=c(0,1), yaxt='n')
          axis(2, at = seq(0, 1, 0.2), labels = seq(0, 100, 20), las=1)
          layout(mat = matrix(1))
        }
      } else {
        stop("to_plot must be one of 'all', 'maps', 'relative_held'")
      }
    }
  }
}

Compare_2D_3D <- function(biodiv_raster,
                          depth_raster,
                          breaks,
                          biodiv_df,
                          val_depth_range=TRUE,
                          priority_weights=NULL,
                          budget_percents=seq(0,1,0.1),
                          budget_weights="equal",
                          penalty = 0,
                          edge_factor = 0.5,
                          gap=0.1,
                          threads=1L,
                          sep_priority_weights=",",
                          portfolio="gap",
                          portfolio_opts=list(number_solutions = 10,
                                              pool_gap = 0.1),
                          sep_biodiv_df=",",
                          locked_in_raster=NULL,
                          locked_out_raster=NULL,
                          verbose = FALSE){

  ## Ensure that biodiv_raster is rast file with depths or path of its file
  if (is.character(biodiv_raster)){
    if (tools::file_ext(biodiv_raster) != ""){
      biodiv_raster <- terra::rast(biodiv_raster) * 1
    } else {
      biodiv_raster <- get_rast(biodiv_raster) * 1
    }
  } else if (!is(biodiv_raster, "SpatRaster")) {
    stop("biodiv_raster can be only either character or SpatRaster class")
  }

  ## Ensure that depth_raster is rast file with depths or path of its file
  if (is.character(depth_raster)){
    depth_raster <- terra::rast(depth_raster)  * 1
  } else if (!is(depth_raster, "SpatRaster")) {
    stop("depth_raster can be only either character or SpatRaster class")
  }
  breaks <- sort(breaks, decreasing = TRUE)
  out_of_analysis <-
    depth_raster > breaks[1] | depth_raster < breaks[length(breaks)]
  depth_raster[out_of_analysis] <- NA


  ## Ensure that locked_in_raster is rast file with depths or path of its file
  if (!is.null(locked_in_raster)){
    if (is.character(locked_in_raster)){
      locked_in_raster <- terra::rast(locked_in_raster)  * 1
    } else if (!is(locked_in_raster, "SpatRaster")) {
      stop("locked_in_raster can be only either character or SpatRaster class")
    }
  }

  ## Ensure that locked_in_raster is rast file with depths or path of its file
  if (!is.null(locked_out_raster)){
    if (is.character(locked_out_raster)){
      locked_out_raster <- terra::rast(locked_out_raster)  * 1
    } else if (!is(locked_out_raster, "SpatRaster")) {
      stop("locked_out_raster can be only either character or SpatRaster class")
    }
  }

  if (!is.null(priority_weights)){
    if (is.character(priority_weights)){
      if (file.exists(priority_weights))
        priority_weights <- switch(tools::file_ext(priority_weights),
                                   csv = read.csv(priority_weights,
                                                  header = TRUE,
                                                  sep=sep_priority_weights),
                                   xls = readxl::read_xls(priority_weights),
                                   xlsx = readxl::read_xlsx(priority_weights)
        )
    } else if (!is(priority_weights, "data.frame")) {
      stop(
        "priority_weights can be only either character path or data.frame class"
        )
    }
  }

  if (is.character(biodiv_df)){
    if (file.exists(biodiv_df))
      biodiv_df <- switch(tools::file_ext(biodiv_df),
                          csv = read.csv(biodiv_df, header = TRUE,
                                         sep=sep_biodiv_df),
                          xls = readxl::read_xls(biodiv_df),
                          xlsx = readxl::read_xlsx(biodiv_df)
      )
  } else if (!is(biodiv_df, "data.frame")) {
    stop("biodiv_df can be only either character path or data.frame class")
  }

  #################################
  proj_depth_raster <- any(res(depth_raster) != res(biodiv_raster)) ||
                       (crs(depth_raster) != crs(biodiv_raster))

  if ( proj_depth_raster ){
    warning(paste0("Different resolution/coordinates among biodiv_raster and ",
                   "depth_raster:\nproject biodiv_raster on depth_raster"))
    depth_raster <- terra::project(depth_raster,biodiv_raster, method="average")
  }

  ext_depth_raster <- ext(depth_raster)
  ext_biodiv_raster <- ext(biodiv_raster)
  if (any(ext_depth_raster != ext_biodiv_raster)){
    warning("Different extents in input SpatRasters")
    min_x <- max( c(ext_depth_raster[1], ext_biodiv_raster[1] ) )
    max_x <- min( c(ext_depth_raster[2], ext_biodiv_raster[2] ) )
    min_y <- max( c(ext_depth_raster[3], ext_biodiv_raster[3] ) )
    max_y <- min( c(ext_depth_raster[4], ext_biodiv_raster[4] ) )
    exts <- terra::ext(min_x, max_x, min_y, max_y)
    if (exts != ext_depth_raster)
      depth_raster <- crop(depth_raster, exts)
    if (exts != ext_biodiv_raster){
      biodiv_raster <- crop(biodiv_raster, exts)
    }
  }
  ## Transform 2D features distribution into a 3D one
  names(biodiv_raster) <- tolower(names(biodiv_raster))
  possible_names <- c("species_name", "pelagic", "min_z", "max_z", "group")
  possible_names <- possible_names[possible_names %in% names(biodiv_df)]
  if ( !any(c("species_name", "pelagic") %in%  possible_names) ){
    stop("Please provide species_name and pelagic in biodiv_df")
  } else if ( !("species_name" %in%  possible_names) ){
    stop("Please provide species_name in biodiv_df")
  } else if ( !("pelagic" %in%  possible_names) ){
    stop("Please provide pelagic in biodiv_df")
  }
  biodiv_df <- biodiv_df[,possible_names]

  biodiv_df[,names(biodiv_df) == "species_name"] <-
    tolower(biodiv_df[,names(biodiv_df) == "species_name"])
  biodiv_raster <-
    biodiv_raster[[names(biodiv_raster) %in% biodiv_df$species_name]]
  biodiv_df <- merge(data.frame(species_name=names(biodiv_raster)), biodiv_df)
  biodiv_df <-
    biodiv_df[order(match(biodiv_df$species_name,names(biodiv_raster))),]

  split_features <- split_rast(biodiv_raster, depth_raster, breaks, biodiv_df,
                               val_depth_range, sep_biodiv_df)

  ## Adjust target, penalty and edge_factor if names of split_features are less
  # from biodiv_raster
  to_keep <- names(split_features[[1]]) %in% biodiv_df$species_name

  if (length(penalty)>1 && !is.matrix(penalty) && !all(to_keep)) {
    penalty <- penalty[which(to_keep)]
  } else if (is.matrix(penalty)) {
    penalty <- penalty[which(to_keep),]
  }

  if (length(edge_factor)>1 && !is.matrix(edge_factor) && !all(to_keep)) {
    edge_factor <- edge_factor[which(to_keep)]
  } else if (is.matrix(edge_factor)) {
    edge_factor <- edge_factor[which(to_keep),]
  }

  ###############
  rev_split_features <- rev(split_features)
  n_depths <- length(split_features)
  pu_rast <- subst(depth_raster >= min(breaks) & depth_raster <= max(breaks),
                   from=FALSE, to=NA)
  if (!is.null(locked_in_raster)){
    pu_rast[locked_in_raster] <- 0
  }
  if (!is.null(locked_out_raster)){
    pu_rast[locked_out_raster] <- 0
  }
  pu_rast <- terra::mask(pu_rast, depth_raster)
  length_budget_percents <- length(budget_percents)
  absolute_held3D <- absolute_held2D <- NULL
  relative_helds3D <- relative_helds2D <- depth_overall_available3D <-
    depth_overall_available2D <-
    matrix(nrow=length_budget_percents, ncol=n_depths)
  mean_overall_helds3D <- mean_overall_helds2D <- sd_overall_helds3D <-
    sd_overall_helds2D <- numeric(length_budget_percents)
  solution3D <- solution2D <- list()
  overall_held3D <- overall_held2D <-
    matrix(nrow=length_budget_percents,ncol=nlyr(split_features[[1]]))
  jaccard_coef <- rep(0, length_budget_percents)
  for (i in 1:length_budget_percents){
    message(paste0("Budget: ", budget_percents[i]))
    solution3D[[i]] <- .single_3D(split_features=split_features,
                                  depth_raster=depth_raster,
                                  breaks=breaks,
                                  biodiv_df=biodiv_df,
                                  priority_weights=priority_weights,
                                  budget_percent=budget_percents[i],
                                  budget_weights=budget_weights,
                                  penalty = penalty,
                                  edge_factor = edge_factor,
                                  gap=gap,
                                  threads = threads,
                                  portfolio=portfolio,
                                  portfolio_opts=portfolio_opts,
                                  locked_in_raster=locked_in_raster,
                                  locked_out_raster=locked_out_raster,
                                  verbose = verbose)

    bar3D <- evaluate_3D(solution3D[[i]], rev_split_features)
    relative_helds3D[i,] <- bar3D$relative_held
    absolute_held3D[[i]] <- bar3D$absolute_held
    overall_held3D[i,] <- bar3D$overall_held
    mean_overall_helds3D[i] <- mean(bar3D$overall_held, na.rm=TRUE)
    sd_overall_helds3D[i] <- sd(bar3D$overall_held, na.rm=TRUE)
    depth_overall_available3D[i,] <- bar3D$depth_overall_available

    solution2D[[i]] <- .prioritize_2D(biodiv_raster=biodiv_raster[[to_keep]],
                                      pu_rast=pu_rast,
                                      biodiv_df=biodiv_df,
                                      priority_weights=priority_weights,
                                      budget_percent=budget_percents[i],
                                      penalty = penalty,
                                      edge_factor = edge_factor,
                                      gap=gap,
                                      threads = threads,
                                      portfolio=portfolio,
                                      portfolio_opts=portfolio_opts,
                                      locked_in_raster=locked_in_raster,
                                      locked_out_raster=locked_out_raster,
                                      verbose = verbose)

    bar2D <- evaluate_3D(solution2D[[i]], rev_split_features)
    relative_helds2D[i,] <- bar2D$relative_held
    absolute_held2D[[i]] <- bar2D$absolute_held
    overall_held2D[i,] <- bar2D$overall_held
    mean_overall_helds2D[i] <- mean(bar2D$overall_held, na.rm=TRUE)
    sd_overall_helds2D[i] <- sd(bar2D$overall_held, na.rm=TRUE)
    depth_overall_available2D[i,] <- bar2D$depth_overall_available

    jaccard_coef[i] <- terra_jaccard(solution2D[[i]], solution3D[[i]])

  }
  relative_helds2D[is.na(relative_helds2D)] <- 0
  relative_helds3D[is.na(relative_helds3D)] <- 0
  depth_overall_available2D[is.na(depth_overall_available2D)] <- 0
  depth_overall_available3D[is.na(depth_overall_available3D)] <- 0
  depth_overall_available2D <- depth_overall_available2D[,n_depths:1]
  depth_overall_available3D <- depth_overall_available3D[,n_depths:1]
  relative_helds2D <- relative_helds2D[,n_depths:1]
  relative_helds3D <- relative_helds3D[,n_depths:1]
  absolute_held2D <- lapply(absolute_held2D, function(x)x[,n_depths:1])
  absolute_held3D <- lapply(absolute_held3D, function(x)x[,n_depths:1])

  total_amount <- bar3D$total_amount[,n_depths:1]
  overall_total_amount <- base::rowSums(total_amount, na.rm=TRUE)

  names_features <- names(split_features[[1]])
  depth_levels_names <- levels(cut(0, breaks = breaks, include.lowest = TRUE))
  rev_depth_levels_names <- rev(depth_levels_names)

  colnames(overall_held3D) <- colnames(overall_held2D) <-
    rownames(total_amount) <-
    names(overall_total_amount) <-
    names_features

  names(solution3D) <- names(solution2D) <-
    names(absolute_held2D) <- names(absolute_held3D) <-
    rownames(overall_held3D) <- rownames(overall_held2D) <-
    names(mean_overall_helds3D) <- names(mean_overall_helds2D) <-
    names(jaccard_coef) <-
    paste0("budget", budget_percents)

  colnames(total_amount) <-
    names(split_features) <-
    rev_depth_levels_names

  if (length_budget_percents > 1){
    names(sd_overall_helds3D) <- names(sd_overall_helds2D) <-
      rownames(depth_overall_available3D) <-
      rownames(depth_overall_available2D) <-
      rownames(relative_helds3D) <- rownames(relative_helds2D) <-
      paste0("budget", budget_percents)

    colnames(depth_overall_available3D) <-
      colnames(depth_overall_available2D) <-
      colnames(relative_helds3D) <- colnames(relative_helds2D) <-
      rev_depth_levels_names

  } else {
    names(depth_overall_available3D) <- names(depth_overall_available2D) <-
      names(relative_helds3D) <- names(relative_helds2D) <-
    rev_depth_levels_names
  }

  absolute_held2D <- lapply(absolute_held2D,
                            function(x, rev_depth_levels_names, names_features){
                              colnames(x) <- rev_depth_levels_names
                              rownames(x) <- names_features
                              return(x)
                            },
  rev_depth_levels_names=rev_depth_levels_names,
  names_features = names_features)

  absolute_held3D <- lapply(absolute_held3D,
                            function(x, rev_depth_levels_names, names_features){
                              colnames(x) <- rev_depth_levels_names
                              rownames(x) <- names_features
                              return(x)
                            },
  rev_depth_levels_names=rev_depth_levels_names,
  names_features = names_features)


  output <- list(
    split_features=split_features,
    solution3D=solution3D,
    absolute_held3D=absolute_held3D,
    overall_held3D=overall_held3D,
    relative_helds3D=relative_helds3D,
    mean_overall_helds3D=mean_overall_helds3D,
    sd_overall_helds3D=sd_overall_helds3D,
    depth_overall_available3D = depth_overall_available3D,
    solution2D=solution2D,
    absolute_held2D=absolute_held2D,
    overall_held2D=overall_held2D,
    relative_helds2D=relative_helds2D,
    mean_overall_helds2D=mean_overall_helds2D,
    sd_overall_helds2D=sd_overall_helds2D,
    depth_overall_available2D = depth_overall_available2D,
    names_features = names_features,
    total_amount = total_amount,
    overall_total_amount = overall_total_amount,
    jaccard_coef = jaccard_coef,
    depth_levels_names=depth_levels_names,
    biodiv_raster=biodiv_raster[[to_keep]],
    depth_raster=depth_raster,
    breaks=breaks,
    biodiv_df=biodiv_df,
    val_depth_range=val_depth_range,
    priority_weights=priority_weights,
    budget_percents=budget_percents,
    budget_weights=budget_weights,
    penalty = penalty,
    edge_factor = edge_factor,
    gap=gap,
    threads=threads,
    sep_priority_weights=sep_priority_weights,
    portfolio=portfolio,
    portfolio_opts=portfolio_opts,
    sep_biodiv_df=sep_biodiv_df,
    locked_in_raster=locked_in_raster,
    locked_out_raster=locked_out_raster
  )
  attr(output, "from.function") <- "Compare_2D_3D"
  return(output)
}



prioritize_3D <- function(split_features,
                          depth_raster,
                          breaks,
                          biodiv_df,
                          priority_weights=NULL,
                          budget_percents=seq(0,1,0.1),
                          budget_weights="equal",
                          penalty = 0,
                          edge_factor = 0.5,
                          gap=0.1,
                          threads=1L,
                          sep_priority_weights=",",
                          portfolio="gap",
                          portfolio_opts=list(number_solutions = 10,
                                              pool_gap = 0.1),
                          sep_biodiv_df=",",
                          locked_in_raster=NULL,
                          locked_out_raster=NULL,
                          verbose = FALSE){

  ## Ensure that depth_raster is rast file with depths or path of its file
  if (is.character(depth_raster)){
    depth_raster <- terra::rast(depth_raster)  * 1
  } else if (!is(depth_raster, "SpatRaster")) {
    stop("depth_raster can be only either character or SpatRaster class")
  }
  breaks <- sort(breaks, decreasing = TRUE)
  out_of_analysis <- depth_raster > breaks[1] |
                     depth_raster < breaks[length(breaks)]
  depth_raster[out_of_analysis] <- NA


  ## Ensure that locked_in_raster is rast file with depths or path of its file
  if (!is.null(locked_in_raster)){
    if (is.character(locked_in_raster)){
      locked_in_raster <- terra::rast(locked_in_raster)  * 1
    } else if (!is(locked_in_raster, "SpatRaster")) {
      stop("locked_in_raster can be only either character or SpatRaster class")
    }
  }

  ## Ensure that locked_in_raster is rast file with depths or path of its file
  if (!is.null(locked_out_raster)){
    if (is.character(locked_out_raster)){
      locked_out_raster <- terra::rast(locked_out_raster)  * 1
    } else if (!is(locked_out_raster, "SpatRaster")) {
      stop("locked_out_raster can be only either character or SpatRaster class")
    }
  }

  if (!is.null(priority_weights)){
    if (is.character(priority_weights)){
      if (file.exists(priority_weights))
        priority_weights <- switch(tools::file_ext(priority_weights),
                                   csv = read.csv(priority_weights,
                                                  header = TRUE,
                                                  sep=sep_priority_weights),
                                   xls = readxl::read_xls(priority_weights),
                                   xlsx = readxl::read_xlsx(priority_weights)
        )
    } else if (!is(priority_weights, "data.frame")) {
      stop(
        "priority_weights can be only either character path or data.frame class"
        )
    }
  }

  if (is.character(biodiv_df)){
    if (file.exists(biodiv_df))
      biodiv_df <- switch(tools::file_ext(biodiv_df),
                          csv = read.csv(biodiv_df, header = TRUE,
                                         sep=sep_biodiv_df),
                          xls = readxl::read_xls(biodiv_df),
                          xlsx = readxl::read_xlsx(biodiv_df)
      )
  } else if (!is(biodiv_df, "data.frame")) {
    stop("biodiv_df can be only either character path or data.frame class")
  }

  #################################
  proj_depth_raster <- any(res(depth_raster) != res(split_features[[1]])) ||
    (crs(depth_raster) != crs(split_features[[1]]))

  if ( proj_depth_raster ){
    warning(paste0("Different resolution/coordinates among split_features[[1]]",
            " and depth_raster:\nproject split_features[[1]] on depth_raster"))
    depth_raster <- terra::project(depth_raster,split_features[[1]],
                                   method="average")
  }

  ext_depth_raster <- ext(depth_raster)
  ext_biodiv_raster <- ext(split_features[[1]])
  if (any(ext_depth_raster != ext_biodiv_raster)){
    warning("Different extents in input SpatRasters")
    min_x <- max( c(ext_depth_raster[1], ext_biodiv_raster[1] ) )
    max_x <- min( c(ext_depth_raster[2], ext_biodiv_raster[2] ) )
    min_y <- max( c(ext_depth_raster[3], ext_biodiv_raster[3] ) )
    max_y <- min( c(ext_depth_raster[4], ext_biodiv_raster[4] ) )
    exts <- terra::ext(min_x, max_x, min_y, max_y)
    if (exts != ext_depth_raster)
      depth_raster <- crop(depth_raster, exts)
    if (exts != ext_biodiv_raster){
      for (i in 1:length(split_features)){
        split_features[[i]] <- crop(split_features[[i]], exts)
      }
    }
  }

  for (i in 1:length(split_features)){
    names(split_features[[i]]) <- tolower(names(split_features[[i]]))
  }

  possible_names <- c("species_name", "pelagic", "min_z", "max_z", "group")
  possible_names <- possible_names[possible_names %in% names(biodiv_df)]
  if ( !any(c("species_name", "pelagic") %in%  possible_names) ){
    stop("Please provide species_name and pelagic in biodiv_df")
  } else if ( !("species_name" %in%  possible_names) ){
    stop("Please provide species_name in biodiv_df")
  } else if ( !("pelagic" %in%  possible_names) ){
    stop("Please provide pelagic in biodiv_df")
  }
  biodiv_df <- biodiv_df[,possible_names]

  biodiv_df[,names(biodiv_df) == "species_name"] <-
    tolower(biodiv_df[,names(biodiv_df) == "species_name"])
  for (i in 1:length(split_features)){
    split_features[[i]] <- split_features[[i]][[names(split_features[[i]]) %in%
                                                  biodiv_df$species_name]]
  }
  biodiv_df <-
    merge(data.frame(species_name=names(split_features[[1]])), biodiv_df)
  biodiv_df <-
    biodiv_df[order(match(biodiv_df$species_name,names(split_features[[1]]))),]

  ## Adjust target, penalty and edge_factor if names of split_features are less
  # from biodiv_raster
  to_keep <- names(split_features[[1]]) %in% biodiv_df$species_name

  if (length(penalty)>1 && !is.matrix(penalty) && !all(to_keep)) {
    penalty <- penalty[which(to_keep)]
  } else if (is.matrix(penalty)) {
    penalty <- penalty[which(to_keep),]
  }

  if (length(edge_factor)>1 && !is.matrix(edge_factor) && !all(to_keep)) {
    edge_factor <- edge_factor[which(to_keep)]
  } else if (is.matrix(edge_factor)) {
    edge_factor <- edge_factor[which(to_keep),]
  }

  ###############
  rev_split_features <- rev(split_features)
  n_depths <- length(split_features)
  pu_rast <- subst(depth_raster >= min(breaks) & depth_raster <= max(breaks),
                   from=FALSE, to=NA)
  if (!is.null(locked_in_raster)){
    pu_rast[locked_in_raster] <- 0
  }
  if (!is.null(locked_out_raster)){
    pu_rast[locked_out_raster] <- 0
  }
  pu_rast <- terra::mask(pu_rast, depth_raster)
  length_budget_percents <- length(budget_percents)
  absolute_held3D <- NULL
  relative_helds3D <- depth_overall_available3D <-
    matrix(nrow=length_budget_percents, ncol=n_depths)
  mean_overall_helds3D <- sd_overall_helds3D <- numeric(length_budget_percents)
  solution3D <- list()
  overall_held3D <-
    matrix(nrow=length_budget_percents,ncol=nlyr(split_features[[1]]))
  for (i in 1:length_budget_percents){
    message(paste0("Budget: ", budget_percents[i]))
    solution3D[[i]] <- .single_3D(split_features=split_features,
                                  depth_raster=depth_raster,
                                  breaks=breaks,
                                  biodiv_df=biodiv_df,
                                  priority_weights=priority_weights,
                                  budget_percent=budget_percents[i],
                                  budget_weights=budget_weights,
                                  penalty = penalty,
                                  edge_factor = edge_factor,
                                  gap=gap,
                                  threads = threads,
                                  portfolio=portfolio,
                                  portfolio_opts=portfolio_opts,
                                  locked_in_raster=locked_in_raster,
                                  locked_out_raster=locked_out_raster,
                                  verbose = verbose)

    bar3D <- evaluate_3D(solution3D[[i]], rev_split_features)
    relative_helds3D[i,] <- bar3D$relative_held
    absolute_held3D[[i]] <- bar3D$absolute_held
    overall_held3D[i,] <- bar3D$overall_held
    mean_overall_helds3D[i] <- mean(bar3D$overall_held, na.rm=TRUE)
    sd_overall_helds3D[i] <- sd(bar3D$overall_held, na.rm=TRUE)
    depth_overall_available3D[i,] <- bar3D$depth_overall_available

  }
  relative_helds3D[is.na(relative_helds3D)] <- 0
  depth_overall_available3D[is.na(depth_overall_available3D)] <- 0
  depth_overall_available3D <- depth_overall_available3D[,n_depths:1]
  relative_helds3D <- relative_helds3D[,n_depths:1]
  absolute_held3D <- lapply(absolute_held3D, function(x)x[,n_depths:1])

  total_amount <- bar3D$total_amount[,n_depths:1]
  overall_total_amount <- base::rowSums(total_amount, na.rm=TRUE)

  names_features <- names(split_features[[1]])
  depth_levels_names <- levels(cut(0, breaks = breaks, include.lowest = TRUE))
  rev_depth_levels_names <- rev(depth_levels_names)

  colnames(overall_held3D) <-
    rownames(total_amount) <-
    names(overall_total_amount) <-
    names_features

  names(solution3D) <-
    rownames(overall_held3D) <-
    names(mean_overall_helds3D) <-
    paste0("budget", budget_percents)

  colnames(total_amount) <-
    names(split_features) <-
    rev_depth_levels_names

  if (length_budget_percents > 1){
    names(sd_overall_helds3D) <-
      rownames(depth_overall_available3D) <-
      rownames(relative_helds3D) <-
      paste0("budget", budget_percents)

    colnames(depth_overall_available3D) <-
      colnames(relative_helds3D) <-
      rev_depth_levels_names

  } else {
    names(depth_overall_available3D) <-
      names(relative_helds3D) <-
      rev_depth_levels_names
  }

  absolute_held3D <- lapply(absolute_held3D,
                            function(x, rev_depth_levels_names, names_features){
                              colnames(x) <- rev_depth_levels_names
                              rownames(x) <- names_features
                              return(x)
                            },
  rev_depth_levels_names=rev_depth_levels_names,
  names_features = names_features)


  output <- list(
    split_features=split_features,
    solution3D=solution3D,
    absolute_held3D=absolute_held3D,
    overall_held3D=overall_held3D,
    relative_helds3D=relative_helds3D,
    mean_overall_helds3D=mean_overall_helds3D,
    sd_overall_helds3D=sd_overall_helds3D,
    depth_overall_available3D = depth_overall_available3D,
    names_features = names_features,
    total_amount = total_amount,
    overall_total_amount = overall_total_amount,
    depth_levels_names=depth_levels_names,
    depth_raster=depth_raster,
    breaks=breaks,
    biodiv_df=biodiv_df,
    priority_weights=priority_weights,
    budget_percents=budget_percents,
    budget_weights=budget_weights,
    penalty = penalty,
    edge_factor = edge_factor,
    gap=gap,
    threads=threads,
    sep_priority_weights=sep_priority_weights,
    portfolio=portfolio,
    portfolio_opts=portfolio_opts,
    sep_biodiv_df=sep_biodiv_df,
    locked_in_raster=locked_in_raster,
    locked_out_raster=locked_out_raster
  )
  attr(output, "from.function") <- "prioritize_3D"
  return(output)
}


coherence <- function(x, w, metric="sa", normalize = TRUE, plot = TRUE,
                      addlines = TRUE, ...){
  if (!is.null(attr(x, "from.function"))){
    if (attr(x, "from.function") == "Compare_2D_3D"){
      raster_data_3D <- sumrast(x$solution3D)
      raster_data_2D <- sumrast(x$solution2D)

      biodiv_raster_copy <- x$biodiv_raster * 1
      biodiv_df_copy <- x$biodiv_df
      names(biodiv_raster_copy) <- tolower(names(biodiv_raster_copy))
      biodiv_df_copy$species_name <- tolower(biodiv_df_copy$species_name)
      ## Remove species for which we have the 2D distribution but not
      ## their depth behavior
      biodiv_raster_copy <-biodiv_raster_copy[[names(biodiv_raster_copy) %in%
                                                 biodiv_df_copy$species_name]]
      raster_data_3D <- raster_data_3D * sum(biodiv_raster_copy, na.rm=TRUE)
      raster_data_2D <- raster_data_2D * sum(biodiv_raster_copy, na.rm=TRUE)
      if (metric == "sa"){
        w = matrix(1, nrow = w, ncol =w)

        # SA 2D
        ff_raster_2D_sa <- geodiv::focal_metrics(raster_data_2D, window=w,
                                                 metrics=list("sa"),
                                                 progress = TRUE, ...)
        val2d <- values(ff_raster_2D_sa$sa)
        val2d[!is.finite(val2d)] <- NA
        mean_val2d <- mean(val2d, na.rm=TRUE)

        # SA 3D
        ff_raster_3D_sa <- geodiv::focal_metrics(raster_data_3D, window=w,
                                                 metrics=list("sa"),
                                                 progress = TRUE, ...)
        val3d <- values(ff_raster_3D_sa$sa)
        val3d[!is.finite(val3d)] <- NA
        mean_val3d <- mean(val3d, na.rm=TRUE)

        if (plot){
          layout(mat = matrix(c(1,2),nrow = 1,ncol = 2,byrow = TRUE))
          terra::plot(ff_raster_2D_sa$sa,
                      main=paste0("2D SA=", round(mean_val2d, 3),
                                  " (window=", nrow(w), "x", ncol(w), ")"),
                      col = rev(terrain.colors(255)))
          if (addlines) maps::map(add=TRUE)
          terra::plot(ff_raster_3D_sa$sa,
                      main=paste0("3D SA=", round(mean_val3d, 3),
                                  " (window=", nrow(w), "x", ncol(w), ")"),
                      col = rev(terrain.colors(255)))
          if (addlines) maps::map(add=TRUE)
          layout(mat = matrix(1))
        }

        return(
          c("sa2D"=round(geodiv::sa(raster_data_2D), 3),
            "sa3D"=round(geodiv::sa(raster_data_3D), 3),
            "sa2Dw"=round(mean_val2d, 3),
            "sa3Dw"=round(mean_val3d, 3)
          ))

      } else if (metric == "sku"){

        w = matrix(1, nrow = w, ncol =w)

        # SKU 2D
        ff_raster_2D_sku <- geodiv::focal_metrics(raster_data_2D, window=w,
                                                  metrics=list("sku"),
                                                  progress = TRUE, ...)
        val2d <- values(ff_raster_2D_sku$sku)
        val2d[!is.finite(val2d)] <- NA
        mean_val2d <- mean(val2d, na.rm=TRUE)

        # SKU 3D
        ff_raster_3D_sku <- geodiv::focal_metrics(raster_data_3D, window=w,
                                                  metrics=list("sku"),
                                                  progress = TRUE, ...)
        val3d <- values(ff_raster_3D_sku$sku)
        val3d[!is.finite(val3d)] <- NA
        mean_val3d <- mean(val3d, na.rm=TRUE)

        if (plot){
          layout(mat = matrix(c(1,2),nrow = 1,ncol = 2,byrow = TRUE))
          terra::plot(ff_raster_2D_sku$sku,
                      main=paste0("2D SKU=", round(mean_val2d, 3),
                                  " (window=", nrow(w), "x", ncol(w), ")"),
                      col = rev(terrain.colors(255)))
          maps::map(add=TRUE)
          terra::plot(ff_raster_3D_sku$sku,
                      main=paste0("3D SKU=", round(mean_val3d, 3),
                                  " (window=", nrow(w), "x", ncol(w), ")"),
                      col = rev(terrain.colors(255)))
          maps::map(add=TRUE)
          layout(mat = matrix(1))
        }

        return(
          c("sku2D"=round(geodiv::sku(raster_data_2D), 3),
            "sku3D"=round(geodiv::sku(raster_data_3D), 3),
            "sku2Dw"=round(mean(val2d, na.rm=TRUE), 3),
            "sku3Dw"=round(mean(val3d, na.rm=TRUE), 3)
          ))

      } else if (metric == "rao"){

        message("2D RAO")
        RAO_2D <- rasterdiv::paRao(x=raster_data_2D,window=w,
                                   na.tolerance=1.0,method="classic", ...)
        mean_2D <- global(RAO_2D[[1]]$alpha.1,"mean",na.rm=T)

        message("3D RAO")
        RAO_3D <- rasterdiv::paRao(x=raster_data_3D,window=w,
                                   na.tolerance=1.0,method="classic", ...)
        mean_3D <- global(RAO_3D[[1]]$alpha.1,"mean",na.rm=T)

        if (plot){
          layout(mat = matrix(c(1,2),nrow = 1,ncol = 2,byrow = TRUE))
          terra::plot(RAO_2D[[1]]$alpha.1,
                      main=paste0("2D RAO=",round(mean_2D,3),
                                  " (window=", w, "x", w, ")"),
                      col = rev(terrain.colors(255)))
          maps::map(add=TRUE)
          terra::plot(RAO_3D[[1]]$alpha.1,
                      main=paste0("3D RAO=",round(mean_3D,3),
                                  " (window=", w, "x", w, ")"),
                      col = rev(terrain.colors(255)))
          maps::map(add=TRUE)
          layout(mat = matrix(1))
        }

        return(
          c("rao2D"=round(mean_2D,3),
            "rao3D"=round(mean_3D,3)
          ))
      }
    }
  }
}
