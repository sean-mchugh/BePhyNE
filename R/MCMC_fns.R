


#####Estimator###############################################################################################################################################################################

Write_to_File = function (dir, outname, ID, iteration, current_vals, append = F)
{
  files <- list(file(file.path(dir, paste(outname, ".", ID,
                                          ".pars.log", sep = "")), open = "a"), file(file.path(dir,
                                                                                               paste(outname, ".", ID, ".log", sep = "")), open = "a"))

  iteration
  if (append == F) {
    names <- lapply(1:length(current_vals[[1]][[1]]), function(x) names(unlist(current_vals[[1]][[1]][[x]])))
    header <- unlist(lapply(1:length(names), function(trait) lapply(1:length(names[[1]]),
                                                                    function(x) paste("pred_", trait, "_", names[[trait]][[x]],
                                                                                      sep = ""))))
    header = c("Iteration", header)
    cat(header, sep = "\t", file = files[[1]], append = TRUE)
    cat("\n", file = files[[1]], append = TRUE)
    header <- names(unlist(current_vals[[1]][[2]]))
    cat(header, sep = "\t", file = files[[2]], append = TRUE)
    cat("\n", file = files[[2]], append = TRUE)
  }
  cat(c(iteration,unname(unlist(current_vals[[1]][[1]]))), sep = "\t",
      file = files[[1]], append = TRUE)
  cat(" ", sep = "", file = files[[1]], append = TRUE)
  cat("\n", file = files[[1]], append = TRUE)
  cat(c(iteration,unname(unlist(current_vals[[1]][[2]]))), sep = "\t",
      file = files[[2]], append = TRUE)
  cat(" ", sep = "", file = files[[2]], append = TRUE)
  cat("\n", file = files[[2]], append = TRUE)
}


metro_haste_full_MV = function (R_corr_start, R_sd_start, A_start, Prior, tree, tibble_data,
                                pa_data, iterations, burnin, move_prob = c(2/5, 2/5, 1/15,
                                                                           1/15, 1/15), n, print.i.freq = 100, print.ac.freq = 10,
                                printing = TRUE, trim = T, trim_freq = 1, H_fixed = T, tuning,
                                k, all_props = F, center_fixed = F, write_file = T, IDlen = 5,
                                dir, outname, prior_only = F, glm_only = F, plot = F, plot_freq,
                                plot_file, True_pars = NULL)
{
  {
    dat_ft = lapply(1:length(tibble_data), function(pred) t(apply(tibble_data[[pred]],
                                                                  1, forwardTransform1)))
    d <- tuning$niche_prop
    w_mu <- tuning$w_mu
    w_sd <- tuning$w_sd
    v_cor <- tuning$v_cor
    for (pred in 1:length(tibble_data)) {
      rownames(dat_ft[[pred]]) <- tree$tip.label
    }
    tibble_data_ft <- lapply(1:length(tibble_data), function(x) make.treedata(tree,
                                                                              dat_ft[[x]])$dat)
    A_start_ft <- lapply(1:length(tibble_data), function(i) forwardTransform1(A_start[[i]])[1:length(A_start[[i]])])
    R <- lapply(1:length(R_sd_start), function(i) rebuild.cov(r = R_corr_start[[i]],
                                                              v = (R_sd_start[[i]])^2))
    starting.jacobian <- lapply(1:length(R_sd_start), function(i) sum(sapply(1:length(R_sd_start[i]),
                                                                             function(x) log(((R_sd_start[[i]][x])^2)))) * log((length(R_sd_start[[i]]) -
                                                                                                                                  1)/2))
    current_vals = list()
    current_vals[[1]] = list()
    current_vals[[1]][[1]] = list()
    current_vals[[1]][[2]] = list()
    for (i in 1:length(tibble_data)) {
      current_vals[[1]][[1]][[i]] = list()
      current_vals[[1]][[1]][[i]]$R = R[[i]]
      current_vals[[1]][[1]][[i]]$R_cor = R_corr_start[[i]]
      current_vals[[1]][[1]][[i]]$R_sd = R_sd_start[[i]]
      current_vals[[1]][[1]][[i]]$A = A_start_ft[[i]]
      current_vals[[1]][[1]][[i]]$dat = tibble_data_ft[[i]]
      current_vals[[1]][[1]][[i]]$curr.jac = starting.jacobian[[i]]
    }
    current_vals[[1]][[2]]$move = NULL
    current_vals[[1]][[2]]$accept = NULL
    current_vals[[1]][[2]]$a.ratio = NULL
    current_vals[[1]][[2]]$lik_sum = NULL
    current_vals[[1]][[2]]$cur.liks = NULL
    current_vals[[1]][[2]]$prop_accept_ratios = NULL
    prior <- Prior
    dt = lapply(1:length(current_vals[[1]][[1]]), function(x) data.frame(species = tree$tip.label,
                                                                         as.data.frame(current_vals[[1]][[1]][[x]]$dat)))
    td = lapply(1:length(dt), function(x) make.treedata(tree,
                                                        dt[[x]]))
    dat.cur = lapply(1:length(current_vals[[1]][[1]]), function(x) current_vals[[1]][[1]][[x]]$dat)
    dat.cur.bt = lapply(1:length(dat.cur), function(x) t(apply(dat.cur[[x]],
                                                               1, backTransform1)))
    A.cur = lapply(1:length(current_vals[[1]][[1]]), function(x) current_vals[[1]][[1]][[x]]$A)
    A.cur.bt = lapply(1:length(A.cur), function(x) backTransform1(A.cur[[x]])[1:length(A.cur[[i]])])
    cur.glm.lik = GLM_Likelihood_MV_fn(dat.cur.bt, pa_data)
    cur.glm.sum.lik = sum(unlist(cur.glm.lik))
    cur.bm.lik = sum(unlist(lapply(1:length(dat.cur), function(x) lnL_ratematrix(dat.cur[[x]],
                                                                                 tree, current_vals[[1]][[1]][[x]]$A, current_vals[[1]][[1]][[x]]$R))))
    prior.root.lik = sum(unlist(lapply(1:length(A.cur), function(x) prior[[x]]$mean.prior(A.cur.bt[[x]]))))
    prior.corr.lik = sum(unlist(lapply(1:length(current_vals[[1]][[1]]),
                                       function(x) prior[[x]]$corr.prior(current_vals[[1]][[1]][[x]]$R_cor))))
    prior.sd.lik = sum(unlist(lapply(1:length(current_vals[[1]][[1]]),
                                     function(x) prior[[x]]$sd.prior(current_vals[[1]][[1]][[x]]$R_sd))))
    prior.heights.lik = lapply(1:length(current_vals[[1]][[1]]),
                               function(pred) {
                                 lapply(1:length(pa_data), function(sp) {
                                   if (is.null(Prior[[pred]]$heights.prior) ==
                                       F) {
                                     Prior[[pred]]$heights.prior[[sp]](dat.cur.bt[[pred]][sp,
                                                                                          3])
                                   }
                                   else {
                                     0
                                   }
                                 })
                               })
    current_vals[[1]][[2]]$cur.liks = list(glm.lik = cur.glm.lik,
                                           glm.sum.lik = cur.glm.sum.lik, bm.lik = cur.bm.lik,
                                           prior.lik = list(prior.root.lik = prior.root.lik,
                                                            prior.corr.lik = prior.corr.lik, prior.sd.lik = prior.sd.lik,
                                                            prior.heights.lik = prior.heights.lik))
    acceptances = 0
    A.center.bm.moves = 0
    A.center.slide.moves = 0
    A.center.mult.moves = 0
    A.width.bm.moves = 0
    A.width.slide.moves = 0
    A.width.mult.moves = 0
    A.theta.moves = 0
    A.R.corr.moves = 0
    A.R.sd.moves = 0
    A.R.sd.moves = 0
    NA_moves = 0
    NA.center.bm.moves = 0
    NA.center.slide.moves = 0
    NA.center.mult.moves = 0
    NA.width.bm.moves = 0
    NA.width.slide.moves = 0
    NA.width.mult.moves = 0
    NA.theta.moves = 0
    NA.R.corr.moves = 0
    NA.R.sd.moves = 0
    Rej.center.bm.moves = 0
    Rej.center.slide.moves = 0
    Rej.center.mult.moves = 0
    Rej.width.bm.moves = 0
    Rej.width.slide.moves = 0
    Rej.width.mult.moves = 0
    Rej.theta.moves = 0
    Rej.R.corr.moves = 0
    Rej.R.sd.moves = 0
    chain = list()
    trimmed_chain = list()
    td = lapply(1:length(dt), function(x) make.treedata(tree,
                                                        dt[[x]]))
  }
  if (all_props == T) {
    props <- list()
  }
  ID <- paste(sample(x = 1:9, size = IDlen, replace = TRUE),
              collapse = "")
  if (write_file == T) {
    Write_to_File(dir = dir, outname = outname, ID = ID, i,
                  current_vals = current_vals)
  }
  for (i in 1:iterations) {
    k = k
    for (x in (1:length(dt))) {
      td[[x]]$dat <- current_vals[[1]][[1]][[x]]$dat
    }
    prop <- sample(c(0, 1, 2, 3, 4, 5), size = 1, prob = move_prob)
    if (prop < 3) {
      node = sample(1:(length(tree$tip.label) + tree$Nnode),
                    1)
      if (node <= length(tree$tip.label)) {
        clade <- list(tip.label = tree$tip.label[node])
      }
      else {
        clade <- extract.clade(tree, node = node)
      }
      clade
      missing <- (1:nrow(td[[x]]$dat))[tree$tip.label %in%
                                         clade$tip.label]
    }
    if (prop == 0) {
      niche_move = sample(c(1, 2), size = 1, prob = c(1/2,
                                                      1/2))
      if (niche_move == 1) {
        proposal = lapply(1:length(td), function(x) Slide_Proposal_byClade_byTrait_fn(tree,
                                                                                      td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor,
                                                                                      current_vals[[1]][[1]][[x]]$R_sd, current_vals[[1]][[1]][[x]]$A,
                                                                                      n, current_vals[[1]][[1]][[x]]$curr.jac, d[[x]]$slide[missing,
                                                                                                                                            3], niche_move = "height", clade))
        move = "Height_Slide"
      }
      else if (niche_move == 2) {
        proposal = lapply(1:length(td), function(x) Multiplier_Proposal_byClade_byTrait_fn(tree,
                                                                                           td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor,
                                                                                           current_vals[[1]][[1]][[x]]$R_sd, current_vals[[1]][[1]][[x]]$A,
                                                                                           n, current_vals[[1]][[1]][[x]]$curr.jac, d[[x]]$multi[missing,
                                                                                                                                                 3], niche_move = "height", clade))
        move = "Height_Multiplier"
      }
      for (z in 2:length(proposal)) {
        proposal[[z]]$td$dat[, 3] <- proposal[[1]]$td$dat[,
                                                          3]
        proposal[[z]]$hr = 0
      }
    }
    else if (prop == 1) {
      niche_move = sample(c(1, 2, 3), size = 1, prob = c(0,
                                                         1/2, 1/2))
      if (niche_move == 1) {
        move = "Center_BM"
      }
      else if (niche_move == 2) {
        proposal = lapply(1:length(td), function(x) Slide_Proposal_byClade_byTrait_fn(tree,
                                                                                      td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor,
                                                                                      current_vals[[1]][[1]][[x]]$R_sd, current_vals[[1]][[1]][[x]]$A,
                                                                                      n, current_vals[[1]][[1]][[x]]$curr.jac, d[[x]]$slide[missing,
                                                                                                                                            1], niche_move = "center", clade))
        if (center_fixed != FALSE) {
          for (pred in center_fixed) {
            proposal[[pred]]$td <- td[[pred]]
          }
        }
        move = "Center_Slide"
      }
      else if (niche_move == 3) {
        proposal = lapply(1:length(td), function(x) Multiplier_Proposal_byClade_byTrait_fn(tree,
                                                                                           td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor,
                                                                                           current_vals[[1]][[1]][[x]]$R_sd, current_vals[[1]][[1]][[x]]$A,
                                                                                           n, current_vals[[1]][[1]][[x]]$curr.jac, d[[x]]$multi[missing,
                                                                                                                                                 1], niche_move = "center", clade))
        if (center_fixed != FALSE) {
          for (pred in center_fixed) {
            proposal[[pred]]$td <- td[[pred]]
          }
        }
        move = "Center_Multiplier"
      }
    }
    else if (prop == 2) {
      niche_move = sample(c(1, 2, 3), size = 1, prob = c(0,
                                                         1/2, 1/2))
      if (niche_move == 1) {
        move = "Width_BM"
      }
      else if (niche_move == 2) {
        proposal = lapply(1:length(td), function(x) Slide_Proposal_byClade_byTrait_fn(tree,
                                                                                      td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor,
                                                                                      current_vals[[1]][[1]][[x]]$R_sd, current_vals[[1]][[1]][[x]]$A,
                                                                                      n, current_vals[[1]][[1]][[x]]$curr.jac, d[[x]]$slide[missing,
                                                                                                                                            2], niche_move = "width", clade))
        move = "Width_Slide"
      }
      else if (niche_move == 3) {
        proposal = lapply(1:length(td), function(x) Multiplier_Proposal_byClade_byTrait_fn(tree,
                                                                                           td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor,
                                                                                           current_vals[[1]][[1]][[x]]$R_sd, current_vals[[1]][[1]][[x]]$A,
                                                                                           n, current_vals[[1]][[1]][[x]]$curr.jac, d[[x]]$multi[missing,
                                                                                                                                                 2], niche_move = "width", clade))
        move = "Width_Multiplier"
      }
    }
    else if (prop == 3) {
      proposal <- lapply(1:length(td), function(x) slideWindow_theta(w_mu[[x]]$slide,
                                                                     tree, td[[x]], current_vals[[1]][[1]][[x]]$R,
                                                                     current_vals[[1]][[1]][[x]]$R_cor, current_vals[[1]][[1]][[x]]$R_sd,
                                                                     current_vals[[1]][[1]][[x]]$A, n, current_vals[[1]][[1]][[x]]$curr.jac))
      move = "theta"
    }
    else if (prop == 4) {
      proposal <- lapply(1:length(td), function(x) Cor_matrix_prop(tree,
                                                                   td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor,
                                                                   current_vals[[1]][[1]][[x]]$R_sd, current_vals[[1]][[1]][[x]]$A,
                                                                   n, current_vals[[1]][[1]][[x]]$curr.jac, v = v_cor[[x]]))
      move = "R_corr"
      if (is.null(proposal) == T) {
        print("Your matrix is fucked yo, this is the signularity error")
        if (trim == TRUE) {
          return(list(acceptances = acceptances, accepted.moves = list(A.dat.bm.moves = A.dat.bm.moves,
                                                                       A.dat.slide.moves = A.dat.slide.moves, A.dat.mult.moves = A.dat.mult.moves,
                                                                       A.theta.moves = A.theta.moves, A.R.corr.moves = A.R.corr.moves,
                                                                       A.R.sd.moves = A.R.sd.moves), rejected.moves = list(Rej.dat.bm.moves = Rej.dat.bm.moves,
                                                                                                                           Rej.dat.slide.moves = Rej.dat.slide.moves,
                                                                                                                           Rej.dat.mult.moves = Rej.dat.mult.moves,
                                                                                                                           Rej.theta.moves = Rej.theta.moves, Rej.R.corr.moves = Rej.R.corr.moves,
                                                                                                                           Rej.R.sd.moves = Rej.R.sd.moves), trimmed_chain = trimmed_chain,
                      NAs = NA_moves, NA.moves = list(NA.dat.bm.moves = NA.dat.bm.moves,
                                                      NA.dat.slide.moves = NA.dat.slide.moves,
                                                      NA.dat.mult.moves = NA.dat.mult.moves,
                                                      NA.theta.moves = NA.theta.moves, NA.R.corr.moves = NA.R.corr.moves,
                                                      NA.R.sd.moves = NA.R.sd.moves), "full"))
        }
        else {
          return(list(acceptances = acceptances, accepted.moves = list(A.dat.bm.moves = A.dat.bm.moves,
                                                                       A.dat.slide.moves = A.dat.slide.moves, A.dat.mult.moves = A.dat.mult.moves,
                                                                       A.theta.moves = A.theta.moves, A.R.corr.moves = A.R.corr.moves,
                                                                       A.R.sd.moves = A.R.sd.moves), rejected.moves = list(Rej.dat.bm.moves = Rej.dat.bm.moves,
                                                                                                                           Rej.dat.slide.moves = Rej.dat.slide.moves,
                                                                                                                           Rej.dat.mult.moves = Rej.dat.mult.moves,
                                                                                                                           Rej.theta.moves = Rej.theta.moves, Rej.R.corr.moves = Rej.R.corr.moves,
                                                                                                                           Rej.R.sd.moves = Rej.R.sd.moves), chain = chain,
                      NAs = NA_moves, NA.moves = list(A.dat.bm.moves = NA.dat.bm.moves,
                                                      NA.dat.slide.moves = NA.dat.slide.moves,
                                                      NA.dat.mult.moves = NA.dat.mult.moves,
                                                      NA.theta.moves = NA.theta.moves, NA.R.corr.moves = NA.R.corr.moves,
                                                      NA.R.sd.moves = NA.R.sd.moves), "full"))
        }
      }
    }
    else if (prop == 5) {
      proposal <- lapply(1:length(td), function(x) SD_prop(tree,
                                                           td[[x]], current_vals[[1]][[1]][[x]]$R, current_vals[[1]][[1]][[x]]$R_cor,
                                                           current_vals[[1]][[1]][[x]]$R_sd, current_vals[[1]][[1]][[x]]$A,
                                                           n, current_vals[[1]][[1]][[x]]$curr.jac, w_sd = w_sd[[x]]$slide))
      move = "R_sd"
      if (is.null(proposal) == T) {
        print("Your matrix is fucked yo, this is the signularity error")
        if (trim == TRUE) {
          return(list(acceptances = acceptances, accepted.moves = list(A.dat.bm.moves = A.dat.bm.moves,
                                                                       A.dat.slide.moves = A.dat.slide.moves, A.dat.mult.moves = A.dat.mult.moves,
                                                                       A.theta.moves = A.theta.moves, A.R.corr.moves = A.R.corr.moves,
                                                                       A.R.sd.moves = A.R.sd.moves), rejected.moves = list(Rej.dat.bm.moves = Rej.dat.bm.moves,
                                                                                                                           Rej.dat.slide.moves = Rej.dat.slide.moves,
                                                                                                                           Rej.dat.mult.moves = Rej.dat.mult.moves,
                                                                                                                           Rej.theta.moves = Rej.theta.moves, Rej.R.corr.moves = Rej.R.corr.moves,
                                                                                                                           Rej.R.sd.moves = Rej.R.sd.moves), trimmed_chain = trimmed_chain,
                      NAs = NA_moves, NA.moves = list(NA.dat.bm.moves = NA.dat.bm.moves,
                                                      NA.dat.slide.moves = NA.dat.slide.moves,
                                                      NA.dat.mult.moves = NA.dat.mult.moves,
                                                      NA.theta.moves = NA.theta.moves, NA.R.corr.moves = NA.R.corr.moves,
                                                      NA.R.sd.moves = NA.R.sd.moves), "full"))
        }
        else {
          return(list(acceptances = acceptances, accepted.moves = list(A.dat.bm.moves = A.dat.bm.moves,
                                                                       A.dat.slide.moves = A.dat.slide.moves, A.dat.mult.moves = A.dat.mult.moves,
                                                                       A.theta.moves = A.theta.moves, A.R.corr.moves = A.R.corr.moves,
                                                                       A.R.sd.moves = A.R.sd.moves), rejected.moves = list(Rej.dat.bm.moves = Rej.dat.bm.moves,
                                                                                                                           Rej.dat.slide.moves = Rej.dat.slide.moves,
                                                                                                                           Rej.dat.mult.moves = Rej.dat.mult.moves,
                                                                                                                           Rej.theta.moves = Rej.theta.moves, Rej.R.corr.moves = Rej.R.corr.moves,
                                                                                                                           Rej.R.sd.moves = Rej.R.sd.moves), chain = chain,
                      NAs = NA_moves, NA.moves = list(A.dat.bm.moves = NA.dat.bm.moves,
                                                      NA.dat.slide.moves = NA.dat.slide.moves,
                                                      NA.dat.mult.moves = NA.dat.mult.moves,
                                                      NA.theta.moves = NA.theta.moves, NA.R.corr.moves = NA.R.corr.moves,
                                                      NA.R.sd.moves = NA.R.sd.moves), "full"))
        }
      }
    }
    probab = acceptance_ratio_MV_fun(proposal, current_vals,
                                     pa_data, tree, prior, prop, k, prior_only = prior_only,
                                     glm_only = glm_only)
    if (is.na(probab$a.ratio)) {
      cat(red("NA acceptance ratio ", i + 1, " ", current_vals[[1]][[2]]$move,
              "\n"))
      NA_moves = NA_moves + 1
      if (move == "Center_Multiplier") {
        NA.center.mult.moves = NA.center.mult.moves +
          1
      }
      else if (move == "Center_Slide") {
        NA.center.slide.moves = NA.center.slide.moves +
          1
      }
      else if (move == "Center_BM") {
        NA.center.bm.moves = NA.center.bm.moves + 1
      }
      else if (move == "Width_Multiplier") {
        NA.width.mult.moves = NA.width.mult.moves + 1
      }
      else if (move == "Width_Slide") {
        NA.width.slide.moves = NA.width.slide.moves +
          1
      }
      else if (move == "Width_BM") {
        NA.width.bm.moves = NA.width.bm.moves + 1
      }
      else if (move == "theta") {
        NA.theta.moves = NA.theta.moves + 1
      }
      else if (move == "R_corr") {
        NA.R.corr.moves = NA.R.corr.moves + 1
      }
      else if (move == "R_sd") {
        NA.R.sd.moves = NA.R.sd.moves + 1
      }
      probab$a.ratio <- -66666666666666
    }
    if (runif(1) < probab$a.ratio) {
      for (y in 1:length(current_vals[[1]][[1]])) {
        current_vals[[1]][[1]][[y]]$R = proposal[[y]]$R
        current_vals[[1]][[1]][[y]]$R_cor = proposal[[y]]$R_cor
        current_vals[[1]][[1]][[y]]$R_sd = proposal[[y]]$R_sd
        current_vals[[1]][[1]][[y]]$A = proposal[[y]]$A
        current_vals[[1]][[1]][[y]]$dat = proposal[[y]]$td$dat
      }
      current_vals[[1]][[2]]$move = move
      current_vals[[1]][[2]]$accept = TRUE
      current_vals[[1]][[2]]$a.ratio = probab$a.ratio
      current_vals[[1]][[2]]$cur.liks = probab$prop.liks
      current_vals[[1]][[2]]$lik_sum = probab
      current_vals[[1]][[2]]$curr.jac = proposal$prop_jac
      current_vals[[1]][[2]]$prop_accept_ratios = list(center.bm = (A.center.bm.moves/(A.center.bm.moves +
                                                                                         Rej.center.bm.moves)), center.slide = (A.center.slide.moves/(A.center.slide.moves +
                                                                                                                                                        Rej.center.slide.moves)), center.mult = (A.center.mult.moves/(A.center.mult.moves +
                                                                                                                                                                                                                        Rej.center.mult.moves)), width.bm = (A.width.bm.moves/(A.width.bm.moves +
                                                                                                                                                                                                                                                                                 Rej.width.bm.moves)), width.slide = (A.width.slide.moves/(A.width.slide.moves +
                                                                                                                                                                                                                                                                                                                                             Rej.width.slide.moves)), width.mult = (A.width.mult.moves/(A.width.mult.moves +
                                                                                                                                                                                                                                                                                                                                                                                                          Rej.width.mult.moves)), theta.slide = (A.theta.moves/(A.theta.moves +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                  Rej.theta.moves)), R.corr.slide = (A.R.corr.moves/(A.R.corr.moves +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       Rej.R.corr.moves)), R.sd.slide = (A.R.sd.moves/(A.R.sd.moves +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         Rej.R.sd.moves)))
      acceptances = acceptances + 1
      if (move == "Center_Multiplier") {
        A.center.mult.moves = A.center.mult.moves + 1
      }
      else if (move == "Center_Slide") {
        A.center.slide.moves = A.center.slide.moves +
          1
      }
      else if (move == "Center_BM") {
        A.center.bm.moves = A.center.bm.moves + 1
      }
      else if (move == "Width_Multiplier") {
        A.width.mult.moves = A.width.mult.moves + 1
      }
      else if (move == "Width_Slide") {
        A.width.slide.moves = A.width.slide.moves + 1
      }
      else if (move == "Width_BM") {
        A.width.bm.moves = A.width.bm.moves + 1
      }
      else if (move == "theta") {
        A.theta.moves = A.theta.moves + 1
      }
      else if (move == "R_corr") {
        A.R.corr.moves = A.R.corr.moves + 1
      }
      else if (move == "R_sd") {
        A.R.sd.moves = A.R.sd.moves + 1
      }
      if (acceptances%%print.ac.freq == 0 & printing ==
          TRUE) {
        cat("iterations:", i + 1, "\n")
        cat("acceptances:", acceptances, "\n")
        cat("acceptance ratio:", probab$a.ratio, "\n")
        cat("move utilized:", move, "\n")
      }
    }
    else {
      current_vals[[1]][[2]]$accept = FALSE
      current_vals[[1]][[2]]$a.ratio = probab$a.ratio
      current_vals[[1]][[2]]$lik_sum = probab
      current_vals[[1]][[2]]$move = move
      current_vals[[1]][[2]]$prop_accept_ratios = list(center.bm = (A.center.bm.moves/(A.center.bm.moves +
                                                                                         Rej.center.bm.moves)), center.slide = (A.center.slide.moves/(A.center.slide.moves +
                                                                                                                                                        Rej.center.slide.moves)), center.mult = (A.center.mult.moves/(A.center.mult.moves +
                                                                                                                                                                                                                        Rej.center.mult.moves)), width.bm = (A.width.bm.moves/(A.width.bm.moves +
                                                                                                                                                                                                                                                                                 Rej.width.bm.moves)), width.slide = (A.width.slide.moves/(A.width.slide.moves +
                                                                                                                                                                                                                                                                                                                                             Rej.width.slide.moves)), width.mult = (A.width.mult.moves/(A.width.mult.moves +
                                                                                                                                                                                                                                                                                                                                                                                                          Rej.width.mult.moves)), theta.slide = (A.theta.moves/(A.theta.moves +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                  Rej.theta.moves)), R.corr.slide = (A.R.corr.moves/(A.R.corr.moves +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       Rej.R.corr.moves)), R.sd.slide = (A.R.sd.moves/(A.R.sd.moves +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         Rej.R.sd.moves)))
      if (move == "Center_Multiplier") {
        Rej.center.mult.moves = Rej.center.mult.moves +
          1
      }
      else if (move == "Center_Slide") {
        Rej.center.slide.moves = Rej.center.slide.moves +
          1
      }
      else if (move == "Center_BM") {
        Rej.center.bm.moves = Rej.center.bm.moves + 1
      }
      else if (move == "Width_Multiplier") {
        Rej.width.mult.moves = Rej.width.mult.moves +
          1
      }
      else if (move == "Width_Slide") {
        Rej.width.slide.moves = Rej.width.slide.moves +
          1
      }
      else if (move == "Width_BM") {
        Rej.width.bm.moves = Rej.width.bm.moves + 1
      }
      else if (move == "theta") {
        Rej.theta.moves = Rej.theta.moves + 1
      }
      else if (move == "R_corr") {
        Rej.R.corr.moves = Rej.R.corr.moves + 1
      }
      else if (move == "R_sd") {
        Rej.R.sd.moves = Rej.R.sd.moves + 1
      }
      if (i%%print.i.freq == 0 & printing == TRUE) {
        cat("iterations:", i + 1, "\n")
        cat("acceptances:", acceptances, "\n")
        cat("acceptance ratio (not accepted):", probab$a.ratio,
            "\n")
        cat("move utilized:", move, "\n")
      }
    }
    if (i > burnin) {
      if (write_file == T) {
        if (trim == TRUE) {
          if (i%%trim_freq == 0) {
            Write_to_File(dir = dir, outname = outname,
                          ID = ID, i,current_vals = current_vals, append = T)
          }
        }
        else {
          Write_to_File(dir = dir, outname = outname,
                        ID = ID, i, current_vals = current_vals, append = T)
        }
      }
      else if (trim == TRUE) {
        if (i%%trim_freq == 0) {
          trimmed_chain = append(trimmed_chain, current_vals[1])
          if (plot == T & i%%plot_freq == 0) {
            print(length(trimmed_chain))
            pdf(file = paste(plot_file))
            if (is.null(True_pars) == T) {
              plotTraces(trimmed_chain = trimmed_chain,
                         plot.true = F)
            }
            else {
              plotTraces(trimmed_chain = trimmed_chain,
                         plot.true = T, True_pars = True_pars)
            }
            dev.off()
          }
        }
      }
      else {
        chain = append(chain, current_vals[1])
      }
    }
    if (all_props == T) {
      props[[i]] <- list(proposal, move)
    }
  }
  if (all_props == T) {
    return(list(accept_ratios = list(center.bm = (A.center.bm.moves/(A.center.bm.moves +
                                                                       Rej.center.bm.moves)), center.slide = (A.center.slide.moves/(A.center.slide.moves +
                                                                                                                                      Rej.center.slide.moves)), center.mult = (A.center.mult.moves/(A.center.mult.moves +
                                                                                                                                                                                                      Rej.center.mult.moves)), width.bm = (A.width.bm.moves/(A.width.bm.moves +
                                                                                                                                                                                                                                                               Rej.width.bm.moves)), width.slide = (A.width.slide.moves/(A.width.slide.moves +
                                                                                                                                                                                                                                                                                                                           Rej.width.slide.moves)), width.mult = (A.width.mult.moves/(A.width.mult.moves +
                                                                                                                                                                                                                                                                                                                                                                                        Rej.width.mult.moves)), theta.slide = (A.theta.moves/(A.theta.moves +
                                                                                                                                                                                                                                                                                                                                                                                                                                                Rej.theta.moves)), R.corr.slide = (A.R.corr.moves/(A.R.corr.moves +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     Rej.R.corr.moves)), R.sd.slide = (A.R.sd.moves/(A.R.sd.moves +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       Rej.R.sd.moves))), accepted.moves = list(A.center.bm.moves = A.center.bm.moves,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                A.center.slide.moves = A.center.slide.moves, A.center.mult.moves = A.center.mult.moves,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                A.width.bm.moves = A.width.bm.moves, A.width.slide.moves = A.width.slide.moves,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                A.width.mult.moves = A.width.mult.moves, A.theta.moves = A.theta.moves,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                A.R.corr.moves = A.R.corr.moves, A.R.sd.moves = A.R.sd.moves),
                rejected.moves = list(Rej.center.bm.moves = Rej.center.bm.moves,
                                      Rej.center.slide.moves = Rej.center.slide.moves,
                                      Rej.center.mult.moves = Rej.center.mult.moves,
                                      Rej.width.bm.moves = Rej.width.bm.moves, Rej.width.slide.moves = Rej.width.slide.moves,
                                      Rej.width.mult.moves = Rej.width.mult.moves,
                                      Rej.theta.moves = Rej.theta.moves, Rej.R.corr.moves = Rej.R.corr.moves,
                                      Rej.R.sd.moves = Rej.R.sd.moves), chain = trimmed_chain,
                NAs = NA_moves, NA.moves = list(NA.center.bm.moves = NA.center.bm.moves,
                                                NA.center.slide.moves = NA.center.slide.moves,
                                                NA.center.mult.moves = NA.center.mult.moves,
                                                NA.width.bm.moves = NA.width.bm.moves, NA.width.slide.moves = NA.width.slide.moves,
                                                NA.width.mult.moves = NA.width.mult.moves, NA.theta.moves = NA.theta.moves,
                                                NA.R.corr.moves = NA.R.corr.moves, NA.R.sd.moves = NA.R.sd.moves),
                all_props = props, "full"))
  }
  if (trim == TRUE) {
    return(list(acceptances = acceptances, accept_ratios = list(center.bm = (A.center.bm.moves/(A.center.bm.moves +
                                                                                                  Rej.center.bm.moves)), center.slide = (A.center.slide.moves/(A.center.slide.moves +
                                                                                                                                                                 Rej.center.slide.moves)), center.mult = (A.center.mult.moves/(A.center.mult.moves +
                                                                                                                                                                                                                                 Rej.center.mult.moves)), width.bm = (A.width.bm.moves/(A.width.bm.moves +
                                                                                                                                                                                                                                                                                          Rej.width.bm.moves)), width.slide = (A.width.slide.moves/(A.width.slide.moves +
                                                                                                                                                                                                                                                                                                                                                      Rej.width.slide.moves)), width.mult = (A.width.mult.moves/(A.width.mult.moves +
                                                                                                                                                                                                                                                                                                                                                                                                                   Rej.width.mult.moves)), theta.slide = (A.theta.moves/(A.theta.moves +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                           Rej.theta.moves)), R.corr.slide = (A.R.corr.moves/(A.R.corr.moves +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                Rej.R.corr.moves)), R.sd.slide = (A.R.sd.moves/(A.R.sd.moves +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  Rej.R.sd.moves))), accepted.moves = list(A.center.bm.moves = A.center.bm.moves,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           A.center.slide.moves = A.center.slide.moves, A.center.mult.moves = A.center.mult.moves,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           A.width.bm.moves = A.width.bm.moves, A.width.slide.moves = A.width.slide.moves,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           A.width.mult.moves = A.width.mult.moves, A.theta.moves = A.theta.moves,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           A.R.corr.moves = A.R.corr.moves, A.R.sd.moves = A.R.sd.moves),
                rejected.moves = list(Rej.center.bm.moves = Rej.center.bm.moves,
                                      Rej.center.slide.moves = Rej.center.slide.moves,
                                      Rej.center.mult.moves = Rej.center.mult.moves,
                                      Rej.width.bm.moves = Rej.width.bm.moves, Rej.width.slide.moves = Rej.width.slide.moves,
                                      Rej.width.mult.moves = Rej.width.mult.moves,
                                      Rej.theta.moves = Rej.theta.moves, Rej.R.corr.moves = Rej.R.corr.moves,
                                      Rej.R.sd.moves = Rej.R.sd.moves), chain = trimmed_chain,
                NAs = NA_moves, NA.moves = list(NA.center.bm.moves = NA.center.bm.moves,
                                                NA.center.slide.moves = NA.center.slide.moves,
                                                NA.center.mult.moves = NA.center.mult.moves,
                                                NA.width.bm.moves = NA.width.bm.moves, NA.width.slide.moves = NA.width.slide.moves,
                                                NA.width.mult.moves = NA.width.mult.moves, NA.theta.moves = NA.theta.moves,
                                                NA.R.corr.moves = NA.R.corr.moves, NA.R.sd.moves = NA.R.sd.moves),
                "full"))
  }
  else {
    return(list(acceptances = acceptances, accept_ratios = list(center.bm = (A.center.bm.moves/(A.center.bm.moves +
                                                                                                  Rej.center.bm.moves)), center.slide = (A.center.slide.moves/(A.center.slide.moves +
                                                                                                                                                                 Rej.center.slide.moves)), center.mult = (A.center.mult.moves/(A.center.mult.moves +
                                                                                                                                                                                                                                 Rej.center.mult.moves)), width.bm = (A.width.bm.moves/(A.width.bm.moves +
                                                                                                                                                                                                                                                                                          Rej.width.bm.moves)), width.slide = (A.width.slide.moves/(A.width.slide.moves +
                                                                                                                                                                                                                                                                                                                                                      Rej.width.slide.moves)), width.mult = (A.width.mult.moves/(A.width.mult.moves +
                                                                                                                                                                                                                                                                                                                                                                                                                   Rej.width.mult.moves)), theta.slide = (A.theta.moves/(A.theta.moves +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                           Rej.theta.moves)), R.corr.slide = (A.R.corr.moves/(A.R.corr.moves +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                Rej.R.corr.moves)), R.sd.slide = (A.R.sd.moves/(A.R.sd.moves +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  Rej.R.sd.moves))), accepted.moves = list(A.center.bm.moves = A.center.bm.moves,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           A.center.slide.moves = A.center.slide.moves, A.center.mult.moves = A.center.mult.moves,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           A.width.bm.moves = A.width.bm.moves, A.width.slide.moves = A.width.slide.moves,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           A.width.mult.moves = A.width.mult.moves, A.theta.moves = A.theta.moves,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           A.R.corr.moves = A.R.corr.moves, A.R.sd.moves = A.R.sd.moves),
                rejected.moves = list(Rej.center.bm.moves = Rej.center.bm.moves,
                                      Rej.center.slide.moves = Rej.center.slide.moves,
                                      Rej.center.mult.moves = Rej.center.mult.moves,
                                      Rej.width.bm.moves = Rej.width.bm.moves, Rej.width.slide.moves = Rej.width.slide.moves,
                                      Rej.width.mult.moves = Rej.width.mult.moves,
                                      Rej.theta.moves = Rej.theta.moves, Rej.R.corr.moves = Rej.R.corr.moves,
                                      Rej.R.sd.moves = Rej.R.sd.moves), chain = chain,
                NAs = NA_moves, NA.moves = list(NA.center.bm.moves = NA.center.bm.moves,
                                                NA.center.slide.moves = NA.center.slide.moves,
                                                NA.center.mult.moves = NA.center.mult.moves,
                                                NA.width.bm.moves = NA.width.bm.moves, NA.width.slide.moves = NA.width.slide.moves,
                                                NA.width.mult.moves = NA.width.mult.moves, NA.theta.moves = NA.theta.moves,
                                                NA.R.corr.moves = NA.R.corr.moves, NA.R.sd.moves = NA.R.sd.moves),
                "full"))
  }
}



