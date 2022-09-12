#' @export
sim_then_infer = 
  function(n_rep, dgf, n_par, sim_n){   
    
    # SIMULATE
    inf        <- simulate_multi(n_par = n_par, n_sim = 100000, n_rep = n_rep, dgf  = dgf)
    samps_flat <- pluck(inf, "fit_flat") %>% posterior_samples("b_id.*") %>% set_names(1:n_par)
    samps_hier <- pluck(inf, "fit_hier") %>% posterior_samples("r_id.*") %>% set_names(1:n_par)
    
    # INFER
    # browser()
    res  = 
      tibble(sim_n = sim_n,
             n_rep = n_rep, 
             n_par = n_par,
             dgf   = dgf,
             tru   = list(pluck(inf, "tru")),
             samps = list(samps_flat, samps_hier),
             prior = c("flat", "hier"),
             labl  = str_c( "dgf" , dgf, "n_rep", n_rep, "prior", prior, sep = "-"),
             res   = pmap(list(samps, labl, tru, n_par, sim_n), my_sum),
             tabs  = map(res, "tabs"), 
             fco   = map_depth(tabs, 2, get_fco))  %>% 
      mutate(algor = map(fco, names)) %>% 
      unnest(tabs, fco, res, algor, .drop = F) %>%
      mutate(fco  = unlist(fco)) %>%
      # mutate(viol = pmap(list(fco, tru), violated)) %>% 
      # # mutate(viol = unlist(viol)) %>% 
      # mutate(fco  = unlist(fco),
      #        viol  = map2_lgl(fco, tru, violated), # redundant I think
      #        qual  = (str_count(fco, "\\|"))/(n_par-1),
      #        fco_universe = fco == "universe", 
      #        tru_family_null = tru == "family_null") %>% 
      select(tru, prior, tabs, fco, algor, viol, qual, fco_universe, tru_family_null)
  } 

#' @export
my_sum = function(samps, labl, tru, n_par, sim_n){
  
  ########################################################
  # RANGE AND EXHAUSTIVE (ONLY FOR LOW DIMENSIONAL PROBLEMS)
  ########################################################
  
  if(n_par <= 3){
    out_all   = get_discrete_posterior(samps, type = "all")
    out_range = get_discrete_posterior(samps, type = "range", n_terms = c(3))
  } else {
    out_all = NULL
    out_range = NULL
  }
  
  ########################################################
  # RECURSION
  ########################################################
  
  out_greedy = get_discrete_posterior(samps, type = "greedy")
  
  if(perty <- T) get_plot(samps, out_greedy, labl = labl, tru = tru, n_par = n_par, sim_n = sim_n)
  
  res            = NULL
  res$tabs       = out_greedy
  res$out_greedy = out_greedy
  res$out_range  = out_range
  res$out_all    = out_all
  return(res = res)
}



#' @export
summarize_mass = function(mass){
  flt_mass <- unlist(flatten(mass)) # view FCO
  
  if(verbose <- F) print(my_summary <- list(credible95 = flt_mass[which(flt_mass >= 0.95)], 
                                            credible90 = flt_mass[which(flt_mass >= 0.90)], 
                                            credible80 = flt_mass[which(flt_mass >= 0.80)], 
                                            credible70 = flt_mass[which(flt_mass >= 0.70)]))
  return(flt_mass)
}

#' @export
get_discrete_posterior = 
  function(V, keep_names = FALSE, type = "greedy", pre = list(list("V2", "V1", "V3"), list("V2", c("V1", "V3"))), domain = "tsx", n_terms = NULL){ 
    # V is a tidy dataframe of samples
    # keepnames = F then replace column names with natural number indices
    # type controls which hypotheses to search. all (prespecified, exhaustive)
    # pre is the specific, a priori hypotheses (only revelant if type == "pre")
    # it is a list of lists of ordered subsets i.e. weak orders
    # domain = the search space for greedy maximization, takes "tsx" or "wux"
    # if type == "range", then n_terms is a vector of integers giving the cardinality of the subset of orders
    
    if(keep_names == FALSE) {
      old_names = colnames(V)
      new_names = 1:length(V)
      # new_names = abbreviate(colnames(V), 1)
      V         = `colnames<-`(V, new_names)
      if(verbose <- F) print(code_book <- tibble(old_names, new_names))
    } 
    
    if(!(type %in% c("pre", "all", "greedy", "range"))) {
      print("The argument `type` should equal one of: `all`, ... ")
      break
    }
    if(!(domain %in% c("wux", "tsx"))) {
      print("The argument `type` should equal one of:`wux`, `tsx`, ... ")
      break
    }
    if(!is.numeric(n_terms) & !is.null(n_terms)) {
      print("The argument `n_terms` should be a vector of integers ")
    }
    
    # browser()
    if(type == "all"){
      
      if(compositions_of_all_subsets <- F){
        filtered =  get_predicates(V, keep_names = keep_names, type = type, pre = pre, n_terms =  2:dim(V)[2])  
      } else {
        # just compositions of whole integer
        filtered =  get_predicates(V, keep_names = keep_names, type = type, pre = pre, n_terms =  dim(V)[2]) 
      }
    } else if(type == "range"){
      
      filtered =  get_predicates(V, keep_names = keep_names, type = type, pre = pre, n_terms = n_terms)  
      
    } else if(type == "pre"){
      
      filtered =  get_predicates(V, keep_names = keep_names, type = type, pre = pre) 
      
    } else if(type == "greedy"){
      
      med = order(matrixStats::colMedians(as.matrix(V))) # initial value for algorithm
      if(keep_names) med = names(V)[med]
      
      # INITIALIZE AT THE MEDIAN, WHICH 
      init   = get_discrete_posterior(V, keep_names = keep_names, type = "pre", pre = list(as.list(med)))$mass
      if(init[[1]]==0) print("no points in the order of the median")
      
      
      # POSSIBLE  FURTHER GROPING AROUND FOR A GOOD ORDER TO INITIALIZE
      dists   <- fields::rdist(V, t(as.matrix(matrixStats::colMeans2(as.matrix(V))))  ) %>% unlist
      special <- which(dists == min(dists))
      special <- special[1] # just in case there are multiple with exactly the same distance!
      special <- order(unlist(V[special, ])) 
      special <- names(V)[special]
      init    <- get_discrete_posterior(V, keep_names = keep_names, type = "pre", pre = list(as.list(as.character(special))))$mass
      init_bf <- init[[1]]/(1/factorial(dim(V)[2]))
      names(init_bf) = names(init)
      
      # MAIN RECURSION
      # browser()
      out = NULL
      out$out_wux  = list(list(max = init, max_bf = init_bf), get_max(V, init, domain = "wux", keep_names = keep_names)) # RECURSION THROUGH WEAKER ORDERS OF THE SAME (UNIVERSAL) SET
      # browser()
      out$out_tsx  = list(list(max = init, max_bf = init_bf), get_max(V, init, domain = "tsx", keep_names = keep_names)) # RECURSION THROUGH TOTAL ORDERS OF INCREASINGLY SMALL, PROPER SUBSETS
      out$out_wux  = out$out_wux %>% unlist %>% enframe %>% separate(name, c("what", "hyp"), sep = "\\.")
      out$out_tsx  = out$out_tsx %>% unlist %>% enframe %>% separate(name, c("what", "hyp"), sep = "\\.")
      
      return(out = out)
      break
      
    }
    
    filtered = map(filtered, ~V[.x, ])
    if(0) check_partition(filtered, V) # for interest, no side effects
    mass     = map_depth(filtered, 1, ~dim(.x)[1]/dim(V)[1]) # relative to all points
    
    return(out = list(filtered = filtered, mass = mass))
  }


#' @export
get_max = 
  function(V, outmax, domain = domain, keep_names, maximize_bf = F){ 
    # maximize_bf means divid by what i think is the prior probability of the constraint 
    if(domain == "wux"){
      next_level =  get_coarse_wux(names(outmax))
    } else if(domain == "tsx"){ 
      next_level =  get_coarse_tsx(names(outmax))
    }
    
    outme    = get_discrete_posterior(V, keep_names = keep_names, type = "pre", pre = next_level)$mass
    outmee   = unlist(outme)
    
    if(maximize_bf){
      # NO!
      # JC EXPERIMENTAL ************ POSSIBLY  adjust for the coset size
      # This is the ratio of prior to posterior mass in the encompassing prior that is
      # satisfied by the contstraint, which equals the Bayes Factor between the
      # constraint and the encompassing prior.
      number_of_elements_of_coset = next_level %>% map_depth(2, length) %>% map(unlist) %>% map(factorial) %>% map_dbl(reduce, `*`) 
      total_permutations          = next_level %>% map(unlist) %>% map(length) %>% map_dbl(factorial)
      prior_probability_of_coset  = number_of_elements_of_coset / total_permutations # e.g. p_1!p_2!p_3!=d!
      print("before size correction:******************")
      if(verbose <- T) print(names(outmee[(outmee) == max(outmee )]))
      outmee                      = outmee / number_of_elements_of_coset
      print("with size correction:******************")
      if(verbose) print(names(outmee[(outmee) == max(outmee )]))
    } else {
      max_bf = NULL
    }
    
    # What to maximize at each iteration. 1 the probability 2. the bayes factor
    # what to return at each iteration. 1 the probability 2. the bayes factor
    
    outmaxer = which(outmee == max(outmee))
    # if(length(outme[outmaxer]) != 1) print(str_c("multiple maxima : ", c(unlist(outme[outmaxer])), sep = "")) # shoudl I print???
    outmax = outme[outmaxer][1] # if many values attain maxima, choose the first arbitrarily to tie break
    
    if(maximize_bf) max_bf = outmee[outmaxer][1] # MAXIMUM BAYES FACTOR VERSUS MAXIMUM POSTERIOR PROBABILITY ABOVE
    
    if(outmax ==1){
      # base case: the trivial order always has probability 1, but may also stop before this
      # return(out = 1) # JC DANGER
      return(out = list(max = outmax))
    } else { 
      out = list(all = outmee,
                 max = outmax, 
                 max_bf = max_bf, 
                 get_max(V, outmax, domain, keep_names))
    }
    return(out = out)
  }



#' @export
get_coarse_wux = function(obj){
  # gets each element of the 1 step coarser/weaker order
  
  for(ii in str_locate_all(obj, "\\|")[[1]][,1]){
    
    init = str_split(obj, "\\|")[[1]]
    aa = head(unlist(str_remove(init,"\\|") %>% map(.,~c(.x, "|"))), -1)
    ll = 1
    qqq = NULL
    naaa = NULL
    for(ii in  which(aa == "|")){
      aa = head(unlist(str_remove(init,"\\|") %>% map(.,~c(.x, "|"))), -1)
      aa[ii] = ","
      naaa[ll] = str_c(aa, collapse = "")
      qqq[ll] = aa %>% str_c(collapse  = "") %>% str_split("\\|") %>% map(as.list) %>% map(~str_split(.x, ","))
      ll = ll + 1
    }
  }
  # browser()
  qqq = set_names(qqq, naaa)
  return(qqq = qqq)
}



#' @export
get_coarse_tsx = function(obj){
  # gets each element of the 1 step coarser/weaker order
  
  init =  str_split(obj, "\\|")[[1]]
  
  ll = 1
  qqq = NULL
  naaa = NULL
  
  for(ii in 1:length(init)){
    
    qqq[[ll]]  = as.list(init[-ii])
    naaa[[ll]] = str_c(init[-ii], collapse = "|")
    ll = ll + 1
  }
  
  qqq = set_names(qqq, naaa)
  return(qqq = qqq)
}


#' @export
get_predicates = function(V, keep_names, type, pre, n_terms){
  
  if(type == "all"){
    
    # browser()
    events = get_orders(V, keep_names = keep_names, type = type, n_terms = n_terms)
    events = flatten(flatten(events)) %>% keep(~length(.x)>1)
    
  } else if(type == "range"){
    
    events = get_orders(V, keep_names = keep_names,  type = type, n_terms = n_terms)
    events = events %>%  map(as.list)
    
  } else if(type == "pre"){
    
    events = pre
    
  }
  out = vector("list", length(events))
  partition_names = NULL
  part_names = NULL
  for(jj in 1:length(events)){ 
    one_event  = events[[jj]] # one partition of paramspace
    n_ranks_in_event = length(one_event)
    
    qq = NULL
    k = 1
    for(ii in 1:(n_ranks_in_event-1)){ 
      
      if(n_ranks_in_event == 1){
        # handle degenerate comparison with only one part, which is tautological 1 <= 1
        qq[[k]] = rep(TRUE, dim(V)[1])
        break
      }
      
      # FOR CONSECUTIVE PAIRS OF VARIABLE SETS IN DEFINITION OF THE EVENT, ASK
      # WHETHER EACH POINT INDEED SATIFYIES THE EVENT PREDICATE. FOR A POINT TO
      # QUALIFY, THE MAXIMUM OF THE LOWER VARIABLE SET SHOULD BE LOWER THAN THE
      # MINIMUM AMOUNG THE HIGHER VARIABLE SET.
      
      xx =   as.matrix(V[, one_event[[ii]],   drop = FALSE])
      yy =   as.matrix(V[, one_event[[ii+1]], drop = FALSE])  
      qq[[k]] = apply(xx, 1, max) <= apply(yy, 1, min) 
      
      if(diagnose <- F) plot(V[qq[[k]], unlist(one_event[k:(k+1)])]) # check pairwise inequalities
      k <- k + 1
    }
    
    out[[jj]] =  as.vector(reduce(qq, `&`))
    names(out)[jj]  =  str_c(map(one_event ,  ~str_c(.x, collapse = ",")), collapse = "|")
    
    if(verbose <- F) {
      print(names(out)[jj])
      print(accumulate(qq, `&`) %>% map_int(sum))  # diagnose number satisfying the increasingly stringent joint constraint
    }
  }
  
  return(out = out)
}



#' @export
get_orders = function(V, keep_names, type, n_terms){
  
  # browser()
  # GET ALL COMBINATIONS OF M FROM INDEX SET
  weak_orders_of_all_subsets = 
    map(n_terms, ~combn(x = 1:dim(V)[2], m = .x, simplify = FALSE)) %>% 
    flatten
  
  if(type == "range") {
    
    # TAKE all subsets, then consider all ordered partitions of each
    all_weak_orders  = weak_orders_of_all_subsets %>% map(combinat::permn) %>% flatten
    
    # CONTRAST WITH: TAKE all 3 part partitions, and consider all orderings: much, much bigger
    # jj <- restrictedparts(15,3, include.zero = F)
    # weak_orders_of_all_subsets = partitions::listParts(jj)
    
    if(keep_names){
      # keep colnames in the original sample matrix, versus numerically index columns with natural numbers (which works fine)
      all_weak_orders = map_depth(all_weak_orders, 2, ~names(V)[.x])
    }
    return(all_weak_orders)
  } 
  
  # A WORK AROUND TO EXTRACT ALL ORDERED PARTITIONS FOR EACH OF THE ABOVE SUBSETS
  # each of these subsets has a size i: find all ordered partitions of {1,...,i}
  events_of_i = 
    weak_orders_of_all_subsets %>% 
    map(length) %>% 
    map(partitions::listParts) %>% 
    map_depth(2, as.character) %>%  
    map_depth(2, ~combinat::permn(.x)) %>% 
    map_depth(4, ~eval(parse(text=.x)))
  
  # the ordered parts of all weak orders of every subset
  all_weak_orders = vector("list")
  for(ii in 1:length(events_of_i)){
    all_weak_orders[[ii]] = 
      map_depth(events_of_i[[ii]], 
                3,
                ~weak_orders_of_all_subsets[[ii]][.x])  
  }
  
  if(keep_names){
    # keep colnames in the original sample matrix, versus numerically index columns with natural numbers (which works fine)
    all_weak_orders = map_depth(all_weak_orders, 4, ~colnames(V)[.x])
  }
  return(out = all_weak_orders)
}

########################################################
# utils
#######################################################


#' @export
get_plot = function(samps, x, labl = NULL, tru, n_par, sim_n){
  # samps as argument to get_discrete_posterior()
  # x as return value of get_discrete_posterior()
  # labl an arbitrary string to plot in title
  
  # Visuals
  p2 = samps %>% gather(parameter, value) %>% mutate(parameter = factor(parameter)                                  ) %>% gf_violin(value~parameter) %>% gf_labs(title = "Natural (arbitary) order")
  p3 = samps %>% gather(parameter, value) %>% mutate(parameter = factor(parameter, levels =names(samps)[order(tru)])) %>% gf_violin(value~parameter) %>% gf_labs(title = "Reordered by true magnitude")
  p4 = map(x[1], plt_max)[[1]]
  p5 = map(x[2], plt_max)[[1]]
  # g = gridExtra::arrangeGrob(p2,p3,p4,p5, top = str_c(labl, " priors: ", str_c(str_c(order(tru), collapse = "|")), collapse = " :"))
  # possibly add bayes factors
  p4b = map(x[1], plt_max_no_axis_limits)[[1]]
  p5b = map(x[2], plt_max_no_axis_limits)[[1]]
  # g = gridExtra::arrangeGrob(p2,p3,p4,p5, p4b, p5b, top = str_c(labl, " priors: ", str_c(str_c(order(tru), collapse = "|")), collapse = " :"))
  g = gridExtra::arrangeGrob(p4, nrow = 1)
  
  # gridExtra::grid.arrange(g)
  flname = as.character(str_c("simules", "-", sim_n,  "-", as.character(n_par), "-", labl, "-", "bfY.pdf", sep = "", collapse = T))
  flname = str_c("r/range/fig/", flname)
  ggsave(file=flname, g) 
  
}


#' @export
get_tab = function(x){
  # x as return value of get_discrete_posterior()
  out = NULL
  (out$tab_tsx  = unlist(extract_tree(x$out_tsx,  "max")) %>% tibble::enframe(name = "order", value = "credibility"))
  (out$tab_wux = unlist(extract_tree(x$out_wux, "max")) %>% tibble::enframe(name = "order", value = "credibility"))
  return(out = out)
}


#' @export
get_fco = function(tab, alpha = 0.9){
  # tab is a field in the output from get_tab(), i.e. it is a cdf over a chain in the partial ranking poset
  
  # get the first value higher than alpha, the desired credibility
  out = tab %>% filter(what == "max") %>% filter(value > alpha)
  out = out$hyp[1]
  
}


#' @export
extract_tree = function(obj, field){ 
  # take a tree and recursively peel the layers away, extracting a named field at each point
  internal_node = NULL
  internal_node[[1]] = flatten(obj[[1]]) 
  for(ii in 1:(vec_depth(obj)-1)){
    obj = flatten(obj[-1]) # move one layer deeper in tree
    internal_node[[ii+1]] = obj[[field]]  # extract desired
  }
  internal_node
}



#' @export
plt_max = function(inn){
  
  if(str_length(inn$hyp[1]) < 25){
    inn %>% 
      filter(what == "max") %>%  
      mutate(hyp = factor(hyp , levels = hyp[1:dim(.)[1]]))  %>%
      rbind(tibble(what = NA, hyp = "any", value = 1)) %>% 
      mutate(Ranking = hyp, Probability = value) %>% 
      gf_point(Probability~Ranking, size = 5) %>% 
      gf_lims(y = c(0,1)) %>% 
      # gf_hline(yintercept = 0.9) %>%
      # gf_labs(title = "CDF of MAP chain") %>% 
      gf_theme(
        # theme = theme_light(),
        panel.background = element_rect(fill = "white",
                                        colour = "lightblue",
                                        size = 0.5, linetype = "solid"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 24),
        axis.text.y = element_text(angle = 45, hjust = 1, size = 24),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30)
      ) + 
      geom_hline(yintercept = 0.9)
    # axis.text=element_text(size=12)
    #+ coord_flip() 
  } else{ 
    inn %>% 
      filter(what == "max") %>%  
      mutate(hyp = factor(hyp , levels = hyp[1:dim(.)[1]]))  %>%
      rbind(tibble(what = NA, hyp = "any", value = 1)) %>% 
      mutate(hyp = 1:dim(.)[1]) %>% 
      mutate(ranking = hyp, probability = value) %>% 
      gf_point(probability~ranking)  
  }
  # gf_theme(axis.text.x = element_text(angle = 90, hjust = 1))
}



#' @export
plt_max_no_axis_limits = function(inn){
  
  # The tick labels get too long
  if(str_length(inn$hyp[1]) < 25){
    inn %>% 
      filter(what == "max_bf") %>%  
      mutate(hyp = factor(hyp , levels = hyp[1:dim(.)[1]]))  %>%
      rbind(tibble(what = NA, hyp = "any", value = 1)) %>% 
      gf_point(value~hyp) + 
      # gf_lims(y = c(0,1)) +
      coord_flip() 
  } else{ 
    inn %>% 
      filter(what == "max_bf") %>%  
      mutate(hyp = factor(hyp , levels = hyp[1:dim(.)[1]]))  %>%
      rbind(tibble(what = NA, hyp = "any", value = 1)) %>% 
      mutate(hyp = 1:dim(.)[1]) %>% 
      gf_point(value~hyp) +
      # gf_lims(y = c(0,1)) +
      coord_flip() 
  }
  # gf_theme(axis.text.x = element_text(angle = 90, hjust = 1))
}



#' @export
check_partition = function(x, V){
  # check x partitions V
  # eats a list of tibbles x
  # and a tibble V
  # Are they a partition? 
  inter    = reduce(x, intersect) 
  unione   = reduce(x, union)
  # browser()
  return(list(union = all_equal(V, unione), disjoint = dim(inter)[1]==0, names = names(x))) # disjoint? union?
}


#' @export
simulate_multi <- function(n_par = 5,  n_sims = 2000, n_rep = 50, dgf = 1){
  
  ########################################################
  # simulate
  ########################################################
  
  id = factor(rep(1:n_par, each = n_rep))
  
  # design
  if(ref_group <- F){
    # reference group parameterization
    X = model.matrix(~id)
  } else {
    X = kronecker(diag(n_par), rep(1, n_rep))
  }
  # outcome
  eps = rnorm(n_par*n_rep)
  tru = rt(n_par, dgf, 0) 
  if(dgf == 0) tru = rep(0, n_par) # the family null
  y = X %*% tru + eps
  D = tibble(y =c(y), id = id)
  
  ########################################################
  # infer:
  ########################################################
  
  out     <- NULL
  out$tru <- tru
  out$D   <- X
  
  if(ref_group) {
    out$fit_flat <- brm(y ~    id,  data = D, sample_prior = T, iter = n_sims) # 2) flat
  } else {
    out$fit_flat <- brm(y ~  -1 +  id,  data = D, sample_prior = T, iter = n_sims) # 2) flat
  } 
  
  out$fit_hier <- brm(y ~ (1|id), data = D, sample_prior = T, iter = n_sims) # 1) multilevel
  
  return(out = out)
}



#' @export
violated <- function(fco, tru){
  # example: violated("2|3,1","2|3|1" )
  
  if(is.character(tru)) {if(tru == "family_null")  return(inconsistent = FALSE)} # universal order always contains - therefore cannot contradict - truth
  if(fco == "universe") {
    return(inconsistent = FALSE) # universal order always contains - therefore cannot contradict - truth
  }
  
  # IF YOU KNOW THE GROUND TRUTH, YOU CAN CHECK WHETHER AN INFERED FCO VIOLATES THE TRUE TOTAL ORDER
  subseet = as.numeric(unlist(unlist(str_split(fco, "\\|")) %>% str_split(",")))
  
  # handle the case of character argument, e.g. 3|1,2|4
  if(is.character(tru)) tru = as.numeric(unlist(unlist(str_split(tru, "\\|")) %>% str_split(",")))
  
  # look for inconsistency between rank in fco and true underlying rank 
  ranker = NULL
  ll = 1
  for(ii in intersect(tru, subseet)){ 
    ranker[ll] = str_split(fco, "\\|") %>% map(str_detect, str_c( ",", as.character(ii), ",|,", as.character(ii), "$|^", as.character(ii), ",|^", as.character(ii), "$")) %>% map(which) %>% unlist
    ll = ll + 1
  }
  # returns the number of rank inconsistencies
  inconsistent = sum(diff(ranker) < 0) > 0
  return(inconsistent = inconsistent)
}



########################################################
# WRITTEN as lmhyp package wrapper
########################################################
# SUM


#' @export
find_local_fco <- function(partial_rank, dat, thresh = .9){
    
    # Check if the most probable full order passes critereon 
    # This condition applies only iteration 1 of recursive call 
    # (which starts at the max probable full order)
    init  = str_replace_all(partial_rank,"\\|", "<")
    init_is_full_order = str_count(init,"w") == (str_count(init, "<") + 1)
    if(init_is_full_order){
      initp = filter(dat, H == init)$Hp
      if(initp > thresh) return(c(partial_rank, initp))
    }
    
    coarser    = 
      get_coarse_wux(partial_rank) %>%
      map_depth(2, combinat::permn) %>%
      map_depth(3, str_c, collapse = "<") %>% 
      map(cross) %>% 
      map_depth(2, str_c, collapse = "<") %>% 
      map(unlist)
    probs      = map(coarser, ~filter(dat, H %in% .x)) %>% map(summarise, sum(Hp)) 
    max_prob   = max(unlist(probs))
    max_event  = probs %>% keep(~.x == max_prob) %>% names %>% `[[`(1) # the last step in case of multiple max
    # print(c(max_event, max_prob))  
    
    if(max_prob > thresh | !str_detect(max_event, "\\|")){ 
      # print("********************")
      return(c(max_event, max_prob)) 
    } else {  
      find_local_fco(max_event, dat)
    } 
  }

# MEAN
#' @export
find_local_fco_mean = 
  function(partial_rank, dat, thresh = .9){
    
    # Check if the most probable full order passes critereon 
    # This condition applies only iteration 1 of recursive call 
    # (which starts at the max probable full order)
    init  = str_replace_all(partial_rank,"\\|", "<")
    init_is_full_order = str_count(init,"w") == (str_count(init, "<") + 1)
    if(init_is_full_order){ 
      initp = filter(dat, H == init)$Hp
      if(initp > thresh) return(c(partial_rank, initp)) 
    }
    
    # Otherwise coarsen
    coarser  = get_coarse_wux(partial_rank) %>% map_depth(2, combinat::permn) %>% map_depth(3, str_c, collapse = "<") %>% map(cross) %>% map_depth(2, str_c, collapse = "<") %>% map(unlist)
    probs    = map(coarser, ~filter(dat, H %in% .x)) %>% map(summarise, sum(Hp))
    # max av prob
    av_probs = map(coarser, ~filter(dat, H %in% .x)) %>% map(summarise, mean(Hp))
    max_av_prob = max(unlist(av_probs))
    max_event = av_probs %>% keep(~.x == max_av_prob) %>% names %>% `[[`(1) # the last step in case of multiple max
    max_prob = unlist(probs[[max_event]][[1]])
    print(c(max_event, max_prob))
    
    if(max_prob > thresh | !str_detect(max_event, "\\|")){ 
      print("********************")
      return(c(max_event, max_prob)) 
    } else {  
      find_local_fco_mean(max_event, dat)
    } 
  }


#' @export
dysfunction = 
  function(partial_rank, Hp, obj, thresh = .9){
    
    init  = str_replace_all(partial_rank,"\\|", "<")
    coarser = get_coarse_wux(partial_rank) %>% map_depth(2, combinat::permn) %>% map_depth(3, str_c, collapse = "<") %>% map(cross) %>% map_depth(2, str_c, collapse = "<") %>% map(unlist)
    hyps = unique(unlist(coarser))
    hyps = str_c(hyps, collapse = ";")
    post = test_hyp(object = obj, hyp = hyps)
    
    # post$post_prob[init]
    pst = post$post_prob
    names(pst) = post$hypotheses
    if(max(pst) == pst[init]){
      reference_index = which(post$hypotheses %in% init)
      correction_factor = Hp / pst[reference_index]
      pst = pst * correction_factor
      prb_coarsening = sum(pst)
      
      return(prb_coarsening = prb_coarsening)
    } else {
      
      print("full order of map parameter is not the most probable full order ")
      return(prb_coarsening = "path")
    }
    
  }