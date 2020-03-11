##################################################################################################################################################################
# PROTOCOL    : WO39392 (IMpassion031)
# OBJECTIVE(s): To simulate the proposed adaptive design under a variety of conditions in order to evaluate its operating characteristics and its effectiveness 
# AUTHOR(s)   : Marcel Wolbers, Dominik Heinzmann
# NOTE(s)     : Decision rules @ end of stage 1 are based on the treatment effect (difference in proportion between experimental group and control group) in S & T
# Last modified: 10 Sep 2018 by Anh Nguyen Duc
##################################################################################################################################################################
require(magrittr)
require(gsDesign)
##################################################################################################################################################################

# IMPORTANT: all sample size is PER ARM with RR=1:1!!!

AdaptiveTrialDesign <- function(p0_S=0.48,              # Probability of response in subgroup S in the control arm
                                        p0_T=0.48,              # Probability of response in complement subgroup T in the control arm
                                        delta_S=0.2,            # Treatment effect (difference in proportion between experimental group and control group) in S 
                                        delta_T=0.2,            # Treatment effect (difference in proportion between experimental group and control group) in T
                                        prev_S=0.47,            # Prevalence of subgroup S
                                        n_F_1=102,              # Sample size per arm for stage 1 related to full population F
                                        n_2=60,                 # Sample size per arm for stage 2 (either for full population F or, if enriched, for S depending on stage 1 decision)
                                        d_T=0.10,               # Threshold in complement group T for decision rule at end of stage 1
                                        d_S=0.15,               # Threshold in subgroup S for decision rule at end of stage 1
                                        sig_level=0.025,        # Overall type I error (one-sided testing)
                                        prop_alpha_spent_1=0.5, # Proportion of type I error to be spent at stage 1 (for early efficacy stop)
                                        drop.out.rate=0.05,     # Drop out proportion similar in both arms (patients "dropping out" prior to surgery are counted as non-responder)
                                        n_2.scale=1,            # scaling factor for n_2 in the calculation of info_1, hence w_1 and w_2
                                        n_1.scale=1,            # scaling factor for n_1 in the calculation of info_1, hence w_1 and w_2
                                        cli_rev_S=0.15,         # Threshold for final clinical rev in S
                                        cli_rev_F=0.15,         # Threshold for fianl cincal rev in ITT
                                        nSimul=1e4,             # number of simulations
                                        stage1_data=NULL,
                                        w_1=NULL, w_2=NULL, sig_comb=NULL
){
  
##---------------- Simulation of a simple adaptive enrichment design with fixed sample sizes for both stages
## !!!! Preliminary version, use at your own risk !!!!
  
# require(gsDesign)

##---------------- Parameter settings
tryCatch({
    
## These must be based on the correct & real prevalence
# response in full population F
p0_F <- prev_S*p0_S+(1-prev_S)*p0_T

# # delta in full population F
delta_F <- prev_S*delta_S+(1-prev_S)*delta_T

## One-sided significance level at stage 1 based on the alpha spent at this stage

sig_1 <- sig_level*prop_alpha_spent_1


foo <- function(tau, p_F, p_S) { # function to calculate Spiessens Dubois p-val for inters. hypo.

  foo1 <- function(x, z, tau) {
    pnorm( q = ( z - sqrt(tau)*x ) / sqrt(1-tau) ) * dnorm(x)
  }
  
  gam <- min(p_F, p_S)*2
  z   <- qnorm(1-gam/2)
  
  1 - integrate(f=foo1, lower=-Inf, upper=z, z=z, tau=tau)$value # assuming equal "weight" for F and S
  # min(2*min(p_F,p_S),max(p_F,p_S))
}


##---------------- Simulate data for stage 1 and 2

#-----------------------------------------------------------------------------------------
#--- Stage 1
#-----------------------------------------------------------------------------------------

## These must be based on the observed prevalence
# Sample size stage 1 per population
n_S_1 <- round(prev_S*n_F_1)
n_T_1 <- n_F_1-n_S_1

# Weights p-value combination for stage 1 and stage 2
# info_1 <- n_F_1/(n_F_1+n_2) # "information fraction stage 1"
## info_1 <- n_S_1/(n_S_1+n_2) # "information fraction stage 1"

# n.all.0 <- n_F_1+n_2
# n_F_1_effective <- n.all.0 - n_2.effective
# info_1 <- n_F_1_effective / n.all.0

n_2.effective  <- ceiling(n_2*n_2.scale)
n_1.effective  <- ceiling(n_F_1*n_1.scale)

n.all.effective <- n_1.effective + n_2.effective
info_1 <- n_1.effective / n.all.effective


w_1 <- sqrt(info_1)
w_2 <- sqrt(1-info_1)




# Calculate nominal one-significance level for stage 1 & 2 combination test using the group sequential design methodology (gsDesign package)
rho <- log(prop_alpha_spent_1)/log(info_1) 
g <- gsDesign(k = 2,timing=c(info_1,1),test.type=1,alpha =sig_level,sfu=sfPower,sfupar=rho) # sfPower=Kim-DeMets (power) spending function (Jennison and Turnbull 2000)
sig_comb <- pnorm(-g$upper$bound)[2]  

if (is.null(stage1_data)) {

  # number of responders (cr0 is control, cr1 is intervention)
  # drop-outs incorporated assuming given p0_s/T are proportions in target population; drop-out assumed prior 
  # surgery (hypothetical responders drop-outs are non-responder by definition) 
  n_cr0_S_1 <- rbinom(nSimul,size=n_S_1,prob=p0_S*(1-drop.out.rate)) # Anh: is this correct, esp. the denom.? => correct as all sample size is per-arm (assuming 1:1 RR)
  n_cr1_S_1 <- rbinom(nSimul,size=n_S_1,prob=(p0_S+delta_S)*(1-drop.out.rate))
  n_cr0_T_1 <- rbinom(nSimul,size=n_T_1,prob=p0_T*(1-drop.out.rate))
  n_cr1_T_1 <- rbinom(nSimul,size=n_T_1,prob=(p0_T+delta_T)*(1-drop.out.rate))
  
  
  n_cr0_F_1 <- n_cr0_S_1+n_cr0_T_1 # Number of pCRs for full pop after stage 1 in control arm
  n_cr1_F_1 <- n_cr1_S_1+n_cr1_T_1 # Number of pCRs for full pop after stage 1 in experimental arm
  n_cr_F_1  <- n_cr0_F_1+n_cr1_F_1 # Number of pCRs for full pop after stage 1 in both arms
  
  # estimated delta and p-value in full pop and subgroup from stage 1
  stage1_data <- data.frame(est_delta_F_1=(n_cr1_F_1-n_cr0_F_1)/n_F_1, 
                            est_delta_S_1=(n_cr1_S_1-n_cr0_S_1)/n_S_1,
                            est_delta_T_1=(n_cr1_T_1-n_cr0_T_1)/n_T_1,
                            n_cr1_S_1=n_cr1_S_1, n_cr0_S_1=n_cr0_S_1, n_S_1=n_S_1,
                            n_cr1_F_1=n_cr1_F_1, n_cr0_F_1=n_cr0_F_1, n_F_1=n_F_1, n_cr_F_1=n_cr_F_1,
                            pval_F_1=apply(cbind(n_cr1_F_1,n_cr0_F_1),1,function(x){ prop.test(x,c(n_F_1,n_F_1),alternative="greater", correct = F)$p.value }),  # prop.test sufficient precise approximation given N
                            pval_S_1=apply(cbind(n_cr1_S_1,n_cr0_S_1),1,function(x){ prop.test(x,c(n_S_1,n_S_1),alternative="greater", correct = F)$p.value }),
                            z_F_1=apply(cbind(n_cr1_F_1,n_cr0_F_1),1,function(x){ prop.test(x,c(n_F_1,n_F_1),alternative="greater", correct = F)$statistic %>% sqrt }),
                            z_S_1=apply(cbind(n_cr1_S_1,n_cr0_S_1),1,function(x){ prop.test(x,c(n_S_1,n_S_1),alternative="greater", correct = F)$statistic %>% sqrt })                            
                            ) # Anh: pertaining to the comment on n_cro_S_1, are the denom. here correct? => correct as all sample size is per-arm (assuming 1:1 RR)
  
  # Simes p-value intersection test
  stage1_data$pval_intersect_1 <- pmin(2*pmin(stage1_data$pval_F_1,stage1_data$pval_S_1),
                                       pmax(stage1_data$pval_F_1,stage1_data$pval_S_1))

} #end of if (is.null(stage1_data))

# Spiessens Dubois p-value intersection test
stage1_data$pval_intersect_1_SD <- sapply(1:nrow(stage1_data), FUN = function(id) {
  foo(tau=prev_S, p_F=stage1_data$pval_F_1[id], p_S=stage1_data$pval_S_1[id])
})
# Exlanation: Given ordered p-values to test H1 and H2, p1 and p2, Sime's rejects H1 intersect H2 if p1<=alpha/2 OR p2<=alpha
# browser()
stage1_data$cp_stage2_S_half_alp_co <- sapply(1:nrow(stage1_data), FUN = function(id) cp_stage2(alp2 = sig_comb/2, z1 = stage1_data$z_S_1[id], n1 = n_S_1, n2 = n_2*prev_S, w1 = w_1, w2 = w_2))

stage1_data$cp_stage2_S_full_alp_en <- sapply(1:nrow(stage1_data), FUN = function(id) cp_stage2(alp2 = sig_comb, z1 = stage1_data$z_S_1[id], n1 = n_S_1, n2 = n_2, w1 = w_1, w2 = w_2))

stage1_data$cp_stage2_F_half_alp_co <- sapply(1:nrow(stage1_data), FUN = function(id) cp_stage2(alp2 = sig_comb/2, z1 = stage1_data$z_F_1[id], n1 = n_F_1, n2 = n_2, w1 = w_1, w2 = w_2))

stage1_data$cp_stage2_F_full_alp_ac <- sapply(1:nrow(stage1_data), FUN = function(id) cp_stage2(alp2 = sig_comb, z1 = stage1_data$z_F_1[id], n1 = n_F_1, n2 = n_2, w1 = w_1, w2 = w_2))

#-----------------------------------------------------------------------------------------
#--- Stage 2 data in case trial continues in overall population F
#-----------------------------------------------------------------------------------------

# browser()
# Sample sizes stage 2
n_F_2 <- n_2
n_S_2 <- round(prev_S*n_F_2) # assuming same prevalence as for stage 1 (take estimate from stage 1 know prior to interim analysis)
n_T_2 <- n_F_2-n_S_2

# number of responders (cr0 is control, cr1 is intervention)
n_cr0_S_2 <- rbinom(nSimul,size=n_S_2,prob=p0_S*(1-drop.out.rate))
n_cr1_S_2 <- rbinom(nSimul,size=n_S_2,prob=(p0_S+delta_S)*(1-drop.out.rate))
n_cr0_T_2 <- rbinom(nSimul,size=n_T_2,prob=p0_T*(1-drop.out.rate))
n_cr1_T_2 <- rbinom(nSimul,size=n_T_2,prob=(p0_T+delta_T)*(1-drop.out.rate))

n_cr0_F_2 <- n_cr0_S_2+n_cr0_T_2
n_cr1_F_2 <- n_cr1_S_2+n_cr1_T_2

# estimated delta and p-value in full pop and subgroup from stage 1
stage2_data_no_enrichment <- data.frame(est_delta_F_2=(n_cr1_F_2-n_cr0_F_2)/n_F_2, 
                                        est_delta_S_2=(n_cr1_S_2-n_cr0_S_2)/n_S_2,
                                        n_cr1_S_2=n_cr1_S_2, n_cr0_S_2=n_cr0_S_2, n_S_2=n_S_2,
                                        n_cr1_F_2=n_cr1_F_2, n_cr0_F_2=n_cr0_F_2, n_F_2=n_F_2,
                                        pval_F_2=apply(cbind(n_cr1_F_2,n_cr0_F_2),1,function(x){ prop.test(x,c(n_F_2,n_F_2),alternative="greater", correct = F)$p.value }),
                                        pval_S_2=apply(cbind(n_cr1_S_2,n_cr0_S_2),1,function(x){ prop.test(x,c(n_S_2,n_S_2),alternative="greater", correct = F)$p.value }))

# Simes p-value intersection test
stage2_data_no_enrichment$pval_intersect_2 <- pmin(2*pmin(stage2_data_no_enrichment$pval_F_2,stage2_data_no_enrichment$pval_S_2),
                                                   pmax(stage2_data_no_enrichment$pval_F_2,stage2_data_no_enrichment$pval_S_2))


stage2_data_no_enrichment$pval_intersect_2_SD <- sapply(1:nrow(stage2_data_no_enrichment), FUN = function(id) {
foo(tau=prev_S, p_F=stage2_data_no_enrichment$pval_F_2[id], 
    p_S=stage2_data_no_enrichment$pval_S_2[id])
})

#-----------------------------------------------------------------------------------------
#--- Stage 2 data in case trial in case only subgroup continues (including enrichment)
#-----------------------------------------------------------------------------------------

# Sample sizes stage 2
n_S_2 <- n_2 # all patients reecruited are patients with subgroup characteristics

# number of responders (cr0 is control, cr1 is intervention)
n_cr0_S_2 <- rbinom(nSimul,size=n_S_2,prob=p0_S*(1-drop.out.rate))
n_cr1_S_2 <- rbinom(nSimul,size=n_S_2,prob=(p0_S+delta_S)*(1-drop.out.rate))

# estimated delta and p-value from stage 2
stage2_enrichment_data <- data.frame(est_delta_S_2enr=(n_cr1_S_2-n_cr0_S_2)/n_S_2,
                                     n_cr1_S_2=n_cr1_S_2, n_cr0_S_2=n_cr0_S_2, n_S_2=n_S_2,
                                     pval_S_2enr=apply(cbind(n_cr1_S_2,n_cr0_S_2),1,function(x){ prop.test(x,c(n_S_2,n_S_2),alternative="greater", correct = F)$p.value }))


#-----------------------------------------------------------------------------------------
#--- Derive final decisions for each simulation
#-----------------------------------------------------------------------------------------

trial_result <- data.frame(stage1_decision=rep(NA,nSimul), stage1_decision_detail=rep(NA,nSimul), 
                           signif_F=NA, signif_S=NA, effect.1=NA,
                           effect_F_overall=NA, effect_S_overall=NA,
                           p_comb_intersect=NA, p_comb_F=NA, p_comb_S=NA,
                           effect_S_1=stage1_data$est_delta_S_1,
                           effect_F_1=stage1_data$est_delta_F_1,
                           effect_T_1=stage1_data$est_delta_T_1, n_cr_F_1=stage1_data$n_cr_F_1,
                           pval_F_1=stage1_data$pval_F_1,
                           pval_S_1=stage1_data$pval_S_1,
                           pval_intersect_1=stage1_data$pval_intersect_1,
                           
                           stage1_decision_SD=rep(NA,nSimul), stage1_decision_detail_SD=rep(NA,nSimul), 
                           signif_F_SD=NA, signif_S_SD=NA, 
                           effect_F_overall_SD=NA, effect_S_overall_SD=NA,                           
                           pval_intersect_1_SD=stage1_data$pval_intersect_1_SD,
                           est_ctr_rate_F_1 = stage1_data$n_cr0_F_1 / n_F_1,
                           est_ctr_rate_S_1 = stage1_data$n_cr0_S_1 / n_S_1)

for (i in 1:nSimul) {
  # Decision to stop for efficacy at the end of stage 1

  if (stage1_data$pval_intersect_1[i]<=sig_1) { # Sime
    trial_result$stage1_decision[i] <- "1 - Stop after stage 1: efficacy"
    trial_result$signif_F[i] <- (stage1_data$pval_F_1[i]<=sig_1)
    trial_result$signif_S[i] <- (stage1_data$pval_S_1[i]<=sig_1)

    if (stage1_data$pval_F_1[i]<=sig_1) {
      trial_result$effect.1[i] <- stage1_data$est_delta_F_1[i]
      trial_result$stage1_decision_detail[i] <- '1a - Stop after stage 1: efficacy in F only'
    }

    if (stage1_data$pval_S_1[i]<=sig_1) {
      trial_result$stage1_decision_detail[i] <- '1b - Stop after stage 1: efficacy in S only'
    }
    
    if (stage1_data$pval_S_1[i]<=sig_1 & stage1_data$pval_F_1[i]<=sig_1) {
      trial_result$stage1_decision_detail[i] <- '1c - Stop after stage 1: efficacy in F and S'
    }
    
    trial_result$effect_F_overall[i] <- stage1_data$est_delta_F_1[i]
    trial_result$effect_S_overall[i] <- stage1_data$est_delta_S_1[i]
  }
  
  if (stage1_data$pval_intersect_1_SD[i]<=sig_1) {
    trial_result$stage1_decision_SD[i] <- "1 - Stop after stage 1: efficacy"
    trial_result$signif_F_SD[i] <- (stage1_data$pval_F_1[i]<=sig_1)
    trial_result$signif_S_SD[i] <- (stage1_data$pval_S_1[i]<=sig_1)
    
    if (stage1_data$pval_F_1[i]<=sig_1) {
      # trial_result$effect.1_SD[i] <- stage1_data$est_delta_F_1[i]
      trial_result$stage1_decision_detail_SD[i] <- '1a - Stop after stage 1: efficacy in F only'
    }
    
    if (stage1_data$pval_S_1[i]<=sig_1) {
      trial_result$stage1_decision_detail_SD[i] <- '1b - Stop after stage 1: efficacy in S only'
    }
    
    if (stage1_data$pval_S_1[i]<=sig_1 & stage1_data$pval_F_1[i]<=sig_1) {
      trial_result$stage1_decision_detail_SD[i] <- '1c - Stop after stage 1: efficacy in F and S'
    }
    
    trial_result$effect_F_overall_SD[i] <- stage1_data$est_delta_F_1[i]
    trial_result$effect_S_overall_SD[i] <- stage1_data$est_delta_S_1[i]
  }

  # Decision to stop for futility at the end of stage 1

  delta_T_interim <- stage1_data$est_delta_T_1[i]
  delta_S_interim <- stage1_data$est_delta_S_1[i]
  delta_F_interim <- stage1_data$est_delta_F_1[i]

  # Sime
  if ((stage1_data$pval_intersect_1[i]>sig_1) & (delta_T_interim<d_T) & (delta_S_interim<d_S)) {
    trial_result$stage1_decision[i] <- trial_result$stage1_decision_detail[i] <- "2 - Stop after stage 1: futility"
    trial_result$signif_F[i] <- F
    trial_result$signif_S[i] <- F
    
    trial_result$effect_F_overall[i] <- stage1_data$est_delta_F_1[i]
    trial_result$effect_S_overall[i] <- stage1_data$est_delta_S_1[i]
  }
  
  # SD
  if ((stage1_data$pval_intersect_1_SD[i]>sig_1) & (delta_T_interim<d_T) & (delta_S_interim<d_S)) {
    trial_result$stage1_decision_SD[i] <- trial_result$stage1_decision_detail_SD[i] <- "2 - Stop after stage 1: futility"
    trial_result$signif_F_SD[i] <- F
    trial_result$signif_S_SD[i] <- F
    
    trial_result$effect_F_overall_SD[i] <- stage1_data$est_delta_F_1[i]
    trial_result$effect_S_overall_SD[i] <- stage1_data$est_delta_S_1[i]
  }

  # Maximum observed treatement effect in S or F
  max_eff <- max(delta_S_interim,delta_F_interim)

  # Decision to continue to stage 2 with S only
  # Sime
  if ((stage1_data$pval_intersect_1[i]>sig_1) & (delta_T_interim<d_T) & (delta_S_interim>=d_S)) {
    trial_result$stage1_decision[i] <- trial_result$stage1_decision_detail[i] <- "3- Continue to stage 2: S only"

    # calculate combination test p-values
    z_comb_intersect <- w_1*qnorm(1-stage1_data$pval_intersect_1[i])+w_2*qnorm(1-stage2_enrichment_data$pval_S_2enr[i])
    z_comb_S <- w_1*qnorm(1-stage1_data$pval_S_1[i])+w_2*qnorm(1-stage2_enrichment_data$pval_S_2enr[i])
    trial_result$p_comb_intersect[i] <- p_comb_intersect <- pnorm(-z_comb_intersect) # 1-p=pnorm(Z), hence p=1-pnorm(Z)=pnorm(-Z)
    trial_result$p_comb_S[i] <- p_comb_S <- pnorm(-z_comb_S)

    # Test decisions
    if ((p_comb_intersect<sig_comb)&(p_comb_S<sig_comb)) {
      trial_result$signif_S[i] <- T
    } else {
      trial_result$signif_S[i] <- F
    }
    trial_result$signif_F[i] <- F

    trial_result$effect_S_overall[i] <- 
      ( (stage1_data$n_cr1_S_1[i] + stage2_enrichment_data$n_cr1_S_2[i]) - (stage1_data$n_cr0_S_1[i] + stage2_enrichment_data$n_cr0_S_2[i]) ) / (stage1_data$n_S_1[i] + stage2_enrichment_data$n_S_2[i])
    trial_result$effect_F_overall[i] <- stage1_data$est_delta_F_1[i]
  }
  # SD
  if ((stage1_data$pval_intersect_1_SD[i]>sig_1) & (delta_T_interim<d_T) & (delta_S_interim>=d_S)) {
    trial_result$stage1_decision_SD[i] <- trial_result$stage1_decision_detail_SD[i] <- "3- Continue to stage 2: S only"
    
    # calculate combination test p-values
    z_comb_intersect <- w_1*qnorm(1-stage1_data$pval_intersect_1_SD[i])+w_2*qnorm(1-stage2_enrichment_data$pval_S_2enr[i])
    z_comb_S <- w_1*qnorm(1-stage1_data$pval_S_1[i])+w_2*qnorm(1-stage2_enrichment_data$pval_S_2enr[i])
    p_comb_intersect <- pnorm(-z_comb_intersect) # 1-p=pnorm(Z), hence p=1-pnorm(Z)=pnorm(-Z)
    p_comb_S <- pnorm(-z_comb_S)
    
    # Test decisions
    if ((p_comb_intersect<sig_comb)&(p_comb_S<sig_comb)) {
      trial_result$signif_S_SD[i] <- T
    } else {
      trial_result$signif_S_SD[i] <- F
    }
    trial_result$signif_F_SD[i] <- F
    
    trial_result$effect_S_overall_SD[i] <- 
      ( (stage1_data$n_cr1_S_1[i] + stage2_enrichment_data$n_cr1_S_2[i]) - (stage1_data$n_cr0_S_1[i] + stage2_enrichment_data$n_cr0_S_2[i]) ) / (stage1_data$n_S_1[i] + stage2_enrichment_data$n_S_2[i])
    trial_result$effect_F_overall_SD[i] <- stage1_data$est_delta_F_1[i]
  }

  # Decision to continue to stage 2 with F only
  # Sime
  if ((stage1_data$pval_intersect_1[i]>sig_1) & (delta_T_interim>=d_T)&(delta_S_interim<d_S)) {
    trial_result$stage1_decision[i] <- trial_result$stage1_decision_detail[i] <- "4- Continue to stage 2: F only"

    # calculate combination test p-values
    z_comb_intersect <- w_1*qnorm(1-stage1_data$pval_intersect_1[i])+w_2*qnorm(1-stage2_data_no_enrichment$pval_F_2[i])
    z_comb_F <- w_1*qnorm(1-stage1_data$pval_F_1[i])+w_2*qnorm(1-stage2_data_no_enrichment$pval_F_2[i])
    trial_result$p_comb_intersect[i] <- p_comb_intersect <- pnorm(-z_comb_intersect)
    trial_result$p_comb_F[i] <- p_comb_F <- pnorm(-z_comb_F)

    # Test decisions
    if ((p_comb_intersect<=sig_comb)&(p_comb_F<=sig_comb)) {
      trial_result$signif_F[i] <- T
    } else {
      trial_result$signif_F[i] <- F
    }
    trial_result$signif_S[i] <- F
    
    trial_result$effect_F_overall[i] <- 
      ( (stage1_data$n_cr1_F_1[i] + stage2_data_no_enrichment$n_cr1_F_2[i]) - (stage1_data$n_cr0_F_1[i] + stage2_data_no_enrichment$n_cr0_F_2[i]) ) / (stage1_data$n_F_1[i] + stage2_data_no_enrichment$n_F_2[i])
    trial_result$effect_S_overall[i] <- stage1_data$est_delta_S_1[i]
  }
  # SD
  if ((stage1_data$pval_intersect_1_SD[i]>sig_1) & (delta_T_interim>=d_T)&(delta_S_interim<d_S)) {
    trial_result$stage1_decision_SD[i] <- trial_result$stage1_decision_detail_SD[i] <- "4- Continue to stage 2: F only"
    
    # calculate combination test p-values
    z_comb_intersect <- w_1*qnorm(1-stage1_data$pval_intersect_1_SD[i])+w_2*qnorm(1-stage2_data_no_enrichment$pval_F_2[i])
    z_comb_F <- w_1*qnorm(1-stage1_data$pval_F_1[i])+w_2*qnorm(1-stage2_data_no_enrichment$pval_F_2[i])
    p_comb_intersect <- pnorm(-z_comb_intersect)
    p_comb_F <- pnorm(-z_comb_F)
    
    # Test decisions
    if ((p_comb_intersect<=sig_comb)&(p_comb_F<=sig_comb)) {
      trial_result$signif_F_SD[i] <- T
    } else {
      trial_result$signif_F_SD[i] <- F
    }
    trial_result$signif_S_SD[i] <- F
    
    trial_result$effect_F_overall_SD[i] <- 
      ( (stage1_data$n_cr1_F_1[i] + stage2_data_no_enrichment$n_cr1_F_2[i]) - (stage1_data$n_cr0_F_1[i] + stage2_data_no_enrichment$n_cr0_F_2[i]) ) / (stage1_data$n_F_1[i] + stage2_data_no_enrichment$n_F_2[i])
    trial_result$effect_S_overall_SD[i] <- stage1_data$est_delta_S_1[i]
  }
  
  # Decision to continue to stage 2 with F and S
  # Sime
  if ((stage1_data$pval_intersect_1[i]>sig_1)&(delta_T_interim>=d_T)&(delta_S_interim>=d_S)) {
    trial_result$stage1_decision[i] <- trial_result$stage1_decision_detail[i] <- "5- Continue to stage 2: F and S"

    # calculate combination test p-values
    z_comb_intersect <- w_1*qnorm(1-stage1_data$pval_intersect_1[i])+w_2*qnorm(1-stage2_data_no_enrichment$pval_intersect_2[i])
    z_comb_F <- w_1*qnorm(1-stage1_data$pval_F_1[i])+w_2*qnorm(1-stage2_data_no_enrichment$pval_F_2[i])
    z_comb_S <- w_1*qnorm(1-stage1_data$pval_S_1[i])+w_2*qnorm(1-stage2_data_no_enrichment$pval_S_2[i])
    trial_result$p_comb_intersect[i] <- p_comb_intersect <- pnorm(-z_comb_intersect)
    trial_result$p_comb_F[i] <- p_comb_F <- pnorm(-z_comb_F)
    trial_result$p_comb_S[i] <- p_comb_S <- pnorm(-z_comb_S)

    # Test decisions
    if (p_comb_intersect<=sig_comb) {
      if (p_comb_F<=sig_comb) {
        trial_result$signif_F[i] <- T
      }
      if (p_comb_F>sig_comb){
        trial_result$signif_F[i] <- F
      }
      if (p_comb_S<=sig_comb) {
        trial_result$signif_S[i] <- T
      }
      if (p_comb_S>sig_comb){
        trial_result$signif_S[i] <- F
      }
    }
    
    if (p_comb_intersect>sig_comb) {
      trial_result$signif_F[i] <- F
      trial_result$signif_S[i] <- F
    }

    trial_result$effect_S_overall[i] <- 
      ( (stage1_data$n_cr1_S_1[i] + stage2_data_no_enrichment$n_cr1_S_2[i]) - (stage1_data$n_cr0_S_1[i] + stage2_data_no_enrichment$n_cr0_S_2[i]) ) / (stage1_data$n_S_1[i] + stage2_data_no_enrichment$n_S_2[i])
    
    trial_result$effect_F_overall[i] <- 
      ( (stage1_data$n_cr1_F_1[i] + stage2_data_no_enrichment$n_cr1_F_2[i]) - (stage1_data$n_cr0_F_1[i] + stage2_data_no_enrichment$n_cr0_F_2[i]) ) / (stage1_data$n_F_1[i] + stage2_data_no_enrichment$n_F_2[i])    
  } # end of if ((stage1_data$pval_intersect_1[i]>sig_1)&(delta_T_interim>=d_T)&(delta_S_interim>=d_S))
  
  # SD
  if ((stage1_data$pval_intersect_1_SD[i]>sig_1)&(delta_T_interim>=d_T)&(delta_S_interim>=d_S)) {
    trial_result$stage1_decision_SD[i] <- trial_result$stage1_decision_detail_SD[i] <- "5- Continue to stage 2: F and S"
    
    # calculate combination test p-values
    z_comb_intersect <- w_1*qnorm(1-stage1_data$pval_intersect_1_SD[i])+w_2*qnorm(1-stage2_data_no_enrichment$pval_intersect_2[i])
    z_comb_F <- w_1*qnorm(1-stage1_data$pval_F_1[i])+w_2*qnorm(1-stage2_data_no_enrichment$pval_F_2[i])
    z_comb_S <- w_1*qnorm(1-stage1_data$pval_S_1[i])+w_2*qnorm(1-stage2_data_no_enrichment$pval_S_2[i])
    p_comb_intersect <- pnorm(-z_comb_intersect)
    p_comb_F <- pnorm(-z_comb_F)
    p_comb_S <- pnorm(-z_comb_S)
    
    # Test decisions
    if (p_comb_intersect<=sig_comb) {
      if (p_comb_F<=sig_comb) {
        trial_result$signif_F_SD[i] <- T
      }
      if (p_comb_F>sig_comb){
        trial_result$signif_F_SD[i] <- F
      }
      if (p_comb_S<=sig_comb) {
        trial_result$signif_S_SD[i] <- T
      }
      if (p_comb_S>sig_comb){
        trial_result$signif_S_SD[i] <- F
      }
    }
    
    if (p_comb_intersect>sig_comb) {
      trial_result$signif_F_SD[i] <- F
      trial_result$signif_S_SD[i] <- F
    }
    
    trial_result$effect_S_overall_SD[i] <- 
      ( (stage1_data$n_cr1_S_1[i] + stage2_data_no_enrichment$n_cr1_S_2[i]) - (stage1_data$n_cr0_S_1[i] + stage2_data_no_enrichment$n_cr0_S_2[i]) ) / (stage1_data$n_S_1[i] + stage2_data_no_enrichment$n_S_2[i])
    
    trial_result$effect_F_overall_SD[i] <- 
      ( (stage1_data$n_cr1_F_1[i] + stage2_data_no_enrichment$n_cr1_F_2[i]) - (stage1_data$n_cr0_F_1[i] + stage2_data_no_enrichment$n_cr0_F_2[i]) ) / (stage1_data$n_F_1[i] + stage2_data_no_enrichment$n_F_2[i])    
  } # end of if ((stage1_data$pval_intersect_1_SD[i]>sig_1)&(delta_T_interim>=d_T)&(delta_S_interim>=d_S))
  
} # end of for (i in 1:nSimul)

# Sime
trial_result$signif_F_or_S  <- trial_result$signif_F|trial_result$signif_S
trial_result$signif_F_and_S <- trial_result$signif_F&trial_result$signif_S

trial_result$signif_clinrev_F <- trial_result$signif_F & !is.na(trial_result$effect_F_overall) & (trial_result$effect_F_overall >= cli_rev_F)

trial_result$signif_clinrev_S <- trial_result$signif_S & !is.na(trial_result$effect_S_overall) & (trial_result$effect_S_overall >= cli_rev_S)

trial_result$signif_clinrev_F_or_S <- trial_result$signif_clinrev_F | trial_result$signif_clinrev_S
trial_result$signif_clinrev_F_and_S <- trial_result$signif_clinrev_F & trial_result$signif_clinrev_S

trial_result$stage1_decision <- factor(trial_result$stage1_decision)

# SD
trial_result$signif_F_or_S_SD  <- trial_result$signif_F_SD|trial_result$signif_S_SD
trial_result$signif_F_and_S_SD <- trial_result$signif_F_SD&trial_result$signif_S_SD

trial_result$signif_clinrev_F_SD <- trial_result$signif_F_SD & !is.na(trial_result$effect_F_overall_SD) & (trial_result$effect_F_overall_SD >= cli_rev_F)

trial_result$signif_clinrev_S_SD <- trial_result$signif_S_SD & !is.na(trial_result$effect_S_overall_SD) & (trial_result$effect_S_overall_SD >= cli_rev_S)

trial_result$signif_clinrev_F_or_S_SD <- trial_result$signif_clinrev_F_SD | trial_result$signif_clinrev_S_SD
trial_result$signif_clinrev_F_and_S_SD <- trial_result$signif_clinrev_F_SD & trial_result$signif_clinrev_S_SD

trial_result$stage1_decision_SD <- factor(trial_result$stage1_decision_SD)


trial_result$cp_stage2_S_half_alp_co <- stage1_data$cp_stage2_S_half_alp_co 

trial_result$cp_stage2_S_full_alp_en <- stage1_data$cp_stage2_S_full_alp_en

trial_result$cp_stage2_F_half_alp_co <- stage1_data$cp_stage2_F_half_alp_co 

trial_result$cp_stage2_F_full_alp_ac <- stage1_data$cp_stage2_F_full_alp_ac 

trial_result$sig_comb <- sig_comb


trial_result$est_ctr_rate_F_2 <- n_cr0_F_2 / n_F_2
trial_result$est_ctr_rate_S_2 <- n_cr0_S_2 / n_S_2

return(trial_result)
  }, error=function(cond) browser())
}
#===============================================================================================================

cp_stage2 <- function(alp2=.025, n1, n2, w1, w2, z1) {
  # conditional prob of success at stage 2 given stage 1 result
  # assuming w1 = sqrt(n1/(n1+n2)) and w2 = sqrt(n2/(n1+n2))
  
  # alp2    one-sided alpha level for stage 2
  # n1      sample size at stage 1
  # n2      sample size at stage 2 (the addition)
  # z1      observed z-score at stage 1
  
  # ref formula 6 in https://onlinelibrary.wiley.com/doi/full/10.1002/sim.4102
  
  # # test code
  # alp2 <- 0.0184; n1=205*.443; n2=120*.443; z1=obs_rate_diff_2_z_score(delta=.15, pa=.3, na=102*.443, nb=103*.443)
  
  # alp2 <- 0.0184; n1=205*.443; n2=120*.443; z1=obs_rate_diff_2_z_score(delta=.12, pa=.35, na=102*.443, nb=103*.443)
  # alp2 <- 0.0184; n1=205; n2=120; z1=obs_rate_diff_2_z_score(delta=.12*.443+.1*(1-.443), pa=.48, na=102*.443, nb=103*.443)
  
  n2sq  <- sqrt(n2)
  n1sq  <- sqrt(n1)
  # N2sq  <- sqrt(n2 + n1)
  # 1 - pnorm( qnorm(1-alp2) * N2sq / n2sq - z1 * (n1sq/n2sq + n2sq/n1sq)) # last term mean assuming alternative exactly the same as observed z1 (stage 1 data) 
  
  1 - pnorm( qnorm(1-alp2)/w2 - z1 * (w1/w2 + n2sq/n1sq)) # from Wassmer's book to account for flexible w1,w2
  
} # end of cp_stage2
#===============================================================================================================

