##################################################################################################################################################################
# OBJECTIVE(s): To simulate the proposed adaptive design under a variety of conditions in order to evaluate its operating characteristics and its effectiveness 
# AUTHOR(s)   : Anh Nguyen Duc based on simpler prototype by Marcel Wolbers and Dominik Heinzmann
# NOTE(s)     : Decision rules @ end of stage 1 are based on the treatment effect (difference in proportion between experimental group and control group) in S & T
# Last modified: 12 Mar 2020 by Anh Nguyen Duc
##################################################################################################################################################################
require(magrittr)
require(rpact)
##################################################################################################################################################################

# IMPORTANT: all sample size is PER ARM with RR=1:1!!!

AdaptiveTrialDesign <- function(p0_S=0.48,              # Probability of response in subgroup S in the control arm
                                p0_C=0.48,              # Probability of response in complement subgroup T in the control arm
                                delta_S=0.2,            # Treatment effect (difference in proportion between experimental group and control group) in S 
                                delta_C=0.2,            # Treatment effect (difference in proportion between experimental group and control group) in T
                                prev_S=0.47,            # Prevalence of subgroup S
                                n_F_1=102,              # Sample size per arm for stage 1 related to full population F
                                n_2=60,                 # Sample size per arm for stage 2 (either for full population F or, if enriched, for S depending on stage 1 decision)
                                d_C=0.10,               # Threshold in complement group T for decision rule at end of stage 1
                                d_S=0.15,               # Threshold in subgroup S for decision rule at end of stage 1
                                sig_level=0.025,        # Overall type I error (one-sided testing)
                                prop_alpha_spent_1=0.5, # Proportion of type I error to be spent at stage 1 (for early efficacy stop)
                                drop.out.rate=0.05,     # Drop out proportion similar in both arms (patients "dropping out" prior to surgery are counted as non-responder)
                                cli_rev_S=0.15,         # Threshold for final clinical rev in S
                                cli_rev_F=0.15,         # Threshold for final cincal rev in ITT
                                nSimul=1e4,             # Number of simulations
                                
                                stage1_data=NULL,       # In case stage 1 data have been calculated
                                w_1=NULL, w_2=NULL      # Weights for stage 1 and 2
){
  
##---------------- Simulation of a simple adaptive enrichment design with fixed sample sizes for both stages


##---------------- Parameter settings
tryCatch({
    
## These must be based on the correct & real prevalence
# response in full population F
p0_F <- prev_S*p0_S+(1-prev_S)*p0_C

# # delta in full population F
delta_F <- prev_S*delta_S+(1-prev_S)*delta_C

## One-sided significance level at stage 1 based on the alpha spent at this stage

sig_1 <- sig_level*prop_alpha_spent_1

##---------------- Simulate data for stage 1 and 2

#-----------------------------------------------------------------------------------------
#--- Stage 1
#-----------------------------------------------------------------------------------------

## These must be based on the observed prevalence
# Sample size stage 1 per population
n_S_1 <- round(prev_S*n_F_1)
n_C_1 <- n_F_1-n_S_1

# Weights p-value combination for stage 1 and stage 2
n_all <- n_F_1 + n_2
info_1<- n_F_1 / n_all

w_1 <- sqrt(info_1)
w_2 <- sqrt(1-info_1)



# Calculate nominal one-significance level for stage 1 & 2 combination test using the group sequential design methodology (gsDesign package)

sig_comb <- getDesignInverseNormal(informationRates = c(w_1^2,1), typeOfDesign ="asUser",
                                   userAlphaSpending = c(sig_1,sig_level),alpha=sig_level)$stageLevels[2] 

if (is.null(stage1_data)) {

  # number of responders (cr0 is control, cr1 is intervention)
  # drop-outs incorporated assuming given p0_s/T are proportions in target population; drop-out assumed prior 
  # surgery (hypothetical responders drop-outs are non-responder by definition) 
  n_cr0_S_1 <- rbinom(nSimul,size=n_S_1,prob=p0_S*(1-drop.out.rate)) 
  n_cr1_S_1 <- rbinom(nSimul,size=n_S_1,prob=(p0_S+delta_S)*(1-drop.out.rate))
  n_cr0_C_1 <- rbinom(nSimul,size=n_C_1,prob=p0_C*(1-drop.out.rate))
  n_cr1_C_1 <- rbinom(nSimul,size=n_C_1,prob=(p0_C+delta_C)*(1-drop.out.rate))
  
  
  n_cr0_F_1 <- n_cr0_S_1+n_cr0_C_1 # Number of pCRs for full pop after stage 1 in control arm
  n_cr1_F_1 <- n_cr1_S_1+n_cr1_C_1 # Number of pCRs for full pop after stage 1 in experimental arm
  n_cr_F_1  <- n_cr0_F_1+n_cr1_F_1 # Number of pCRs for full pop after stage 1 in both arms
  
  # estimated delta and p-value in full pop and subgroup from stage 1
  stage1_data <- data.frame(est_delta_F_1=(n_cr1_F_1-n_cr0_F_1)/n_F_1, 
                            est_delta_S_1=(n_cr1_S_1-n_cr0_S_1)/n_S_1,
                            est_delta_C_1=(n_cr1_C_1-n_cr0_C_1)/n_C_1,
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


#-----------------------------------------------------------------------------------------
#--- Stage 2 data in case trial continues in overall population F
#-----------------------------------------------------------------------------------------

# browser()
# Sample sizes stage 2
n_F_2 <- n_2
n_S_2 <- round(prev_S*n_F_2) # assuming same prevalence as for stage 1 (take estimate from stage 1 know prior to interim analysis)
n_C_2 <- n_F_2-n_S_2

# number of responders (cr0 is control, cr1 is intervention)
n_cr0_S_2 <- rbinom(nSimul,size=n_S_2,prob=p0_S*(1-drop.out.rate))
n_cr1_S_2 <- rbinom(nSimul,size=n_S_2,prob=(p0_S+delta_S)*(1-drop.out.rate))
n_cr0_C_2 <- rbinom(nSimul,size=n_C_2,prob=p0_C*(1-drop.out.rate))
n_cr1_C_2 <- rbinom(nSimul,size=n_C_2,prob=(p0_C+delta_C)*(1-drop.out.rate))

n_cr0_F_2 <- n_cr0_S_2+n_cr0_C_2
n_cr1_F_2 <- n_cr1_S_2+n_cr1_C_2

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
                           effect_C_1=stage1_data$est_delta_C_1, n_cr_F_1=stage1_data$n_cr_F_1,
                           pval_F_1=stage1_data$pval_F_1,
                           pval_S_1=stage1_data$pval_S_1,
                           pval_intersect_1=stage1_data$pval_intersect_1,
                           
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
  

  # Decision to stop for futility at the end of stage 1

  delta_C_interim <- stage1_data$est_delta_C_1[i]
  delta_S_interim <- stage1_data$est_delta_S_1[i]
  delta_F_interim <- stage1_data$est_delta_F_1[i]

  # Sime
  if ((stage1_data$pval_intersect_1[i]>sig_1) & (delta_C_interim<d_C) & (delta_S_interim<d_S)) {
    trial_result$stage1_decision[i] <- trial_result$stage1_decision_detail[i] <- "2 - Stop after stage 1: futility"
    trial_result$signif_F[i] <- F
    trial_result$signif_S[i] <- F
    
    trial_result$effect_F_overall[i] <- stage1_data$est_delta_F_1[i]
    trial_result$effect_S_overall[i] <- stage1_data$est_delta_S_1[i]
  }
  

  # Maximum observed treatement effect in S or F
  max_eff <- max(delta_S_interim,delta_F_interim)

  # Decision to continue to stage 2 with S only
  # Sime
  if ((stage1_data$pval_intersect_1[i]>sig_1) & (delta_C_interim<d_C) & (delta_S_interim>=d_S)) {
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


  # Decision to continue to stage 2 with F only
  # Sime
  if ((stage1_data$pval_intersect_1[i]>sig_1) & (delta_C_interim>=d_C)&(delta_S_interim<d_S)) {
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

  
  # Decision to continue to stage 2 with F and S
  # Sime
  if ((stage1_data$pval_intersect_1[i]>sig_1)&(delta_C_interim>=d_C)&(delta_S_interim>=d_S)) {
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
  } # end of if ((stage1_data$pval_intersect_1[i]>sig_1)&(delta_C_interim>=d_C)&(delta_S_interim>=d_S))
  

  
} # end of for (i in 1:nSimul)

# Sime
trial_result$signif_F_or_S  <- trial_result$signif_F|trial_result$signif_S
trial_result$signif_F_and_S <- trial_result$signif_F&trial_result$signif_S

trial_result$signif_clinrev_F <- trial_result$signif_F & !is.na(trial_result$effect_F_overall) & (trial_result$effect_F_overall >= cli_rev_F)

trial_result$signif_clinrev_S <- trial_result$signif_S & !is.na(trial_result$effect_S_overall) & (trial_result$effect_S_overall >= cli_rev_S)

trial_result$signif_clinrev_F_or_S <- trial_result$signif_clinrev_F | trial_result$signif_clinrev_S
trial_result$signif_clinrev_F_and_S <- trial_result$signif_clinrev_F & trial_result$signif_clinrev_S

trial_result$stage1_decision <- factor(trial_result$stage1_decision)


trial_result$sig_comb <- sig_comb


trial_result$est_ctr_rate_F_2 <- n_cr0_F_2 / n_F_2
trial_result$est_ctr_rate_S_2 <- n_cr0_S_2 / n_S_2

return(trial_result)
  }, error=function(cond) browser())
}
#===============================================================================================================


