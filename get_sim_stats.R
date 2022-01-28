load("/home/kweave23/gamete_data/gen_model_results_noDNM/g1000_s30000_c0.01_se0.005_r1/assess_out_rs_1848.Rdata")
ao1848 <- assess_out
load("/home/kweave23/gamete_data/gen_model_results_noDNM/g1000_s30000_c0.01_se0.005_r1/assess_out_rs_357.Rdata")
ao357 <- assess_out
load("/home/kweave23/gamete_data/gen_model_results_noDNM/g1000_s30000_c0.01_se0.005_r1/assess_out_rs_42.Rdata")
ao42 <- assess_out

# Phasing
## accuracy
mean(c(ao1848$phasing$acc, ao42$phasing$acc, ao357$phasing$acc)) #>99.99325
sd(c(ao1848$phasing$acc, ao42$phasing$acc, ao357$phasing$acc)) #>0.003375583

## completeness
mean(c(ao1848$phasing$com, ao42$phasing$com, ao357$phasing$com)) #>0.9995835
sd(c(ao1848$phasing$com, ao42$phasing$com, ao357$phasing$com)) #>0.0002533294

## switch error rate
mean(c(ao1848$phasing$ser, ao42$phasing$ser, ao357$phasing$ser)) #>6.753981e-05
sd(c(ao1848$phasing$ser, ao42$phasing$ser, ao357$phasing$ser)) #>3.375583e-05

## largest haplotype segment
mean(c(ao1848$phasing$lhs, ao42$phasing$lhs, ao357$phasing$lhs)) #> 10093.33
sd(c(ao1848$phasing$lhs, ao42$phasing$lhs, ao357$phasing$lhs)) #> 4871.999

# Imputation
## accuracy
mean(c(mean(ao1848$gam_imputation$acc), mean(ao42$gam_imputation$acc), mean(ao357$gam_imputation$acc))) #> 99.96159
sd(c(mean(ao1848$gam_imputation$acc), mean(ao42$gam_imputation$acc), mean(ao357$gam_imputation$acc))) #> 0.002440232

## completeness
mean(c(mean(ao1848$gam_imputation$com), mean(ao42$gam_imputation$com), mean(ao357$gam_imputation$com))) #> 0.9933872
sd(c(mean(ao1848$gam_imputation$com), mean(ao42$gam_imputation$com), mean(ao357$gam_imputation$com))) #> 0.0003731026

## switch error rate
mean(c(mean(ao1848$gam_imputation$ser), mean(ao42$gam_imputation$ser), mean(ao357$gam_imputation$ser))) #> 0.0001190001
sd(c(mean(ao1848$gam_imputation$ser), mean(ao42$gam_imputation$ser), mean(ao357$gam_imputation$ser))) #> 3.311384e-05

## largest haplotype segement
mean(c(mean(ao1848$gam_imputation$lhs), mean(ao42$gam_imputation$lhs), mean(ao357$gam_imputation$lhs))) #> 8247.589
sd(c(mean(ao1848$gam_imputation$lhs), mean(ao42$gam_imputation$lhs), mean(ao357$gam_imputation$lhs))) #> 2965.562

# Discovery
## precision
mean(c(ao1848$recomb$precision, ao42$recomb$precision, ao357$recomb$precision)) #> 0.9823685
sd(c(ao1848$recomb$precision, ao42$recomb$precision, ao357$recomb$precision)) #> 0.002915218

## recall
mean(c(ao1848$recomb$recall, ao42$recomb$recall, ao357$recomb$recall)) #> 0.9369934
sd(c(ao1848$recomb$recall, ao42$recomb$recall, ao357$recomb$recall)) #> 0.003158956

## accuracy
mean(c(ao1848$recomb$accuracy, ao42$recomb$accuracy, ao357$recomb$accuracy)) #> 0.9429396
sd(c(ao1848$recomb$accuracy, ao42$recomb$accuracy, ao357$recomb$accuracy)) #> 0.003574476

## specificity
mean(c(ao1848$recomb$specificity, ao42$recomb$specificity, ao357$recomb$specificity)) #> 0.9579396
sd(c(ao1848$recomb$specificity, ao42$recomb$specificity, ao357$recomb$specificity)) #> 0.00526997

## f1 score
mean(c(ao1848$recomb$f1, ao42$recomb$f1, ao357$recomb$f1)) #> 0.9591434
sd(c(ao1848$recomb$f1, ao42$recomb$f1, ao357$recomb$f1)) #> 0.002744348

## false discovery rate
mean(c(ao1848$recomb$fdr, ao42$recomb$fdr, ao357$recomb$fdr)) #> 0.01763148
sd(c(ao1848$recomb$fdr, ao42$recomb$fdr, ao357$recomb$fdr)) #> 0.002915218

## false positive rate
mean(c(ao1848$recomb$fpr, ao42$recomb$fpr, ao357$recomb$fpr)) #> 0.04206039
sd(c(ao1848$recomb$fpr, ao42$recomb$fpr, ao357$recomb$fpr)) #> 0.00526997