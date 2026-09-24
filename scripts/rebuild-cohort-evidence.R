#!/usr/bin/env Rscript
# Rebuild aggregate cohort evidence from controlled local inputs.
# No participant-level data is written or included in the web repository.
suppressPackageStartupMessages({library(survival); library(jsonlite); library(openxlsx)})
options(stringsAsFactors=FALSE)
args <- commandArgs(FALSE)
entry <- sub('^--file=', '', args[grepl('^--file=', args)][1])
repo <- normalizePath(file.path(dirname(entry), '..'))
root <- Sys.getenv('DISCO_DATA_ROOT')
questionnaire <- Sys.getenv('DISCO_RULAS_QUESTIONNAIRE')
if (!nzchar(root) || !nzchar(questionnaire)) stop('Set DISCO_DATA_ROOT and DISCO_RULAS_QUESTIONNAIRE.')
out <- Sys.getenv('DISCO_EVIDENCE_OUTPUT', file.path(repo, 'data', 'cohort-evidence.json'))
line <- grep('^const STANDARD_REF = ', readLines(file.path(repo,'index.html'),warn=FALSE), value=TRUE)
ref <- fromJSON(sub(';$','',sub('^const STANDARD_REF = ','',line)))
markers <- ref$markers
source_paths <- character()
read_rds <- function(p) {source_paths <<- unique(c(source_paths,p));readRDS(file.path(root,p))}
read_csv <- function(p,...) {source_paths <<- unique(c(source_paths,p));read.csv(file.path(root,p),...)}
num <- function(x) as.numeric(as.character(x))
factor_codes <- function(x, valid) {v<-num(x);factor(ifelse(v %in% valid,v,NA))}
binary <- function(x, yes, no) factor(ifelse(x %in% yes,1,ifelse(x %in% no,0,NA)),levels=c(0,1))
checked_match <- function(ids, reference_ids) {
 if(anyDuplicated(reference_ids)) stop('Duplicate join keys in reference input.')
 match(as.character(ids),as.character(reference_ids))
}

# The same aggregate-reference equations as the calculator; subset the matrices
# for cohorts lacking panel markers. Do not renormalize the subset weights.
score_panel <- function(values, panel) {
 ii <- match(panel,markers)
 y <- as.matrix(values[,panel,drop=FALSE]);storage.mode(y)<-'double'
 valid <- rowSums(!is.finite(y) | y<0)==0
 out <- data.frame(DISCO=rep(NA_real_,nrow(y)),DM=rep(NA_real_,nrow(y)))
 if(!any(valid))stop('No valid biomarker rows.')
 y <- y[valid,,drop=FALSE];y[,'CRP']<-log(y[,'CRP']+1)
 z <- sweep(sweep(y,2,ref$ref_mean[ii],'-'),2,ref$ref_sd[ii],'/')
 delta <- sweep(z,2,ref$mu[ii],'-');ss<-ref$ss[ii,ii,drop=FALSE]
 fac <- ref$n_ref/(ref$n_ref+1)
 dg <- sweep(fac*delta^2,2,diag(ss),'+')
 acc <- numeric(nrow(delta))
 for(j in seq_len(ncol(delta)-1))for(k in (j+1):ncol(delta)) {
  rnew<-(ss[j,k]+fac*delta[,j]*delta[,k])/sqrt(dg[,j]*dg[,k])
  acc<-acc+2*ref$weight[ii[j],ii[k]]*(ref$Rref[ii[j],ii[k]]-rnew)^2
 }
 disco<-log(acc*ref$n_ref^2)
 S<-ref$S[ii,ii,drop=FALSE]
 dm<-sqrt(rowSums((delta %*% solve(S))*delta))
 out$DISCO[valid]<-disco;out$DM[valid]<-dm
 # Verify the vectorized update against literal per-person matrix updates.
 for(i in unique(round(seq(1,nrow(delta),length.out=min(5,nrow(delta)))))) {
  newss<-ss+fac*tcrossprod(delta[i,]);sn<-sqrt(diag(newss))
  direct<-log(sum(ref$weight[ii,ii]*(ref$Rref[ii,ii]-newss/tcrossprod(sn))^2)*ref$n_ref^2)
  stopifnot(isTRUE(all.equal(disco[i],direct,tolerance=1e-9)))
 }
 out
}

cohorts <- list()
# UK Biobank: fixed UKB reference, matching all 10 calculator inputs.
bio<-read_rds('Result1/Clean.UKB.Clinical.rds')
a<-read_rds('UKB_Data/extended_covariates_clinical.rds')
bio<-bio[checked_match(a$ID,bio$sampleID),]
rawcov<-read_rds('UKB_Data/covariates_raw_full.rds')
edu<-num(rawcov$edu_q0[checked_match(a$ID,rawcov$eid)])
race<-read_csv('UKB_Data/ukb_race_21000.csv',colClasses='character')
race<-race$race[checked_match(a$ID,race$eid)];race[!race %in% c('White','Non-White')]<-NA
values<-data.frame(CRP=bio$CRP,GLU=bio$GLU,CHOL=bio$TC,HDL=bio$HDL.C,CREA=bio$CREA,UA=bio$GML,WBC=bio$wbc,UREA=bio$U,RBC=bio$RBC,ALB=bio$ALB)
d<-data.frame(time=num(a$time),status=num(a$status),age=num(a$age),sex=factor(a$sex),bmi=num(a$bmi_raw),race=factor(race),education=binary(edu,c(1,2),c(-7,3,4,5,6)),smoking=factor_codes(a$smoke_status,0:2),alcohol=factor_codes(a$alcohol_freq,1:6),physical_activity=factor_codes(a$ipaq_group,0:2))
d$bmi[d$bmi<=0]<-NA
cohorts$UKB<-list(label='UK Biobank',wave='Baseline 2006–2010',values=values,d=d,panel=markers,covariates=c('age','sex','bmi','race','education','smoking','alcohol','physical_activity'))
rm(bio,a,rawcov,values,d);invisible(gc(FALSE))

# NHANES: recover raw CRP mg/L before applying ln(CRP + 1).
source_paths<-c(source_paths,'Package/DISCO4HD/data/NHANES4.rda')
load(file.path(root,'Package/DISCO4HD/data/NHANES4.rda'))
a<-NHANES4
cov<-read_rds('UKB_Data/NHANES_covariates.rds');cov<-cov[checked_match(a$sampleID,cov$sampleID),]
values<-data.frame(CRP=(exp(a$lncrp)-1)*10,GLU=a$glucose_mmol,CHOL=a$totchol*.02586,HDL=a$hdl*.02586,CREA=a$creat*88.4,UA=a$uap*59.5,WBC=a$wbc,UREA=a$bun*.357,RBC=a$rbc,ALB=a$albumin*10)
d<-data.frame(time=num(a$time),status=num(a$status),age=num(a$age),sex=factor_codes(cov$sex,c(1,2)),bmi=num(a$bmi),race=factor(cov$race,levels=c('White','Non-White')),education=factor(cov$education,levels=c('low','high')),smoking=factor_codes(cov$smoke,0:1),alcohol=factor_codes(cov$drink,0:1))
d$bmi[d$bmi<=0]<-NA
cohorts$NHANES<-list(label='NHANES',wave='Harmonized NHANES4 sample',values=values,d=d,panel=markers,covariates=c('age','sex','bmi','race','education','smoking','alcohol'))
rm(NHANES4,a,cov,values,d);invisible(gc(FALSE))

# CHARLS baseline: RBC and ALB absent from the analyzed panel.
a<-read_rds('Result4/CHARLS.DISCO.rds')
source_paths<-c(source_paths,'Result4/CHARLS/Blood_Biomarkers.txt')
bio<-read.delim(file.path(root,'Result4/CHARLS/Blood_Biomarkers.txt'))
bio<-bio[checked_match(a$id,bio$ID),]
ckm<-read_csv('Result4/CHARLS/CKM.Charls.stat2011.csv')
bmi<-num(ckm$BMI[checked_match(a$id,ckm$id)]);bmi[bmi<10|bmi>60]<-NA
values<-data.frame(CRP=bio$CRP,GLU=bio$Glu/18.0182,CHOL=bio$CHOL*.02586,HDL=bio$HDL*.02586,CREA=bio$Crea*88.4,UA=bio$UA*59.5,WBC=bio$WBC,UREA=bio$Bun*.357)
edu<-suppressWarnings(num(a$edu))
d<-data.frame(time=num(a$time),status=num(a$status),age=num(a$age),sex=factor_codes(a$sex,c(1,2)),bmi=bmi,education=binary(edu,4:11,1:3),smoking=binary(a$smoke,1,2),alcohol=binary(a$drink,1,2))
cohorts$CHARLS<-list(label='CHARLS',wave='Baseline 2011',values=values,d=d,panel=setdiff(markers,c('RBC','ALB')),covariates=c('age','sex','bmi','education','smoking','alcohol'))
rm(a,bio,ckm,values,d);invisible(gc(FALSE))

# CLHLS: join covariates and outcomes by ID to the biomarker rows. Verify the
# historical survival object matches the reconstructed age/sex/status ordering.
bio<-read_csv('Result4/CLHLS/New/biomarker_dataset_CLHLS_2014-1.tab.csv')
a<-read_rds('Result4/CLHLS.DISCO.rds')
lon<-read_rds('Result4/CLHLS/New/CLHLS.2014_2018.longitudinal.datasets.rds')
idx<-checked_match(bio$id,lon$id);stopifnot(!anyNA(idx));lon<-lon[idx,]
stopifnot(isTRUE(all.equal(num(a$age),num(lon$trueage))),isTRUE(all.equal(num(a$sex),num(lon$a1))),isTRUE(all.equal(num(a$status),num(lon$dth14_18))))
edu<-num(lon$f1);edu[edu %in% c(88,99)|edu<0]<-NA
weight<-num(lon$g101);height<-num(lon$g1021);weight[weight %in% c(999,888)]<-NA;height[height %in% c(999,888)]<-NA
bmi<-weight/(height/100)^2;bmi[bmi<15|bmi>50]<-NA
values<-data.frame(CRP=bio$crphs,GLU=bio$glu,CHOL=bio$cho,HDL=bio$hdlc,CREA=bio$crea,UA=bio$ua,WBC=bio$wbc,UREA=bio$bun,RBC=bio$rbc,ALB=bio$alb)
d<-data.frame(time=num(a$time),status=num(a$status),age=num(a$age),sex=factor_codes(a$sex,c(1,2)),bmi=bmi,education=factor(ifelse(is.na(edu),NA,ifelse(edu==0,0,1))),smoking=binary(lon$d71,1,2),alcohol=binary(lon$d81,1,2))
cohorts$CLHLS<-list(label='CLHLS',wave='Biomarker wave 2014; follow-up to 2018',values=values,d=d,panel=markers,covariates=c('age','sex','bmi','education','smoking','alcohol'))
rm(a,bio,lon,values,d);invisible(gc(FALSE))

# RuLAS wave 2 (2016): UREA unavailable; questionnaire covariates matched by ID.
a<-read_rds('Result4/RuLAS.DISCO.rds')
bio<-read_csv('Result4/Rugao/2016Rugao.data.x.csv');bio<-bio[checked_match(a$ID,bio$X2016ID),]
sh<-read.xlsx(questionnaire,sheet='Sheet2');sh<-sh[checked_match(a$ID,sh[[1]]),]
stopifnot(all(c('A61','J1','K1') %in% names(sh)))
values<-data.frame(CRP=bio$CRP2016,GLU=bio$GLU2016,CHOL=bio$CHOL2016,HDL=bio$HDLC2016,CREA=bio$CREA2016,UA=bio$UA2016,WBC=bio$WBC2016,RBC=bio$RBC2016,ALB=bio$ALB2016)
bmi<-num(bio$BMI2016);bmi[bmi<=0]<-NA
edu<-suppressWarnings(num(sh$A61))
d<-data.frame(time=num(a$time),status=num(a$status),age=num(a$age),sex=factor_codes(a$sex,c(1,2)),bmi=bmi,education=binary(edu,2:7,1),smoking=binary(sh$J1,c('2','3','水烟'),'1'),alcohol=binary(sh$K1,c('2','3'),'1'))
cohorts$RuLAS<-list(label='RuLAS',wave='Wave 2, 2016',values=values,d=d,panel=setdiff(markers,'UREA'),covariates=c('age','sex','bmi','education','smoking','alcohol'))
rm(a,bio,sh,values,d);invisible(gc(FALSE))

# Full available adjustment follows the local expanded-covariate specifications.
# Use identical complete-case participants for DISCO and raw DM; no new imputation.
results<-list()
for(nm in names(cohorts)) {
 c<-cohorts[[nm]];scores<-score_panel(c$values,c$panel);d<-cbind(c$d,scores)
 fields<-c('time','status',c$covariates,'DISCO','DM')
 ok<-complete.cases(d[,fields]);numeric_fields<-fields[vapply(d[,fields],is.numeric,logical(1))]
 ok<-ok & rowSums(!is.finite(as.matrix(d[,numeric_fields])))==0 & d$status %in% c(0,1) & d$time>0
 ok[is.na(ok)]<-FALSE
 counts<-list(source_n=nrow(d),invalid_outcome_n=sum(!is.finite(d$time)|d$time<=0|!d$status %in% c(0,1),na.rm=TRUE),invalid_score_n=sum(!is.finite(d$DISCO)|!is.finite(d$DM)),complete_case_n=sum(ok))
 d<-droplevels(d[ok,,drop=FALSE]);if(nrow(d)<100 || sum(d$status)<20)stop('Insufficient complete cases for ',nm)
 estimates<-list()
 for(metric in c('DISCO','DM')) {
  metric_sd<-sd(d[[metric]]);d$metric_z<-as.numeric(scale(d[[metric]]));warnings<-character()
  fit<-withCallingHandlers(coxph(reformulate(c('metric_z',c$covariates),response='Surv(time,status)'),data=d,ties='efron'),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  if(length(warnings))stop(nm,' ',metric,' model warning: ',paste(warnings,collapse=';'))
  s<-summary(fit);b<-unname(coef(fit)['metric_z']);se<-sqrt(vcov(fit)['metric_z','metric_z']);ci<-unname(s$concordance[1]);cse<-unname(s$concordance[2])
  stopifnot(fit$n==nrow(d),fit$nevent==sum(d$status),all(is.finite(c(b,se,ci,cse))))
  ph<-cox.zph(fit,terms=TRUE);ph_p<-unname(ph$table['metric_z','p'])
  estimates[[metric]]<-list(hr=exp(b),hr_low=exp(b-1.96*se),hr_high=exp(b+1.96*se),c_index=ci,c_low=max(0,ci-1.96*cse),c_high=min(1,ci+1.96*cse),score_sd=metric_sd,ph_test_p=ph_p)
 }
 results[[nm]]<-list(cohort=nm,label=c$label,wave=c$wave,markers=c$panel,marker_count=length(c$panel),n=nrow(d),events=sum(d$status),covariates=c$covariates,selection=counts,estimates=estimates)
 cat(nm,'markers=',length(c$panel),'N=',nrow(d),'events=',sum(d$status),'DISCO C=',round(estimates$DISCO$c_index,4),'HR=',round(estimates$DISCO$hr,3),'DM C=',round(estimates$DM$c_index,4),'HR=',round(estimates$DM$hr,3),'\n')
}
result<-list(schema_version=1,generated_utc=format(Sys.time(),tz='UTC',format='%Y-%m-%dT%H:%M:%SZ'),analysis=list(outcome='All-cause mortality',model='Cox proportional hazards, Efron ties',hr_scale='Per 1 SD within the complete-case analysis sample',disco_scale='Natural log DISCO',dm_scale='Untransformed Mahalanobis distance',c_index='Harrell concordance of the full adjusted model; apparent in-sample estimate',confidence_intervals='95% Wald intervals; C-index bounds truncated to [0,1]',sampling_weights='Unweighted; survey design not modeled',missing_data='Complete cases using supplied cleaned biomarker inputs; no additional imputation. Valid binary mortality status and positive follow-up required.',reference='Fixed UK Biobank standard matrix embedded in index.html; subset matrices for 8/9-marker panels, without re-normalizing weights.',covariate_policy='Age, sex, BMI, education, smoking and alcohol in all cohorts; ethnicity/race in UKB and NHANES; physical activity in UKB.',source_notes='Availability is specific to these analyzed waves and source files; it is not a claim about all waves of a cohort.'),markers=markers,cohorts=unname(results),sources=list(controlled_data_files=source_paths,rulas_questionnaire='Wave 2 questionnaire, Sheet2: A61, J1, K1',upstream_methods=c('19_rebuild_logcrp_plus1.R','01_covariate_models.R','06_crosscohort_covariates.R','run_DISCO_DM_table.R'),software=list(R=as.character(getRversion()),survival=as.character(packageVersion('survival')))))
dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
write_json(result,out,pretty=TRUE,auto_unbox=TRUE,digits=12,na='null')
cat('Aggregate evidence saved:',out,'\n')
