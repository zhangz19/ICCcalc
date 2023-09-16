
/*
Walter Reinisch, Vivek Pradhan, Saira Ahmad, Zhen Zhang, Jeremy D Gale, 
Alternative endoscopy reading paradigms determine score reliability and effect size in ulcerative colitis, 
Journal of Crohn's and Colitis, 2023;, jjad134, https://doi.org/10.1093/ecco-jcc/jjad134
*/

proc import datafile='yourfilepath/foo.csv' out=dsin;
run;

/* Normal  */
ods output AdditionalEstimates=icc1 FitStatistics=fit1;  /*ParameterEstimates*/
proc nlmixed data=dsin; 
  by window; 
  parms mu=0 verror=1 vpat=1;
  pred=mu+beta;  
  model mcsflxn~ normal(pred, verror); 
  random beta ~normal(0, vpat) subject=subjid; 
  estimate 'icc' vpat/(vpat + verror); 
run; 
proc export data=icc1  replace
outfile='yourfilepath/icc1.csv'  
dbms=csv;
run;
proc export data=fit1  replace
outfile='yourfilepath/fit1.csv'  
dbms=csv;
run;


/* Normal with covariates  */
ods output AdditionalEstimates=icc2 FitStatistics=fit2; 
proc nlmixed data=dsin; 
  by window; 
  parms mu=0 b1=0 verror=1 vpat=1;
  pred=mu + b1*ANTITNFN + beta;  
  model mcsflxn~ normal(pred, verror); 
  random beta ~normal(0, vpat) subject=subjid; 
  estimate 'icc' vpat/(vpat + verror); 
run; 
proc export data=icc2  replace
outfile='yourfilepath/icc2.csv'  
dbms=csv;
run;
proc export data=fit2  replace
outfile='yourfilepath/fit2.csv'  
dbms=csv;
run;


/*  Ordinal cumulative logit */
ods output AdditionalEstimates=icc3 FitStatistics=fit3;  /*ParameterEstimates*/
PROC NLMIXED DATA=dsin QPOINTS=21;
  by window; 
PARMS vpat=1 mu=0 g1=-6 g2=-3 g3=-.7;
pred=mu+beta;
IF (mcsflxn=0) THEN p = 1 / (1 + EXP(-(g1-pred)));
ELSE IF (mcsflxn=1) THEN p = (1/(1 + EXP(-(g2-pred)))) - (1/(1 + EXP(-(g1-pred))));
ELSE IF (mcsflxn=2) THEN p = (1/(1 + EXP(-(g3-pred)))) - (1/(1 + EXP(-(g2-pred))));
ELSE IF (mcsflxn=3) THEN p = 1 - (1 / (1 + EXP(-(g3-pred))));
loglike = LOG(p);
MODEL mcsflxn ~ GENERAL(loglike);
random beta ~normal(0, vpat) subject=subjid; 
ESTIMATE 'icc' vpat/((((ATAN(1)*4)**2)/3)+vpat); 
RUN;
proc export data=icc3  replace
outfile='yourfilepath/icc3.csv'  
dbms=csv;
run;
proc export data=fit3  replace
outfile='yourfilepath/fit3.csv'  
dbms=csv;
run;


/*  Ordinal cumulative logit with covariates    */
ods output AdditionalEstimates=icc4 FitStatistics=fit4;  /*ParameterEstimates*/
PROC NLMIXED DATA=dsin QPOINTS=21;
  by window; 
PARMS vpat=1 mu=0 b1=0 g1=-6 g2=-3 g3=-.7;
pred=mu + b1*ANTITNFN  + beta;
IF (mcsflxn=0) THEN p = 1 / (1 + EXP(-(g1-pred)));
ELSE IF (mcsflxn=1) THEN p = (1/(1 + EXP(-(g2-pred)))) - (1/(1 + EXP(-(g1-pred))));
ELSE IF (mcsflxn=2) THEN p = (1/(1 + EXP(-(g3-pred)))) - (1/(1 + EXP(-(g2-pred))));
ELSE IF (mcsflxn=3) THEN p = 1 - (1 / (1 + EXP(-(g3-pred))));
loglike = LOG(p);
MODEL mcsflxn ~ GENERAL(loglike);
random beta ~normal(0, vpat) subject=subjid; 
ESTIMATE 'icc' vpat/((((ATAN(1)*4)**2)/3)+vpat); 
RUN;
proc export data=icc4  replace
outfile='yourfilepath/icc4.csv'  
dbms=csv;
run;
proc export data=fit4  replace
outfile='yourfilepath/fit4.csv'  
dbms=csv;
run;


/*   Ordinal cumulative logit with scale parameter: maybe identifiability issue    */
/*
ods output AdditionalEstimates=icc5 FitStatistics=fit5;  /*ParameterEstimates*/
/*
PROC NLMIXED DATA=dsin QPOINTS=21;
  by window; 
PARMS vpat=1 s=1 mu=0 g1=-6 g2=-3 g3=-.7;
pred=mu+beta;
IF (mcsflxn=0) THEN p = 1 / (1 + EXP(-(g1-pred)/s));
ELSE IF (mcsflxn=1) THEN p = (1/(1 + EXP(-(g2-pred)/s))) - (1/(1 + EXP(-(g1-pred)/s)));
ELSE IF (mcsflxn=2) THEN p = (1/(1 + EXP(-(g3-pred)/s))) - (1/(1 + EXP(-(g2-pred)/s)));
ELSE IF (mcsflxn=3) THEN p = 1 - (1 / (1 + EXP(-(g3-pred)/s)));
loglike = LOG(p);
MODEL mcsflxn ~ GENERAL(loglike);
random beta ~normal(0, vpat) subject=subjid; 
ESTIMATE 'icc' vpat/((((ATAN(1)*4*s)**2)/3)+vpat); 
RUN;
proc export data=icc5  replace
outfile='yourfilepath/icc5.csv'  
dbms=csv;
run;
proc export data=fit5  replace
outfile='yourfilepath/fit5.csv'  
dbms=csv;
run;
*/

/*   Ordinal cumulative logit with scale parameter and covariates   */
/*
ods output AdditionalEstimates=icc6 FitStatistics=fit6;  /*ParameterEstimates*/
/*
PROC NLMIXED DATA=dsin QPOINTS=21;
  by window; 
PARMS vpat=1 s=1 mu=0 b1=0 g1=-6 g2=-3 g3=-.7;
pred=mu + b1*ANTITNFN + beta;
IF (mcsflxn=0) THEN p = 1 / (1 + EXP(-(g1-pred)/s));
ELSE IF (mcsflxn=1) THEN p = (1/(1 + EXP(-(g2-pred)/s))) - (1/(1 + EXP(-(g1-pred)/s)));
ELSE IF (mcsflxn=2) THEN p = (1/(1 + EXP(-(g3-pred)/s))) - (1/(1 + EXP(-(g2-pred)/s)));
ELSE IF (mcsflxn=3) THEN p = 1 - (1 / (1 + EXP(-(g3-pred)/s)));
loglike = LOG(p);
MODEL mcsflxn ~ GENERAL(loglike);
random beta ~normal(0, vpat) subject=subjid; 
ESTIMATE 'icc' vpat/((((ATAN(1)*4*s)**2)/3)+vpat); 
RUN;
proc export data=icc6  replace
outfile='yourfilepath/icc6.csv'  
dbms=csv;
run;
proc export data=fit6  replace
outfile='yourfilepath/fit6.csv'  
dbms=csv;
run;
*/


/*  truncated Poisson/NB  */
/*
PROC NLMIXED data=dsin QPOINTS=21;
  by window; 
parms vpat=1 mu=0;
loglambda = mu + beta; 
lambda = exp(loglambda);
lambda2 = lambda*lambda; 
lambda3 = lambda2*lambda; 
ll = mcsflxn*loglambda - lambda - lgamma(mcsflxn + 1) - log( exp(-lambda) * (1 + lambda + lambda2/2 + lambda3/6)  );
model mcsflxn ~ general(ll);
 random beta ~ normal(0, vpat) subject=subjid;
RUN;

PROC NLMIXED data=dsin QPOINTS=11;
  by window; 
parms vpat=1 mu=0 alpha=0.1;
pred = exp(mu + beta);
m = 1/alpha;
apred = alpha*pred; 
ll = lgamma(mcsflxn+m)-lgamma(mcsflxn+1)-lgamma(m) + mcsflxn*log(apred)-(mcsflxn+m)*log(1+apred);
model mcsflxn ~ general(ll);
 random beta ~ normal(0, vpat) subject=subjid;
RUN;


PROC NLMIXED data=dsin QPOINTS=11;
  by window; 
parms vpat=1 mu=0 alpha=0.1;
pred = exp(mu + beta);
m = 1/alpha;
apred = alpha*pred; 
p0 = ( 1 + apred)**(-m); 
p1 = m*apred* ( (1+apred)**(-(m+1)) ); 
p2 = lgamma(2+m)-lgamma(2+1)-lgamma(m) * (apred**2) * ( (1+apred)**(-(m+2)) ); 
p3 = lgamma(3+m)-lgamma(3+1)-lgamma(m) * (apred**3) * ( (1+apred)**(-(m+3)) ); 
ll = lgamma(mcsflxn+m)-lgamma(mcsflxn+1)-lgamma(m) + mcsflxn*log(apred)-(mcsflxn+m)*log(1+apred) - log( p0+p1+p2+p3 );
   model mcsflxn ~ general(ll);
 random beta ~ normal(0, vpat) subject=subjid;
RUN;
*/


/*
proc nlmixed data=dsin; 
  by window; 
  parms mu=0 verror=1 vpat=1;
  pred=mu+beta;  
  model mcsflxn~ poisson(pred, verror); 
  random beta ~normal(0, vpat) subject=subjid; 
  estimate 'icc' vpat/(vpat + verror); 
run; 

PROC NLMIXED data=dsin;
parms vpat=1 beta0=0;
 lambda = exp(beta0 + u);
 model mcsflxn ~ poisson(lambda);
 random u ~ normal(0,vpat) subject=subjid;
RUN;


PROC NLMIXED DATA=dsin QPOINTS=11;
PARMS vpat=1 mu=0 g1=-6 g2=-3 g3=-.7;
pred=mu+beta;
IF (mcsflxn=0) THEN p = 1 / (1 + EXP(-(g1-pred)));
ELSE IF (mcsflxn=1) THEN p = (1/(1 + EXP(-(g2-pred)))) - (1/(1 + EXP(-(g1-pred))));
ELSE IF (mcsflxn=2) THEN p = (1/(1 + EXP(-(g3-pred)))) - (1/(1 + EXP(-(g2-pred))));
ELSE IF (mcsflxn=3) THEN p = 1 - (1 / (1 + EXP(-(g3-pred))));
loglike = LOG(p);
MODEL mcsflxn ~ GENERAL(loglike);
random beta ~normal(0, vpat) subject=subjid; 
RUN;

*/
