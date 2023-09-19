// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <gsl/gsl_integration.h>


using namespace Rcpp;


double A2simple(double u, void * params){
 double* ptr = (double *) params;  // change the parameter vector from type void
                                   // to type double
 double a = *ptr;     // ptr[0]
 double b = *(ptr+1); // ptr[1]
 double bprime = *(ptr+2); // ptr[2]
 double c = *(ptr+3); // ptr[3]
 double x = *(ptr+4); // ptr[4]
 double y = *(ptr+5); // ptr[5]
 double M = *(ptr+6); // ptr[6] 
                      // exp(-M+log(f(u)))
                      // M has been added to avoid over/under-flow issues
 return exp(-M+(a-1)*log(u)+(c-a-1)*log(1-u)-b*log(1-u*x)-bprime*log(1-u*y));
} 


// [[Rcpp::export]]
List pepmarginallikelihoodc(arma::mat x, NumericVector y, arma::mat xnew, bool pred =false, 
                            bool intrinsic = false, double d0=0, double d1=0){
 double res2;
 double R1;
 int p = x.n_cols;    // number of explanatory variables
 double n = y.size(); // sample size
 p = p+1;             // increase p by 1 (dimension of model under consideration)
 int k0 = 1, k1 = p, nstar = n; // parameters k0 (dimension of null), k1 and n*
 double delta = n;    // delta = n for PEP
 if (intrinsic){      // if prior is intrinsic
   nstar = k1+1;      // change n*
   delta = n/nstar;   // and delta parameters
 }
 // Additional parameters
 double a1 = (nstar+d0-k1)/2;     // a1
 double b1 = (nstar+d1-d0-k1)/2;  // b1
 double a0 = (n+d0-k0)/2;         // a0
 double an1 = (n+d0-k1)/2;        // an1
 int ke = (k1-k0);   // models' dimension difference
 double res;
 double R0 = 0;      // coefficient of determination for null model
 double R;
 arma::colvec y2 = y; // response vector - as 1-column matrix
 arma::mat z(n,p);    // matrix z of size nxp - will be the design matrix
 z.ones();            // fill z with 1's
                      // add input matrix
 z(arma::span(0,n-1),arma::span(1,p-1)) = x;
 arma::colvec coef = arma::solve(z, y2); // beta coefficients
 arma::colvec resid = y2 - z*coef;       // residuals
 R1 = var(y2-resid)/var(y2); // coefficient of determination for model under consideration
 R  = (1-R1)/(1-R0);       // compute R
                           // allocate workspace for the integration
 gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);  
    
 double result, error;
 double c = 0.5*ke+a1+b1; // c
 // sequence of points to look the max (M) of the function over
 // M will be needed to avoid over/under-flow issues
 NumericVector sq = {pow(10,-10),0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99,0.999,1-pow(10,-10)};
 double u = 0;
 if (intrinsic){ // if intrinsic prior then add an extra point to the sequence
  // because for small R the integrating function becomes nearly singular
  double zc = 1/(1+delta*R);
  // point arisen by maximizing part of the integrand
  u = ((n-1)*zc-k1)/((n-1)*zc-k1*zc); 
  if (0<u && u<1){ // if this point is between 0 and 1
                   // add it to the sequence
   sq = {pow(10,-10),0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99,0.999,1-pow(10,-10),u};
  }
 }                      // end if (intrinsic)
 double s = sq.size();  // sequence length
 NumericVector vsq(s);  // function values at the sequence points initialization
 int k = 0;
 for (int i = 0; i < s; ++i){ // for point (i+1)
                              // compute function value
  double temp = (b1-1)*log(sq[i])+(c-b1-1)*log(1-sq[i])-a0*log(1-sq[i]*1/(1+delta*R))+an1*log(1-sq[i]*1/(1+delta));  
  if (arma::is_finite(temp)){ // if this is finite
   vsq[k] = temp;    // add it to the vector
   k = k+1;
  }                  // end if
 }                   // end for
 double M = max(vsq[Rcpp::Range(0,k-1)]); // compute max across the sequence
                     // vector of parameters needed for the integration
 NumericVector parameterv = Rcpp::NumericVector::create(b1,a0,
                               -an1,c,1/(1+delta*R),
                              1/(1+delta),M);
  
 gsl_function F;
 F.function = &A2simple;    // function to integrate
 F.params = &parameterv[0]; // with parameters
   
 // Gaussian-Kronrod quadrature routine
 gsl_integration_qags(&F, 0,1,0, 1e-7, 1000, w, &result, &error); 

 // if  the integration result is infinite or NAN then
 // stop and output an error message
 if (isinf(result)||isnan(result)) {
       Rcpp::Rcout << "Error with the integral!"
  << std::endl;
  stop("Error with the integral!");
 }

 double resulto = result;

 // For small R, the integrand becomes nearly singular in the case 
 // of intrinsic-- thus the previous integration result might be wrong. 
 // For that, if prior is intrinsic check if the result is below the relative 
 // error limit: if so try yet another type of integration (adaptive integration 
 // with known singular points) where the point u is given as a singular point.
 if ((result<1e-7) &&(intrinsic) && (0<u) && (u<1)){
  double pts[3];
  pts[0] = 0;
  pts[1] = u;
  pts[2] = 1;
  gsl_integration_qagp(&F, pts,3,0, 1e-7, 1000, w, &result, &error);
 }
 
 gsl_integration_workspace_free (w); // free the memory associated with the workspace

 List results;
 if (!pred){    // compute remaining terms needed for the result
  res = Rf_lbeta(0.5*ke+a1, b1) - Rf_lbeta(a1, b1) + 
        an1*log(1+delta)-a0*log(1+delta*R);
                                          // marginal likelihood
  res2 = log(result)+lgamma(c)-(lgamma(b1)+lgamma(c-b1))+M+res; 
  results["marglikel"] = res2;            // store it
  results["Rsquar"] = R1;                 // store Rsq
 }else{
  // allocate workspace for the integration
  gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);
  parameterv[2] = -an1+1; // slightly change this parameter for the second integral
  double result2;
  // use the same approach for the 2nd integral since it is not expected to 
  // differ a lot from the 1st
  if ((resulto<1e-7) &&(intrinsic) && (0<u) && (u<1)){
    double pts[3];
    pts[0] = 0;
    pts[1] = u; // give u as a known (nearly) singular point
    pts[2] = 1;
    gsl_integration_qagp(&F, pts,3,0, 1e-7, 1000, w, &result2, &error);     
  }else{        // compute second integral (needed for prediction)
    gsl_integration_qags(&F, 0,1,0, 1e-7, 1000, w, &result2, &error);
  }
  // if the integration result is infinite or NAN then stop and output an error message
  if (isinf(result2)||isnan(result2)){
       Rcpp::Rcout << "Error with the integral!"
   << std::endl;
   stop("Error with the integral!");
  } 
  gsl_integration_workspace_free (w); // free the memory associated with the workspace
  double Ew = log(delta/(1+delta)) + log(result2)-log(result); // eqn (10) (log) from Econometrics paper
  arma::colvec ypred = coef(0) + (xnew*coef(arma::span(1,p-1)))* exp(Ew); // eqn (9)
  results["predict"] = ypred;           // predicted/fitted values
 }                                      // end if
 return results;
}                                       // end function pepmarginallikelihoodc


// [[Rcpp::export]]
List full_enumeration_pepc(arma::mat x, arma::mat gamma, NumericVector y, 
                           bool intrinsic = false, bool reference_prior = true){
 int nmodels = gamma.n_rows;   // number of models
 // vectors for marginal likelihood and R-squared initialization
 NumericVector MargLikel(nmodels), Rsq(nmodels); 
 double d0 = 0, d1 = 0;
 if (!reference_prior)
   d0 = 1;
 // set the gsl error handler off because otherwise the R console
 // will automatically close (by default) in the presence of an error
 gsl_set_error_handler_off();

 for (int i = 1; i < nmodels; ++i){ // i starts from 1 (and not 0) because for the null model 
  // both the log marginal likelihood and Rsq are 0 (already set in initialization)
  arma::rowvec armag = gamma.row(i); // get the model representation
  int nvar = sum(armag); // number of input variables included
  if (nvar==1) {         // if only one
    arma::colvec xg = x.cols(find(armag==1)); //corresponding 1-column matrix
    if (!reference_prior) 
       d1 = 2;
                                   // compute marginal likelihood and Rsq
    List modelresults = pepmarginallikelihoodc(xg,y,xg,false,intrinsic,d0,d1);
    MargLikel[i] = modelresults["marglikel"]; // store marginal likelihood
    Rsq[i] = modelresults["Rsquar"];          // and Rsq
  }
  else {arma::mat xg = x.cols(find(armag==1)); //corresponding input matrix
    if (!reference_prior)
       d1 = xg.n_cols+1;
                                   // compute marginal likelihood and Rsq    
    List modelresults = pepmarginallikelihoodc(xg,y,xg,false,intrinsic,d0,d1);
    MargLikel[i] = modelresults["marglikel"]; // store marginal likelihood
    Rsq[i] = modelresults["Rsquar"];          // and Rsq
  }                                     // end if
 }                                      // end for
 gsl_set_error_handler(NULL);           // set the gsl error handler on again
 List resultspep;                       // list with the results
 resultspep["marglikel"] = MargLikel;   // 1st element: log marginal likelihoods
 resultspep["Rsquar"] = Rsq;            // 2nd element: Rsq's
 return resultspep;
}                                       // end function full_enumeration_pepc


// [[Rcpp::export]]
arma::mat mc3_pepc(arma::mat x, NumericVector y, bool betabinom, int itermc3, 
                   double current_Rsq, double current_marginal, double current_logpost, arma::rowvec current_gamma,
                   bool intrinsic = false, bool reference_prior = true){
 int p = x.n_cols;   // number of explanatory variables
 arma::mat resultM(itermc3,p+2); // matrix with the results, will contain the selected 
 // model after each iteration together with its marginal likelihood and Rsq
 NumericVector betabinomv(p+1); // prior probs vector initialization
 
 if (betabinom){     // if model prior is beta-binomial
  for (int j = 0; j < (p+1); ++j)
   betabinomv[j] = -Rf_lchoose(p,j); // compute log prior probs (shifted)
 }
 double d0 = 0, d1 = 0;
 if (!reference_prior)
   d0 = 1;
 // set the gsl error handler off because otherwise the R console
 // will automatically close (by default) in the presence of an error
 gsl_set_error_handler_off();

 for (int i = 0; i < itermc3; ++i){ // for iteration (i+1)
  IntegerVector permut = sample(p,p)-1; // permuted vector of (1,2,...p)
  // for element i in the permuted vector
  for(IntegerVector::iterator it = permut.begin(); it != permut.end(); ++it) {
   arma::rowvec proposed_gamma = current_gamma; // current model becomes proposed
   proposed_gamma[*it] = 1-proposed_gamma[*it]; // change status of the corresponding
                                                // input variable (1->0, 0->1)
   int pl = sum(proposed_gamma); // number of input variables included
   double proposed_marginal;
   double proposed_Rsq;
   if (pl==0){                   // if none
     proposed_marginal = 0;      // proposed marginal likelihood is set to 0
     proposed_Rsq = 0;           // same for Rsq
   }else if (pl==1){
     arma::colvec xg = x.cols(find(proposed_gamma==1)); //corresponding 1-column matrix
     if (!reference_prior) 
       d1 = 2;
                                 // compute marginal likelihood and Rsq
     List modelresults = pepmarginallikelihoodc(xg,y,xg,false,intrinsic,d0,d1);
     proposed_marginal = modelresults["marglikel"]; // store marginal likelihood
     proposed_Rsq = modelresults["Rsquar"];         // and Rsq
   }else{
     arma::mat xg = x.cols(find(proposed_gamma==1)); //corresponding input matrix
     if (!reference_prior)
       d1 = xg.n_cols+1;
                                 // compute marginal likelihood and Rsq
     List modelresults = pepmarginallikelihoodc(xg,y,xg,false,intrinsic,d0,d1);
     proposed_marginal = modelresults["marglikel"]; // store marginal likelihood
     proposed_Rsq = modelresults["Rsquar"];         // and Rsq
   }                                            // end if
   double proposed_logpost = proposed_marginal; // log posterior (proposed)
   if (betabinom){     // if model prior is beta-binomial
                       // add log prior to log marginal likelihood
    proposed_logpost = proposed_marginal+betabinomv[pl];
   }                    // difference between proposed and current log posterior
   double logalpha = proposed_logpost-current_logpost;
   NumericVector u = runif(1); // sample one value from uniform distribution
   if (log(u)[0]<logalpha){    // if that value is smaller than the difference
     current_gamma = proposed_gamma; // proposed model becomes current
     current_logpost = proposed_logpost; // with its log posterior prob
     current_marginal = proposed_marginal; // log marginal likelihood
     current_Rsq = proposed_Rsq;        // and Rsq
   }                                    // end if  
  }                                     // end inner for
  // save the final model in the resultsM matrix
  resultM(arma::span(i,i),arma::span(0,p-1)) = current_gamma;
  resultM(i,p) = current_marginal;      // with its marginal likelihood 
  resultM(i,p+1) = current_Rsq;         // and Rsq
 }                                      // end outer for
 gsl_set_error_handler(NULL);           // reset the gsl error handler
 return resultM;
}                                       // end function mc3_pepc


// [[Rcpp::export]]
List predict_pepc(arma::mat x, arma::mat gamma, NumericVector y,arma::mat xnew, 
                           bool intrinsic = false, bool reference_prior = true){
 int nmodels = gamma.n_rows, n = x.n_rows;    // number of models & sample size
 List predictionsspep(nmodels); // list of length number of models
 double d0 = 0, d1 = 0;
 if (!reference_prior)
   d0 = 1;
 // set the gsl error handler off because otherwise the R console
 // will automatically close (by default) in the presence of an error
 gsl_set_error_handler_off();

 for (int i = 0; i < nmodels; ++i){     // for model (i+1)
  arma::rowvec armag = gamma.row(i);    // get model representation
  int nvar = sum(armag); // number of input variables it includes
  if (nvar==0){                         // if none
    predictionsspep[i] = rep(mean(y),n); // predicted values will be the mean response
  }else if (nvar==1) {                  // if only one
    arma::colvec xg = x.cols(find(armag==1)); // corresponding 1-column matrix
    arma::colvec xgnew = xnew.cols(find(armag==1)); // corresponding 1-column new data matrix
    if (!reference_prior) 
       d1 = 2;
                                        // compute fitted/predicted values
    List predresults = pepmarginallikelihoodc(xg,y,xgnew,true,intrinsic,d0,d1);
    predictionsspep[i] = predresults["predict"]; // and store them
  }
  else {arma::mat xg = x.cols(find(armag==1)); // corresponding input matrix    
    arma::mat xgnew = xnew.cols(find(armag==1)); // corresponding new data matrix
    if (!reference_prior)
       d1 = xg.n_cols+1;
                                        // compute fitted/predicted values
    List predresults = pepmarginallikelihoodc(xg,y,xgnew,true,intrinsic,d0,d1);
    predictionsspep[i] = predresults["predict"]; // and store them
  }                                     // end if  
 }                                      // end for 
 gsl_set_error_handler(NULL);           // set the gsl error handler on again
 return predictionsspep;
}                                       // end function predict_pepc
