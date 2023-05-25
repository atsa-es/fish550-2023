// Dynamic Factor Analysis for multivariate time series
#include <TMB.hpp>

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(obs); /*  timeSteps x stateDim*/
  DATA_MATRIX(Covar);
  PARAMETER_VECTOR(logsdObs);
  PARAMETER_VECTOR(cholCorr);
  PARAMETER_MATRIX(covState); /* x[t] - x[t-1] */
  PARAMETER_MATRIX(covinitState); /* x[1] */
  PARAMETER_MATRIX(D);
  PARAMETER_MATRIX(Z);
  PARAMETER_MATRIX(u); /* State */
  
  
  int timeSteps=obs.col(0).size();
  int obsDim=obs.row(0).size();
  
  vector<Type> sdObs=exp(logsdObs);
  
  using namespace density;
  UNSTRUCTURED_CORR_t<Type> corMatGen(cholCorr);// This is the full Cormat
  matrix<Type> FullCorrMat=corMatGen.cov();
  
 
  MVNORM_t<Type> initialState(covinitState);
  MVNORM_t<Type> neg_log_density_process(covState);
  /* Define likelihood */
  Type ans=0;
  //ans -= dnorm(vector<Type>(u.row(0)),Type(0),Type(1),1).sum();
  ans += initialState(u.row(0));
  for(int i=1;i<timeSteps;i++){ 
    ans+= neg_log_density_process(u.row(i)-u.row(i-1)); // Process likelihood 
  }
  
  matrix<Type> pred(timeSteps,obsDim);  
  pred = (Z * u.transpose()) + (D * Covar);
  
  for(int i=0;i<timeSteps;i++){ //move one time step at a time
     int nonNAcount = 0; //start at zero NA values
	 vector<int> GoodVals(obs.row(i).size());
	 for(int j=0;j<obs.row(i).size();j++){//loop over all time series for this time step
	    if(!isNA(obs.row(i)(j))){//if value is not NA
			GoodVals(nonNAcount) = j; //add position to good values (good values only stored in beginning of vector)
			nonNAcount++; //increment the values of
		}
	 }
	 if(nonNAcount<obs.row(i).size()){ //if NA values present
		matrix<Type> subCorr(nonNAcount,nonNAcount);
		vector<Type> subSds(nonNAcount);
	 	vector<Type> subData(nonNAcount);
	 	vector<Type> subPred(nonNAcount);
	 	
	 	for(int j=0;j<nonNAcount;j++){
	 		subData(j) = obs.row(i)(GoodVals(j));
			subPred(j) = pred.transpose().row(i)(GoodVals(j));
			subSds(j) = sdObs(GoodVals(j));
			for(int k=0;k<nonNAcount;k++){
				subCorr(j,k) = FullCorrMat(GoodVals(j),GoodVals(k));
			}//end of loop through for truncated cormat
		}//end of removal of NA's from sets
		vector<Type> subDiffer = subData-subPred;
		ans += VECSCALE(MVNORM(subCorr),subSds)(subDiffer);
	 }else{
	   	vector<Type> differ = obs.row(i)-pred.transpose().row(i);
	 	ans += VECSCALE(corMatGen,sdObs)(differ);
	 }//end of data likelihood for this time step
  }//end of loop over time steps
  
  matrix<Type> FullCovMat(obsDim,obsDim);
  matrix<Type> dSD(obsDim,1);
  dSD = sdObs;
  FullCovMat = dSD.asDiagonal() * FullCorrMat * dSD.asDiagonal();
  ADREPORT(Z);
  ADREPORT(D);
  ADREPORT(u);
  ADREPORT(FullCovMat);
  
  REPORT(FullCorrMat);
  
  return ans;
}
