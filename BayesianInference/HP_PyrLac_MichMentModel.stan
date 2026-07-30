
functions{

    real[] PyruvateHP_NMR_ODEs4(real t, real[] y, real[] p, real[] x_r, int[] x_i){

      // Parameters
      real kin = p[1];
      real kf = p[2];
      real kr = p[3];
      real kpl = p[4];
      real T1_P = p[5];
      real T1_L = p[6];

      // ODEs
      real dInd_dt[5];

      real Pout = y[1];
      real Pin = y[2]; 
      real LDH = y[3]; 
      real LDHPyr = y[4]; 
      real Lac = y[5]; 

      dInd_dt[1] =  -Pout*kin - Pout/T1_P;
      dInd_dt[2] =  Pout*kin + kr*LDHPyr - kf*Pin*LDH - Pin/T1_P; 
      dInd_dt[3] =  kr*LDHPyr + kpl*LDHPyr - kf*Pin*LDH; 
      dInd_dt[4] =  kf*Pin*LDH - kr*LDHPyr - kpl*LDHPyr - LDHPyr/T1_P;
      dInd_dt[5] =  kpl*LDHPyr - (Lac/T1_L);


      // Results
      return dInd_dt;

    }


    vector SteadyState(vector init, vector p, real[] x_r, int[] x_i){

      return init; 

    }

}


data {
    // Observables
    int m; // Total number of data series
    int stslm; // Maximum number of rows for all the observable matrices
    array[1,m] int stsl; // Number of elements at each time series for each series m
    array[stslm,m] int sts; // Sampling times for each series m

    int obser;//-> Introduce this so we have all the data in one same array (easier generalisation). Work on generalisation in case different experiments have different obsevables?
    array[1,obser] int obSta; // -> This variable will be to know which are the observable states
    array[1,m] int ncells;
    array[1,m] int nts;

   //  real itp[m];
    array[stslm,m,obser] real Means; // ---> General arrays of means and errors
    array[stslm,m,obser] real Erros;

    // Inputs
    int tml; // Maximum length of the rows for the sampling times matrix
    array[tml, m] real ts; // Time series for each serie m
    array[5,m] real Y0us; // Y0 vectors
}


transformed data {
    int nParms = 6; // Number of parameters of the model //--------> Introduce number in generation of script
    int Neq = 5; // Total number of equations of the model //-----> Introduce number in generation of script
    array[0] int x_i; // Empty x_i object (needs to be defined)
    array[0] real x_r; // Input values for each event ordered as IPTG, aTc, IPTG, aTc, ...
    array[Neq,m] real ivss = Y0us; // Initial experimental values for the calculation of the steady state ordered as LacI+RFP, TetR+GFP -------------------------> Careful to how I define this
    real itp = 0;   
}


parameters {
    real<lower=-2,upper=2> kin;
    real<lower=-2,upper=2> kf;
    real<lower=-2,upper=2> kr;
    real<lower=-2,upper=2> kpl;
    real T1_P;
    real T1_X;

    real<lower=-2, upper=3.82> LDH;
}


transformed parameters {

   array[nParms] real theta; 
   array[1] real inits;

   // theta[1] = 0+(((kin-(-2))*(1-0))/(2-(-2)));
   // theta[2] = 0+(((kf-(-2))*(1-0))/(2-(-2)));
   // theta[3] = 0+(((kr-(-2))*(1-0))/(2-(-2)));
   // theta[4] = 0+(((kpl-(-2))*(1-0))/(2-(-2)));
   theta[1] = 0+(((kin-(-2))*(0.2-0.00000000001))/(2-(-2)));
   theta[2] = 0+(((kf-(-2))*(0.5-0.00000000001))/(2-(-2)));
   theta[3] = 0+(((kr-(-2))*(1.1-0.00000000001))/(2-(-2)));
   theta[4] = 0+(((kpl-(-2))*(0.4-0.00000000001))/(2-(-2)));
   theta[5] = 51;
   theta[6] = 41;

   inits[1] = (LDH*275)+550;


//    theta[1] = 0.0009999+(((k_f1-(-2))*(0.09999-0.0009999))/(2-(-2))); 
//    theta[2] = 0.000789+(((k_r1-(-2))*(0.0789-0.000789))/(2-(-2))); 
//    theta[3] = (T1_P*5.12)+51.21; 
//    theta[4] = (T1_X*8.34)+41.73; 
//    theta[5] = 4.896e-5+(((kin-(-2))*(0.004896-4.896e-5))/(2-(-2)));  

   // theta[1] = 0.00009999+(((k_f1-(-2))*(0.9999-0.00009999))/(2-(-2))); 
   // theta[2] = 0.000000789+(((k_r1-(-2))*(0.0789-0.000000789))/(2-(-2))); 
   // theta[3] = (T1_P*5.12)+51.21; 
   // theta[4] = (T1_X*8.34)+41.73; 
   // theta[5] = 4.896e-6+(((kin-(-2))*(0.04896-4.896e-6))/(2-(-2)));  


}


model {
   // Intermediate parameters
  vector[Neq] ing; // Vector that will include the solution of the algebraic solution for the steady state of the model
  
  array[Neq,m] real Y0; // Initial values for the ODEs variables at the first event
  array[stslm,m,Neq] real yhat; // ---> Generall array to include all the observables (easier generalisation)


   kin ~ uniform(-2,2); 
   kf ~ uniform(-2,2); 
   kr ~ uniform(-2,2); 
   kpl ~ uniform(-2,2); 
   LDH ~ normal(0,1); 


   // print(theta);
   // Likelihood
   for (j in 1:m){

      array[Neq] real ivst; // Initial value of the states
      array[tml+1,Neq] real y_hat; // Object to include the ODEs solutions for each state
      int lts = num_elements(ts[1:nts[1,j]+1,j]);  // Length of the time series for each event
      array[lts,Neq] real part1; // Temporary object that will include the solution of the ODE for each event at each loop
      vector[stsl[1,j]] yhatPyr;
      // Loop (over the number of events/inputs) to solve the ODE for each event stopping the solver and add them to the final object y_hat

      int q=1;
      
      // Calculation of the solution for the ODEs where for events that are not the first one. The time series starts one minute before the original point of the time serie overlaping with the last point of the previous event with same state values at the time
      ivst = ivss[,j];

      ivst[3] = inits[1];
      
      // print(itp);
      // print(ts[,j]);

 
         part1[1:nts[1,j]+1,] = integrate_ode_bdf(PyruvateHP_NMR_ODEs4,ivst,itp,ts[1:nts[1,j]+1,j],theta,x_r, x_i, 1.0e-9, 1.0e-9, 1e7);



      // Introduction of the result of part1 into the object y_hat
      y_hat[(1),] = ivst;
      for (y in 1:lts){
         y_hat[(y)+1,]=part1[(y),];
      };

      for (t in 1:stsl[1,j]){ //----> General form
            yhat[t,j,1] = part1[(sts[t,j]+1),1];
            yhat[t,j,2] = part1[(sts[t,j]+1),2];
            yhat[t,j,3] = part1[(sts[t,j]+1),3];
            yhat[t,j,4] = part1[(sts[t,j]+1),4];
            yhat[t,j,5] = part1[(sts[t,j]+1),5];

            
      }


      
      yhatPyr = to_vector(yhat[1:stsl[1,j],j,1]) + to_vector(yhat[1:stsl[1,j],j,2]) + to_vector(yhat[1:stsl[1,j],j,4]);

      // Means[1:stsl[1,j],j,1] ~ normal(yhatPyr[1:stsl[1,j]],Erros[1:stsl[1,j],j,1]);
      target += normal_lpdf(Means[1:stsl[1,j],j,1] | yhatPyr[1:stsl[1,j]], Erros[1:stsl[1,j],j,1]);


      // Means[1:stsl[1,j],j,2] ~ normal(yhat[1:stsl[1,j],j,3],Erros[1:stsl[1,j],j,2]);
      target += normal_lpdf(Means[1:stsl[1,j],j,2] | yhat[1:stsl[1,j],j,5], Erros[1:stsl[1,j],j,2]);
      // for (t in 1:stsl[1,j]){
      //    Means[t,j,1] ~ normal(yhat[t,j,1]+yhat[t,j,2],Erros[t,j,1]);
      //    // Means[t,j,2] ~ normal(yhat[t,j,3],Erros[t,j,2]);
      // }
   }
}




































