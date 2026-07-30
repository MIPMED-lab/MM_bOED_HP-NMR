
functions{

    real[] PyruvateHP_NMR_ODEs2(real t, real[] y, real[] p, real[] x_r, int[] x_i){

      // Parameters
      real k_f1 = 2*p[1]; 
      real k_r1 = 2*p[2]; 
      real T1_P = p[3]; 
      real T1_X = p[4]; 

      // ODEs
      real dInd_dt[2];

      real Php = y[1]; 
      real Xhp = y[2]; 

      dInd_dt[1] =  k_r1*Xhp - ((1/T1_P)+k_f1)*Php; 
      dInd_dt[2] =  k_f1*Php - ((1/T1_X)+k_r1)*Xhp; 

      // Results
      return dInd_dt;

    }

   real[] PyruvateHP_NMR_ODEs4(real t, real[] y, real[] p, real[] x_r, int[] x_i){

      // Parameters
      real k_f1 = 4*p[1]; 
      real k_r1 = 4*p[2]; 
      real T1_P = p[3]; 
      real T1_X = p[4]; 

      // ODEs
      real dInd_dt[2];

      real Php = y[1]; 
      real Xhp = y[2]; 

      dInd_dt[1] =  k_r1*Xhp - ((1/T1_P)+k_f1)*Php; 
      dInd_dt[2] =  k_f1*Php - ((1/T1_X)+k_r1)*Xhp; 

      // Results
      return dInd_dt;

    }

    real[] PyruvateHP_NMR_ODEs8(real t, real[] y, real[] p, real[] x_r, int[] x_i){

      // Parameters
      real k_f1 = 8*p[1]; 
      real k_r1 = 8*p[2]; 
      real T1_P = p[3]; 
      real T1_X = p[4]; 

      // ODEs
      real dInd_dt[2];

      real Php = y[1]; 
      real Xhp = y[2]; 

      dInd_dt[1] =  k_r1*Xhp - ((1/T1_P)+k_f1)*Php; 
      dInd_dt[2] =  k_f1*Php - ((1/T1_X)+k_r1)*Xhp; 

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
    array[2,m] real Y0us; // Y0 vectors
}


transformed data {
    int nParms = 4; // Number of parameters of the model //--------> Introduce number in generation of script
    int Neq = 2; // Total number of equations of the model //-----> Introduce number in generation of script
    array[0] int x_i; // Empty x_i object (needs to be defined)
    array[0] real x_r; // Input values for each event ordered as IPTG, aTc, IPTG, aTc, ...
    array[Neq,m] real ivss = Y0us; // Initial experimental values for the calculation of the steady state ordered as LacI+RFP, TetR+GFP -------------------------> Careful to how I define this
    real itp = 0;   
}


parameters {
    real k_f1; 
    real k_r1; 
    real<lower=-2, upper=2> T1_P; 
    real<lower=-2, upper=2> T1_X; 
}


transformed parameters {

   array[nParms] real theta; 
// theta[1] = 0+(((k_f1-(-2))*(10^-2-0))/(2-(-2)));
   theta[1] = 1.68e-5+(((k_f1-(-2))*(0.00168-1.68e-5))/(2-(-2)));
   theta[2] = 2.23e-8+(((k_r1-(-2))*((2.23e-4)-2.23e-8))/(2-(-2))); 
   theta[3] = (T1_P*5)+50; // sd 10% of estimated
   theta[4] = (T1_X*10.5)+55; // sd 20 % of estimated

}


model {
   // Intermediate parameters
  vector[Neq] ing; // Vector that will include the solution of the algebraic solution for the steady state of the model
  
  array[Neq,m] real Y0; // Initial values for the ODEs variables at the first event
  array[stslm,m,Neq] real yhat; // ---> Generall array to include all the observables (easier generalisation)


   k_f1 ~ uniform(-2,2); 
   k_r1 ~ uniform(-2,2); 
   T1_X ~ normal(0,1); 
   T1_P ~ normal(0,1); 


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
      
      // print(itp);
      // print(ts[,j]);

      

      if (ncells[1,j]==2){
         part1[1:nts[1,j]+1,] = integrate_ode_bdf(PyruvateHP_NMR_ODEs2,ivst,itp,ts[1:nts[1,j]+1,j],theta,x_r, x_i, 1.0e-9, 1.0e-9, 1e7);
      };
      if (ncells[1,j]==4){
         part1[1:nts[1,j]+1,] = integrate_ode_bdf(PyruvateHP_NMR_ODEs4,ivst,itp,ts[1:nts[1,j]+1,j],theta,x_r, x_i, 1.0e-9, 1.0e-9, 1e7);
      };
      if (ncells[1,j]==8){
         part1[1:nts[1,j]+1,] = integrate_ode_bdf(PyruvateHP_NMR_ODEs8,ivst,itp,ts[1:nts[1,j]+1,j],theta,x_r, x_i, 1.0e-9, 1.0e-9, 1e7);
      };


      // Introduction of the result of part1 into the object y_hat
      y_hat[(1),] = ivst;
      for (y in 1:lts){
         y_hat[(y)+1,]=part1[(y),];
      };

      for (t in 1:stsl[1,j]){ //----> General form
            yhat[t,j,1] = part1[(sts[t,j]+1),1];
            yhat[t,j,2] = part1[(sts[t,j]+1),2];
      }


      

      // Means[1:stsl[1,j],j,1] ~ normal(yhatPyr[1:stsl[1,j]],Erros[1:stsl[1,j],j,1]);
      target += normal_lpdf(Means[1:stsl[1,j],j,1] | yhat[1:stsl[1,j],j,1], Erros[1:stsl[1,j],j,1]);


      // Means[1:stsl[1,j],j,2] ~ normal(yhat[1:stsl[1,j],j,3],Erros[1:stsl[1,j],j,2]);
      target += normal_lpdf(Means[1:stsl[1,j],j,2] | yhat[1:stsl[1,j],j,2], Erros[1:stsl[1,j],j,2]);
      // for (t in 1:stsl[1,j]){
      //    Means[t,j,1] ~ normal(yhat[t,j,1]+yhat[t,j,2],Erros[t,j,1]);
      //    // Means[t,j,2] ~ normal(yhat[t,j,3],Erros[t,j,2]);
      // }
   }
}




































