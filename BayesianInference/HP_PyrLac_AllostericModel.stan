
functions{

   real[] PyruvateHP_NMR_ODEs4(real t, real[] y, real[] p, real[] x_r, int[] x_i){

      // Parameters

      real T1_X = p[1];
      real T1_P = p[2];
      real kin = p[3];
      real kpl = p[4];
      real kf = p[5];
      real kr = p[6];
      real a1 = p[7];
      real a2 = p[8];
      real a3 = p[9];
      real b1 = p[10];
      real b2 = p[11];
      real b3 = p[12];
      real ScFm = p[13];

      // ODEs
      real dInd_dt[12];

      real Pout = y[1];
      real Php = y[2];
      real Xhp = y[3];
      real LDH = y[4];
      real LDHP1 = y[5];
      real LDHP2 = y[6];
      real LDHP3 = y[7];
      real LDHP4 = y[8];
      real XhpHyper = y[9];
      real PhpHyper = y[10];
      real Obs_XhpHyper = y[11];
      real Obs_PhpHyper = y[12];

      dInd_dt[1] = - (Pout*kin);
      dInd_dt[2] = (Pout*kin) + (kr*LDHP1) + (a1*kr*LDHP2) + (a2*kr*LDHP3) + (a3*kr*LDHP4) - (kf*LDH*Php) - (a1*kf*LDHP1*Php) - (a2*kf*LDHP2*Php) - (a3*kf*LDHP3*Php);
      dInd_dt[3] = (kpl*LDHP1)*3 + (b1*kpl*LDHP2)*2 + (b2*kpl*LDHP3)*2 + (b3*kpl*LDHP4);
      dInd_dt[4] = (kr*LDHP1) + (kpl*LDHP1)*3 - (kf*LDH*Php);
      dInd_dt[5] = (kf*LDH*Php) + (a1*kr*LDHP2) + (b1*kpl*LDHP2)*2 - (kr*LDHP1) - (a1*kf*LDHP1*Php) - (kpl*LDHP1)*3;
      dInd_dt[6] = (a1*kf*LDHP1*Php) + (a2*kr*LDHP3) + (b2*kpl*LDHP3)*2 - (a1*kr*LDHP2) - (a2*kf*LDHP2*Php) - (b1*kpl*LDHP2)*2;
      dInd_dt[7] = (a2*kf*LDHP2*Php) + (a3*kr*LDHP4) + (b3*kpl*LDHP4) - (a2*kr*LDHP3) - (a3*kf*LDHP3*Php) - (b2*kpl*LDHP3)*2;
      dInd_dt[8] = (a3*kf*LDHP3*Php) - (a3*kr*LDHP4) - (b3*kpl*LDHP4);
      dInd_dt[9] = dInd_dt[3] - (XhpHyper/T1_X);
      dInd_dt[10] = dInd_dt[1] + dInd_dt[2]  + dInd_dt[5] + dInd_dt[6]*2 + dInd_dt[7]*3 + dInd_dt[8]*4  - (PhpHyper/T1_P);
      dInd_dt[11] = dInd_dt[9]*(ScFm*0.12);
      dInd_dt[12] = dInd_dt[10]*(ScFm*0.12);

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
    array[12,m] real Y0us; // Y0 vectors
}


transformed data {
    int nParms = 13; // Number of parameters of the model //--------> Introduce number in generation of script
    int Neq = 12; // Total number of equations of the model //-----> Introduce number in generation of script
    array[0] int x_i; // Empty x_i object (needs to be defined)
    array[0] real x_r; // Input values for each event ordered as IPTG, aTc, IPTG, aTc, ...
    array[Neq,m] real ivss = Y0us; // Initial experimental values for the calculation of the steady state ordered as LacI+RFP, TetR+GFP -------------------------> Careful to how I define this
    real itp = 0;   
}


parameters {

    real<lower=-2, upper=2> T1_P; 
    real<lower=-2, upper=2> T1_X; 
    real kin; 
    real kpl;
    real kf;
    real kr;
    real a1;
    real a2;
    real a3;
    real b1;
    real b2;
    real b3;
    real ScFm;

    real LDH;

}


transformed parameters {

   array[nParms] real theta; 
   array[1] real inits;

   theta[1] = (T1_X*10)+40.92;
   theta[2] = (T1_P*10)+51.23;
   theta[3] = 0.0001+(((kin-(-2))*(0.01-0.0001))/(2-(-2)));
   theta[4] = 0.007+(((kpl-(-2))*(0.7-0.007))/(2-(-2)));
   theta[5] = 0.001+(((kf-(-2))*(0.5-0.001))/(2-(-2)));
   theta[6] = 0.00001+(((kr-(-2))*(0.1-0.00001))/(2-(-2)));
   theta[7] = 0.00001+(((a1-(-2))*(0.001-0.00001))/(2-(-2)));
   theta[8] = 0.001+(((a2-(-2))*(1-0.001))/(2-(-2)));
   theta[9] = 0.00001+(((a3-(-2))*(0.001-0.00001))/(2-(-2)));
   theta[10] = 0.001+(((b1-(-2))*(1-0.001))/(2-(-2)));
   theta[11] = 0.00001+(((b2-(-2))*(0.001-0.00001))/(2-(-2)));
   theta[12] = 0.005+(((b3-(-2))*(0.6-0.005))/(2-(-2)));
   theta[13] = (ScFm*2.29)+18.5;
   
   inits[1] = 1+(((LDH-(-2))*(700-1))/(2-(-2)));

}


model {
   // Intermediate parameters
  vector[Neq] ing; // Vector that will include the solution of the algebraic solution for the steady state of the model
  
  array[Neq,m] real Y0; // Initial values for the ODEs variables at the first event
  array[stslm,m,Neq] real yhat; // ---> Generall array to include all the observables (easier generalisation)

   T1_X ~ normal(0,1); 
   T1_P ~ normal(0,1); 
   kin ~ uniform(-2,2);
   kpl ~ uniform(-2,2);
   kf ~ uniform(-2,2);
   kr ~ uniform(-2,2);
   a1 ~ uniform(-2,2);
   a2 ~ uniform(-2,2);
   a3 ~ uniform(-2,2);
   b1 ~ uniform(-2,2);
   b2 ~ uniform(-2,2);
   b3 ~ uniform(-2,2);
   ScFm ~ uniform(-2,2);

   LDH ~ uniform(-2,2);

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
      
      ivst[4] = inits[1];

      // print(itp);
      // print(ts[,j]);

      

      
      part1[1:nts[1,j]+1,] = integrate_ode_bdf(PyruvateHP_NMR_ODEs4,ivst,itp,ts[1:nts[1,j]+1,j],theta,x_r, x_i, 1.0e-12, 1.0e-12, 1e7);
      


      // Introduction of the result of part1 into the object y_hat
      y_hat[(1),] = ivst;
      for (y in 1:lts){
         y_hat[(y)+1,]=part1[(y),];
      };

      for (t in 1:stsl[1,j]){ //----> General form
         for (ob in 1:Neq){
            yhat[t,j,ob] = part1[(sts[t,j]+1),ob];
         }
      }


      
      yhatPyr = to_vector(yhat[1:stsl[1,j],j,1]);

      // Means[1:stsl[1,j],j,1] ~ normal(yhatPyr[1:stsl[1,j]],Erros[1:stsl[1,j],j,1]);
      target += normal_lpdf(Means[1:stsl[1,j],j,1] | yhat[1:stsl[1,j],j,12], Erros[1:stsl[1,j],j,1]);


      // Means[1:stsl[1,j],j,2] ~ normal(yhat[1:stsl[1,j],j,3],Erros[1:stsl[1,j],j,2]);
      target += normal_lpdf(Means[1:stsl[1,j],j,2] | yhat[1:stsl[1,j],j,11], Erros[1:stsl[1,j],j,2]);
      // for (t in 1:stsl[1,j]){
      //    Means[t,j,1] ~ normal(yhat[t,j,1]+yhat[t,j,2],Erros[t,j,1]);
      //    // Means[t,j,2] ~ normal(yhat[t,j,3],Erros[t,j,2]);
      // }
   }
}




































