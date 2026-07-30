
functions{

   real[] PyruvateHP_NMR_ODEs4(real t, real[] y, real[] p, real[] x_r, int[] x_i){

      // Parameters
      real T1_X = p[1]; 
      real T1_P = p[2];
      real kin = p[3];
      real kpl = p[4];
      real kbn = p[5];
      real kun = p[6];
      real kbp = p[7];
      real kup = p[8];
      real kunE = p[9];
      real kbnE = p[10];
      real ki = p[11];
      real kr = p[12];
      real ScFm = p[13];

      // ODEs
      real dInd_dt[14];

      real Pout = y[1];
      real Php = y[2]; 
      real Xhp = y[3]; 
      real NADH = y[4];
      real NAD = y[5];
      real LDH = y[6];
      real LDHac = y[7];
      real LDHacE = y[8];
      real LDHna = y[9];
      real LDHP = y[10];
      real XhpHyper = y[11];
      real PhpHyper = y[12];
      real Obs_XhpHyper = y[13];
      real Obs_PhpHyper = y[14];

      dInd_dt[1] = - (Pout*kin);
      dInd_dt[2] = (Pout*kin) + kup*LDHP + kr*LDHna - kbp*LDHac*Php - ki*LDH*Php; 
      dInd_dt[3] = kpl*LDHP; 
      dInd_dt[4] = kun*LDHac - kbn*LDH*NADH;
      dInd_dt[5] = kunE*LDHacE - kbnE*LDH*NAD;
      dInd_dt[6] = kun*LDHac + kr*LDHna + kunE*LDHacE - kbn*LDH*NADH - ki*LDH*Php - kbnE*LDH*NAD;
      dInd_dt[7] = kbn*LDH*NADH + kup*LDHP - kun*LDHac - kbp*LDHac*Php;
      dInd_dt[8] = kpl*LDHP + kbnE*LDH*NAD - kunE*LDHacE;
      dInd_dt[9] = ki*LDH*Php - kr*LDHna;
      dInd_dt[10] = kbp*LDHac*Php - kup*LDHP - kpl*LDHP;
      dInd_dt[11] = dInd_dt[3] - (XhpHyper/T1_X);
      dInd_dt[12] = (dInd_dt[1] + dInd_dt[2] + dInd_dt[9] + dInd_dt[10]) - (PhpHyper/T1_P);
      dInd_dt[13] = dInd_dt[11] * (ScFm*0.12);
      dInd_dt[14] = dInd_dt[12] * (ScFm*0.12);

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
    array[14,m] real Y0us; // Y0 vectors
}


transformed data {
    int nParms = 13; // Number of parameters of the model //--------> Introduce number in generation of script
    int Neq = 14; // Total number of equations of the model //-----> Introduce number in generation of script
    array[0] int x_i; // Empty x_i object (needs to be defined)
    array[0] real x_r; // Input values for each event ordered as IPTG, aTc, IPTG, aTc, ...
    array[Neq,m] real ivss = Y0us; // Initial experimental values for the calculation of the steady state ordered as LacI+RFP, TetR+GFP -------------------------> Careful to how I define this
    real itp = 0;   
}


parameters {

    // real<lower=-2, upper=2> T1_P; 
    // real<lower=-2, upper=2> T1_X; 
    real kin; 
    real kpl;
    real kbn;
    real kun;
    real kbp;
    real kup;
    real kunE;
    real kbnE;
    real ki;
    real kr;
    // real ScFm;

    real NADH;
    real NAD;
    real LDH;
}


transformed parameters {

   array[nParms] real theta; 
   array[3] real inits;

   theta[1] = 48;
   theta[2] = 55;
   theta[3] = (kin*0.00275)+0.0055;
   theta[4] = (kpl*0.01375)+0.0275;
   theta[5] = (kbn*0.0225)+0.045;
   theta[6] = (kun*0.225)+0.45;
   theta[7] = (kbp*0.475)+0.95;
   theta[8] = (kup*0.475)+0.95;
   theta[9] = (kunE*0.475)+0.95;
   theta[10] = (kbnE*0.00425)+0.0085;
   theta[11] = (ki*0.475)+0.95;
   theta[12] = (kr*0.00145)+0.0029;
   theta[13] = 15;
   
   inits[1] = (NADH*260)+650;
   inits[2] = (NAD*6400)+16000;
   inits[3] = (LDH*275)+550;

}




model {
   // Intermediate parameters
  vector[Neq] ing; // Vector that will include the solution of the algebraic solution for the steady state of the model
  
  array[Neq,m] real Y0; // Initial values for the ODEs variables at the first event
  array[stslm,m,Neq] real yhat; // ---> Generall array to include all the observables (easier generalisation)

//    T1_X ~ normal(0,1); 
//    T1_P ~ normal(0,1); 
   kin ~ normal(0,1);
   kpl ~ normal(0,1);
   kbn ~ normal(0,1);
   kun ~ normal(0,1);
   kbp ~ normal(0,1);
   kup ~ normal(0,1);
   kunE ~ normal(0,1);
   kbnE ~ normal(0,1);
   ki ~ normal(0,1);
   kr ~ normal(0,1);
//    ScFm ~ uniform(-2,2);

   NADH ~ normal(0,1);
   NAD ~ normal(0,1);
   LDH ~ normal(0,1);


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
      ivst[5] = inits[2];
      ivst[6] = inits[3];

      // print(itp);
      
      part1[1:nts[1,j]+1,] = integrate_ode_bdf(PyruvateHP_NMR_ODEs4,ivst,itp,ts[1:nts[1,j]+1,j],theta,x_r, x_i, 1.0e-9, 1.0e-9, 1e7);
      
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


      

      // Means[1:stsl[1,j],j,1] ~ normal(yhatPyr[1:stsl[1,j]],Erros[1:stsl[1,j],j,1]);
      // target += normal_lpdf(Means[1:stsl[1,j],j,1] | yhat[1:stsl[1,j],j,14], Erros[1:stsl[1,j],j,1]);


      // Means[1:stsl[1,j],j,2] ~ normal(yhat[1:stsl[1,j],j,3],Erros[1:stsl[1,j],j,2]);
      target += normal_lpdf(Means[1:stsl[1,j],j,2] | yhat[1:stsl[1,j],j,13], Erros[1:stsl[1,j],j,2]);


      // for (t in 1:stsl[1,j]){
      //    Means[t,j,1] ~ normal(yhat[t,j,1]+yhat[t,j,2],Erros[t,j,1]);
      //    // Means[t,j,2] ~ normal(yhat[t,j,3],Erros[t,j,2]);
      // }
   }
}




































