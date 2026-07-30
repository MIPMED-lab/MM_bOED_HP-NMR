
functions{

   real[] PyruvateHP_NMR_ODEs4(real t, real[] y, real[] p, real[] x_r, int[] x_i){

      // Parameters
      real T1_X = p[1]; 
      real T1_P = p[2];
      real kpl = p[3];
      real kbn = p[4];
      real kun = p[5];
      real kbp = p[6];
      real kup = p[7];
      real kunE = p[8];
      real kbnE = p[9];
      real ki = p[10];
      real kr = p[11];
      real nmCl = p[12];
      real ScFm = p[13];

      // ODEs
      real dInd_dt[13];

      real Php = y[1]; 
      real Xhp = y[2]; 
      real NADH = y[3];
      real NAD = y[4];
      real LDH = y[5];
      real LDHac = y[6];
      real LDHacE = y[7];
      real LDHna = y[8];
      real LDHP = y[9];
      real XhpHyper = y[10];
      real PhpHyper = y[11];
      real Obs_XhpHyper = y[12];
      real Obs_PhpHyper = y[13];

      dInd_dt[1] = kup*LDHP + kr*LDHna - kbp*LDHac*Php - ki*LDH*Php; 
      dInd_dt[2] = kpl*LDHP; 
      dInd_dt[3] = kun*LDHac - kbn*LDH*NADH;
      dInd_dt[4] = kunE*LDHacE - kbnE*LDH*NAD;
      dInd_dt[5] = kun*LDHac + kr*LDHna + kunE*LDHacE - kbn*LDH*NADH - ki*LDH*Php - kbnE*LDH*NAD;
      dInd_dt[6] = kbn*LDH*NADH + kup*LDHP - kun*LDHac - kbp*LDHac*Php;
      dInd_dt[7] = kpl*LDHP + kbnE*LDH*NAD - kunE*LDHacE;
      dInd_dt[8] = ki*LDH*Php - kr*LDHna;
      dInd_dt[9] = kbp*LDHac*Php - kup*LDHP - kpl*LDHP;
      dInd_dt[10] = dInd_dt[2]*nmCl - (XhpHyper/T1_X);
      dInd_dt[11] = (dInd_dt[1] + dInd_dt[8] + dInd_dt[9]) - (PhpHyper/T1_P);
      dInd_dt[12] = dInd_dt[10] * (ScFm*0.12);
      dInd_dt[13] = dInd_dt[11] * (ScFm*0.12);

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
    array[13,m] real Y0us; // Y0 vectors

    # scalling factor
    array[1,m] real ScFms;
    array[1,m] real nmCls;
}


transformed data {
    int nParms = 13; // Number of parameters of the model //--------> Introduce number in generation of script
    int Neq = 13; // Total number of equations of the model //-----> Introduce number in generation of script
    array[0] int x_i; // Empty x_i object (needs to be defined)
    array[0] real x_r; // Input values for each event ordered as IPTG, aTc, IPTG, aTc, ...
    array[Neq,m] real ivss = Y0us; // Initial experimental values for the calculation of the steady state ordered as LacI+RFP, TetR+GFP -------------------------> Careful to how I define this
    real itp = 0;   
}


parameters {

    // real<lower=-2, upper=2> T1_P; 
   //  real<lower=-2, upper=2> T1_X; 
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

   //  real<lower= -1.9201680672268908> NADH;
   //  real<lower=-2.4771653543307086> NAD;
    real LDH;
}


transformed parameters {

   array[nParms,m] real theta; 
   array[3] real inits;

   for (j in 1:m){
      theta[1,j] = 51.9806429875;
      theta[2,j] = 51.23000000000003;
      theta[3,j] = (kpl*247)+2470;
      theta[4,j] = (kbn*523)+5229;
      theta[5,j] = (kun*15)+152;
      theta[6,j] = (kbp*1)+9.5;
      theta[7,j] = (kup*91)+915;
      theta[8,j] = (kunE*459)+4589;
      theta[9,j] = (kbnE*87)+879;
      theta[10,j] = (ki*1)+11;
      theta[11,j] = (kr*81)+816;
      theta[12,j] = 3e6;
      theta[13,j] = 1;
      
      
      inits[1] = 220;
      inits[2] = 0;
      inits[3] = (LDH*0.000268)+0.00268;
      
      
   }

}





model {
   // Intermediate parameters
  vector[Neq] ing; // Vector that will include the solution of the algebraic solution for the steady state of the model
  
  array[Neq,m] real Y0; // Initial values for the ODEs variables at the first event
  array[stslm,m,Neq] real yhat; // ---> Generall array to include all the observables (easier generalisation)

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

   LDH ~ normal(0,1);

   
   // Likelihood
   for (j in 1:m){

      array[Neq] real ivst; // Initial value of the states
      array[tml,Neq] real y_hat; // Object to include the ODEs solutions for each state
      int lts = num_elements(ts[1:nts[1,j]+1,j]);  // Length of the time series for each event
      array[lts,Neq] real part1; // Temporary object that will include the solution of the ODE for each event at each loop
      vector[stsl[1,j]] yhatPyr;

      real nadhCn_muM;
      real nadh_fmcmin;

      // Loop (over the number of events/inputs) to solve the ODE for each event stopping the solver and add them to the final object y_hat

      int q=1;
      

      // Calculation of the solution for the ODEs where for events that are not the first one. The time series starts one minute before the original point of the time serie overlaping with the last point of the previous event with same state values at the time
      ivst = ivss[,j];
      
      ivst[3] = inits[1];
      ivst[4] = inits[2];
      ivst[5] = inits[3];

      // print(ts[1:nts[1,j]+1,j]);
      
      
      part1[1:nts[1,j]+1,] = integrate_ode_bdf(PyruvateHP_NMR_ODEs4,ivst,itp,ts[1:nts[1,j]+1,j],theta[,j],x_r, x_i, 1.0e-9, 1.0e-9, 1e7);
      

      // Introduction of the result of part1 into the object y_hat
      y_hat[(1),] = ivst;
      for (y in 1:lts){
         y_hat[(y),]=part1[(y),];
      };


      for (t in 1:2){ //----> General form
         for (ob in 1:Neq){
            yhat[t,j,ob] = part1[t,ob];
         }
      }

      nadhCn_muM = 220 - (yhat[2,j,3]+yhat[2,j,6]+yhat[2,j,9]);
      
      nadh_fmcmin = ((nadhCn_muM*2e-4)/1e5)*1e9;

      target += normal_lpdf(Means[2,j,1] | nadh_fmcmin, Erros[2,j,1]);


   }
}




































