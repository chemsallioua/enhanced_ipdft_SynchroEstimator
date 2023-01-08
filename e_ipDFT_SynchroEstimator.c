#define S_FUNCTION_NAME e_ipDFT_SynchroEstimator /* Defines and Includes */
#define S_FUNCTION_LEVEL 2 

#include "simstruc.h"
#include <math.h>
#define PI 3.141592653589793238

typedef struct {
    double r;
    double i; 
}complex;



complex hann_frt(double k, double N);
double hann(double* out_ptr, unsigned int out_len);
int e_ipDFT(complex* v_dft, int v_dft_len, double* delta, complex* phasor, int P, int n);
int dft_r(double* in_ptr, complex* out_ptr , unsigned int out_len, int n_bins, double norm_factor );



#define SAMPLE_RATE_IDX 0                                      //Signal Sample Rate or Sampling Frequency Index
#define SAMPLE_RATE(S) ssGetSFcnParam(S,SAMPLE_RATE_IDX)       //Signal Sample Rate MACRO getter
 
#define NOMINAL_FREQ_IDX 1                                         //Signal Nominal Frequency Index
#define NOMINAL_FREQ_PARAM(S) ssGetSFcnParam(S,NOMINAL_FREQ_IDX)   //Signal Nominal frequuency MACRO getter
 
#define WINDOW_LENGTH_IDX 2                                         //Window length Index
#define WINDOW_LENGTH_PARAM(S) ssGetSFcnParam(S,WINDOW_LENGTH_IDX)  //Window length expressed in number of samples, MACRO getter

#define P_IDX 3                                       //Number of iterations in the e-IPDT part
#define P(S) ssGetSFcnParam(S,P_IDX)                  //Number of iterations in the e-IPDT part MACRO getter

#define FRAME_PER_SEC_IDX 4                                      //Output Frames per Second
#define FRAME_PER_SEC(S) ssGetSFcnParam(S,FRAME_PER_SEC_IDX)    //Output Frames per Second MACRO getter

#define PHASE_MODE_IDX 5                                  //Frame Phase Mode : 1 for Power Frequency as Reference, 2 for any frequecny as reference 
#define PHASE_MODE(S) ssGetSFcnParam(S,PHASE_MODE_IDX)    //Frame Phase Mode MACRO getter

#define THRS_IDX 6                                  //Thresholds to trigger change in the states S1,S2 for rocof estimation
#define THRS(S) ssGetSFcnParam(S,THRS_IDX)          //Thresholds MACRO getter

#define LP_IDX 7                                  //Lowpass Filter Coefficients a1 , b0, b1 respectively (for rocof estimation)
#define LP(S) ssGetSFcnParam(S,LP_IDX)            //Lowpass Filter Coefficients MACRO getter
 
#define NPARAMS 8   //Total number of Parameters


#define HANN_IDX 0
#define HANN_PTR(S) ((double*)ssGetDWork(S,HANN_IDX))     //Hann Window Coefficients DWork Vector

#define NDWORK 1   //total number of DWork Vectorss


#define OLD_FREQ_IDX 0
#define OLD_ROCOF_IDX 1
#define STATE_IDX 2

#define NDSCSTAT 3   //total number of Descrete States

double norm_fact;

/*====================*
 * S-function methods *
 *====================*/
/* Function: mdlInitializeSizes ===============================================
 *
 *    The sizes information is used by Simulink to determine the S-function
 *    block's characteristics (number of inputs, outputs, states, etc.).  
 */

static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumSFcnParams(S, NPARAMS);  //Number of expected parameters

    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        return; /* Parameter mismatch reported by the Simulink engine*/
    }

    int nInputPorts = 2;   //Number of Import Ports
    int nOutputPorts = 5;   //Set number of Output Ports

    if (!ssSetNumInputPorts(S, nInputPorts)) return;
    ssSetInputPortWidth(S, 0, *mxGetPr(WINDOW_LENGTH_PARAM(S)));
    ssSetInputPortDirectFeedThrough(S, 0, 1);
    ssSetInputPortWidth(S, 1, *mxGetPr(WINDOW_LENGTH_PARAM(S)));
    ssSetInputPortDirectFeedThrough(S, 1, 1);

    //Setting Number of Outputs
    if (!ssSetNumOutputPorts(S,nOutputPorts)) return;
    ssSetOutputPortWidth(S, 0, 1);
    ssSetOutputPortWidth(S, 1, 1);
    ssSetOutputPortWidth(S, 2, 1);
    ssSetOutputPortWidth(S, 3, 1);
    ssSetOutputPortWidth(S, 4, 1);

    //Setting Number of Sample Times to 1
    ssSetNumSampleTimes(S, 1);

    //Setting Number of Continuous and descrete States
    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, NDSCSTAT);

    //Setting Number of different Work Vectors
    ssSetNumDWork(S, NDWORK);
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 0);
    ssSetNumPWork(S, 0);
    ssSetNumModes(S, 0);

    //Setting type and Width of Dwork Vectors
    //int N = (int*)mxGetPr(WINDOW_LENGTH_PARAM(S));
    //int M = 3*(N);


    ssSetDWorkWidth(S, HANN_IDX, *mxGetPr(WINDOW_LENGTH_PARAM(S)));
    ssSetDWorkDataType(S, HANN_IDX, SS_DOUBLE);

    /* operating point save/restore compliance to be same as a built-in block */
    ssSetOperatingPointCompliance(S, USE_DEFAULT_OPERATING_POINT);


    /* Set this S-function as runtime thread-safe for multicore execution */
    ssSetRuntimeThreadSafetyCompliance(S, RUNTIME_THREAD_SAFETY_COMPLIANCE_TRUE);
    /* Take care when specifying exception free code - see sfuntmpl.doc */
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);

    
}

/* Function: mdlInitializeSampleTimes =========================================
 *
 *    This function is used to specify the sample time(s) for your S-function.
 *    It must be the same number of sample times as specified in
 *    ssSetNumSampleTimes. 
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, 1/(*mxGetPr(FRAME_PER_SEC(S))) );
    ssSetOffsetTime(S, 0, 0.0);
}

/* Function: mdlInitializeSampleTimes =========================================
 *   This function initializes the Dwork Vectors relative to the DFTVector and the 
 *   WLS Marix used during the execution of the function mdlOutputs.   
 */
#define MDL_START
static void mdlStart(SimStruct *S)
{
    
    const double *n      = mxGetPr(WINDOW_LENGTH_PARAM(S));
    const double *f0     = mxGetPr(NOMINAL_FREQ_PARAM(S));
    const double *fs     = mxGetPr(SAMPLE_RATE(S));
    const double *frms   = mxGetPr(FRAME_PER_SEC(S)); 
    const double *p      = mxGetPr(P(S));
    double *x = ssGetRealDiscStates(S);

    double *hann_coeff = HANN_PTR(S);

    /*Check parameters values and stops the simulation by reporting an error in 
    *case of not meeting the parameters conditions.
    */

    if( *f0 != 50.0 && *f0 != 60.0 ){
        ssSetErrorStatus(S,"[Synchorestimator]Parameter Error: Nominal Frequency value could exclusively be be 50 or 60");
    }
    /*if( ((int)*fs%1.0!=0.0) || (*n%1.0!=0.0) || (*f0%1.0!=0.0) || (*frms%1.0!=0.0) || (*p%1.0!=0.0) ){
        ssSetErrorStatus(S,"[Synchorestimator]Parameter Error: Window Length, Nominal Frequency, Sample Rate, Frames per Second and Number of Iter (P) must be an integer.");
    }*/
    if ( (*fs < 500) || ((*fs/(*f0)) - (long)(*fs/(*f0)) != 0) || (*fs < 2*(*f0)) ||(*fs - (long)(*fs) != 0) ){
         ssSetErrorStatus(S,"[Synchorestimator]Parameter Error: one of the Sampling Frequency conditions not met: fs >= 500, fs must be a multiple of f0, fs> 2*f0, fs must be an integer.");
    }
    if((*n - (long)(*n) != 0) || (*n < (*fs)/(*f0)) || (((*n)*(*f0)/(*fs)) - (long)((*n)*(*f0)/(*fs)) != 0) ){
          ssSetErrorStatus(S,"[Synchorestimator]Parameter Error: one of the Window length conditions not met: WL must be an integer, WL must be a multiple of (fs/f0), WL >= fs/f00\n Check also Nominal Frequency ");
    }
    if(*frms <= 0){
          ssSetErrorStatus(S,"[Synchorestimator]Parameter Error: Frames per second must be greater than zero");
    }
    if(*p <= 0){
        ssSetErrorStatus(S,"[Synchorestimator]Parameter Error: Number of Iter (P) must be greater than zero");
    }

    /*  Initialize DWork vectors */
    norm_fact = hann(hann_coeff, *n);

    /*  Initialize states */
    x[STATE_IDX] = 0.0;
    x[OLD_ROCOF_IDX] = 0.0;
    x[OLD_FREQ_IDX] = *f0;

}

/* Function: mdlOutputs =======================================================
 * Abstract:
 *    In this function, you compute the outputs of your S-function
 *    block. Generally outputs are placed in the output vector(s),
 *    ssGetOutputPortSignal.
 */


static void mdlOutputs(SimStruct *S, int_T tid)
{
    //Getting the pointers to the DWork Vectors
    double *hann_coeff = HANN_PTR(S);

    int i,j,p;
    
    //Getting poiters to access the S Function Parameters
    const double *n      = mxGetPr(WINDOW_LENGTH_PARAM(S));
    const double *f0     = mxGetPr(NOMINAL_FREQ_PARAM(S));
    const double *fs     = mxGetPr(SAMPLE_RATE(S));
    const double *frms   = mxGetPr(FRAME_PER_SEC(S)); 
    const double *P      = mxGetPr(P(S));
    const double *th     = mxGetPr(THRS(S)); 
    const double *lpf    = mxGetPr(LP(S));
    const double *phase_mode    = mxGetPr(PHASE_MODE(S));

    //Getting pointers to access the S Function Satates
    double *x = ssGetRealDiscStates(S);

    //Getting poiter to access the Input Signal Port
    InputRealPtrsType SamplePtrs = ssGetInputPortRealSignalPtrs(S,0);
    InputRealPtrsType FRACSECPtrs = ssGetInputPortRealSignalPtrs(S,1);

    //Getting poiters to access the Output Signal Ports
    double *Mag = ssGetOutputPortRealSignal(S,0);
    double *Phase = ssGetOutputPortRealSignal(S,1);
    double *Freq = ssGetOutputPortRealSignal(S,2);
    double *rocof = ssGetOutputPortRealSignal(S,3);
    double *FrameFRACSEC = ssGetOutputPortRealSignal(S,4);
    
    //Calculated quantities declaration
    double e_ipdft_amp;
    double e_ipdft_ph;
    double e_ipdft_freq;
    double e_ipdft_rocof;

    int k1 = (int)(*f0)*(*n)/(*fs); //index of the highest dft magnitude
    
    //Input windowing
    //printf("%lf ,", *SamplePtrs[0]);
    double *v_win = (double*) malloc(*n * sizeof(double));  //windowed input signal array
    for(i=0; i<(*n); i++){
        v_win[i] = *SamplePtrs[i] * hann_coeff[i];
    }

    //computing the three DFT bins for the windowed input  
    complex v_dft[(int)(*n)]; //dft of the windowed input array
    dft_r(v_win, v_dft , *n, 12 , norm_fact );

    //interpolating the three bins to get the fractional correction term "delta_corr"
    complex v1;
    double delta_corr;
    k1 = e_ipDFT(v_dft, 12, &delta_corr, &v1, *P, *n);

    //extracting the amplitude and phase ------------------------------------------------------------------------------------
    e_ipdft_amp = 2*sqrt(v1.r*v1.r+v1.i*v1.i);
    
    if(v1.r < 0){
    if(v1.i < 0){e_ipdft_ph = atan(v1.i/v1.r) - PI ;}
    else{e_ipdft_ph = atan(v1.i/v1.r) + PI;}}
    else if( v1.i == 0 ){ e_ipdft_ph = 0;}
    else{ e_ipdft_ph = atan(v1.i/v1.r);}
    
    e_ipdft_freq =(delta_corr == 0.0) ? *f0 : (k1+delta_corr)*(*fs)/(*n);
    

    //Two state ROCOF Estimation---------------------------------------------------------------------------------------------
    double rcf = (e_ipdft_freq - x[OLD_FREQ_IDX])*(*frms); 
    double rcf_der = (rcf - x[OLD_ROCOF_IDX])*(*frms);

    x[STATE_IDX] = (!x[STATE_IDX] && ( rcf > th[0] || rcf_der > th[1] )) ? 1.0 : x[STATE_IDX];
    x[STATE_IDX] = (x[STATE_IDX] && rcf < th[2]) ? 0.0 : x[STATE_IDX];

    if(!x[STATE_IDX]){
         e_ipdft_rocof = lpf[1]*rcf + (lpf[2]- lpf[0]) * x[OLD_ROCOF_IDX];  
    }else {e_ipdft_rocof = rcf;}

    x[OLD_FREQ_IDX] = e_ipdft_freq;
    x[OLD_ROCOF_IDX] = e_ipdft_rocof;

    //Phase Correction -----------------------------------------------------------------------------------------------------

    double fracsec= *FRACSECPtrs[(int)((*n)/2)];   //the frame timestamp correponds to the mid point of the estimation window (int)((*n)/2)
    
//     for(i=0; i<*n ; i++){
//        p
//     }
    //printf("--------------\n");
    fracsec = fracsec - ((long)fracsec);  //extracting only the fractional part of timestamp
    //printf("%lf\n", e_ipdft_ph );


    //e_ipdft_ph = -(((long)(e_ipdft_ph/PI))%2)*PI + (e_ipdft_ph/PI - (long)(e_ipdft_ph/PI))*PI;
    if (*phase_mode == 1){ 
        e_ipdft_ph = e_ipdft_ph -2*PI*(*f0)*fracsec; 
        }
    else if(*phase_mode == 2){
        e_ipdft_ph = e_ipdft_ph -2*PI*(e_ipdft_freq)*fracsec;
        }
    
    e_ipdft_ph = -(((long)(e_ipdft_ph/PI))%2)*PI + (e_ipdft_ph/PI - (long)(e_ipdft_ph/PI))*PI;
    //e_ipdft_ph =(e_ipdft_ph/(2*PI) - (long)(e_ipdft_ph/(2*PI)))*2*PI;
    //printf("c = %lf , f = %lf, ph = %lf\n", e_ipdft_ph/e_ipdft_freq, e_ipdft_freq,  e_ipdft_ph );
    
     
    //End--------------------------------------------------------------------------------------------------------------------
    
    //OUTPUT FRAME
    *Mag = e_ipdft_amp;
    *Phase = e_ipdft_ph;
    *Freq = e_ipdft_freq;
    *rocof = e_ipdft_rocof;   
    *FrameFRACSEC = fracsec;

}
static void mdlTerminate(SimStruct *S){}


int e_ipDFT(complex* v_dft, int v_dft_len, double* delta, complex* phasor, int P, int n){


    //computing the magnitude of the DFT and extracting the largest magnitude and its relative index

    int i, p;
    double v_dft_mag[v_dft_len]; //magnitude of dft
    int k1; double dft_max = -INFINITY; 

    for(i=0; i<(v_dft_len); i++){
      v_dft_mag[i] = sqrt((v_dft[i].r*v_dft[i].r) + (v_dft[i].i*v_dft[i].i)); 
      if (v_dft_mag[i] > dft_max){

            dft_max = v_dft_mag[i];
            k1 = i;
        }

    }

    int sigma = ( v_dft_mag[k1+1] > v_dft_mag[k1-1] ) ? 1:-1; //sign of the delta correction
    double delta_corr =2*sigma*(fabs(v_dft_mag[k1+1]-v_dft_mag[k1-1])/(v_dft_mag[k1]*2 +v_dft_mag[k1-1]+ v_dft_mag[k1+1])); 
    
    
    if(delta_corr == 0){

			*phasor = v_dft[k1];
			*delta = 0;
		}
    else{
        
        double a = fabs((delta_corr*delta_corr-1)*(PI*delta_corr)/sin(PI*delta_corr)); 
        complex v_ipdft;  //the interpolated dft phasor estimation
        v_ipdft.r =a*(v_dft[1].r*cos(PI*delta_corr)+v_dft[1].i*sin(PI*delta_corr));
        v_ipdft.i =a*(-v_dft[1].r*sin(PI*delta_corr)+v_dft[1].i*cos(PI*delta_corr));


        //This part is iterated P times
        complex v1 ={v_ipdft.r, v_ipdft.i};
        double e_delta_corr = delta_corr;
        complex e_ipdft_3max[3]={v_dft[k1-1], v_dft[k1] , v_dft[k1+1]};

        double e_ipdft_3mag[3];
        complex v_e_ipdft_max; 
        double e_a; 
        complex hann_ft; 
        complex e_ipdft_3max_new[3];

        for(p=0 ; p<P ; p++){ //e-ipdft iterations-----------------------------------------------------------------------------

            for(i=0; i<3; i++){

                hann_ft = hann_frt((i-1) + e_delta_corr + 2*(double)k1 , n); 

                e_ipdft_3max_new[i].r = e_ipdft_3max[i].r-(v1.r*hann_ft.r + v1.i*hann_ft.i)/norm_fact; 
                e_ipdft_3max_new[i].i = e_ipdft_3max[i].i+(v1.i*hann_ft.r - v1.r*hann_ft.i)/norm_fact;
                e_ipdft_3mag[i] = sqrt((e_ipdft_3max_new[i].r*e_ipdft_3max_new[i].r)+(e_ipdft_3max_new[i].i*e_ipdft_3max_new[i].i));


            }

            //interpolating the three bins to get the fractional correction term "delta_corr"
            e_delta_corr =2*sigma*(fabs(e_ipdft_3mag[2]-e_ipdft_3mag[0])/(e_ipdft_3mag[1]*2 +e_ipdft_3mag[0]+ e_ipdft_3mag[2])); 

            v_e_ipdft_max = e_ipdft_3max_new[1];
            e_a = (1-e_delta_corr*e_delta_corr)*fabs((PI*e_delta_corr)/sin(PI*e_delta_corr)); 
            v1.r =e_a*(v_e_ipdft_max.r*cos(PI*e_delta_corr)+v_e_ipdft_max.i*sin(PI*e_delta_corr));
            v1.i =e_a*(-v_e_ipdft_max.r*sin(PI*e_delta_corr)+v_e_ipdft_max.i*cos(PI*e_delta_corr));

            *phasor = v1;
            *delta = e_delta_corr;

        }
    }

    return k1;
}

complex hann_frt(double k, double N){

    complex w;
    w.r =cos(PI*k*(N-1)/N)*0.5*sin(PI*k)/sin(PI*k/N)   
        -cos(PI*(k+1)*(N-1)/N)*0.25*sin(PI*(k+1))/sin(PI*(k+1)/N)
        -cos(PI*(k-1)*(N-1)/N)*0.25*sin(PI*(k-1))/sin(PI*(k-1)/N);

    w.i = -sin(PI*k*(N-1)/N)*0.5*sin(PI*k)/sin(PI*k/N)   
        +sin(PI*(k+1)*(N-1)/N)*0.25*sin(PI*(k+1))/sin(PI*(k+1)/N)
        +sin(PI*(k-1)*(N-1)/N)*0.25*sin(PI*(k-1))/sin(PI*(k-1)/N);

    return w;

}
double hann(double* out_ptr, unsigned int out_len){
    
    double norm_fact =0;
    int i=0;
    for (i=0; i < out_len; i++){
 	   out_ptr[i] = 0.5*(1-cos(2*PI*i/out_len));
       norm_fact += out_ptr[i]; 
    }

    return norm_fact;

}

int dft_r(double* in_ptr, complex* out_ptr , unsigned int out_len, int n_bins, double norm_factor ){

    int k,n;
    for (k = 0 ; k < n_bins ; ++k)
    {
        // Real part of X[k]
        out_ptr[k].r = 0;
        for (n=0 ; n<out_len ; ++n) out_ptr[k].r += (in_ptr[n] * cos(n * k * PI*2 / out_len))/norm_factor;
         
        // Imaginary part of X[k]
        out_ptr[k].i = 0;
        for (n=0 ; n<out_len ; ++n) out_ptr[k].i -= (in_ptr[n] * sin(n * k * PI*2 / out_len))/norm_factor;

    }

    return 0;

}



#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif

