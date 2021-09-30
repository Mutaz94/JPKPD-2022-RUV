Wed Sep 29 20:47:13 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	Two-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/7/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/spa1/B/dat3.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER NSIG=2
$PK
ET1 = EXP(ETA(1)*THETA(6))
ET2 = EXP(ETA(2)*THETA(7))
ET3 = EXP(ETA(3)*THETA(8))
ET4 = EXP(ETA(4)*THETA(9))
ET5 = EXP(ETA(5)*THETA(10))

CL = 5.0 * THETA(1) * ET1
V2 = 35  * THETA(2) * ET2
Q  = 50  * THETA(3) * ET3
V3 = 50  * THETA(4) * ET4
KA = 0.7 * THETA(5) * ET5
SC = V2
$ERROR
CVERR = 0.05
W = THETA(11)*F*CVERR

Y 	= F + W*ERR(1)

$THETA
(0,1) ; CL
(0,1) ; V2
(0,1) ; Q
(0,1) ; V3
(0,1) ; KA
(0,1) ; IIVCL
(0,1) ; IIVV2
(0,1) ; IIVQ
(0,1) ; IIVV3
(0,1) ; IIVKA
(0,1) ; CVPropErr

$OMEGA  (0.09 FIX)x5
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       29 SEP 2021
Days until program expires : 200
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 template control stream
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:      600
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      500
 TOT. NO. OF INDIVIDUALS:      100
0LENGTH OF THETA:  11
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   5
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
0INITIAL ESTIMATE OF OMEGA:
 0.9000E-01
 0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
0OMEGA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0INITIAL ESTIMATE OF SIGMA:
 0.1000E+01
0SIGMA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:     NOSLOW
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 Cholesky Transposition of R Matrix (CHOLROFF):0
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0
 SIMPLE POS. DEFINITE SETTING: (POSDEF):-1
1DOUBLE PRECISION PREDPP VERSION 7.5.0

 TWO COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN4)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   5
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   BASIC PK PARAMETER NO.  1: ELIMINATION RATE (K)
   BASIC PK PARAMETER NO.  2: CENTRAL-TO-PERIPH. RATE (K23)
   BASIC PK PARAMETER NO.  3: PERIPH.-TO-CENTRAL RATE (K32)
   BASIC PK PARAMETER NO.  5: ABSORPTION RATE (KA)
 TRANSLATOR WILL CONVERT PARAMETERS
 CL, V2, Q, V3 TO K, K23, K32 (TRANS4)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         PERIPH.      ON         NO         YES        NO         NO
    4         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            6           *           *           *           *
    3            *           *           *           *           *
    4            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      6
   TIME DATA ITEM IS DATA ITEM NO.:          2
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   4

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.
1


 #TBLN:      1
 #METH: First Order Conditional Estimation with Interaction

 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            10000
 NO. OF SIG. FIGURES REQUIRED:            2
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      100
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     100
 NOPRIOR SETTING (NOPRIOR):                 0
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          1
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      0
 RAW OUTPUT FILE (FILE): m3.ext
 EXCLUDE TITLE (NOTITLE):                   NO
 EXCLUDE COLUMN LABELS (NOLABEL):           NO
 FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL
 KNUTHSUMOFF:                               0
 INCLUDE LNTWOPI:                           NO
 INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO
 INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO
 ADDITIONAL CONVERGENCE TEST (CTYPE=4)?:    NO
 EM OR BAYESIAN METHOD USED:                 NONE


 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI

 MONITORING OF SEARCH:


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2149.46938035431        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0504E+02 -3.7281E+01 -4.1960E+01  1.4974E+01  4.6923E+01  6.9736E+01  7.4091E+00  1.4512E+01  1.7744E+01  2.8529E+01
            -3.9453E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2157.43065858128        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      192
 NPARAMETR:  1.0201E+00  1.0838E+00  1.1909E+00  1.0101E+00  1.0619E+00  1.0082E+00  9.4490E-01  8.9857E-01  9.2612E-01  7.7015E-01
             1.0877E+00
 PARAMETER:  1.1987E-01  1.8047E-01  2.7473E-01  1.1006E-01  1.6010E-01  1.0821E-01  4.3324E-02 -6.9472E-03  2.3253E-02 -1.6117E-01
             1.8408E-01
 GRADIENT:  -8.3555E+00  4.3603E+00  2.2787E+01 -1.5151E+01  9.3479E+00  1.6003E+00 -3.7437E+00 -8.6171E+00 -1.3419E+01 -1.3138E+01
             1.0073E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2158.25371037371        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  1.0236E+00  1.0798E+00  9.2796E-01  1.0017E+00  9.3657E-01  9.9390E-01  1.1331E+00  7.8765E-01  8.9953E-01  6.0073E-01
             1.0810E+00
 PARAMETER:  1.2333E-01  1.7673E-01  2.5232E-02  1.0174E-01  3.4464E-02  9.3877E-02  2.2500E-01 -1.3870E-01 -5.8819E-03 -4.0962E-01
             1.7792E-01
 GRADIENT:  -3.2122E+00  1.3584E+01  1.3227E+01 -6.9521E+00 -3.4887E+00 -4.6240E+00  6.1807E+00 -1.4861E+00 -1.0208E+01 -1.1787E+01
             7.9273E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2160.71341875727        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      546
 NPARAMETR:  1.0253E+00  1.1782E+00  8.3177E-01  9.3401E-01  9.4945E-01  1.0068E+00  9.4851E-01  5.5097E-01  1.0075E+00  7.5909E-01
             1.0570E+00
 PARAMETER:  1.2499E-01  2.6396E-01 -8.4202E-02  3.1735E-02  4.8131E-02  1.0673E-01  4.7139E-02 -4.9608E-01  1.0747E-01 -1.7563E-01
             1.5542E-01
 GRADIENT:  -8.4190E-01 -2.5326E+00 -2.6794E-02 -1.7079E+00 -3.0658E-01  3.0115E-02 -2.0492E-01  3.0116E-01  5.6497E-01  9.1869E-01
            -1.1546E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2160.80981767016        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      722
 NPARAMETR:  1.0272E+00  1.3289E+00  7.3326E-01  8.3838E-01  9.7530E-01  1.0074E+00  8.8879E-01  3.9065E-01  1.0807E+00  7.5068E-01
             1.0620E+00
 PARAMETER:  1.2684E-01  3.8435E-01 -2.1026E-01 -7.6286E-02  7.4995E-02  1.0740E-01 -1.7891E-02 -8.3993E-01  1.7760E-01 -1.8678E-01
             1.6016E-01
 GRADIENT:   6.7404E-01  2.1950E+00 -1.4972E+00  2.1148E+00 -5.3316E-01 -3.1944E-01  6.2417E-01  4.2116E-01 -4.4927E-01  4.6503E-01
             8.6135E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2160.84401424703        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      901
 NPARAMETR:  1.0289E+00  1.3856E+00  7.0960E-01  7.9901E-01  9.9439E-01  1.0101E+00  8.5706E-01  2.1875E-01  1.1220E+00  7.5497E-01
             1.0619E+00
 PARAMETER:  1.2849E-01  4.2617E-01 -2.4305E-01 -1.2439E-01  9.4371E-02  1.1009E-01 -5.4251E-02 -1.4198E+00  2.1512E-01 -1.8108E-01
             1.6002E-01
 GRADIENT:   3.8138E+00  1.5069E+00  5.9202E+00 -3.0635E+00 -3.4050E+00  6.3289E-01 -1.1147E+00 -5.7126E-02 -8.5846E-01 -1.6033E+00
            -7.5222E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2160.90657015541        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1076
 NPARAMETR:  1.0266E+00  1.3955E+00  6.8508E-01  7.9301E-01  9.8869E-01  1.0081E+00  8.5994E-01  1.3756E-01  1.1277E+00  7.5740E-01
             1.0614E+00
 PARAMETER:  1.2628E-01  4.3323E-01 -2.7822E-01 -1.3192E-01  8.8629E-02  1.0808E-01 -5.0889E-02 -1.8837E+00  2.2021E-01 -1.7787E-01
             1.5961E-01
 GRADIENT:  -1.4367E+00 -9.7723E-01 -2.0504E+00  1.7331E+00  1.7431E+00 -2.7043E-01  3.2178E-01  5.5475E-02  1.4377E-01  6.1926E-01
             4.9353E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2160.92480126796        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1252
 NPARAMETR:  1.0260E+00  1.3987E+00  6.8242E-01  7.8967E-01  9.8878E-01  1.0077E+00  8.5862E-01  4.5257E-02  1.1307E+00  7.5474E-01
             1.0616E+00
 PARAMETER:  1.2562E-01  4.3551E-01 -2.8211E-01 -1.3614E-01  8.8718E-02  1.0768E-01 -5.2431E-02 -2.9954E+00  2.2285E-01 -1.8138E-01
             1.5975E-01
 GRADIENT:  -2.9358E+00 -1.3978E+00 -2.5019E-01 -1.4165E-01  3.9956E-01 -4.4170E-01  5.5825E-02  4.6574E-03  1.9541E-02  8.8321E-02
            -2.9688E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2160.92913148555        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1428
 NPARAMETR:  1.0263E+00  1.3984E+00  6.8167E-01  7.8974E-01  9.8804E-01  1.0086E+00  8.6009E-01  1.0000E-02  1.1297E+00  7.5311E-01
             1.0615E+00
 PARAMETER:  1.2598E-01  4.3533E-01 -2.8321E-01 -1.3605E-01  8.7970E-02  1.0859E-01 -5.0718E-02 -5.6996E+00  2.2199E-01 -1.8355E-01
             1.5972E-01
 GRADIENT:  -2.1476E+00 -1.1281E+00 -9.7876E-02 -2.7336E-01  1.5701E-01 -7.6709E-02  1.7660E-01  0.0000E+00 -8.8691E-02 -1.6157E-02
            -9.3973E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2160.93435114334        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     1588
 NPARAMETR:  1.0284E+00  1.3970E+00  6.8138E-01  7.8997E-01  9.8739E-01  1.0093E+00  8.5877E-01  1.0000E-02  1.1298E+00  7.5267E-01
             1.0617E+00
 PARAMETER:  1.2796E-01  4.3435E-01 -2.8363E-01 -1.3576E-01  8.7305E-02  1.0929E-01 -5.2252E-02 -7.7637E+00  2.2205E-01 -1.8413E-01
             1.5983E-01
 GRADIENT:   5.2047E+02  2.8977E+02  4.0530E+00  6.5195E+01  7.3203E+00  6.7963E+01  5.0428E+00  0.0000E+00  1.3018E+01  6.5099E-01
             1.5894E+00

0ITERATION NO.:   49    OBJECTIVE VALUE:  -2160.93450048481        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:     1725
 NPARAMETR:  1.0286E+00  1.3967E+00  6.8141E-01  7.9032E-01  9.8765E-01  1.0096E+00  8.5877E-01  1.0000E-02  1.1298E+00  7.5287E-01
             1.0617E+00
 PARAMETER:  1.2819E-01  4.3496E-01 -2.8310E-01 -1.3506E-01  8.7070E-02  1.0961E-01 -5.1252E-02 -7.7637E+00  2.2228E-01 -1.8367E-01
             1.5984E-01
 GRADIENT:   5.3242E-02  1.1528E+00  2.3633E-01  2.3184E-01 -7.9032E-01  1.5190E-02  5.4339E-02  0.0000E+00  3.8958E-02  1.7207E-02
            -2.3255E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1725
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.8608E-04 -1.7390E-02 -3.0999E-04  1.1967E-02 -2.2851E-02
 SE:             2.9876E-02  2.2461E-02  1.5838E-04  2.4715E-02  2.1466E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9236E-01  4.3879E-01  5.0314E-02  6.2824E-01  2.8710E-01

 ETASHRINKSD(%)  1.0000E-10  2.4752E+01  9.9469E+01  1.7201E+01  2.8087E+01
 ETASHRINKVR(%)  1.0000E-10  4.3378E+01  9.9997E+01  3.1444E+01  4.8285E+01
 EBVSHRINKSD(%)  3.5958E-01  2.4440E+01  9.9522E+01  1.7405E+01  2.8151E+01
 EBVSHRINKVR(%)  7.1788E-01  4.2907E+01  9.9998E+01  3.1780E+01  4.8377E+01
 RELATIVEINF(%)  9.9167E+01  3.9441E+00  3.8903E-04  5.7130E+00  6.3997E+00
 EPSSHRINKSD(%)  3.2718E+01
 EPSSHRINKVR(%)  5.4732E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2160.9345004848142     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1241.9959672801415     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.85
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.86
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2160.935       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.40E+00  6.82E-01  7.91E-01  9.87E-01  1.01E+00  8.60E-01  1.00E-02  1.13E+00  7.53E-01  1.06E+00
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        9.00E-02
 
 ETA2
+        0.00E+00  9.00E-02
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-02
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  9.00E-02
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  9.00E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        3.00E-01
 
 ETA2
+        0.00E+00  3.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  3.00E-01
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  3.00E-01
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  3.00E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.02E+03
 
 TH 2
+       -6.13E+00  4.57E+02
 
 TH 3
+        6.35E+00  2.45E+02  6.55E+02
 
 TH 4
+       -5.21E+00  3.26E+02 -3.15E+02  9.57E+02
 
 TH 5
+        1.22E+00 -3.66E+02 -6.13E+02  3.57E+02  1.03E+03
 
 TH 6
+        1.88E-01 -1.37E+00  2.23E+00 -1.69E+00  3.04E-01  1.93E+02
 
 TH 7
+        1.06E+00  2.17E+01 -1.62E+01 -8.55E+00 -9.72E+00  6.31E-01  9.38E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.65E-01 -2.53E+01 -2.32E+01  4.64E+01 -6.35E+00 -3.03E-01  1.68E+01  0.00E+00  8.33E+01
 
 TH10
+        9.63E-01 -7.37E+00 -5.13E+01 -1.82E+01 -5.73E+01  5.59E-01  3.01E+01  0.00E+00  7.63E+00  9.69E+01
 
 TH11
+       -7.53E+00 -1.99E+01 -2.57E+01 -5.17E+00 -5.86E+00  9.15E-01  7.65E+00  0.00E+00  6.01E+00  2.78E+01  3.67E+02
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       32.776
Stop Time:
Wed Sep 29 20:47:47 CDT 2021
