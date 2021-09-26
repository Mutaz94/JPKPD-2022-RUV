Sat Sep 25 11:42:46 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat43.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
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

$OMEGA  0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       25 SEP 2021
Days until program expires : 204
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
 NO. OF DATA RECS IN DATA SET:      500
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

 TOT. NO. OF OBS RECS:      400
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
 NO. OF SIG. FIGURES REQUIRED:            3
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
 RAW OUTPUT FILE (FILE): m43.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1611.26316388170        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.7169E+01 -5.1211E+01 -1.9394E+01 -2.9305E+01  7.8161E+01  6.1500E+00  5.6192E+00 -2.2905E+00  2.9175E+01 -3.5864E+01
            -6.7129E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1622.25409001884        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.7841E-01  1.0249E+00  9.7409E-01  1.0252E+00  9.3649E-01  9.8125E-01  9.3199E-01  1.0023E+00  8.0103E-01  1.1558E+00
             1.1672E+00
 PARAMETER:  7.8173E-02  1.2457E-01  7.3752E-02  1.2490E-01  3.4388E-02  8.1076E-02  2.9567E-02  1.0230E-01 -1.2185E-01  2.4482E-01
             2.5463E-01
 GRADIENT:  -7.2090E+00  1.8199E+01  1.6013E+00  1.4376E+01 -1.5475E+01 -1.4140E+00 -3.4361E+00  1.7035E+00 -8.7860E+00  5.9418E+00
             6.1010E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1622.45744838364        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      268
 NPARAMETR:  9.8214E-01  8.9639E-01  1.1043E+00  1.1204E+00  9.3117E-01  9.8327E-01  9.8445E-01  9.6605E-01  7.8774E-01  1.2291E+00
             1.1661E+00
 PARAMETER:  8.1976E-02 -9.3813E-03  1.9923E-01  2.1367E-01  2.8685E-02  8.3130E-02  8.4328E-02  6.5461E-02 -1.3858E-01  3.0632E-01
             2.5370E-01
 GRADIENT:  -2.4319E+01  2.9095E+01  1.1264E+01  2.8293E+01 -3.0475E+01 -3.1216E+00 -1.1333E+00 -2.0320E+00 -6.3077E+00  1.0337E+01
             5.7020E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1624.36128071443        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      444
 NPARAMETR:  9.9142E-01  8.2031E-01  1.4599E+00  1.1681E+00  1.0584E+00  9.9005E-01  8.0541E-01  1.4337E+00  8.6226E-01  1.2923E+00
             1.1510E+00
 PARAMETER:  9.1382E-02 -9.8070E-02  4.7837E-01  2.5539E-01  1.5677E-01  9.0003E-02 -1.1640E-01  4.6023E-01 -4.8199E-02  3.5643E-01
             2.4066E-01
 GRADIENT:   2.7494E+00  6.8577E+00 -1.9871E-01  9.9562E+00 -4.2361E+00  6.5803E-01  3.2119E+00  2.5057E+00  7.1069E+00  1.7080E+00
             1.3497E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1626.36271148035        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      620
 NPARAMETR:  9.8684E-01  5.4220E-01  1.5451E+00  1.3397E+00  9.9266E-01  9.8553E-01  2.9807E-01  1.3715E+00  7.7749E-01  1.2563E+00
             1.1542E+00
 PARAMETER:  8.6752E-02 -5.1211E-01  5.3509E-01  3.9247E-01  9.2630E-02  8.5425E-02 -1.1104E+00  4.1590E-01 -1.5169E-01  3.2817E-01
             2.4339E-01
 GRADIENT:  -1.3751E+00  6.2165E+00  1.1538E+00  1.2577E+01 -3.7434E+00  1.4328E-01  1.3043E-01  2.8714E-01 -1.9388E+00 -8.3736E-01
             1.9580E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1626.95750138136        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      795
 NPARAMETR:  9.8438E-01  3.4721E-01  1.6095E+00  1.4605E+00  9.5746E-01  9.8309E-01  6.1272E-02  1.3983E+00  7.1362E-01  1.2481E+00
             1.1474E+00
 PARAMETER:  8.4261E-02 -9.5782E-01  5.7592E-01  4.7875E-01  5.6529E-02  8.2947E-02 -2.6924E+00  4.3527E-01 -2.3741E-01  3.2161E-01
             2.3754E-01
 GRADIENT:  -8.6000E-01  1.5030E+00  5.6785E-01  6.0206E+00 -6.9533E-01  4.1095E-02  3.8635E-03 -4.3580E-01 -2.1284E-01 -7.2944E-02
             6.9575E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1627.05724823106        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      970
 NPARAMETR:  9.8290E-01  2.2805E-01  1.6653E+00  1.5343E+00  9.4040E-01  9.8153E-01  1.0000E-02  1.4618E+00  6.7549E-01  1.2406E+00
             1.1444E+00
 PARAMETER:  8.2755E-02 -1.3782E+00  6.1000E-01  5.2810E-01  3.8550E-02  8.1353E-02 -4.5565E+00  4.7967E-01 -2.9231E-01  3.1559E-01
             2.3486E-01
 GRADIENT:  -4.5689E-02  1.5497E-01 -6.2484E-02  6.6900E-01 -1.8265E-01 -2.0532E-02  0.0000E+00 -5.9129E-03 -3.2243E-01  2.7926E-01
             9.7807E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1627.06401746322        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1145
 NPARAMETR:  9.8201E-01  1.6673E-01  1.7195E+00  1.5744E+00  9.3862E-01  9.8085E-01  1.0000E-02  1.5152E+00  6.5769E-01  1.2378E+00
             1.1440E+00
 PARAMETER:  8.1844E-02 -1.6914E+00  6.4203E-01  5.5387E-01  3.6654E-02  8.0660E-02 -6.0927E+00  5.1557E-01 -3.1902E-01  3.1334E-01
             2.3455E-01
 GRADIENT:   9.5071E-02  6.9526E-02  4.7517E-02  7.0690E-01 -3.0974E-02  9.8190E-04  0.0000E+00 -2.7972E-02 -9.0195E-02 -5.1232E-02
            -3.9268E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1627.06442882622        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1321
 NPARAMETR:  9.8162E-01  1.4556E-01  1.7307E+00  1.5878E+00  9.3602E-01  9.8060E-01  1.0000E-02  1.5270E+00  6.5178E-01  1.2368E+00
             1.1440E+00
 PARAMETER:  8.1451E-02 -1.8272E+00  6.4854E-01  5.6236E-01  3.3883E-02  8.0408E-02 -6.7901E+00  5.2332E-01 -3.2804E-01  3.1252E-01
             2.3451E-01
 GRADIENT:  -2.1776E-02  4.6791E-02 -1.1357E-02  5.9480E-01 -7.3877E-02  1.1542E-03  0.0000E+00  5.5010E-03 -4.1819E-02  1.5236E-02
            -8.5298E-03

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1627.06448183217        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1503
 NPARAMETR:  9.8153E-01  1.3879E-01  1.7357E+00  1.5919E+00  9.3570E-01  9.8052E-01  1.0000E-02  1.5318E+00  6.5000E-01  1.2366E+00
             1.1440E+00
 PARAMETER:  8.1355E-02 -1.8748E+00  6.5143E-01  5.6490E-01  3.3535E-02  8.0330E-02 -7.0281E+00  5.2643E-01 -3.3078E-01  3.1233E-01
             2.3452E-01
 GRADIENT:   1.9260E-02 -2.0399E-02 -2.5785E-02 -4.2435E-01  8.2441E-02  3.7700E-03  0.0000E+00  4.3996E-04  2.6760E-02 -6.1499E-03
             9.9983E-03

0ITERATION NO.:   49    OBJECTIVE VALUE:  -1627.06449073677        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     1664
 NPARAMETR:  9.8153E-01  1.3894E-01  1.7359E+00  1.5919E+00  9.3560E-01  9.8051E-01  1.0000E-02  1.5318E+00  6.4997E-01  1.2366E+00
             1.1440E+00
 PARAMETER:  8.1355E-02 -1.8743E+00  6.5141E-01  5.6491E-01  3.3535E-02  8.0329E-02 -7.0281E+00  5.2646E-01 -3.3079E-01  3.1230E-01
             2.3451E-01
 GRADIENT:   5.4947E-03 -1.1372E-02 -2.9200E-02 -1.4647E-01  4.6960E-02  3.5343E-03  0.0000E+00  2.7448E-03  1.1894E-02 -9.3980E-03
             4.7925E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1664
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7165E-04 -5.1385E-05 -3.5804E-02 -9.6388E-03 -4.4372E-02
 SE:             2.9818E-02  2.6004E-05  1.7176E-02  2.8689E-02  2.1319E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9541E-01  4.8151E-02  3.7111E-02  7.3689E-01  3.7399E-02

 ETASHRINKSD(%)  1.0457E-01  9.9913E+01  4.2458E+01  3.8884E+00  2.8580E+01
 ETASHRINKVR(%)  2.0903E-01  1.0000E+02  6.6890E+01  7.6257E+00  4.8991E+01
 EBVSHRINKSD(%)  5.2495E-01  9.9914E+01  4.6861E+01  4.1883E+00  2.4385E+01
 EBVSHRINKVR(%)  1.0471E+00  1.0000E+02  7.1763E+01  8.2013E+00  4.2824E+01
 RELATIVEINF(%)  9.7174E+01  3.4207E-06  7.7089E+00  4.7762E+00  1.1428E+01
 EPSSHRINKSD(%)  4.4127E+01
 EPSSHRINKVR(%)  6.8782E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1627.0644907367730     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -891.91366417303482     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    20.75
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.83
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1627.064       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.82E-01  1.39E-01  1.74E+00  1.59E+00  9.36E-01  9.81E-01  1.00E-02  1.53E+00  6.50E-01  1.24E+00  1.14E+00
 


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
+        1.18E+03
 
 TH 2
+       -3.20E+01  4.39E+02
 
 TH 3
+        9.08E-02  1.43E+01  4.86E+01
 
 TH 4
+       -1.72E+01  6.00E+02 -2.34E+01  9.54E+02
 
 TH 5
+        1.23E-02 -1.60E+02 -1.12E+02 -6.46E+01  4.90E+02
 
 TH 6
+        2.24E+00 -2.14E+00  1.16E+00 -3.45E+00  2.11E+00  2.08E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        7.46E-01 -1.46E+00 -1.75E+01 -4.77E+00 -5.41E+00  4.02E-01  0.00E+00  2.07E+01
 
 TH 9
+        2.37E+00 -1.12E+02  5.14E+00 -2.37E+00 -7.21E-02 -5.71E+00  0.00E+00 -3.31E-01  3.98E+02
 
 TH10
+        1.69E-02  5.83E+00 -8.63E-01 -8.60E-01 -5.54E+01  7.56E-02  0.00E+00  6.60E+00  1.81E+00  5.20E+01
 
 TH11
+       -8.39E+00 -1.57E+01 -4.87E+00 -1.56E+01 -2.00E+00  2.14E+00  0.00E+00  4.00E+00  2.12E+01  8.18E+00  1.61E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       26.636
Stop Time:
Sat Sep 25 11:43:14 CDT 2021
