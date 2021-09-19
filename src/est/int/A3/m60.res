Sat Sep 18 01:50:04 CDT 2021
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
$DATA ../../../../data/int/A3/dat60.csv ignore=@
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
Current Date:       18 SEP 2021
Days until program expires : 211
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m60.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   112.872621846424        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.5079E+01  2.4622E+02  3.2965E+02 -9.7421E+01  3.4050E+02  2.6097E+01 -2.7555E+02 -2.8630E+02 -1.9592E+01 -2.7001E+02
            -7.0787E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2376.34479021365        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0368E+00  8.7702E-01  8.0226E-01  1.2221E+00  7.7796E-01  7.8817E-01  1.0647E+00  1.0334E+00  8.0688E-01  9.8106E-01
             5.1585E+00
 PARAMETER:  1.3619E-01 -3.1224E-02 -1.2032E-01  3.0056E-01 -1.5108E-01 -1.3804E-01  1.6273E-01  1.3285E-01 -1.1458E-01  8.0878E-02
             1.7406E+00
 GRADIENT:  -3.0887E+01 -9.5577E+00 -2.3550E+01  2.7741E+01  2.2624E+01 -2.8211E+01  5.0969E+00  1.2793E+01  1.6154E+01  2.5017E+01
             8.0059E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2417.10275888353        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.0150E+00  4.9097E-01  3.4008E-01  1.4676E+00  3.4099E-01  7.7527E-01  1.2496E+00  6.3599E-01  8.5171E-01  4.7334E-01
             4.6447E+00
 PARAMETER:  1.1494E-01 -6.1138E-01 -9.7857E-01  4.8365E-01 -9.7591E-01 -1.5454E-01  3.2285E-01 -3.5257E-01 -6.0507E-02 -6.4794E-01
             1.6357E+00
 GRADIENT:  -9.8458E+01  1.7408E+02  1.1240E+01  3.4921E+02 -2.1241E+02 -5.6578E+01 -2.9553E+01  6.4254E+00 -4.9179E+01 -1.5583E+01
             6.8174E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2729.18371572036        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  9.8580E-01  3.0101E-01  2.1121E-01  1.2067E+00  2.3476E-01  8.1403E-01  1.2549E+00  3.3371E-02  1.2989E+00  7.3517E-01
             2.4365E+00
 PARAMETER:  8.5695E-02 -1.1006E+00 -1.4549E+00  2.8791E-01 -1.3492E+00 -1.0575E-01  3.2702E-01 -3.3001E+00  3.6150E-01 -2.0765E-01
             9.9056E-01
 GRADIENT:  -4.3412E+01  1.1230E+02  3.9373E+01  1.1127E+02 -6.3753E+01 -5.6370E+01 -8.5727E+00 -3.4306E-02  2.9832E+01 -1.3634E+01
            -1.3273E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2754.37254338345        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      296
 NPARAMETR:  1.0123E+00  2.0037E-01  1.2950E-01  9.6043E-01  1.7154E-01  1.0041E+00  1.2417E+00  4.5740E-02  1.2917E+00  7.5854E-01
             2.5405E+00
 PARAMETER:  1.1225E-01 -1.5076E+00 -1.9441E+00  5.9621E-02 -1.6630E+00  1.0406E-01  3.1645E-01 -2.9848E+00  3.5599E-01 -1.7635E-01
             1.0324E+00
 GRADIENT:   2.2803E+01  3.0779E+00 -2.2398E+00  3.4238E+01 -9.2045E+00  2.7410E+01 -3.9158E+00 -7.3926E-02 -3.0725E+01 -1.9002E+00
             1.2207E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2758.65006339159        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      433
 NPARAMETR:  1.0037E+00  2.1195E-01  1.4058E-01  9.6350E-01  1.8302E-01  9.3173E-01  1.2878E+00  5.2964E-02  1.3594E+00  7.6478E-01
             2.5243E+00
 PARAMETER:  1.0371E-01 -1.4514E+00 -1.8620E+00  6.2821E-02 -1.5981E+00  2.9287E-02  3.5296E-01 -2.8382E+00  4.0708E-01 -1.6817E-01
             1.0260E+00
 GRADIENT:   1.0448E+00  4.3656E+00 -3.5062E+00 -3.8713E+00 -4.1905E-01  1.2518E+00  2.4279E+00 -8.3285E-02 -1.6024E+00  2.0106E+00
             1.4340E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2760.18393809074        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      614
 NPARAMETR:  9.9920E-01  2.1104E-01  1.3907E-01  9.6066E-01  1.8054E-01  9.0511E-01  1.2250E+00  3.9716E-01  1.3642E+00  7.4871E-01
             2.5279E+00
 PARAMETER:  9.9204E-02 -1.4557E+00 -1.8728E+00  5.9865E-02 -1.6118E+00  2.9721E-04  3.0295E-01 -8.2341E-01  4.1056E-01 -1.8940E-01
             1.0274E+00
 GRADIENT:  -1.1155E+01  1.0342E+01  4.6043E+00 -5.4379E+00 -1.8914E+01 -9.7200E+00 -6.0253E+00 -2.8497E+00 -1.3807E+00  2.4711E+00
             1.5159E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2762.97211595909        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      792
 NPARAMETR:  1.0012E+00  2.0729E-01  1.3797E-01  9.7449E-01  1.8038E-01  9.2268E-01  1.2958E+00  8.0459E-01  1.3922E+00  6.7594E-01
             2.4738E+00
 PARAMETER:  1.0118E-01 -1.4737E+00 -1.8807E+00  7.4157E-02 -1.6127E+00  1.9529E-02  3.5914E-01 -1.1742E-01  4.3088E-01 -2.9165E-01
             1.0058E+00
 GRADIENT:  -3.0191E+00  1.7816E+00  2.3515E+00 -1.4221E+00 -4.0634E+00 -2.9140E+00 -1.7506E+00  5.5226E-02 -1.8804E-02 -7.9155E-01
             3.1539E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2763.00511668567        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      971            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0023E+00  2.0648E-01  1.3730E-01  9.7478E-01  1.8015E-01  9.3006E-01  1.3092E+00  8.0599E-01  1.3966E+00  6.7763E-01
             2.4689E+00
 PARAMETER:  1.0227E-01 -1.4776E+00 -1.8856E+00  7.4457E-02 -1.6140E+00  2.7498E-02  3.6939E-01 -1.1569E-01  4.3405E-01 -2.8916E-01
             1.0038E+00
 GRADIENT:   5.0016E+00  6.6030E+00  7.4115E+00  8.4570E-01  4.8417E+01  3.4055E-01  2.8806E-01  2.8027E-02  7.8244E-01  3.3701E-01
             1.8167E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2763.00511941792        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1152
 NPARAMETR:  1.0023E+00  2.0647E-01  1.3730E-01  9.7481E-01  1.8014E-01  9.2994E-01  1.3092E+00  8.0580E-01  1.3966E+00  6.7765E-01
             2.4689E+00
 PARAMETER:  1.0226E-01 -1.4776E+00 -1.8856E+00  7.4495E-02 -1.6140E+00  2.7470E-02  3.6940E-01 -1.1580E-01  4.3401E-01 -2.8914E-01
             1.0038E+00
 GRADIENT:   1.5067E-02 -6.7162E-03  9.0366E-03  3.2161E-03  4.6463E-03  1.1835E-02 -2.6813E-03  2.8177E-04 -8.8381E-03 -2.2596E-03
            -2.6136E-03
 NUMSIGDIG:         4.2         5.3         5.3         4.5         6.0         3.0         4.4         2.9         4.1         4.3
                    5.9

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1152
 NO. OF SIG. DIGITS IN FINAL EST.:  3.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4231E-03  1.3139E-02  1.3726E-02 -4.5126E-03  1.5798E-02
 SE:             2.9368E-02  2.4082E-02  1.6241E-02  2.6604E-02  2.3640E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6135E-01  5.8534E-01  3.9803E-01  8.6531E-01  5.0395E-01

 ETASHRINKSD(%)  1.6148E+00  1.9323E+01  4.5590E+01  1.0873E+01  2.0802E+01
 ETASHRINKVR(%)  3.2035E+00  3.4912E+01  7.0396E+01  2.0564E+01  3.7276E+01
 EBVSHRINKSD(%)  1.7629E+00  1.8204E+01  4.4376E+01  8.1734E+00  2.2151E+01
 EBVSHRINKVR(%)  3.4948E+00  3.3094E+01  6.9059E+01  1.5679E+01  3.9396E+01
 RELATIVEINF(%)  9.6494E+01  1.9648E+01  4.0342E+00  3.8131E+01  6.9125E+00
 EPSSHRINKSD(%)  1.9871E+01
 EPSSHRINKVR(%)  3.5794E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2763.0051194179182     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1108.9157596495074     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    30.10
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.38
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2763.005       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  2.06E-01  1.37E-01  9.75E-01  1.80E-01  9.30E-01  1.31E+00  8.06E-01  1.40E+00  6.78E-01  2.47E+00
 


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
+        1.22E+03
 
 TH 2
+        1.11E+01  1.03E+04
 
 TH 3
+        1.32E+01 -2.10E+03  2.69E+04
 
 TH 4
+        7.74E+00 -1.92E+02 -5.87E+02  4.78E+02
 
 TH 5
+        9.91E+00 -7.63E+03 -2.52E+04 -7.66E+02  4.64E+04
 
 TH 6
+        7.06E+00 -4.04E+00  7.55E+00  1.92E+00 -9.28E+00  1.84E+02
 
 TH 7
+        6.17E+00  5.97E+01 -5.09E+00 -5.94E+00 -5.47E+01 -6.44E+00  5.23E+01
 
 TH 8
+        8.64E-01  4.31E+01 -4.80E-01 -1.01E+01  5.04E+01  8.23E+00  3.31E+00 -4.33E+00
 
 TH 9
+        6.29E+00  2.59E+00  2.99E+02 -1.03E+01  1.29E+02 -3.24E+00 -1.26E+00 -1.19E+00  6.31E+01
 
 TH10
+        8.46E+00  2.49E+01  2.54E+01  1.18E+01  1.35E+02 -3.69E+00  1.85E+01  2.94E+01  9.39E+00  1.75E+02
 
 TH11
+       -1.83E+01 -2.69E+01 -9.51E+01  3.18E+00  9.38E+01  2.95E+00  6.95E+00  1.34E+01  7.33E+00  1.15E+00  1.72E+02
 
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
 #CPUT: Total CPU Time in Seconds,       44.612
Stop Time:
Sat Sep 18 01:50:50 CDT 2021
