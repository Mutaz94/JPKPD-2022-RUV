Sat Sep 18 15:19:43 CDT 2021
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
$DATA ../../../../data/spa/D/dat40.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m40.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   12246.6692946634        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.0771E+02  2.0595E+02 -2.4701E+01  1.6857E+02  8.3071E+01 -1.5572E+03 -7.2782E+02 -9.1512E+01 -1.0117E+03 -4.3062E+02
            -2.3973E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -645.798410107413        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.4370E+00  1.2076E+00  9.8100E-01  1.5625E+00  1.2106E+00  1.6661E+00  1.2664E+00  9.8393E-01  1.1861E+00  1.0075E+00
             1.4732E+01
 PARAMETER:  4.6256E-01  2.8867E-01  8.0813E-02  5.4628E-01  2.9114E-01  6.1046E-01  3.3615E-01  8.3798E-02  2.7065E-01  1.0749E-01
             2.7900E+00
 GRADIENT:   9.1685E+00  2.5516E+01 -2.9426E+00  4.3207E+01 -9.6869E+00  3.0464E+01 -4.2458E-01  3.5531E+00  7.0430E+00  2.7114E+00
             1.2486E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -665.884810789852        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.2958E+00  7.0351E-01  1.2890E+00  1.7212E+00  1.4889E+00  1.5684E+00  3.0701E+00  3.3883E-01  1.1927E+00  1.6857E+00
             1.3193E+01
 PARAMETER:  3.5913E-01 -2.5167E-01  3.5385E-01  6.4302E-01  4.9804E-01  5.5005E-01  1.2217E+00 -9.8226E-01  2.7625E-01  6.2217E-01
             2.6797E+00
 GRADIENT:  -1.8221E+01  1.6366E+01  2.3146E+00  1.5029E+01 -1.1153E+01  1.0688E+01  1.1199E+01  3.3734E-01  1.7113E+01  2.3780E+00
             1.1884E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -698.309006831334        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  1.0738E+00  2.5938E-01  9.2818E-01  1.6396E+00  4.1379E+00  1.5086E+00  3.8339E+00  1.6642E-01  6.8841E-01  7.0991E+00
             1.0347E+01
 PARAMETER:  1.7119E-01 -1.2495E+00  2.5469E-02  5.9447E-01  1.5202E+00  5.1120E-01  1.4439E+00 -1.6932E+00 -2.7336E-01  2.0600E+00
             2.4367E+00
 GRADIENT:  -4.2252E+01  1.0962E+01  5.5143E+00  7.3626E+01 -5.6698E+00  1.0904E+01  4.0513E-01  5.6049E-03 -4.6869E+00 -3.8310E+00
            -3.7157E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -713.542319262621        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      296
 NPARAMETR:  9.9596E-01  1.8451E-01  4.3360E-01  1.2441E+00  8.5643E+00  1.2507E+00  1.7971E+00  1.9795E-02  2.2134E-01  6.7468E+00
             1.0856E+01
 PARAMETER:  9.5956E-02 -1.5901E+00 -7.3563E-01  3.1844E-01  2.2476E+00  3.2371E-01  6.8619E-01 -3.8223E+00 -1.4080E+00  2.0091E+00
             2.4847E+00
 GRADIENT:  -1.1084E+00  4.6975E+00 -9.2151E+00  9.3965E+00 -5.8106E+00 -2.3417E+01  4.6206E-01 -5.1624E-05  1.2121E+00  1.2832E+01
            -9.6526E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -747.023839033071        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  7.1937E-01  3.4207E-02  8.6529E-02  6.6518E-01  5.8854E+01  1.4049E+00  3.0002E-01  1.0000E-02  1.0000E-02  3.1744E+00
             1.0417E+01
 PARAMETER: -2.2938E-01 -3.2753E+00 -2.3473E+00 -3.0770E-01  4.1751E+00  4.3997E-01 -1.1039E+00 -1.1055E+01 -4.6117E+00  1.2551E+00
             2.4435E+00
 GRADIENT:   3.7997E+01 -5.3269E-01 -3.2568E+01  6.9559E+01  8.1561E-02 -1.9157E+01  2.4141E-04  0.0000E+00  0.0000E+00 -9.5732E-04
            -4.4301E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -757.284372906075        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      484
 NPARAMETR:  4.8716E-01  1.8984E-02  2.9993E-02  3.1978E-01  2.0887E+02  1.4618E+00  5.2977E-02  1.0000E-02  1.0000E-02  2.0522E+00
             1.0508E+01
 PARAMETER: -6.1916E-01 -3.8642E+00 -3.4068E+00 -1.0401E+00  5.4417E+00  4.7970E-01 -2.8379E+00 -1.5765E+01 -7.2014E+00  8.1890E-01
             2.4521E+00
 GRADIENT:  -1.5581E+00  6.7776E-01  2.0036E+01 -2.5374E+01 -3.8782E-03 -3.0010E+00  2.5748E-04  0.0000E+00  0.0000E+00 -2.7737E-06
             7.1574E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -757.717275773051        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      661
 NPARAMETR:  4.3561E-01  1.1568E-02  2.2429E-02  2.5935E-01  3.2274E+02  1.4860E+00  3.1269E-02  1.0000E-02  1.0000E-02  1.7709E+00
             1.0398E+01
 PARAMETER: -7.3101E-01 -4.3595E+00 -3.6974E+00 -1.2496E+00  5.8768E+00  4.9606E-01 -3.3651E+00 -1.7384E+01 -8.1716E+00  6.7146E-01
             2.4416E+00
 GRADIENT:   1.2620E-01  1.3415E-01  2.2761E+00 -3.6591E+00  7.5899E-04  1.6519E-01  3.9695E-06  0.0000E+00  0.0000E+00 -1.6516E-06
             2.3835E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -757.739254651817        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      837            RESET HESSIAN, TYPE II
 NPARAMETR:  4.4372E-01  1.0000E-02  2.3599E-02  2.7007E-01  3.1687E+02  1.4871E+00  3.3013E-02  1.0000E-02  1.0000E-02  1.7615E+00
             1.0405E+01
 PARAMETER: -7.1257E-01 -4.5464E+00 -3.6465E+00 -1.2091E+00  5.8585E+00  4.9680E-01 -3.3108E+00 -1.7431E+01 -8.2252E+00  6.6614E-01
             2.4423E+00
 GRADIENT:   4.9938E+00  0.0000E+00  5.7384E+00  2.9664E+00  1.3092E-03  5.2191E-01  1.8008E-06  0.0000E+00  0.0000E+00 -1.8323E-06
             1.9356E+00

0ITERATION NO.:   42    OBJECTIVE VALUE:  -757.739254651817        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:      894
 NPARAMETR:  4.4372E-01  1.0000E-02  2.3599E-02  2.7007E-01  3.1687E+02  1.4871E+00  3.3013E-02  1.0000E-02  1.0000E-02  1.7615E+00
             1.0405E+01
 PARAMETER: -7.1257E-01 -4.5464E+00 -3.6465E+00 -1.2091E+00  5.8585E+00  4.9680E-01 -3.3108E+00 -1.7431E+01 -8.2252E+00  6.6614E-01
             2.4423E+00
 GRADIENT:   2.4317E-02  0.0000E+00 -1.3286E-01  1.5526E-01  1.3119E-03 -1.4429E-02  1.6818E-06  0.0000E+00  0.0000E+00 -1.8596E-06
            -1.8296E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      894
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7476E-03  7.8376E-06  1.1586E-04 -2.1325E-04  2.6945E-06
 SE:             2.8891E-02  3.9049E-06  3.2704E-04  4.0862E-04  6.7956E-06
 N:                     100         100         100         100         100

 P VAL.:         9.5177E-01  4.4737E-02  7.2313E-01  6.0175E-01  6.9173E-01

 ETASHRINKSD(%)  3.2110E+00  9.9987E+01  9.8904E+01  9.8631E+01  9.9977E+01
 ETASHRINKVR(%)  6.3189E+00  1.0000E+02  9.9988E+01  9.9981E+01  1.0000E+02
 EBVSHRINKSD(%)  3.4767E+00  9.9984E+01  9.8937E+01  9.8674E+01  9.9980E+01
 EBVSHRINKVR(%)  6.8325E+00  1.0000E+02  9.9989E+01  9.9982E+01  1.0000E+02
 RELATIVEINF(%)  1.2261E-01  3.2702E-08  2.1516E-05  1.8060E-05  1.2528E-08
 EPSSHRINKSD(%)  7.3396E+00
 EPSSHRINKVR(%)  1.4140E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -757.73925465181685     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -22.588428088078672     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.56
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.38
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -757.739       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.44E-01  1.00E-02  2.36E-02  2.70E-01  3.17E+02  1.49E+00  3.30E-02  1.00E-02  1.00E-02  1.76E+00  1.04E+01
 


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
+        2.35E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.08E+04  0.00E+00  1.27E+06
 
 TH 4
+       -6.72E+02  0.00E+00 -1.37E+05  1.63E+04
 
 TH 5
+        4.63E-04  0.00E+00 -8.74E-03  6.23E-04 -5.19E-09
 
 TH 6
+       -6.30E+00  0.00E+00  1.14E+03 -1.29E+02  1.34E-05  7.03E+01
 
 TH 7
+       -2.05E-01  0.00E+00 -6.66E-01 -4.37E-02 -1.84E-05 -7.10E-02  4.29E-01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+       -2.93E-03  0.00E+00 -5.66E-02  1.69E-02  2.04E-06 -1.27E-02 -3.60E-02  0.00E+00  0.00E+00  3.39E-03
 
 TH11
+       -2.25E+01  0.00E+00  3.39E+02 -2.33E+01 -6.27E-06  9.55E-01  3.92E-03  0.00E+00  0.00E+00 -2.27E-04  4.01E+00
 
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
 #CPUT: Total CPU Time in Seconds,       17.004
Stop Time:
Sat Sep 18 15:20:02 CDT 2021
