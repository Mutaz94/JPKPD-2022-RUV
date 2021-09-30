Wed Sep 29 14:34:52 CDT 2021
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
$DATA ../../../../data/spa/S1/dat79.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m79.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1665.72270393682        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0400E+02 -1.0014E+01  2.3840E-01  2.2182E+01 -1.5631E+01  5.9482E+01 -9.7782E+00  5.1297E-02  2.1708E+01  1.0332E+01
            -6.3290E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1672.47401898860        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      192
 NPARAMETR:  1.0135E+00  1.1749E+00  1.0813E+00  9.0954E-01  1.1146E+00  8.8431E-01  1.1444E+00  1.0070E+00  8.7133E-01  9.1696E-01
             1.3185E+00
 PARAMETER:  1.1342E-01  2.6119E-01  1.7819E-01  5.1830E-03  2.0848E-01 -2.2951E-02  2.3489E-01  1.0702E-01 -3.7739E-02  1.3309E-02
             3.7651E-01
 GRADIENT:  -1.2441E+01 -1.4676E+01  2.0192E+01 -4.0525E+01 -2.2297E+01 -2.1207E+01 -3.9629E-01 -6.4260E+00 -2.0677E+00  9.5211E-01
             4.7650E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1674.33906761122        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  1.0280E+00  1.3524E+00  1.0128E+00  8.2417E-01  1.1750E+00  9.0339E-01  1.0834E+00  1.4984E+00  8.7866E-01  8.2510E-01
             1.3040E+00
 PARAMETER:  1.2760E-01  4.0186E-01  1.1272E-01 -9.3379E-02  2.6124E-01 -1.6054E-03  1.8011E-01  5.0441E-01 -2.9357E-02 -9.2251E-02
             3.6546E-01
 GRADIENT:   2.2463E+01  8.8509E+00  1.1708E+00  2.7865E+00  6.0869E+00 -1.3240E+01  5.5453E+00 -2.2611E-01 -3.1277E+00 -6.7151E+00
             4.1860E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1678.19003629319        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      544
 NPARAMETR:  1.0190E+00  1.3077E+00  1.0028E+00  8.4963E-01  1.1543E+00  9.3591E-01  1.0612E+00  1.3379E+00  9.3123E-01  9.2858E-01
             1.1503E+00
 PARAMETER:  1.1881E-01  3.6824E-01  1.0280E-01 -6.2953E-02  2.4353E-01  3.3769E-02  1.5936E-01  3.9112E-01  2.8756E-02  2.5906E-02
             2.3998E-01
 GRADIENT:   2.5527E+00  4.0734E+00  1.4371E+00  3.6281E+00 -1.7515E+00 -9.8795E-03  1.2275E-01 -1.9922E-01 -2.6490E-01 -4.9370E-01
            -7.5109E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1678.38273913155        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      721
 NPARAMETR:  1.0199E+00  1.5383E+00  6.4787E-01  6.9748E-01  1.0997E+00  9.3813E-01  9.4360E-01  9.6189E-01  1.0408E+00  8.6353E-01
             1.1528E+00
 PARAMETER:  1.1970E-01  5.3066E-01 -3.3407E-01 -2.6029E-01  1.9506E-01  3.6138E-02  4.1946E-02  6.1147E-02  1.3998E-01 -4.6721E-02
             2.4220E-01
 GRADIENT:  -1.4522E+00  1.6699E+01 -4.1940E-01  1.3528E+01 -4.9437E+00 -1.5303E-01 -1.1268E-01  2.9826E-01 -1.1484E+00  5.2840E-01
             1.1514E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1678.87177087509        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      900
 NPARAMETR:  1.0200E+00  1.9066E+00  3.9854E-01  4.6007E-01  1.2098E+00  9.3049E-01  8.0851E-01  5.8633E-01  1.3584E+00  9.6010E-01
             1.1721E+00
 PARAMETER:  1.1985E-01  7.4532E-01 -8.1995E-01 -6.7638E-01  2.9050E-01  2.7955E-02 -1.1256E-01 -4.3388E-01  4.0634E-01  5.9284E-02
             2.5884E-01
 GRADIENT:  -3.5930E+00  2.6828E+01 -1.0313E+00  1.3514E+01 -4.4865E+00 -3.5213E+00  1.5150E+00  8.8190E-02 -1.6987E+00  3.9713E+00
             7.4710E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1679.15794096278        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1076
 NPARAMETR:  1.0210E+00  2.0605E+00  3.0839E-01  3.5881E-01  1.2748E+00  9.3721E-01  7.5437E-01  4.2868E-01  1.6179E+00  9.7539E-01
             1.1589E+00
 PARAMETER:  1.2080E-01  8.2296E-01 -1.0764E+00 -9.2497E-01  3.4280E-01  3.5154E-02 -1.8187E-01 -7.4705E-01  5.8112E-01  7.5086E-02
             2.4746E-01
 GRADIENT:  -4.6983E-01  2.8377E+01 -4.4586E-01  1.1791E+01 -3.8117E+00 -5.5420E-01 -2.3356E+00  9.7057E-02  8.7196E-01  1.0136E-01
             1.7095E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1679.55324297825        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1261             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0217E+00  2.0627E+00  2.9212E-01  3.3424E-01  1.2899E+00  9.3826E-01  7.5450E-01  1.0064E-01  1.6637E+00  9.8684E-01
             1.1582E+00
 PARAMETER:  1.2148E-01  8.2400E-01 -1.1306E+00 -9.9589E-01  3.5456E-01  3.6268E-02 -1.8170E-01 -2.1962E+00  6.0905E-01  8.6757E-02
             2.4690E-01
 GRADIENT:   4.0840E+02  8.6861E+02  2.9373E+00  7.0076E+01  1.3132E+01  2.6479E+01  1.2573E+01  1.0750E-02  1.2249E+01  8.3681E-01
             2.3682E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1679.55787199934        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1444             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0217E+00  2.0640E+00  2.9272E-01  3.3463E-01  1.2895E+00  9.3824E-01  7.5450E-01  1.5781E-02  1.6660E+00  9.8482E-01
             1.1575E+00
 PARAMETER:  1.2143E-01  8.2462E-01 -1.1285E+00 -9.9473E-01  3.5428E-01  3.6253E-02 -1.8170E-01 -4.0489E+00  6.1044E-01  8.4706E-02
             2.4622E-01
 GRADIENT:   4.0867E+02  8.7337E+02  3.1939E+00  7.0586E+01  1.2610E+01  2.6495E+01  1.2559E+01  3.6030E-04  1.2399E+01  5.3510E-01
             1.9393E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1679.55944954614        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1627
 NPARAMETR:  1.0217E+00  2.0633E+00  2.9261E-01  3.3511E-01  1.2899E+00  9.3824E-01  7.5460E-01  1.0000E-02  1.6661E+00  9.8409E-01
             1.1574E+00
 PARAMETER:  1.2143E-01  8.2431E-01 -1.1289E+00 -9.9330E-01  3.5455E-01  3.6252E-02 -1.8156E-01 -5.0790E+00  6.1051E-01  8.3964E-02
             2.4620E-01
 GRADIENT:   1.9509E+00 -1.5586E+01  4.6386E-02 -4.4913E-01 -2.7181E-01  6.0534E-02  1.0862E-02  0.0000E+00  3.9397E-03  1.1387E-02
            -2.3737E-02

0ITERATION NO.:   46    OBJECTIVE VALUE:  -1679.55944954614        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1649
 NPARAMETR:  1.0217E+00  2.0633E+00  2.9261E-01  3.3511E-01  1.2899E+00  9.3824E-01  7.5460E-01  1.0000E-02  1.6661E+00  9.8409E-01
             1.1574E+00
 PARAMETER:  1.2143E-01  8.2431E-01 -1.1289E+00 -9.9330E-01  3.5455E-01  3.6252E-02 -1.8156E-01 -5.0790E+00  6.1051E-01  8.3964E-02
             2.4620E-01
 GRADIENT:   1.9509E+00 -1.5586E+01  4.6386E-02 -4.4913E-01 -2.7181E-01  6.0534E-02  1.0862E-02  0.0000E+00  3.9397E-03  1.1387E-02
            -2.3737E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1649
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2170E-04 -2.8050E-02 -1.9981E-04  3.6509E-02 -4.2968E-02
 SE:             2.9774E-02  2.6416E-02  6.7716E-05  2.0192E-02  2.0615E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9674E-01  2.8830E-01  3.1705E-03  7.0597E-02  3.7128E-02

 ETASHRINKSD(%)  2.5330E-01  1.1503E+01  9.9773E+01  3.2354E+01  3.0938E+01
 ETASHRINKVR(%)  5.0595E-01  2.1683E+01  9.9999E+01  5.4240E+01  5.2305E+01
 EBVSHRINKSD(%)  6.5335E-01  1.1801E+01  9.9804E+01  3.5832E+01  2.8823E+01
 EBVSHRINKVR(%)  1.3024E+00  2.2209E+01  1.0000E+02  5.8824E+01  4.9338E+01
 RELATIVEINF(%)  9.8634E+01  7.4813E+00  3.2692E-05  2.6474E+00  1.4074E+01
 EPSSHRINKSD(%)  4.3173E+01
 EPSSHRINKVR(%)  6.7707E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1679.5594495461389     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -944.40862298240074     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.71
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     6.38
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1679.559       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  2.06E+00  2.93E-01  3.35E-01  1.29E+00  9.38E-01  7.55E-01  1.00E-02  1.67E+00  9.84E-01  1.16E+00
 


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
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.29E+03
 
 TH 2
+        5.69E+01  3.78E+02
 
 TH 3
+        6.72E+01  1.54E+02  5.73E+02
 
 TH 4
+       -5.12E+01  3.48E+02 -6.51E+02  1.58E+03
 
 TH 5
+       -6.82E+01 -1.04E+02 -3.22E+02  3.65E+02  3.36E+02
 
 TH 6
+        3.07E+01 -2.54E+01 -2.97E+01 -1.14E+01 -6.31E+01  5.93E+01
 
 TH 7
+        9.25E+01  4.06E+00 -1.89E+01 -3.25E+01 -6.13E+00  1.68E+01  1.90E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.30E+00 -2.74E+00 -5.26E+01  7.31E+01  5.40E+00  1.51E+01  6.64E+00  0.00E+00  8.91E+00
 
 TH10
+        3.90E+00 -2.38E+01 -2.84E+01 -4.59E+00 -5.67E+01  3.94E+01  8.76E+00  0.00E+00  1.33E+01  3.75E+01
 
 TH11
+       -6.42E+01 -2.21E+01 -3.86E+01  1.73E+01 -2.57E+01 -1.09E+01  2.55E+01  0.00E+00  9.19E+00  2.97E+01  1.31E+02
 
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
+        1.19E+03
 
 TH 2
+       -9.38E+00  3.69E+02
 
 TH 3
+        8.13E+00  1.32E+02  5.83E+02
 
 TH 4
+       -2.61E+01  3.21E+02 -6.46E+02  1.55E+03
 
 TH 5
+       -6.28E+00 -1.08E+02 -3.07E+02  3.40E+02  3.59E+02
 
 TH 6
+        2.67E-01 -1.98E+00  8.59E-01 -8.46E+00 -2.00E+00  2.21E+02
 
 TH 7
+        1.69E+00  7.60E+00 -2.63E+01 -2.33E+01 -5.78E-01 -5.24E-01  2.26E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.17E+00 -1.13E+01 -4.54E+01  8.28E+01 -2.97E+00 -3.72E-01  8.76E+00  0.00E+00  2.39E+01
 
 TH10
+       -1.85E-01 -1.20E+01 -2.85E+01  7.82E+00 -5.84E+01  4.31E-01  8.68E+00  0.00E+00  5.38E+00  6.55E+01
 
 TH11
+       -9.42E+00 -1.65E+01 -2.61E+01  1.58E+01 -7.77E+00  3.12E+00  1.27E+01  0.00E+00  4.06E+00  1.70E+01  1.61E+02
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.19E+03
 
 TH 2
+       -9.42E+01  3.61E+02
 
 TH 3
+       -3.41E+01  1.37E+02  5.29E+02
 
 TH 4
+        1.40E+00  2.94E+02 -6.30E+02  1.53E+03
 
 TH 5
+        4.38E+01 -1.17E+02 -3.14E+02  3.45E+02  3.65E+02
 
 TH 6
+       -5.60E+01  2.64E+01 -1.51E+01  7.24E+01 -9.52E+00  1.50E+02
 
 TH 7
+       -1.09E+02  1.13E+01 -2.17E+01 -2.02E+01  9.24E+00  2.03E+01  2.77E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.96E+01 -9.22E+00 -4.11E+01  8.90E+01 -9.05E-01  1.03E+01  1.22E+01  0.00E+00  2.55E+01
 
 TH10
+       -1.79E+01 -1.72E+01 -2.67E+01  1.24E-01 -5.49E+01  1.56E+01 -1.40E+00  0.00E+00  3.38E+00  6.35E+01
 
 TH11
+        6.56E+01 -1.40E+01 -2.07E+01  2.47E+01  1.17E+01 -1.27E+01 -7.73E+00  0.00E+00  3.72E+00  1.10E+01  2.05E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,       29.137
Stop Time:
Wed Sep 29 14:35:23 CDT 2021
