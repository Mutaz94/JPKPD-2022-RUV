Sat Sep 25 07:57:07 CDT 2021
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
$DATA ../../../../data/spa/A1/dat23.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m23.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1475.16674097378        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -8.6365E+00  1.8083E+01  2.6149E+01  1.8566E+01  1.0417E+01  1.8655E+01  1.7837E+00 -2.6509E+01  6.3308E+00  7.9822E+00
            -3.9401E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1565.25692979024        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0069E+00  1.0462E+00  9.5552E-01  9.6510E-01  9.6655E-01  9.1262E-01  9.4270E-01  1.1679E+00  8.4413E-01  7.4714E-01
             1.5790E+00
 PARAMETER:  1.0690E-01  1.4515E-01  5.4499E-02  6.4481E-02  6.5982E-02  8.5660E-03  4.0990E-02  2.5517E-01 -6.9446E-02 -1.9150E-01
             5.5680E-01
 GRADIENT:  -3.8520E+01  1.5215E+01  2.2258E+01 -8.1331E+00 -2.2363E+01 -1.4718E+01 -8.5165E+00 -1.2726E+01 -1.4768E+01  2.3410E+00
            -3.9044E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1569.43303612356        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0427E+00  1.1047E+00  8.3669E-01  9.2676E-01  9.6587E-01  9.4050E-01  9.8172E-01  1.4809E+00  8.9405E-01  3.9021E-01
             1.7122E+00
 PARAMETER:  1.4181E-01  1.9954E-01 -7.8306E-02  2.3936E-02  6.5276E-02  3.8658E-02  8.1553E-02  4.9265E-01 -1.1990E-02 -8.4108E-01
             6.3779E-01
 GRADIENT:   4.7849E+01 -1.6712E+01 -1.8532E+01  1.4613E+01  4.7595E+01 -1.7693E+00 -2.1981E+00 -1.4506E-01  1.3187E+00 -1.9684E+00
            -1.2603E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1571.53537432171        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0149E+00  1.1803E+00  3.6825E-01  8.1368E-01  6.7349E-01  9.4879E-01  9.8032E-01  7.8928E-01  8.5329E-01  3.0714E-01
             1.6940E+00
 PARAMETER:  1.1480E-01  2.6578E-01 -8.9900E-01 -1.0619E-01 -2.9528E-01  4.7427E-02  8.0129E-02 -1.3664E-01 -5.8651E-02 -1.0805E+00
             6.2710E-01
 GRADIENT:  -2.3676E+01  3.4805E+01  2.0593E+01  4.2253E+00 -4.3831E+01  1.7603E+00  4.5516E+00 -1.1687E-01  6.1005E+00  2.6058E+00
             6.7258E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1572.21721594704        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  1.0054E+00  1.1146E+00  2.4438E-01  7.8192E-01  5.5528E-01  9.4249E-01  9.3384E-01  5.6735E-01  7.8056E-01  2.1527E-01
             1.6888E+00
 PARAMETER:  1.0534E-01  2.0852E-01 -1.3090E+00 -1.4600E-01 -4.8828E-01  4.0768E-02  3.1549E-02 -4.6677E-01 -1.4774E-01 -1.4359E+00
             6.2403E-01
 GRADIENT:  -3.1914E+01  2.1433E+01  1.7071E+01 -5.4313E+00 -4.7052E+01  3.4595E-01 -2.6767E+00 -5.0010E-01  6.1999E-01  2.3375E+00
             9.2576E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1595.76715278450        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      377
 NPARAMETR:  9.9980E-01  1.6547E+00  1.4899E-01  4.7802E-01  7.9316E-01  9.4006E-01  7.5437E-01  1.4821E+00  1.1982E+00  1.0000E-02
             1.6115E+00
 PARAMETER:  9.9804E-02  6.0364E-01 -1.8038E+00 -6.3811E-01 -1.3174E-01  3.8190E-02 -1.8187E-01  4.9348E-01  2.8082E-01 -4.6103E+00
             5.7715E-01
 GRADIENT:  -2.8352E+01  1.1242E+01 -4.5477E+00  4.3565E+01  1.1581E+01  5.1729E+00  8.6502E+00 -1.7634E+01  1.4755E+01  0.0000E+00
             4.3370E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1601.90122108273        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      455
 NPARAMETR:  1.0005E+00  1.7005E+00  1.5281E-01  4.5942E-01  8.3331E-01  9.3637E-01  7.4772E-01  1.8608E+00  1.2153E+00  1.0000E-02
             1.5416E+00
 PARAMETER:  1.0053E-01  6.3092E-01 -1.7786E+00 -6.7780E-01 -8.2349E-02  3.4257E-02 -1.9073E-01  7.2101E-01  2.9499E-01 -4.7794E+00
             5.3285E-01
 GRADIENT:  -2.6743E+01  1.3172E+01 -1.1563E+01  4.6322E+01  2.0587E+01  4.0608E+00  8.6326E+00 -8.6681E+00  9.4805E+00  0.0000E+00
             2.9908E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1602.00401713492        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:      622             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0006E+00  1.7007E+00  1.5291E-01  4.5924E-01  8.3272E-01  9.3550E-01  7.4471E-01  1.8637E+00  1.2153E+00  1.0000E-02
             1.5409E+00
 PARAMETER:  1.0059E-01  6.3106E-01 -1.7779E+00 -6.7819E-01 -8.3059E-02  3.3329E-02 -1.9477E-01  7.2255E-01  2.9498E-01 -4.7789E+00
             5.3238E-01
 GRADIENT:  -2.6547E+01  1.4723E+01 -1.1054E+01  4.5431E+01  1.8024E+01  3.7318E+00  7.5679E+00 -8.5014E+00  9.4133E+00  0.0000E+00
             2.9763E+01

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1602.01094277358        NO. OF FUNC. EVALS.: 103
 CUMULATIVE NO. OF FUNC. EVALS.:      725
 NPARAMETR:  1.0007E+00  1.6997E+00  1.5264E-01  4.5955E-01  8.3251E-01  9.3523E-01  7.4444E-01  1.8623E+00  1.2156E+00  1.0000E-02
             1.5417E+00
 PARAMETER:  1.0059E-01  6.3106E-01 -1.7779E+00 -6.7819E-01 -8.3209E-02  3.3132E-02 -1.9532E-01  7.2255E-01  2.9498E-01 -4.7789E+00
             5.3238E-01
 GRADIENT:  -4.7964E+04  1.5260E+04  5.4173E+03 -1.4179E+04  4.8218E+04  4.8205E+04 -4.9351E+04  1.3312E+04 -3.2675E+04  0.0000E+00
            -1.8077E+04

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      725
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.9219E-02 -6.6283E-03  1.1507E-03 -7.8182E-03 -7.0963E-04
 SE:             2.9423E-02  2.6555E-02  1.6556E-02  2.1174E-02  3.0202E-04
 N:                     100         100         100         100         100

 P VAL.:         5.1362E-01  8.0289E-01  9.4459E-01  7.1195E-01  1.8791E-02

 ETASHRINKSD(%)  1.4306E+00  1.1038E+01  4.4535E+01  2.9064E+01  9.8988E+01
 ETASHRINKVR(%)  2.8406E+00  2.0858E+01  6.9237E+01  4.9681E+01  9.9990E+01
 EBVSHRINKSD(%)  1.0363E+00  9.3238E+00  4.8051E+01  2.5492E+01  9.8994E+01
 EBVSHRINKVR(%)  2.0619E+00  1.7778E+01  7.3013E+01  4.4486E+01  9.9990E+01
 RELATIVEINF(%)  9.6531E+01  1.1589E+01  9.0268E+00  7.5923E+00  1.3135E-03
 EPSSHRINKSD(%)  4.3343E+01
 EPSSHRINKVR(%)  6.7900E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1602.0109427735777     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -866.86011620983948     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.38
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.26
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1602.011       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.70E+00  1.53E-01  4.59E-01  8.33E-01  9.35E-01  7.44E-01  1.86E+00  1.22E+00  1.00E-02  1.54E+00
 


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
+        2.38E+08
 
 TH 2
+        3.97E+02  2.09E+06
 
 TH 3
+        1.58E+03  2.21E+04  3.27E+07
 
 TH 4
+       -1.48E+03  7.21E+06 -7.66E+04  2.49E+07
 
 TH 5
+        2.29E+04 -2.74E+03 -1.00E+04  8.29E+03  3.48E+08
 
 TH 6
+       -4.79E+03  4.49E+02  1.78E+03 -1.57E+03  5.83E+03  2.76E+08
 
 TH 7
+       -3.52E+03  3.48E+02  1.26E+03 -1.14E+03 -1.99E+08 -3.35E+03  1.14E+08
 
 TH 8
+        3.33E+02 -1.67E+06  1.72E+04 -2.78E+03 -1.76E+03  3.62E+02  2.65E+02  1.32E+06
 
 TH 9
+       -1.24E+03 -5.05E+02 -1.90E+03  1.70E+03  6.41E+03 -1.36E+03 -9.77E+02 -3.90E+02  1.88E+07
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -5.55E+02  7.12E+02  2.93E+03 -2.51E+03  2.78E+03 -5.92E+02 -4.23E+02  5.81E+02 -2.15E+03  0.00E+00  3.58E+06
 
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
 #CPUT: Total CPU Time in Seconds,       13.706
Stop Time:
Sat Sep 25 07:57:22 CDT 2021
