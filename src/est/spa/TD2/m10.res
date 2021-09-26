Sat Sep 25 13:19:41 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat10.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m10.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1651.93421947408        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.5368E+01 -2.1890E+01 -1.1852E+01 -2.4044E+01 -8.1058E+00 -2.6356E+01  7.3255E+00  6.0813E+00  1.0845E+01  2.3386E+01
            -1.6955E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1655.01782517301        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0088E+00  1.0305E+00  1.0863E+00  9.9086E-01  1.0435E+00  1.1069E+00  9.4307E-01  9.6496E-01  9.3375E-01  8.4183E-01
             1.0841E+00
 PARAMETER:  1.0874E-01  1.3007E-01  1.8274E-01  9.0815E-02  1.4256E-01  2.0154E-01  4.1387E-02  6.4330E-02  3.1450E-02 -7.2183E-02
             1.8079E-01
 GRADIENT:   8.1396E+01 -6.3696E+00  1.3680E+01 -2.5705E+01  3.8087E+00  1.9656E+01 -6.1972E-01 -5.2632E+00 -8.7603E+00 -4.9012E+00
             8.4755E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1656.54552850662        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      161
 NPARAMETR:  9.8971E-01  1.0885E+00  1.0591E+00  9.6186E-01  1.0518E+00  1.0878E+00  8.1186E-01  1.1527E+00  1.0126E+00  8.2203E-01
             1.0613E+00
 PARAMETER:  8.9653E-02  1.8484E-01  1.5745E-01  6.1118E-02  1.5050E-01  1.8414E-01 -1.0843E-01  2.4207E-01  1.1254E-01 -9.5975E-02
             1.5947E-01
 GRADIENT:   4.5161E+01  2.4027E+00  1.0331E+00 -5.0692E+00  3.8579E+00  1.3731E+01 -4.5069E-01 -1.6404E-01 -4.9747E+00 -5.7673E+00
             2.1743E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1657.50765638757        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      277
 NPARAMETR:  9.7838E-01  1.0718E+00  1.2638E+00  9.8515E-01  1.1304E+00  1.0689E+00  5.7860E-01  1.4101E+00  1.1251E+00  9.5554E-01
             1.0477E+00
 PARAMETER:  7.8143E-02  1.6932E-01  3.3414E-01  8.5038E-02  2.2261E-01  1.6665E-01 -4.4714E-01  4.4367E-01  2.1784E-01  5.4525E-02
             1.4657E-01
 GRADIENT:  -1.2135E+01  3.2070E+00 -1.4987E+00  9.8627E+00 -2.2582E-01 -2.3734E+00  2.6518E+00  5.3511E-01  6.4682E+00  1.4537E+00
            -4.2807E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1658.60290376160        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      452
 NPARAMETR:  9.8361E-01  1.2564E+00  1.3445E+00  8.4418E-01  1.2334E+00  1.0722E+00  2.0612E-01  1.7905E+00  1.3419E+00  1.0333E+00
             1.0492E+00
 PARAMETER:  8.3478E-02  3.2827E-01  3.9601E-01 -6.9386E-02  3.0974E-01  1.6969E-01 -1.4793E+00  6.8249E-01  3.9408E-01  1.3272E-01
             1.4807E-01
 GRADIENT:  -3.0240E+00 -5.0842E-01  9.3491E-01 -3.0757E+00 -1.3579E+00 -7.4706E-01  5.0256E-01  4.4338E-01 -1.6489E-01  8.4636E-01
             4.9812E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1658.91697136732        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      628
 NPARAMETR:  9.8363E-01  1.1813E+00  1.3599E+00  8.9752E-01  1.2142E+00  1.0724E+00  5.2508E-02  1.7046E+00  1.2876E+00  1.0340E+00
             1.0477E+00
 PARAMETER:  8.3494E-02  2.6665E-01  4.0744E-01 -8.1185E-03  2.9407E-01  1.6991E-01 -2.8468E+00  6.3332E-01  3.5281E-01  1.3344E-01
             1.4658E-01
 GRADIENT:  -2.2792E+00  1.8850E-01 -2.1611E+00  2.4837E+00  3.0069E+00 -5.9252E-01  3.8531E-02  3.4658E-01  4.8463E-01  1.1028E+00
             3.7693E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1658.95842712846        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      804
 NPARAMETR:  9.8482E-01  1.1536E+00  1.3898E+00  9.1533E-01  1.2067E+00  1.0741E+00  1.9831E-02  1.6983E+00  1.2585E+00  1.0219E+00
             1.0478E+00
 PARAMETER:  8.4706E-02  2.4291E-01  4.2914E-01  1.1529E-02  2.8787E-01  1.7149E-01 -3.8205E+00  6.2964E-01  3.2995E-01  1.2169E-01
             1.4673E-01
 GRADIENT:   4.1325E-01  3.2765E-01 -1.7324E-01  3.4817E-01  2.4571E-01  1.0512E-01  5.3172E-03 -3.6800E-03 -2.8639E-01  8.9077E-02
             1.0919E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1658.96104486666        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      979
 NPARAMETR:  9.8463E-01  1.1591E+00  1.3894E+00  9.1124E-01  1.2083E+00  1.0738E+00  1.0000E-02  1.7058E+00  1.2651E+00  1.0221E+00
             1.0476E+00
 PARAMETER:  8.4511E-02  2.4766E-01  4.2885E-01  7.0540E-03  2.8923E-01  1.7124E-01 -4.6899E+00  6.3405E-01  3.3514E-01  1.2185E-01
             1.4645E-01
 GRADIENT:  -1.9610E-03 -8.1313E-03 -6.4167E-03 -1.1266E-03  1.2041E-02 -6.2163E-04  0.0000E+00 -5.8560E-07 -1.2365E-03  3.5116E-04
             1.1185E-03

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1658.96104486666        NO. OF FUNC. EVALS.:  28
 CUMULATIVE NO. OF FUNC. EVALS.:     1007
 NPARAMETR:  9.8473E-01  1.1592E+00  1.3898E+00  9.1125E-01  1.2082E+00  1.0740E+00  1.0000E-02  1.7058E+00  1.2651E+00  1.0221E+00
             1.0475E+00
 PARAMETER:  8.4511E-02  2.4766E-01  4.2885E-01  7.0540E-03  2.8923E-01  1.7124E-01 -4.6899E+00  6.3405E-01  3.3514E-01  1.2185E-01
             1.4645E-01
 GRADIENT:  -1.7031E-02 -5.0592E-03 -5.8775E-03 -6.1725E-04  1.2827E-02 -6.1154E-03  0.0000E+00 -1.0515E-04 -1.1802E-03 -3.8432E-04
             1.0693E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1007
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0662E-03 -9.8650E-04 -3.3122E-02 -1.6721E-03 -3.6673E-02
 SE:             2.9893E-02  2.7379E-04  1.5865E-02  2.9180E-02  2.1489E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7155E-01  3.1451E-04  3.6815E-02  9.5430E-01  8.7892E-02

 ETASHRINKSD(%)  1.0000E-10  9.9083E+01  4.6852E+01  2.2447E+00  2.8010E+01
 ETASHRINKVR(%)  1.0000E-10  9.9992E+01  7.1753E+01  4.4390E+00  4.8174E+01
 EBVSHRINKSD(%)  4.0708E-01  9.9203E+01  4.9538E+01  2.9643E+00  2.6206E+01
 EBVSHRINKVR(%)  8.1251E-01  9.9994E+01  7.4536E+01  5.8407E+00  4.5545E+01
 RELATIVEINF(%)  9.9068E+01  7.5218E-04  9.6897E+00  1.3107E+01  2.1397E+01
 EPSSHRINKSD(%)  4.3799E+01
 EPSSHRINKVR(%)  6.8415E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1658.9610448666583     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -923.81021830292013     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.85
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.04
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1658.961       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.85E-01  1.16E+00  1.39E+00  9.11E-01  1.21E+00  1.07E+00  1.00E-02  1.71E+00  1.27E+00  1.02E+00  1.05E+00
 


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
+        9.87E+02
 
 TH 2
+       -1.43E+01  5.87E+02
 
 TH 3
+        2.14E+00  6.06E+01  5.66E+01
 
 TH 4
+       -1.02E+01  5.52E+02 -4.57E+00  7.88E+02
 
 TH 5
+       -3.41E-02 -1.78E+02 -9.70E+01 -3.72E+01  4.11E+02
 
 TH 6
+       -4.32E-01 -9.71E-01  5.26E-01 -2.52E-01  1.03E-01  1.70E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -6.08E-02 -1.98E+01 -1.60E+01 -7.59E-01 -1.56E+00  2.05E-01  0.00E+00  1.38E+01
 
 TH 9
+        2.06E+00 -1.04E+02  1.40E+00  7.64E+00  5.04E+00 -7.49E-01  0.00E+00 -8.14E-02  1.10E+02
 
 TH10
+        1.63E+00 -4.91E+00 -9.07E-01 -1.76E+00 -5.86E+01  6.55E-01  0.00E+00  6.76E+00  2.60E+00  7.12E+01
 
 TH11
+       -5.21E+00 -2.42E+01 -5.87E+00 -1.44E+01 -6.78E+00  4.37E+00  0.00E+00  3.80E+00  9.09E+00  1.43E+01  1.95E+02
 
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
 #CPUT: Total CPU Time in Seconds,       17.952
Stop Time:
Sat Sep 25 13:20:01 CDT 2021
