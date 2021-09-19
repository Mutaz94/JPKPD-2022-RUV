Sat Sep 18 12:00:50 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat99.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m99.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1691.22077970142        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.0130E+01 -7.5459E+01 -3.5941E+01 -4.1514E+01  5.4285E+01  3.3731E+01 -1.6660E+01 -2.1742E+00  6.5453E+00  1.5494E+01
             2.0068E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1699.76565496527        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.9936E-01  1.1676E+00  1.0281E+00  9.3106E-01  1.0217E+00  8.6216E-01  1.2141E+00  1.1106E+00  9.3937E-01  8.0292E-01
             9.4003E-01
 PARAMETER:  9.9359E-02  2.5495E-01  1.2770E-01  2.8565E-02  1.2152E-01 -4.8312E-02  2.9401E-01  2.0490E-01  3.7454E-02 -1.1950E-01
             3.8153E-02
 GRADIENT:   2.9289E+01  1.7983E+01  1.6015E+01 -5.9699E+00  2.6786E+00 -2.3914E+01  6.5145E+00 -1.0355E+01 -3.2975E+00 -2.8402E+00
            -1.1114E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1701.15694074474        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.9802E-01  1.0448E+00  1.1891E+00  1.0056E+00  1.0259E+00  8.8107E-01  1.2577E+00  1.5459E+00  8.7751E-01  7.3828E-01
             9.3846E-01
 PARAMETER:  9.8019E-02  1.4382E-01  2.7324E-01  1.0556E-01  1.2561E-01 -2.6619E-02  3.2927E-01  5.3563E-01 -3.0672E-02 -2.0343E-01
             3.6483E-02
 GRADIENT:   2.9442E+01 -8.6211E+00 -9.1789E-01 -9.5450E+00  1.4565E+01 -1.4015E+01  1.2186E+00 -8.2148E-02 -3.5844E+00 -5.9744E+00
            -1.1274E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1702.05813686980        NO. OF FUNC. EVALS.: 125
 CUMULATIVE NO. OF FUNC. EVALS.:      281
 NPARAMETR:  9.9833E-01  1.0847E+00  1.2174E+00  9.9495E-01  1.0370E+00  9.2656E-01  1.2003E+00  1.6327E+00  9.2411E-01  8.0683E-01
             9.6673E-01
 PARAMETER:  9.8333E-02  1.8135E-01  2.9670E-01  9.4937E-02  1.3635E-01  2.3727E-02  2.8258E-01  5.9025E-01  2.1078E-02 -1.1464E-01
             6.6162E-02
 GRADIENT:   2.8629E+01  6.6131E+00  1.9968E+00  3.3912E+00 -1.3734E+01  6.3524E+00  2.1038E+00  3.1907E+00  1.2269E+00  2.6971E+00
             3.9295E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1702.22980444878        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      437
 NPARAMETR:  1.0054E+00  1.0847E+00  1.2324E+00  9.9753E-01  1.0494E+00  9.2302E-01  1.1946E+00  1.6327E+00  9.2071E-01  7.9966E-01
             9.5949E-01
 PARAMETER:  1.0540E-01  1.8135E-01  3.0898E-01  9.7527E-02  1.4825E-01  1.9895E-02  2.7780E-01  5.9025E-01  1.7385E-02 -1.2357E-01
             5.8645E-02
 GRADIENT:   1.3103E+00 -1.6207E+00 -9.4936E-01  4.7990E-01  1.1703E-01  2.5604E-02 -2.2995E-01  1.3972E+00 -4.2123E-01 -4.3812E-01
            -2.5453E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1702.24331390676        NO. OF FUNC. EVALS.: 155
 CUMULATIVE NO. OF FUNC. EVALS.:      592
 NPARAMETR:  1.0047E+00  1.0861E+00  1.2502E+00  9.9821E-01  1.0555E+00  9.2284E-01  1.1918E+00  1.6218E+00  9.2538E-01  8.1163E-01
             9.5999E-01
 PARAMETER:  1.0472E-01  1.8263E-01  3.2327E-01  9.8213E-02  1.5403E-01  1.9705E-02  2.7547E-01  5.8355E-01  2.2453E-02 -1.0870E-01
             5.9167E-02
 GRADIENT:   4.8103E+01  7.8165E+00  1.2098E+00  8.4029E+00  6.7590E-01  4.8586E+00  1.3651E+00  2.5708E-01  4.7711E-01 -1.2414E-01
            -1.1364E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1702.24705045431        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      665
 NPARAMETR:  1.0035E+00  1.0905E+00  1.2368E+00  9.9496E-01  1.0543E+00  9.2180E-01  1.1890E+00  1.5930E+00  9.2745E-01  8.1368E-01
             9.5998E-01
 PARAMETER:  1.0346E-01  1.8664E-01  3.1257E-01  9.4948E-02  1.5285E-01  1.8568E-02  2.7315E-01  5.6561E-01  2.4679E-02 -1.0619E-01
             5.9155E-02
 GRADIENT:   4.4109E+01  8.0827E+00  1.3894E+00  8.3477E+00  1.4836E+00  4.4072E+00  1.2882E+00 -3.7927E-01  4.4814E-01 -3.2724E-02
            -8.4392E-03

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1702.24965069402        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:      859             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0038E+00  1.0906E+00  1.2361E+00  9.9448E-01  1.0541E+00  9.2255E-01  1.1892E+00  1.5934E+00  9.2779E-01  8.1392E-01
             9.6006E-01
 PARAMETER:  1.0382E-01  1.8675E-01  3.1195E-01  9.4464E-02  1.5271E-01  1.9389E-02  2.7330E-01  5.6590E-01  2.5046E-02 -1.0590E-01
             5.9236E-02
 GRADIENT:   4.5197E+01  7.7171E+00  1.3630E+00  7.8269E+00  1.3409E+00  4.7274E+00  1.3469E+00 -3.0687E-01  5.1136E-01  4.9056E-02
             5.5510E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1702.24972559081        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1041
 NPARAMETR:  1.0039E+00  1.0905E+00  1.2361E+00  9.9426E-01  1.0540E+00  9.2296E-01  1.1891E+00  1.5930E+00  9.2842E-01  8.1388E-01
             9.6001E-01
 PARAMETER:  1.0387E-01  1.8666E-01  3.1195E-01  9.4241E-02  1.5256E-01  1.9833E-02  2.7317E-01  5.6562E-01  2.5732E-02 -1.0595E-01
             5.9186E-02
 GRADIENT:  -2.6594E+00 -4.6614E-01  1.0475E+00 -1.0706E-01  2.9265E-02  1.2652E-03  3.1982E-03 -5.5468E-01 -6.2942E-03  1.6448E-02
            -3.2076E-02

0ITERATION NO.:   44    OBJECTIVE VALUE:  -1702.25033667038        NO. OF FUNC. EVALS.: 133
 CUMULATIVE NO. OF FUNC. EVALS.:     1174
 NPARAMETR:  1.0041E+00  1.0907E+00  1.2361E+00  9.9434E-01  1.0539E+00  9.2296E-01  1.1891E+00  1.5930E+00  9.2847E-01  8.1371E-01
             9.6003E-01
 PARAMETER:  1.0413E-01  1.8678E-01  3.1199E-01  9.4328E-02  1.5252E-01  1.9830E-02  2.7317E-01  5.6562E-01  2.5779E-02 -1.0615E-01
             5.9208E-02
 GRADIENT:   5.8266E+05 -1.6242E+05 -1.9445E+05  6.6645E-02  3.9780E+05 -1.3569E-04 -2.2211E+05 -2.1657E+02 -6.0674E+05 -1.0648E-03
            -3.0337E+05

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1174
 NO. OF SIG. DIGITS IN FINAL EST.:  3.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3573E-03 -6.6153E-03 -4.0716E-02  1.5460E-03 -4.0663E-02
 SE:             2.9852E-02  2.1137E-02  1.7067E-02  2.2627E-02  1.9132E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6373E-01  7.5430E-01  1.7051E-02  9.4553E-01  3.3553E-02

 ETASHRINKSD(%)  1.0000E-10  2.9189E+01  4.2823E+01  2.4197E+01  3.5906E+01
 ETASHRINKVR(%)  1.0000E-10  4.9858E+01  6.7307E+01  4.2539E+01  5.8920E+01
 EBVSHRINKSD(%)  4.6264E-01  2.9800E+01  4.5242E+01  2.5012E+01  3.3537E+01
 EBVSHRINKVR(%)  9.2315E-01  5.0720E+01  7.0015E+01  4.3767E+01  5.5826E+01
 RELATIVEINF(%)  9.8285E+01  1.6394E+00  2.9843E+00  1.9595E+00  1.0300E+01
 EPSSHRINKSD(%)  4.5357E+01
 EPSSHRINKVR(%)  7.0142E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1702.2503366703838     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -967.09951010664565     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.74
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.50
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1702.250       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.09E+00  1.24E+00  9.94E-01  1.05E+00  9.23E-01  1.19E+00  1.59E+00  9.28E-01  8.14E-01  9.60E-01
 


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
+        1.39E+09
 
 TH 2
+       -7.12E+08  3.66E+08
 
 TH 3
+        3.76E+08 -1.93E+08  1.02E+08
 
 TH 4
+       -4.92E+04  2.57E+04 -3.95E+08  7.29E+02
 
 TH 5
+        1.37E+04 -7.20E+03  2.45E+08 -3.19E+04  5.87E+08
 
 TH 6
+        7.64E+03 -3.92E+03 -4.26E+08 -3.88E+00  4.97E+03  1.78E+09
 
 TH 7
+       -4.47E+08  2.29E+08 -1.21E+08  1.58E+04 -4.42E+03 -2.46E+03  1.44E+08
 
 TH 8
+       -2.05E-01 -1.29E+01 -8.73E+07  2.05E+01 -4.01E+00 -1.88E-01  2.31E+00  3.74E+07
 
 TH 9
+       -4.52E+04  2.32E+04 -4.24E+08  5.55E+04 -1.02E+09  1.77E+09  1.46E+04  2.49E+01  1.76E+09
 
 TH10
+        3.00E+03 -1.53E+03 -4.55E+08 -1.89E+01  1.85E+03  1.90E+09 -9.57E+02  1.08E+01 -3.37E+03  7.63E+01
 
 TH11
+        4.16E+02 -2.23E+02 -4.10E+08  5.36E+04 -1.50E+04 -8.32E+03 -1.33E+02  1.00E+01  4.92E+04 -3.25E+03  1.65E+09
 
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
 #CPUT: Total CPU Time in Seconds,       22.307
Stop Time:
Sat Sep 18 12:01:14 CDT 2021
