Wed Sep 29 11:23:18 CDT 2021
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
$DATA ../../../../data/spa/B/dat61.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m61.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1712.50682343844        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.8208E+02 -1.4063E+01 -1.7185E+01  2.1451E+01  2.2812E+01  5.9529E+01  1.1317E+01  1.0752E+01  2.8536E+01  1.3485E+01
             5.9346E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1716.61359407965        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      194
 NPARAMETR:  1.0276E+00  1.0513E+00  1.0456E+00  1.0218E+00  1.0192E+00  9.6607E-01  9.3161E-01  9.2168E-01  8.4677E-01  9.1762E-01
             1.0059E+00
 PARAMETER:  1.2721E-01  1.5001E-01  1.4456E-01  1.2154E-01  1.1900E-01  6.5485E-02  2.9157E-02  1.8441E-02 -6.6329E-02  1.4024E-02
             1.0586E-01
 GRADIENT:   7.4468E+00  1.9851E+01  8.1912E+00  1.8063E+01  6.4607E+00 -5.0195E+00 -4.1441E+00  1.1600E+00 -1.1690E+01 -6.6712E+00
            -3.0406E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1719.34388282950        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  1.0296E+00  9.7618E-01  7.9468E-01  1.0500E+00  8.7124E-01  1.0041E+00  1.1260E+00  4.2484E-01  8.0740E-01  8.2003E-01
             1.0001E+00
 PARAMETER:  1.2917E-01  7.5897E-02 -1.2982E-01  1.4878E-01 -3.7837E-02  1.0405E-01  2.1871E-01 -7.5605E-01 -1.1393E-01 -9.8412E-02
             1.0010E-01
 GRADIENT:   7.2709E+00  4.8200E+00 -1.9234E+01  3.2201E+01  3.6356E+01  9.0793E+00  1.7863E+00  1.1596E+00 -2.6190E+00 -3.3525E+00
             1.4411E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1720.51936019634        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      550
 NPARAMETR:  1.0222E+00  8.4519E-01  7.6960E-01  1.1152E+00  7.8814E-01  9.7783E-01  1.2568E+00  2.7937E-01  7.6720E-01  7.9999E-01
             9.8954E-01
 PARAMETER:  1.2199E-01 -6.8195E-02 -1.6189E-01  2.0901E-01 -1.3808E-01  7.7576E-02  3.2861E-01 -1.1752E+00 -1.6501E-01 -1.2316E-01
             8.9486E-02
 GRADIENT:  -8.2004E+00  5.1625E+00 -2.9437E+00  1.0637E+01  3.2289E+00 -8.5363E-01  6.1212E-01  5.3491E-01 -7.4924E-01  4.9816E-02
            -6.6077E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1720.77760126520        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      726
 NPARAMETR:  1.0243E+00  7.2302E-01  7.8551E-01  1.1817E+00  7.5206E-01  9.7777E-01  1.4171E+00  1.3665E-01  7.3591E-01  8.1730E-01
             9.8871E-01
 PARAMETER:  1.2402E-01 -2.2432E-01 -1.4142E-01  2.6692E-01 -1.8494E-01  7.7517E-02  4.4862E-01 -1.8903E+00 -2.0665E-01 -1.0175E-01
             8.8648E-02
 GRADIENT:  -9.0898E-01  5.8489E-01  2.0856E-01  4.7862E-01 -4.4575E-01 -1.4536E-01  9.7030E-02  9.1481E-02  4.1095E-01  3.9290E-02
            -1.9080E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1720.80652600437        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      886
 NPARAMETR:  1.0257E+00  7.1177E-01  7.8799E-01  1.1866E+00  7.5012E-01  9.7812E-01  1.4330E+00  5.8912E-02  7.3106E-01  8.1963E-01
             9.8909E-01
 PARAMETER:  1.2540E-01 -2.4000E-01 -1.3827E-01  2.7109E-01 -1.8753E-01  7.7881E-02  4.5977E-01 -2.7317E+00 -2.1326E-01 -9.8902E-02
             8.9029E-02
 GRADIENT:   2.6735E+00 -6.7467E-01  1.7459E+00 -3.2689E+00 -6.0982E-01  9.4476E-02 -1.9003E-01  1.2664E-02 -2.4547E-01 -7.9430E-01
            -3.1893E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1720.81275681600        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1062
 NPARAMETR:  1.0243E+00  7.1180E-01  7.8643E-01  1.1879E+00  7.4946E-01  9.7782E-01  1.4340E+00  2.9257E-02  7.3140E-01  8.2435E-01
             9.8911E-01
 PARAMETER:  1.2405E-01 -2.3995E-01 -1.4025E-01  2.7217E-01 -1.8840E-01  7.7573E-02  4.6044E-01 -3.4316E+00 -2.1279E-01 -9.3162E-02
             8.9050E-02
 GRADIENT:  -5.9578E-01  1.1223E-01 -2.4360E-01  2.2061E-01  5.1252E-01 -4.3292E-02 -2.3642E-02  4.0175E-03  2.6585E-02  1.8511E-01
             4.1553E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1720.81446030273        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1238
 NPARAMETR:  1.0244E+00  7.0807E-01  7.8403E-01  1.1894E+00  7.4652E-01  9.7794E-01  1.4410E+00  1.0000E-02  7.3000E-01  8.2103E-01
             9.8897E-01
 PARAMETER:  1.2415E-01 -2.4522E-01 -1.4331E-01  2.7345E-01 -1.9234E-01  7.7689E-02  4.6533E-01 -5.2091E+00 -2.1471E-01 -9.7196E-02
             8.8912E-02
 GRADIENT:  -3.4810E-01  3.2438E-02  1.0048E-01 -3.0623E-01  9.7962E-02  9.9806E-03 -1.6670E-02  0.0000E+00 -9.8350E-03 -3.5660E-02
            -8.3930E-03

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1720.81617583597        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:     1336
 NPARAMETR:  1.0262E+00  7.0804E-01  7.8395E-01  1.1887E+00  7.4647E-01  9.7825E-01  1.4421E+00  1.0000E-02  7.3014E-01  8.2121E-01
             9.8899E-01
 PARAMETER:  1.2588E-01 -2.4526E-01 -1.4341E-01  2.7288E-01 -1.9240E-01  7.8011E-02  4.6612E-01 -5.2091E+00 -2.1452E-01 -9.6975E-02
             8.8931E-02
 GRADIENT:   3.7134E+00 -3.6523E-01  4.7878E-01 -2.0311E+00 -2.8692E-01  1.3936E-01  9.9951E-02  0.0000E+00  6.3059E-02  6.3438E-02
             1.3716E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1336
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -5.7457E-04  1.0810E-02 -4.6819E-04 -1.1987E-02 -5.6027E-03
 SE:             2.9849E-02  2.0096E-02  2.1763E-04  2.4655E-02  2.3554E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8464E-01  5.9062E-01  3.1448E-02  6.2682E-01  8.1199E-01

 ETASHRINKSD(%)  2.4862E-03  3.2677E+01  9.9271E+01  1.7403E+01  2.1090E+01
 ETASHRINKVR(%)  4.9724E-03  5.4676E+01  9.9995E+01  3.1777E+01  3.7732E+01
 EBVSHRINKSD(%)  4.3082E-01  3.3634E+01  9.9308E+01  1.7100E+01  1.9214E+01
 EBVSHRINKVR(%)  8.5977E-01  5.5956E+01  9.9995E+01  3.1275E+01  3.4736E+01
 RELATIVEINF(%)  9.8472E+01  3.4216E+00  4.2992E-04  6.6750E+00  4.4905E+00
 EPSSHRINKSD(%)  4.3395E+01
 EPSSHRINKVR(%)  6.7958E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1720.8161758359699     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -985.66534927223177     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.98
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.91
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1720.816       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  7.08E-01  7.84E-01  1.19E+00  7.46E-01  9.78E-01  1.44E+00  1.00E-02  7.30E-01  8.21E-01  9.89E-01
 


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
+        1.09E+03
 
 TH 2
+       -1.22E+01  4.92E+02
 
 TH 3
+        2.26E+01  2.75E+02  9.10E+02
 
 TH 4
+       -1.10E+01  4.26E+02 -3.42E+02  1.02E+03
 
 TH 5
+       -5.68E+00 -5.12E+02 -1.21E+03  3.51E+02  1.92E+03
 
 TH 6
+       -1.55E+00 -3.02E+00  4.52E+00 -2.03E+00 -1.71E+00  2.05E+02
 
 TH 7
+        1.37E+00  3.58E+01  2.35E+00 -7.06E+00 -9.48E+00  1.49E-01  2.57E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.86E+00 -2.35E+01 -4.26E+01  1.13E+01  2.32E+01 -5.40E-01  1.99E+01  0.00E+00  1.82E+02
 
 TH10
+       -1.53E+00  6.41E-01 -7.94E+01 -3.13E+01 -4.23E+01  5.23E-02  8.10E+00  0.00E+00  1.32E+01  1.22E+02
 
 TH11
+       -7.37E+00 -1.38E+01 -4.68E+01 -7.14E+00  1.72E+01  2.06E+00  3.90E+00  0.00E+00  1.38E+01  2.92E+01  2.19E+02
 
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
 #CPUT: Total CPU Time in Seconds,       22.938
Stop Time:
Wed Sep 29 11:23:42 CDT 2021
