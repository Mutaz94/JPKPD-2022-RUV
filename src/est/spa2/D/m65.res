Thu Sep 30 09:27:14 CDT 2021
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
$DATA ../../../../data/spa2/D/dat65.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 NO. OF DATA RECS IN DATA SET:      700
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

 TOT. NO. OF OBS RECS:      600
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
 RAW OUTPUT FILE (FILE): m65.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   19198.4711711548        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.7149E+02  4.1950E+02 -7.8663E+00  1.4473E+02  5.5685E+02 -3.0985E+03 -1.0238E+03 -7.1941E+01 -1.8198E+03 -1.1479E+03
            -3.5867E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -674.728944163489        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.3182E+00  1.3012E+00  9.4088E-01  2.9059E+00  9.5750E-01  3.6385E+00  2.6278E+00  9.4130E-01  4.1681E+00  1.3957E+00
             1.0974E+01
 PARAMETER:  3.7623E-01  3.6332E-01  3.9058E-02  1.1668E+00  5.6575E-02  1.3916E+00  1.0662E+00  3.9511E-02  1.5275E+00  4.3342E-01
             2.4956E+00
 GRADIENT:  -3.5452E+00 -1.8534E+01 -5.3282E+01  9.4095E+01  1.4827E+01  1.2520E+02 -4.9307E+01  3.1090E+00  4.5753E+01  1.2021E+01
             1.5137E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -749.538148358602        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.4043E+00  2.4131E+00  3.5627E+00  2.9084E+00  3.1806E+00  3.0072E+00  4.0315E+00  6.4345E-01  8.2619E+00  2.4604E+00
             1.0249E+01
 PARAMETER:  4.3955E-01  9.8092E-01  1.3705E+00  1.1676E+00  1.2571E+00  1.2010E+00  1.4941E+00 -3.4091E-01  2.2117E+00  1.0003E+00
             2.4272E+00
 GRADIENT:   2.5635E+01  3.7744E+00 -2.0462E+01  5.5232E+01 -1.6500E+01  6.7340E+01  2.3629E+01  6.3391E-01  7.9575E+01  3.4809E+01
             1.5504E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -791.029924054556        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      341
 NPARAMETR:  1.1808E+00  3.3026E+00  3.9166E+00  1.2309E+00  2.6249E+00  2.4766E+00  4.4074E+00  1.5767E-01  9.0677E+00  1.4871E+00
             9.8126E+00
 PARAMETER:  2.6622E-01  1.2947E+00  1.4652E+00  3.0775E-01  1.0650E+00  1.0069E+00  1.5833E+00 -1.7472E+00  2.3047E+00  4.9682E-01
             2.3837E+00
 GRADIENT:  -2.1526E+01 -7.5966E+00 -4.7921E-01  2.0588E+01  5.3416E+00 -5.7918E+00  3.0235E+01  7.3945E-02  6.2177E+00  1.4944E+01
             8.6902E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -834.182045498021        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      519
 NPARAMETR:  1.2448E+00  3.8648E+00  3.1549E+00  3.5512E-01  1.9461E+00  2.5826E+00  3.6451E+00  8.2789E-02  9.9712E+00  5.1384E-01
             8.4247E+00
 PARAMETER:  3.1899E-01  1.4519E+00  1.2490E+00 -9.3529E-01  7.6582E-01  1.0488E+00  1.3934E+00 -2.3915E+00  2.3997E+00 -5.6584E-01
             2.2312E+00
 GRADIENT:  -1.0612E+01  1.0538E+01 -7.2176E+00  8.4196E+00 -2.9300E+01 -3.7825E+01 -2.5593E+01  3.6817E-02  1.0894E+01  2.6864E+00
            -2.7421E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -838.259314299473        NO. OF FUNC. EVALS.: 124
 CUMULATIVE NO. OF FUNC. EVALS.:      643
 NPARAMETR:  1.2532E+00  3.4952E+00  3.3582E+00  3.5646E-01  2.0971E+00  2.5725E+00  3.6257E+00  5.3831E-02  1.0062E+01  3.4685E-01
             8.4681E+00
 PARAMETER:  3.2572E-01  1.3514E+00  1.3114E+00 -9.3154E-01  8.4054E-01  1.0449E+00  1.3880E+00 -2.8219E+00  2.4088E+00 -9.5886E-01
             2.2363E+00
 GRADIENT:  -8.8528E+00 -6.6637E+00 -1.7009E+00  6.0350E+00 -5.1597E+00 -3.9746E+01 -2.5999E+01  1.6377E-02  7.3698E+00  1.1977E+00
            -9.3647E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -844.861332956279        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      804
 NPARAMETR:  1.2842E+00  3.4894E+00  3.4601E+00  1.9327E-01  2.1588E+00  2.8601E+00  3.8641E+00  1.0000E-02  8.5795E+00  2.6597E-02
             8.5327E+00
 PARAMETER:  3.5012E-01  1.3497E+00  1.3413E+00 -1.5437E+00  8.6955E-01  1.1509E+00  1.4517E+00 -4.7917E+00  2.2494E+00 -3.5269E+00
             2.2439E+00
 GRADIENT:   1.5365E+01  5.4973E+01 -1.8556E-02  3.2571E+00  7.3892E+00  5.5641E+01  8.3545E+01  0.0000E+00  7.7914E+00  7.5719E-03
             3.6686E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -855.895376508131        NO. OF FUNC. EVALS.: 121
 CUMULATIVE NO. OF FUNC. EVALS.:      925
 NPARAMETR:  1.2909E+00  3.2708E+00  4.4430E+00  5.4403E-01  2.0836E+00  2.8545E+00  3.5984E+00  1.0000E-02  1.4795E-02  1.6427E-02
             8.5501E+00
 PARAMETER:  3.5535E-01  1.2850E+00  1.5913E+00 -5.0874E-01  8.3410E-01  1.1489E+00  1.3805E+00 -4.7917E+00 -4.1134E+00 -4.0088E+00
             2.2459E+00
 GRADIENT:   1.0435E+01  6.3439E+01 -4.5950E-02  1.4265E+01  2.6128E+00  4.8361E+01 -3.9240E+01  0.0000E+00  1.7951E-03  2.9471E-03
             1.5921E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -877.556749092316        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:     1090
 NPARAMETR:  1.3023E+00  2.3884E+00  6.2773E+00  8.4564E-01  1.9714E+00  2.8548E+00  4.9390E+00  1.0000E-02  1.5611E-02  1.0000E-02
             8.5185E+00
 PARAMETER:  3.6415E-01  9.7061E-01  1.9369E+00 -6.7666E-02  7.7874E-01  1.1490E+00  1.6972E+00 -4.7917E+00 -4.0597E+00 -4.5738E+00
             2.2422E+00
 GRADIENT:  -2.6069E-01  2.1490E+01  7.1479E-02  3.8316E+01 -1.2783E+01  2.1633E+00 -4.4149E+01  0.0000E+00 -4.0058E-03  0.0000E+00
            -4.8784E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -886.941353942640        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1272
 NPARAMETR:  1.3519E+00  1.8246E+00  2.3185E+01  8.7310E-01  2.0843E+00  2.8526E+00  5.8011E+00  1.0000E-02  1.0000E-02  1.0000E-02
             8.5245E+00
 PARAMETER:  4.0155E-01  7.0123E-01  3.2217E+00 -3.6709E-02  8.3459E-01  1.1485E+00  1.8577E+00 -4.7917E+00 -4.5657E+00 -4.7208E+00
             2.2434E+00
 GRADIENT:   1.1360E+03 -1.2989E+03 -4.0590E-02 -6.5348E+00  4.5686E-01  7.3107E+02 -4.8641E+02  0.0000E+00 -3.8827E-04  0.0000E+00
             4.6571E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1272
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.4966E-03  3.0113E-04  1.0662E-06 -7.8742E-04  2.3203E-05
 SE:             3.0873E-02  2.6893E-02  1.4637E-06  1.8439E-04  7.1695E-05
 N:                     100         100         100         100         100

 P VAL.:         8.3333E-01  9.9107E-01  4.6633E-01  1.9516E-05  7.4621E-01

 ETASHRINKSD(%)  1.0000E-10  9.9044E+00  9.9995E+01  9.9382E+01  9.9760E+01
 ETASHRINKVR(%)  1.0000E-10  1.8828E+01  1.0000E+02  9.9996E+01  9.9999E+01
 EBVSHRINKSD(%)  2.6807E+00  5.1693E+00  9.9994E+01  9.9641E+01  9.9694E+01
 EBVSHRINKVR(%)  5.2895E+00  1.0071E+01  1.0000E+02  9.9999E+01  9.9999E+01
 RELATIVEINF(%)  9.4148E+01  5.6932E+01  1.1087E-07  7.9160E-04  2.6057E-04
 EPSSHRINKSD(%)  8.4997E+00
 EPSSHRINKVR(%)  1.6277E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -886.94135394264003     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       215.78488590296706     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    31.83
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    11.30
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -886.941       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.35E+00  1.82E+00  2.27E+01  8.72E-01  2.08E+00  2.85E+00  5.80E+00  1.00E-02  1.00E-02  1.00E-02  8.53E+00
 


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
+        7.71E+04
 
 TH 2
+       -3.28E+04  1.40E+04
 
 TH 3
+        1.41E-04  1.46E-03  1.37E-04
 
 TH 4
+       -5.21E+00  2.21E+01  3.34E-03  3.00E+06
 
 TH 5
+        9.16E+01 -4.18E+01 -3.93E-02 -1.51E+05  7.60E+03
 
 TH 6
+       -7.95E+01  2.09E+00  2.78E-04  2.77E+00  1.48E+01  2.46E+01
 
 TH 7
+        2.49E+01 -1.60E-01 -8.17E-04 -2.08E+01 -3.61E+00  5.96E+01  2.29E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00 -2.34E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        1.16E+02 -4.91E+01  4.25E-04 -1.40E+04  7.03E+02  3.14E+02 -5.27E+00  0.00E+00  0.00E+00  0.00E+00  8.19E+01
 
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
 #CPUT: Total CPU Time in Seconds,       43.242
Stop Time:
Thu Sep 30 09:27:58 CDT 2021
