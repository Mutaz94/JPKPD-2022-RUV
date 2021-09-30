Thu Sep 30 03:37:55 CDT 2021
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
$DATA ../../../../data/spa1/D/dat78.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m78.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   24452.5166557307        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.1656E+02  4.0164E+02 -1.0597E+02  2.0928E+02  2.3981E+02 -2.1997E+03 -1.0485E+03 -4.9465E+01 -1.8875E+03 -4.8648E+02
            -4.6835E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -519.317174765893        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2164E+00  9.0652E-01  9.1044E-01  1.6988E+00  1.3775E+00  2.8619E+00  1.2709E+00  9.2884E-01  1.6798E+00  9.7649E-01
             1.3599E+01
 PARAMETER:  2.9587E-01  1.8613E-03  6.1745E-03  6.2991E-01  4.2024E-01  1.1515E+00  3.3971E-01  2.6178E-02  6.1869E-01  7.6212E-02
             2.7100E+00
 GRADIENT:  -8.0962E+00  1.8818E+01 -1.6185E+01 -3.8264E+00 -3.0593E+00  1.0030E+02 -2.8868E-01  6.9220E+00 -3.1197E+01  1.2782E+00
             1.0410E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -558.280138086960        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      162
 NPARAMETR:  1.2631E+00  9.1770E-01  1.2773E+00  1.9272E+00  6.8419E+00  2.2420E+00  1.5989E+00  9.8013E-02  2.2344E+00  1.6461E+00
             1.3457E+01
 PARAMETER:  3.3355E-01  1.4120E-02  3.4474E-01  7.5606E-01  2.0231E+00  9.0736E-01  5.6931E-01 -2.2227E+00  9.0396E-01  5.9840E-01
             2.6995E+00
 GRADIENT:  -2.4717E+00  2.6341E+01  9.0318E+00  4.2399E+01 -3.4126E+00  1.1342E+01  4.0566E+00  4.3827E-03 -2.6632E+00  9.6119E-02
             1.0955E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -604.373042976264        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      237
 NPARAMETR:  9.1917E-01  9.7980E-02  2.9080E-01  1.5136E+00  1.3193E+01  1.5729E+00  8.0451E-01  2.9206E-02  1.0286E+00  4.5111E+00
             1.1247E+01
 PARAMETER:  1.5713E-02 -2.2230E+00 -1.1351E+00  5.1448E-01  2.6797E+00  5.5290E-01 -1.1752E-01 -3.4334E+00  1.2816E-01  1.6065E+00
             2.5201E+00
 GRADIENT:  -3.4042E+01  5.6997E+00  2.7216E+00  1.8742E+02  6.3730E-01 -9.6315E+01  1.0572E-02 -1.6581E-02 -3.9992E+01 -9.1358E-02
            -6.2962E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -672.030634031579        NO. OF FUNC. EVALS.: 121
 CUMULATIVE NO. OF FUNC. EVALS.:      358
 NPARAMETR:  6.0469E-01  2.3257E-02  6.4001E-02  5.7807E-01  2.1908E+02  1.5656E+00  2.4273E-02  1.0000E-02  7.3768E-01  4.9096E+00
             9.8300E+00
 PARAMETER: -4.0303E-01 -3.6612E+00 -2.6489E+00 -4.4806E-01  5.4894E+00  5.4824E-01 -3.6184E+00 -7.2539E+00 -2.0424E-01  1.6912E+00
             2.3854E+00
 GRADIENT:  -1.1315E+01  1.3524E-01  4.2893E+01 -1.7707E+01  1.6114E-03 -4.8764E+01  1.5406E-06  0.0000E+00 -1.2698E+01 -1.0146E-04
            -1.4149E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -694.853275682831        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      533
 NPARAMETR:  4.4809E-01  1.2083E-02  2.4619E-02  2.9955E-01  2.0517E+03  1.6232E+00  1.0000E-02  1.0000E-02  6.4623E-01  4.9599E+00
             1.1163E+01
 PARAMETER: -7.0275E-01 -4.3160E+00 -3.6042E+00 -1.1055E+00  7.7264E+00  5.8439E-01 -9.0999E+00 -1.0653E+01 -3.3660E-01  1.7014E+00
             2.5126E+00
 GRADIENT:  -4.7134E+00  3.1854E-01  4.8931E+00 -5.3329E+00 -1.1075E-04 -1.5426E+00  0.0000E+00  0.0000E+00  1.7598E+00 -5.4766E-08
            -6.9999E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -695.466847045898        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:      723
 NPARAMETR:  4.4270E-01  1.0000E-02  2.3535E-02  2.9362E-01  1.2129E+03  1.6178E+00  1.0000E-02  1.0000E-02  4.0965E-01  3.4855E+00
             1.1530E+01
 PARAMETER: -7.1486E-01 -4.5514E+00 -3.6492E+00 -1.1255E+00  7.2007E+00  5.8108E-01 -9.0614E+00 -1.0851E+01 -7.9245E-01  1.3486E+00
             2.5449E+00
 GRADIENT:   3.4255E+00  0.0000E+00  1.2601E+00 -4.2136E+00  3.2765E-04  1.0295E+00  0.0000E+00  0.0000E+00 -8.3532E-02 -1.9048E-07
            -5.4992E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -695.528437118643        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      899
 NPARAMETR:  4.3944E-01  1.0000E-02  2.3519E-02  2.9584E-01  1.5701E+01  1.6110E+00  1.0000E-02  1.0000E-02  2.9836E-01  4.4339E+00
             1.1719E+01
 PARAMETER: -7.2226E-01 -4.5514E+00 -3.6499E+00 -1.1179E+00  2.8537E+00  5.7683E-01 -9.0614E+00 -1.0851E+01 -1.1095E+00  1.5893E+00
             2.5612E+00
 GRADIENT:  -3.0110E+00  0.0000E+00  1.6291E+00 -1.4845E+00 -4.2705E-03  5.6713E-01  0.0000E+00  0.0000E+00  7.1152E-02 -1.3732E-03
             1.7871E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -695.541730567467        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1076
 NPARAMETR:  4.3999E-01  1.0000E-02  2.3485E-02  2.9600E-01  1.6954E+01  1.6066E+00  1.0000E-02  1.0000E-02  2.6570E-01  9.8064E+00
             1.1733E+01
 PARAMETER: -7.2099E-01 -4.5514E+00 -3.6514E+00 -1.1174E+00  2.9305E+00  5.7410E-01 -9.0614E+00 -1.0851E+01 -1.2254E+00  2.3830E+00
             2.5624E+00
 GRADIENT:  -3.2916E-01  0.0000E+00 -4.7863E-01 -3.8758E-02  3.9906E-03 -6.3519E-02  0.0000E+00  0.0000E+00  6.4092E-04 -2.4226E-03
            -1.5706E-01

0ITERATION NO.:   43    OBJECTIVE VALUE:  -695.542003349906        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1168
 NPARAMETR:  4.4005E-01  1.0000E-02  2.3485E-02  2.9595E-01  1.7081E+01  1.6064E+00  1.0000E-02  1.0000E-02  2.6562E-01  1.1004E+01
             1.1733E+01
 PARAMETER: -7.2087E-01 -4.5514E+00 -3.6514E+00 -1.1176E+00  2.9380E+00  5.7401E-01 -9.0614E+00 -1.0851E+01 -1.2257E+00  2.4983E+00
             2.5624E+00
 GRADIENT:  -2.2686E-01  0.0000E+00 -2.5716E-01 -3.8454E-01 -1.7459E-03 -6.6437E-02  0.0000E+00  0.0000E+00 -2.4477E-03  1.0274E-03
            -1.2034E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1168
 NO. OF SIG. DIGITS IN FINAL EST.:  2.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.1445E-03  1.8285E-06  1.0739E-04 -6.7500E-03 -4.6193E-04
 SE:             2.8945E-02  1.0306E-06  3.0468E-04  9.7858E-03  6.9845E-04
 N:                     100         100         100         100         100

 P VAL.:         9.1349E-01  7.6039E-02  7.2448E-01  4.9033E-01  5.0838E-01

 ETASHRINKSD(%)  3.0303E+00  9.9997E+01  9.8979E+01  6.7216E+01  9.7660E+01
 ETASHRINKVR(%)  5.9687E+00  1.0000E+02  9.9990E+01  8.9252E+01  9.9945E+01
 EBVSHRINKSD(%)  3.2083E+00  9.9996E+01  9.9017E+01  6.8040E+01  9.7809E+01
 EBVSHRINKVR(%)  6.3136E+00  1.0000E+02  9.9990E+01  8.9786E+01  9.9952E+01
 RELATIVEINF(%)  1.2619E+00  3.3442E-08  3.0979E-05  2.3909E-02  9.0645E-04
 EPSSHRINKSD(%)  5.7014E+00
 EPSSHRINKVR(%)  1.1078E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -695.54200334990583     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       223.39652985476687     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.00
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.99
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -695.542       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.40E-01  1.00E-02  2.35E-02  2.96E-01  1.71E+01  1.61E+00  1.00E-02  1.00E-02  2.66E-01  1.10E+01  1.17E+01
 


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
+        2.05E+03
 
 TH 2
+        0.00E+00  5.09E+03
 
 TH 3
+       -1.50E+04  0.00E+00  1.30E+06
 
 TH 4
+        3.48E+01  0.00E+00 -1.22E+05  1.25E+04
 
 TH 5
+        1.40E-01  0.00E+00 -2.69E+00  2.14E-01  4.40E-04
 
 TH 6
+        2.19E+00  0.00E+00  3.43E+02 -7.02E+01  3.97E-04  6.41E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -4.54E+01  0.00E+00 -7.51E+02  1.10E+02 -1.39E-03 -2.23E+00  0.00E+00  0.00E+00  4.89E+00
 
 TH10
+       -1.94E-03  0.00E+00 -1.50E-01 -8.18E-03 -2.94E-04  2.73E-03  0.00E+00  0.00E+00  4.72E-04  2.70E-04
 
 TH11
+       -2.23E+01  0.00E+00  4.16E+02 -2.89E+01 -2.03E-03  8.75E-01  0.00E+00  0.00E+00  2.74E+00  7.55E-05  3.53E+00
 
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
 #CPUT: Total CPU Time in Seconds,       29.999
Stop Time:
Thu Sep 30 03:38:27 CDT 2021
