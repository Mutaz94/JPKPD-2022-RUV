Sat Sep 25 13:14:35 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat99.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1684.31087806319        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.0053E+01 -6.2841E+01 -2.6529E+01 -4.5342E+01  4.9003E+01  3.4548E+01 -1.2386E+01  2.5953E+00  1.8300E+00  7.9455E+00
             8.1377E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1689.87401416856        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0024E+00  1.0903E+00  9.7477E-01  9.8928E-01  9.8821E-01  8.9868E-01  1.1756E+00  1.0091E+00  9.6558E-01  8.7326E-01
             9.9936E-01
 PARAMETER:  1.0237E-01  1.8650E-01  7.4442E-02  8.9223E-02  8.8136E-02 -6.8279E-03  2.6176E-01  1.0904E-01  6.4969E-02 -3.5523E-02
             9.9362E-02
 GRADIENT:   3.2320E+01  1.0499E+01 -1.1916E+00  1.8068E+01  2.0724E+01 -5.3541E+00  2.5029E+00 -1.3455E+00  2.2351E+00 -8.6881E-01
             4.2329E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1690.58407789809        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0031E+00  1.0365E+00  8.1305E-01  1.0115E+00  8.6057E-01  9.0452E-01  1.2237E+00  1.0017E+00  9.1334E-01  6.6165E-01
             9.8020E-01
 PARAMETER:  1.0306E-01  1.3586E-01 -1.0697E-01  1.1143E-01 -5.0162E-02 -3.5001E-04  3.0185E-01  1.0172E-01  9.3528E-03 -3.1302E-01
             8.0005E-02
 GRADIENT:   3.2621E+01  1.2776E+01  1.0927E+00  2.5022E+01  8.9919E-01 -3.3627E+00 -1.1027E+00  1.2789E-01  8.3417E-01 -3.5571E+00
            -1.6512E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1690.84077699452        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      252
 NPARAMETR:  9.9593E-01  1.0963E+00  6.9208E-01  9.5475E-01  8.3134E-01  9.0903E-01  1.1917E+00  7.5045E-01  9.1705E-01  6.7914E-01
             9.8297E-01
 PARAMETER:  9.5924E-02  1.9197E-01 -2.6805E-01  5.3694E-02 -8.4718E-02  4.6193E-03  2.7540E-01 -1.8708E-01  1.3409E-02 -2.8692E-01
             8.2821E-02
 GRADIENT:  -3.3680E+01 -1.3946E+00 -2.6975E+00  4.4710E+00  1.6429E+00 -6.7334E+00 -1.2604E+00  3.9319E-01 -1.3399E+00 -8.8178E-01
             4.9920E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1691.41397061398        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      429
 NPARAMETR:  1.0077E+00  1.2854E+00  6.0488E-01  8.3382E-01  8.7748E-01  9.2341E-01  1.0558E+00  5.4252E-01  1.0131E+00  7.3007E-01
             9.8355E-01
 PARAMETER:  1.0763E-01  3.5108E-01 -4.0273E-01 -8.1741E-02 -3.0700E-02  2.0319E-02  1.5427E-01 -5.1153E-01  1.1302E-01 -2.1461E-01
             8.3408E-02
 GRADIENT:  -4.0216E+00  1.1921E+00  6.2892E-01  1.5026E+00 -7.0350E-01 -5.4750E-01 -1.1061E+00 -2.9347E-02 -6.3774E-01 -1.9344E-01
             8.6949E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1691.46206731550        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      604
 NPARAMETR:  1.0094E+00  1.3847E+00  5.4038E-01  7.6593E-01  8.9346E-01  9.2516E-01  1.0036E+00  4.0177E-01  1.0705E+00  7.3262E-01
             9.8362E-01
 PARAMETER:  1.0936E-01  4.2550E-01 -5.1548E-01 -1.6667E-01 -1.2654E-02  2.2207E-02  1.0358E-01 -8.1187E-01  1.6808E-01 -2.1112E-01
             8.3482E-02
 GRADIENT:  -2.2825E-02 -1.1335E-01 -1.2766E-01 -2.4266E-01 -3.4893E-01  3.5551E-02  3.2891E-02  8.0444E-02  9.3067E-02  2.3673E-01
             7.9480E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1691.47844373883        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      781
 NPARAMETR:  1.0094E+00  1.4025E+00  5.2097E-01  7.5331E-01  8.9199E-01  9.2520E-01  9.9641E-01  2.1401E-01  1.0809E+00  7.3632E-01
             9.8398E-01
 PARAMETER:  1.0934E-01  4.3826E-01 -5.5207E-01 -1.8328E-01 -1.4302E-02  2.2251E-02  9.6400E-02 -1.4418E+00  1.7777E-01 -2.0609E-01
             8.3850E-02
 GRADIENT:  -1.9730E-01  6.6937E-01 -1.2377E-01  3.4268E-01 -6.1427E-01 -4.2412E-03  2.1035E-01  1.8706E-02  9.4009E-02  2.6138E-01
             1.4454E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1691.48473913177        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      956
 NPARAMETR:  1.0094E+00  1.4092E+00  5.1522E-01  7.4805E-01  8.9312E-01  9.2523E-01  9.9276E-01  5.7127E-02  1.0860E+00  7.3787E-01
             9.8399E-01
 PARAMETER:  1.0937E-01  4.4304E-01 -5.6315E-01 -1.9029E-01 -1.3038E-02  2.2289E-02  9.2733E-02 -2.7625E+00  1.8250E-01 -2.0398E-01
             8.3857E-02
 GRADIENT:  -1.0468E-01 -1.1311E-01 -4.6593E-02 -2.1900E-02  7.3862E-02 -6.3900E-03  5.0001E-02  8.8525E-04 -2.8233E-03  1.5723E-02
             3.1751E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1691.48518502095        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1131
 NPARAMETR:  1.0095E+00  1.4122E+00  5.1378E-01  7.4615E-01  8.9389E-01  9.2525E-01  9.9094E-01  1.0661E-02  1.0880E+00  7.3817E-01
             9.8397E-01
 PARAMETER:  1.0941E-01  4.4518E-01 -5.6596E-01 -1.9283E-01 -1.2175E-02  2.2305E-02  9.0900E-02 -4.4412E+00  1.8435E-01 -2.0358E-01
             8.3836E-02
 GRADIENT:  -1.7855E-02 -2.9678E-03 -2.6358E-02  3.5334E-02  2.3524E-02 -2.5084E-03  1.9824E-03  3.2401E-05  1.3034E-02  1.0310E-02
             9.7197E-03

0ITERATION NO.:   42    OBJECTIVE VALUE:  -1691.48519105494        NO. OF FUNC. EVALS.:  63
 CUMULATIVE NO. OF FUNC. EVALS.:     1194
 NPARAMETR:  1.0095E+00  1.4119E+00  5.1388E-01  7.4633E-01  8.9379E-01  9.2527E-01  9.9102E-01  1.0000E-02  1.0878E+00  7.3802E-01
             9.8394E-01
 PARAMETER:  1.0941E-01  4.4498E-01 -5.6578E-01 -1.9263E-01 -1.2292E-02  2.2310E-02  9.1082E-02 -4.6721E+00  1.8415E-01 -2.0372E-01
             8.3819E-02
 GRADIENT:  -1.7922E-02  2.8686E-03 -4.5562E-04 -9.3450E-03 -5.0738E-03 -1.8468E-03  4.3933E-03  0.0000E+00  1.1141E-03  1.3275E-03
             7.4249E-04

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1194
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2418E-04 -1.3985E-02 -3.5620E-04  1.1913E-02 -2.2638E-02
 SE:             2.9834E-02  2.4871E-02  1.4877E-04  2.3921E-02  2.0537E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9668E-01  5.7391E-01  1.6655E-02  6.1849E-01  2.7032E-01

 ETASHRINKSD(%)  5.1956E-02  1.6680E+01  9.9502E+01  1.9860E+01  3.1200E+01
 ETASHRINKVR(%)  1.0389E-01  3.0577E+01  9.9998E+01  3.5776E+01  5.2666E+01
 EBVSHRINKSD(%)  4.8101E-01  1.6499E+01  9.9566E+01  2.0307E+01  3.1243E+01
 EBVSHRINKVR(%)  9.5972E-01  3.0276E+01  9.9998E+01  3.6490E+01  5.2725E+01
 RELATIVEINF(%)  9.9017E+01  5.1956E+00  1.6693E-04  4.4865E+00  5.3447E+00
 EPSSHRINKSD(%)  4.4653E+01
 EPSSHRINKVR(%)  6.9367E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1691.4851910549376     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -956.33436449119938     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.71
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.65
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1691.485       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.41E+00  5.14E-01  7.46E-01  8.94E-01  9.25E-01  9.91E-01  1.00E-02  1.09E+00  7.38E-01  9.84E-01
 


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
+        1.26E+03
 
 TH 2
+       -5.62E+00  4.07E+02
 
 TH 3
+        8.90E+00  2.59E+02  8.48E+02
 
 TH 4
+       -1.50E+01  2.59E+02 -5.56E+02  1.09E+03
 
 TH 5
+       -4.01E+00 -3.55E+02 -8.52E+02  5.39E+02  1.22E+03
 
 TH 6
+       -1.30E+00 -1.01E+00  1.95E+00 -3.76E+00 -1.33E+00  2.30E+02
 
 TH 7
+       -9.44E-01  2.06E+01 -3.76E+01 -9.95E+00  1.03E+01  1.29E+00  1.02E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.03E+00 -1.99E+01 -4.87E+01  5.74E+01 -1.45E+01 -1.41E+00  1.32E+01  0.00E+00  7.71E+01
 
 TH10
+       -5.86E-01 -1.53E+01 -5.18E+01 -2.66E+01 -8.03E+01  2.83E-01  1.96E+01  0.00E+00  1.67E+01  9.31E+01
 
 TH11
+       -8.76E+00 -1.36E+01 -2.86E+01 -1.95E+00 -6.49E+00  4.53E+00  9.65E+00  0.00E+00  9.31E+00  2.12E+01  2.15E+02
 
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
 #CPUT: Total CPU Time in Seconds,       19.385
Stop Time:
Sat Sep 25 13:14:56 CDT 2021
