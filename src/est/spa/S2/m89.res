Sat Sep 18 13:42:38 CDT 2021
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
$DATA ../../../../data/spa/S2/dat89.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m89.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1674.34969608406        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0753E+01 -4.7961E+01 -2.2018E+01 -2.9190E+01  4.5142E+01  4.0049E+01  1.2727E+00  4.2552E+00  1.8353E+01 -6.6878E+00
             1.0993E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1678.82762905834        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.9497E-01  1.0337E+00  1.0097E+00  9.9993E-01  9.7818E-01  8.8543E-01  9.9058E-01  9.8411E-01  9.2399E-01  1.0215E+00
             9.7475E-01
 PARAMETER:  9.4958E-02  1.3314E-01  1.0964E-01  9.9932E-02  7.7938E-02 -2.1683E-02  9.0539E-02  8.3980E-02  2.0942E-02  1.2127E-01
             7.4421E-02
 GRADIENT:   4.1412E+01 -2.6654E+00  5.3588E+00 -1.2364E+01 -7.5535E+00 -3.5220E+00 -6.2855E-01  1.4443E+00  1.3604E+00 -6.9624E-01
            -1.6178E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1679.91197052013        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.9333E-01  8.6253E-01  8.6551E-01  1.1215E+00  8.4573E-01  8.9258E-01  1.2839E+00  6.2232E-01  7.9002E-01  9.1727E-01
             9.8066E-01
 PARAMETER:  9.3303E-02 -4.7881E-02 -4.4436E-02  2.1463E-01 -6.7553E-02 -1.3634E-02  3.4992E-01 -3.7431E-01 -1.3570E-01  1.3650E-02
             8.0469E-02
 GRADIENT:   3.3195E+01  1.7896E+01 -1.2848E+01  5.5256E+01  2.3237E+01 -1.0488E+00  1.1110E+00  2.6588E-01 -4.8985E+00  1.9068E+00
             5.0645E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1680.98509294505        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.8145E-01  8.4009E-01  5.9920E-01  1.0853E+00  6.7208E-01  8.9790E-01  1.3063E+00  2.9945E-01  7.8968E-01  6.9678E-01
             9.7997E-01
 PARAMETER:  8.1276E-02 -7.4241E-02 -4.1216E-01  1.8188E-01 -2.9737E-01 -7.7013E-03  3.6719E-01 -1.1058E+00 -1.3613E-01 -2.6129E-01
             7.9771E-02
 GRADIENT:  -8.2671E+00  6.4880E+00 -2.1336E+01  3.2936E+01  2.8370E+01 -2.8044E-01  1.1302E+00  1.5160E+00 -1.0791E+00  3.2334E+00
             1.1919E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1681.15474035326        NO. OF FUNC. EVALS.: 110
 CUMULATIVE NO. OF FUNC. EVALS.:      339
 NPARAMETR:  9.7997E-01  7.9157E-01  5.1987E-01  1.0957E+00  5.9842E-01  8.9910E-01  1.3723E+00  2.1842E-01  7.6818E-01  6.0489E-01
             9.8055E-01
 PARAMETER:  7.9764E-02 -1.3374E-01 -5.5418E-01  1.9136E-01 -4.1346E-01 -6.3597E-03  4.1647E-01 -1.4213E+00 -1.6373E-01 -4.0270E-01
             8.0354E-02
 GRADIENT:  -5.3145E+01  8.9286E+00 -2.6035E+01  2.5164E+01  2.1333E+01 -2.9668E+00  1.0615E+00  1.2001E+00 -1.7216E+00  4.4787E+00
             8.5278E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1683.65991062382        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      517
 NPARAMETR:  9.9567E-01  6.3689E-01  5.5650E-01  1.1946E+00  5.6329E-01  8.9262E-01  1.6907E+00  1.0874E-01  7.3467E-01  6.4167E-01
             9.7822E-01
 PARAMETER:  9.5659E-02 -3.5116E-01 -4.8609E-01  2.7781E-01 -4.7395E-01 -1.3591E-02  6.2516E-01 -2.1188E+00 -2.0833E-01 -3.4368E-01
             7.7983E-02
 GRADIENT:  -7.3797E+00  2.2221E+01 -6.2573E+00  3.5705E+01 -1.1108E+01 -4.8259E+00  3.6732E+00  1.8958E-01 -1.8810E-01  7.5268E+00
            -1.3498E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1685.58244035368        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      692
 NPARAMETR:  9.9663E-01  4.6266E-01  5.5695E-01  1.2615E+00  5.2154E-01  9.0091E-01  2.1102E+00  1.0000E-02  7.0752E-01  5.9990E-01
             9.9041E-01
 PARAMETER:  9.6628E-02 -6.7076E-01 -4.8528E-01  3.3233E-01 -5.5096E-01 -4.3503E-03  8.4680E-01 -5.2237E+00 -2.4599E-01 -4.1099E-01
             9.0359E-02
 GRADIENT:   1.0118E+00  3.0830E+00  3.1118E+00  3.9603E+00 -4.8359E+00 -6.4756E-01  7.4441E-01  0.0000E+00 -1.1200E-01 -6.3153E-01
            -1.7313E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1685.62493452649        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      867
 NPARAMETR:  9.9569E-01  4.3313E-01  5.6248E-01  1.2752E+00  5.1907E-01  9.0205E-01  2.1945E+00  1.0000E-02  7.0293E-01  6.1238E-01
             9.9102E-01
 PARAMETER:  9.5685E-02 -7.3673E-01 -4.7541E-01  3.4311E-01 -5.5572E-01 -3.0894E-03  8.8598E-01 -5.8586E+00 -2.5250E-01 -3.9041E-01
             9.0975E-02
 GRADIENT:   6.3602E-03  6.2874E-03 -2.2948E-02  4.3140E-02  1.9743E-02 -2.0999E-03  3.9920E-03  0.0000E+00 -1.4052E-02 -4.8901E-03
            -8.2001E-04

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1685.62493452649        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      889
 NPARAMETR:  9.9569E-01  4.3313E-01  5.6248E-01  1.2752E+00  5.1907E-01  9.0205E-01  2.1945E+00  1.0000E-02  7.0293E-01  6.1238E-01
             9.9102E-01
 PARAMETER:  9.5685E-02 -7.3673E-01 -4.7541E-01  3.4311E-01 -5.5572E-01 -3.0894E-03  8.8598E-01 -5.8586E+00 -2.5250E-01 -3.9041E-01
             9.0975E-02
 GRADIENT:   6.3602E-03  6.2874E-03 -2.2948E-02  4.3140E-02  1.9743E-02 -2.0999E-03  3.9920E-03  0.0000E+00 -1.4052E-02 -4.8901E-03
            -8.2001E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      889
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.0286E-04  3.7448E-02 -6.7821E-04 -2.5139E-02  2.1054E-02
 SE:             2.9838E-02  2.0804E-02  2.9663E-04  2.5890E-02  2.1820E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7586E-01  7.1857E-02  2.2233E-02  3.3156E-01  3.3460E-01

 ETASHRINKSD(%)  3.9639E-02  3.0304E+01  9.9006E+01  1.3265E+01  2.6901E+01
 ETASHRINKVR(%)  7.9263E-02  5.1425E+01  9.9990E+01  2.4771E+01  4.6566E+01
 EBVSHRINKSD(%)  5.0967E-01  3.2008E+01  9.9057E+01  1.2249E+01  2.4607E+01
 EBVSHRINKVR(%)  1.0167E+00  5.3771E+01  9.9991E+01  2.2997E+01  4.3160E+01
 RELATIVEINF(%)  9.8440E+01  9.3151E+00  5.3984E-04  2.2827E+01  2.8810E+00
 EPSSHRINKSD(%)  4.4172E+01
 EPSSHRINKVR(%)  6.8832E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1685.6249345264921     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -950.47410796275392     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.26
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.70
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1685.625       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.96E-01  4.33E-01  5.62E-01  1.28E+00  5.19E-01  9.02E-01  2.19E+00  1.00E-02  7.03E-01  6.12E-01  9.91E-01
 


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
+        1.37E+03
 
 TH 2
+       -2.14E+01  6.09E+02
 
 TH 3
+        3.67E+01  5.93E+02  3.43E+03
 
 TH 4
+       -1.02E+01  3.20E+02 -6.73E+02  1.09E+03
 
 TH 5
+       -1.07E+01 -1.10E+03 -4.44E+03  6.11E+02  6.48E+03
 
 TH 6
+        5.09E+00 -3.76E+00  4.06E+00 -2.30E+00 -1.46E+00  2.44E+02
 
 TH 7
+        1.95E+00  4.80E+01 -2.65E+01 -1.19E+01  2.30E+01  2.02E-03  1.68E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.71E+00 -2.32E+01 -1.40E+01  1.94E+01 -1.97E+01 -2.61E+00  5.03E-01  0.00E+00  2.62E+02
 
 TH10
+       -3.77E+00  7.43E-02 -2.55E+02 -5.94E+01  1.24E+02  1.24E+00  1.04E+01  0.00E+00  2.24E+01  1.75E+02
 
 TH11
+       -1.09E+01 -3.11E+00 -5.49E+01 -1.58E+01  5.13E+01 -7.71E-01  1.85E+00  0.00E+00  1.76E+01  3.13E+01  2.21E+02
 
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
 #CPUT: Total CPU Time in Seconds,       15.034
Stop Time:
Sat Sep 18 13:42:54 CDT 2021
