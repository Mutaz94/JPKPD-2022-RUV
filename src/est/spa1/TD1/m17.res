Thu Sep 30 01:12:27 CDT 2021
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
$DATA ../../../../data/spa1/TD1/dat17.csv ignore=@
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m17.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2105.44468875858        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.5529E+02 -3.0665E+01  2.4380E+01 -3.2898E+01 -1.9625E+01  4.4949E+01 -1.3788E+01  4.0536E+00  1.7161E+01 -3.1786E+01
            -4.5867E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2119.67371096015        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.0424E+00  1.1538E+00  1.0048E+00  9.9456E-01  1.0935E+00  1.0145E+00  1.1811E+00  9.2944E-01  8.4600E-01  1.3344E+00
             1.0429E+00
 PARAMETER:  1.4157E-01  2.4302E-01  1.0477E-01  9.4549E-02  1.8939E-01  1.1439E-01  2.6647E-01  2.6832E-02 -6.7233E-02  3.8850E-01
             1.4200E-01
 GRADIENT:   1.2056E+01  1.9850E+01  5.5827E+00  9.8108E+00 -2.6891E+01  2.9023E+00 -2.8408E+00  3.2986E+00 -4.1685E+00  1.1004E+01
            -9.3736E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2120.37019767032        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  1.0395E+00  9.8348E-01  1.0704E+00  1.0995E+00  1.0580E+00  1.0044E+00  1.3866E+00  7.5058E-01  7.4824E-01  1.2724E+00
             1.0488E+00
 PARAMETER:  1.3875E-01  8.3344E-02  1.6807E-01  1.9488E-01  1.5640E-01  1.0434E-01  4.2687E-01 -1.8691E-01 -1.9004E-01  3.4093E-01
             1.4767E-01
 GRADIENT:   9.7279E+00  2.0709E+01  1.6277E+01  1.2174E+01 -1.4661E+01 -1.2192E-01 -2.1799E+00 -3.3372E+00 -1.1193E+01  7.0620E-01
            -7.6336E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2121.32427184044        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      545
 NPARAMETR:  1.0358E+00  9.6539E-01  9.3926E-01  1.0886E+00  9.8949E-01  1.0063E+00  1.3824E+00  5.9237E-01  8.0468E-01  1.1814E+00
             1.0526E+00
 PARAMETER:  1.3515E-01  6.4779E-02  3.7333E-02  1.8488E-01  8.9437E-02  1.0631E-01  4.2380E-01 -4.2362E-01 -1.1731E-01  2.6671E-01
             1.5131E-01
 GRADIENT:  -3.7535E-01  4.0638E-01 -3.0602E-01  5.2516E-02  6.4317E-01  1.6713E-01 -1.6469E-01  2.0479E-02 -2.7424E-01 -4.2145E-01
            -4.4222E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2121.34183982721        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      720
 NPARAMETR:  1.0339E+00  8.4188E-01  1.0714E+00  1.1744E+00  1.0003E+00  1.0016E+00  1.5468E+00  7.3208E-01  7.6531E-01  1.2176E+00
             1.0524E+00
 PARAMETER:  1.3333E-01 -7.2116E-02  1.6901E-01  2.6075E-01  1.0033E-01  1.0158E-01  5.3618E-01 -2.1186E-01 -1.6748E-01  2.9687E-01
             1.5103E-01
 GRADIENT:   3.8121E-01  6.4331E+00  7.1159E+00  4.2914E+00 -7.5787E+00 -6.2046E-01  2.4263E-01 -1.3811E+00 -3.1842E-01 -8.5163E-01
            -2.1141E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2121.34613395564        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      895
 NPARAMETR:  1.0329E+00  7.9919E-01  1.1303E+00  1.2058E+00  1.0101E+00  9.9998E-01  1.6135E+00  8.0864E-01  7.5151E-01  1.2344E+00
             1.0527E+00
 PARAMETER:  1.3242E-01 -1.2415E-01  2.2248E-01  2.8715E-01  1.1005E-01  9.9979E-02  5.7841E-01 -1.1240E-01 -1.8567E-01  3.1060E-01
             1.5138E-01
 GRADIENT:   1.6329E-01  8.6214E+00  8.2866E+00  8.1212E+00 -1.0086E+01 -8.5249E-01  4.3516E-01 -1.5403E+00 -3.5176E-01 -6.1552E-01
            -2.0496E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2121.35549929412        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1070
 NPARAMETR:  1.0319E+00  7.5146E-01  1.1853E+00  1.2397E+00  1.0160E+00  9.9849E-01  1.6947E+00  8.7842E-01  7.3745E-01  1.2495E+00
             1.0537E+00
 PARAMETER:  1.3139E-01 -1.8574E-01  2.6999E-01  3.1484E-01  1.1590E-01  9.8489E-02  6.2751E-01 -2.9629E-02 -2.0456E-01  3.2272E-01
             1.5233E-01
 GRADIENT:  -2.8572E-01  1.0103E+01  8.1631E+00  1.2622E+01 -1.1763E+01 -1.0201E+00  6.7676E-01 -1.3885E+00 -3.4178E-01 -4.5402E-02
            -1.3012E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2121.48145907603        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1250
 NPARAMETR:  1.0294E+00  6.3009E-01  1.2590E+00  1.3151E+00  1.0108E+00  9.9671E-01  1.9248E+00  9.6972E-01  7.1013E-01  1.2702E+00
             1.0567E+00
 PARAMETER:  1.2900E-01 -3.6190E-01  3.3032E-01  3.7388E-01  1.1071E-01  9.6703E-02  7.5481E-01  6.9250E-02 -2.4230E-01  3.3917E-01
             1.5516E-01
 GRADIENT:  -1.3915E+00  7.4026E+00  1.2515E+00  1.6522E+01 -7.5895E+00 -7.5014E-01  8.0475E-01  1.5246E-02  4.0405E-02  1.6763E+00
             1.5918E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2121.60173589365        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1426
 NPARAMETR:  1.0275E+00  5.1342E-01  1.3862E+00  1.3867E+00  1.0347E+00  9.9557E-01  2.1940E+00  1.1125E+00  6.8700E-01  1.3014E+00
             1.0569E+00
 PARAMETER:  1.2712E-01 -5.6665E-01  4.2654E-01  4.2690E-01  1.3416E-01  9.5556E-02  8.8573E-01  2.0662E-01 -2.7543E-01  3.6346E-01
             1.5538E-01
 GRADIENT:  -9.3482E-01  3.4000E+00 -1.9000E+00  1.0536E+01 -1.7711E+00 -1.2909E-01  3.6051E-01  4.0071E-01 -4.3723E-02  1.3690E+00
             1.5091E+00

0ITERATION NO.:   42    OBJECTIVE VALUE:  -2121.60225941283        NO. OF FUNC. EVALS.:  64
 CUMULATIVE NO. OF FUNC. EVALS.:     1490
 NPARAMETR:  1.0275E+00  5.1167E-01  1.3885E+00  1.3878E+00  1.0352E+00  9.9553E-01  2.1987E+00  1.1149E+00  6.8668E-01  1.3019E+00
             1.0569E+00
 PARAMETER:  1.2710E-01 -5.7007E-01  4.2820E-01  4.2770E-01  1.3456E-01  9.5524E-02  8.8785E-01  2.0880E-01 -2.7589E-01  3.6379E-01
             1.5537E-01
 GRADIENT:  -5.4809E+04  1.2222E+04 -1.6247E+04  1.6287E+04  1.0354E+05 -1.4143E-01  1.5692E+04  3.3254E+04 -5.0514E+04 -3.8310E+04
            -8.9692E+04

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1490
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.0085E-03  1.9368E-02 -3.9409E-02 -2.9514E-02 -3.2728E-02
 SE:             2.9918E-02  1.9261E-02  1.5689E-02  2.3094E-02  2.1557E-02
 N:                     100         100         100         100         100

 P VAL.:         9.1990E-01  3.1463E-01  1.2007E-02  2.0125E-01  1.2896E-01

 ETASHRINKSD(%)  1.0000E-10  3.5474E+01  4.7441E+01  2.2633E+01  2.7780E+01
 ETASHRINKVR(%)  1.0000E-10  5.8364E+01  7.2376E+01  4.0143E+01  4.7843E+01
 EBVSHRINKSD(%)  4.0885E-01  3.8461E+01  5.2137E+01  2.0552E+01  2.2452E+01
 EBVSHRINKVR(%)  8.1604E-01  6.2129E+01  7.7092E+01  3.6881E+01  3.9863E+01
 RELATIVEINF(%)  9.8336E+01  4.4453E+00  6.3629E+00  7.4461E+00  1.7756E+01
 EPSSHRINKSD(%)  3.4071E+01
 EPSSHRINKVR(%)  5.6533E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2121.6022594128262     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1202.6637262081535     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.57
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.03
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2121.602       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  5.12E-01  1.39E+00  1.39E+00  1.04E+00  9.96E-01  2.20E+00  1.11E+00  6.87E-01  1.30E+00  1.06E+00
 


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
+        2.04E+07
 
 TH 2
+       -2.96E+01  4.09E+06
 
 TH 3
+        7.46E+00  5.02E+01  9.83E+05
 
 TH 4
+       -5.90E+00  3.85E+02  2.60E+03  9.88E+05
 
 TH 5
+       -2.64E+01  1.75E+03 -1.05E+03  9.58E+02  1.80E+07
 
 TH 6
+        5.50E+01 -2.69E+01  1.17E+01 -1.16E+01 -4.77E+01  1.98E+02
 
 TH 7
+       -1.43E+02  9.23E+01 -3.03E+01  2.33E+01  1.35E+02 -3.34E+00  9.14E+04
 
 TH 8
+       -1.90E+01 -3.28E+04  1.60E+04 -1.61E+04  2.30E+03 -2.91E+01  8.01E+01  6.39E+06
 
 TH 9
+        2.45E+01 -4.56E+03  2.22E+03 -2.27E+03 -1.32E+07  3.52E+01 -8.65E+01 -5.64E+03  9.71E+06
 
 TH10
+        1.12E+01 -2.11E+03  1.03E+03 -1.05E+03 -1.18E+03  1.49E+01 -3.91E+01 -3.15E+06  2.80E+03  1.55E+06
 
 TH11
+        1.01E+04 -4.51E+03  2.20E+03 -2.23E+03 -9.44E+03  4.24E+01 -1.09E+06 -5.61E+03  6.95E+03  2.78E+03  1.29E+07
 
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
 #CPUT: Total CPU Time in Seconds,       31.671
Stop Time:
Thu Sep 30 01:13:01 CDT 2021
