Sat Sep 25 02:25:47 CDT 2021
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
$DATA ../../../../data/int/SL3/dat56.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      970
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      870
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
 RAW OUTPUT FILE (FILE): m56.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -590.426908623787        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.1843E+01 -1.0799E+02  1.8826E+02  8.1291E+01  1.0434E+02 -1.6979E+01 -1.3137E+02 -2.6209E+02 -1.1408E+02 -2.5785E+01
            -5.8109E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2682.59850754285        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0273E+00  1.5704E+00  9.0328E-01  7.6197E-01  1.2648E+00  1.1060E+00  1.2384E+00  9.4583E-01  9.4359E-01  9.4905E-01
             2.6285E+00
 PARAMETER:  1.2695E-01  5.5131E-01 -1.7180E-03 -1.7185E-01  3.3493E-01  2.0076E-01  3.1381E-01  4.4304E-02  4.1937E-02  4.7709E-02
             1.0664E+00
 GRADIENT:   1.0133E+01  5.7471E+01 -7.3009E+00  3.1747E+01  5.7197E+00  2.1917E+01  3.8491E+01  2.0926E+00 -1.0848E+01 -3.8458E+01
            -1.5313E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2689.80082648274        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0514E+00  1.7265E+00  7.1230E-01  6.9100E-01  1.2975E+00  1.0310E+00  1.0640E+00  3.7838E-01  7.5760E-01  1.3871E+00
             2.6641E+00
 PARAMETER:  1.5016E-01  6.4612E-01 -2.3926E-01 -2.6962E-01  3.6042E-01  1.3048E-01  1.6201E-01 -8.7185E-01 -1.7760E-01  4.2719E-01
             1.0799E+00
 GRADIENT:   5.4974E+01  1.0159E+02 -1.2760E+01  9.1357E+01 -1.2275E+01 -4.7987E+00  1.1651E+01  2.7297E-01 -9.5805E+00  8.7359E+00
             2.2005E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2708.52549912938        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0116E+00  2.0060E+00  7.5318E-01  4.4742E-01  1.6504E+00  1.0415E+00  7.6062E-01  1.0505E-01  1.5422E+00  1.5139E+00
             2.5572E+00
 PARAMETER:  1.1153E-01  7.9614E-01 -1.8345E-01 -7.0426E-01  6.0103E-01  1.4065E-01 -1.7362E-01 -2.1534E+00  5.3319E-01  5.1467E-01
             1.0389E+00
 GRADIENT:  -1.8938E+01  1.3190E+01 -9.9102E-01  2.2348E+01  9.6614E+00  5.2146E-01 -1.0919E+01 -8.0943E-03  1.5309E+00  2.2690E+00
            -2.6990E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2717.96717414591        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  1.0200E+00  2.3375E+00  3.1772E-01  2.0516E-01  1.8205E+00  1.0437E+00  7.4537E-01  1.0000E-02  2.1317E+00  1.6685E+00
             2.5566E+00
 PARAMETER:  1.1982E-01  9.4907E-01 -1.0466E+00 -1.4839E+00  6.9911E-01  1.4276E-01 -1.9388E-01 -4.6369E+00  8.5691E-01  6.1191E-01
             1.0387E+00
 GRADIENT:  -1.5527E+00  2.7952E+00  3.8420E-01  1.4847E+00 -2.1906E+00  1.5339E+00  9.0170E-01  0.0000E+00 -1.4640E+00  4.2520E+00
             8.1991E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2719.72084977804        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  1.0206E+00  2.4222E+00  2.2001E-01  1.4561E-01  1.8894E+00  1.0395E+00  7.2782E-01  1.0000E-02  2.5692E+00  1.6853E+00
             2.5507E+00
 PARAMETER:  1.2042E-01  9.8469E-01 -1.4141E+00 -1.8268E+00  7.3625E-01  1.3871E-01 -2.1770E-01 -6.0318E+00  1.0436E+00  6.2197E-01
             1.0364E+00
 GRADIENT:  -1.7708E-01  2.0755E-01 -6.5836E-01  9.6516E-01  1.7930E+00 -9.8205E-03 -5.1803E-02  0.0000E+00 -2.6320E-01  2.4669E-01
             2.0708E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2721.05202857664        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      536
 NPARAMETR:  1.0262E+00  2.5488E+00  1.4488E-01  8.7249E-02  1.9752E+00  1.0436E+00  7.1454E-01  1.0000E-02  3.2385E+00  1.7497E+00
             2.5541E+00
 PARAMETER:  1.2583E-01  1.0356E+00 -1.8318E+00 -2.3390E+00  7.8069E-01  1.4269E-01 -2.3612E-01 -7.9400E+00  1.2751E+00  6.5946E-01
             1.0377E+00
 GRADIENT:   1.3690E+00  1.5960E+01  2.0872E-01  5.4429E-01  1.2277E+00 -3.5524E-02  2.1477E-01  0.0000E+00 -2.1385E+00  1.3712E+00
             1.7770E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2721.10503858057        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:      721
 NPARAMETR:  1.0254E+00  2.5404E+00  1.4379E-01  8.5863E-02  1.9719E+00  1.0436E+00  7.1377E-01  1.0000E-02  3.2616E+00  1.7411E+00
             2.5521E+00
 PARAMETER:  1.2508E-01  1.0323E+00 -1.8394E+00 -2.3550E+00  7.7898E-01  1.4266E-01 -2.3719E-01 -7.9885E+00  1.2822E+00  6.5454E-01
             1.0369E+00
 GRADIENT:   1.5008E-01  9.1489E-01  1.9190E-01 -2.9769E-01  4.2648E-01  2.0390E-02 -2.4062E-02  0.0000E+00 -2.0216E+00  7.4294E-02
             4.5503E-01

0ITERATION NO.:   36    OBJECTIVE VALUE:  -2721.10503858057        NO. OF FUNC. EVALS.:  29
 CUMULATIVE NO. OF FUNC. EVALS.:      750
 NPARAMETR:  1.0253E+00  2.5399E+00  1.4389E-01  8.5787E-02  1.9703E+00  1.0435E+00  7.1382E-01  1.0000E-02  3.2600E+00  1.7408E+00
             2.5531E+00
 PARAMETER:  1.2508E-01  1.0323E+00 -1.8394E+00 -2.3550E+00  7.7898E-01  1.4266E-01 -2.3719E-01 -7.9885E+00  1.2822E+00  6.5454E-01
             1.0369E+00
 GRADIENT:   8.3187E-02  9.8543E-01 -2.2063E+04  1.7236E+04  3.1179E-01  1.9473E-02 -2.5049E-02  0.0000E+00  3.1638E+04  7.3484E-02
            -3.9214E+04

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      750
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1453E-03 -1.6046E-02 -1.1161E-05  2.4330E-02 -2.0816E-02
 SE:             2.9488E-02  2.7693E-02  1.0579E-05  1.3141E-02  2.6104E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6902E-01  5.6229E-01  2.9142E-01  6.4101E-02  4.2520E-01

 ETASHRINKSD(%)  1.2116E+00  7.2259E+00  9.9965E+01  5.5976E+01  1.2550E+01
 ETASHRINKVR(%)  2.4085E+00  1.3930E+01  1.0000E+02  8.0619E+01  2.3525E+01
 EBVSHRINKSD(%)  1.3648E+00  7.0959E+00  9.9946E+01  6.7135E+01  9.6205E+00
 EBVSHRINKVR(%)  2.7109E+00  1.3688E+01  1.0000E+02  8.9199E+01  1.8315E+01
 RELATIVEINF(%)  9.7243E+01  3.2496E+01  1.7167E-05  3.2375E+00  6.0428E+01
 EPSSHRINKSD(%)  1.6821E+01
 EPSSHRINKVR(%)  3.0812E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          870
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1598.9530477761305     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2721.1050385805652     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1122.1519908044347     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.91
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.84
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2721.105       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  2.54E+00  1.44E-01  8.59E-02  1.97E+00  1.04E+00  7.14E-01  1.00E-02  3.26E+00  1.74E+00  2.55E+00
 


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
+        9.47E+02
 
 TH 2
+        3.02E+07  1.48E+06
 
 TH 3
+       -1.27E+04 -4.53E+02  1.45E+08
 
 TH 4
+        1.67E+04  9.61E+02 -1.90E+08  2.48E+08
 
 TH 5
+       -1.95E+00 -1.30E+01 -8.72E+02  1.20E+03  7.75E+01
 
 TH 6
+        4.29E+00 -3.30E+00 -1.01E+04  1.32E+04 -7.68E-01  1.76E+02
 
 TH 7
+        2.46E+00  2.54E+00  4.52E+02 -6.88E+02  4.20E-01 -2.32E+00  3.06E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.10E+02  2.77E+01 -9.16E+06  1.20E+07  5.32E+01  6.37E+02 -3.07E+01  0.00E+00  5.79E+05
 
 TH10
+        5.11E-01 -8.94E-01 -5.99E+02  8.12E+02 -6.52E+00  1.78E-01  1.07E+00  0.00E+00  3.81E+01  4.17E+01
 
 TH11
+       -1.29E+03  1.46E+06  1.45E+07 -1.90E+07 -8.44E+01 -1.00E+03  5.70E+01  0.00E+00 -9.18E+05 -5.52E+01  1.45E+06
 
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
 #CPUT: Total CPU Time in Seconds,       28.869
Stop Time:
Sat Sep 25 02:26:18 CDT 2021
