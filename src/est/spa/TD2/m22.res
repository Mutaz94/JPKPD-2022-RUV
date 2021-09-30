Wed Sep 29 18:51:15 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat22.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m22.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1673.06591171571        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.9124E+02  8.7780E+00 -2.7776E+01  8.2383E+01  6.7386E+01  8.0464E+01  5.0798E+00 -6.6260E-01  2.6813E+01  1.0274E+01
             5.0752E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1681.75278088933        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  9.6390E-01  1.0353E+00  1.0093E+00  9.7081E-01  9.7091E-01  8.1163E-01  9.8996E-01  1.0142E+00  9.1692E-01  9.3752E-01
             9.8491E-01
 PARAMETER:  6.3230E-02  1.3472E-01  1.0927E-01  7.0371E-02  7.0478E-02 -1.0872E-01  8.9908E-02  1.1412E-01  1.3264E-02  3.5479E-02
             8.4796E-02
 GRADIENT:   1.0231E+00 -4.1902E+00  7.1851E+00 -6.2512E+00  9.4036E-01 -2.0791E+01 -3.0095E+00 -4.8119E+00  2.1218E+00  6.1524E+00
            -4.2146E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1683.19185265115        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      350
 NPARAMETR:  9.6179E-01  1.0625E+00  9.1297E-01  9.5693E-01  9.2640E-01  8.2150E-01  1.0668E+00  1.1027E+00  8.7639E-01  8.0436E-01
             9.8197E-01
 PARAMETER:  6.1040E-02  1.6063E-01  8.9527E-03  5.5980E-02  2.3549E-02 -9.6624E-02  1.6471E-01  1.9780E-01 -3.1946E-02 -1.1771E-01
             8.1810E-02
 GRADIENT:  -8.6190E+00  9.0781E+00  4.0566E+00  7.0388E+00 -5.0732E+00 -1.5981E+01  5.8624E-01 -6.2560E-01  3.2038E-02 -3.3846E-01
            -4.7128E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1683.69290633653        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      528
 NPARAMETR:  9.6583E-01  1.1265E+00  7.6651E-01  9.0475E-01  8.9006E-01  8.5712E-01  1.0312E+00  9.0552E-01  8.9477E-01  7.7606E-01
             9.9265E-01
 PARAMETER:  6.5236E-02  2.1911E-01 -1.6590E-01 -9.3354E-05 -1.6467E-02 -5.4177E-02  1.3075E-01  7.5019E-04 -1.1186E-02 -1.5352E-01
             9.2619E-02
 GRADIENT:   9.5168E-01  3.5941E+00  8.2104E-01  2.9714E+00 -2.6750E+00  8.2511E-01  4.8977E-01  3.8087E-02 -3.8579E-01  2.6989E-01
             2.7407E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1683.76374817779        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      710
 NPARAMETR:  9.6626E-01  1.3024E+00  6.2721E-01  7.8578E-01  9.0659E-01  8.5633E-01  9.2385E-01  7.4997E-01  9.8586E-01  7.7725E-01
             9.9179E-01
 PARAMETER:  6.5675E-02  3.6420E-01 -3.6648E-01 -1.4108E-01  1.9311E-03 -5.5103E-02  2.0795E-02 -1.8772E-01  8.5760E-02 -1.5199E-01
             9.1754E-02
 GRADIENT:  -5.4555E-01  3.2237E+00  9.2401E-01  1.6479E+00 -1.7721E+00  4.4202E-02  4.6682E-01 -7.9422E-02 -3.2672E-01 -1.3323E-01
            -1.6673E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1683.76730567416        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      886
 NPARAMETR:  9.6660E-01  1.3709E+00  5.7847E-01  7.4072E-01  9.1586E-01  8.5729E-01  8.8757E-01  6.8579E-01  1.0282E+00  7.8200E-01
             9.9205E-01
 PARAMETER:  6.6032E-02  4.1550E-01 -4.4737E-01 -2.0013E-01  1.2112E-02 -5.3982E-02 -1.9267E-02 -2.7719E-01  1.2776E-01 -1.4590E-01
             9.2021E-02
 GRADIENT:  -5.5482E-02  5.2677E+00  1.3134E+00  2.4053E+00 -3.4399E+00  3.8302E-01  2.6782E-01 -1.3378E-02 -1.4013E-01  1.1229E-01
            -2.7426E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1683.77848359838        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1067
 NPARAMETR:  9.6694E-01  1.3972E+00  5.5725E-01  7.1973E-01  9.2271E-01  8.5691E-01  8.7214E-01  6.5672E-01  1.0484E+00  7.8398E-01
             9.9214E-01
 PARAMETER:  6.6382E-02  4.3445E-01 -4.8474E-01 -2.2888E-01  1.9558E-02 -5.4417E-02 -3.6808E-02 -3.2050E-01  1.4723E-01 -1.4337E-01
             9.2112E-02
 GRADIENT:   9.0751E-01 -2.8179E+00 -9.3563E-01  2.9334E-01  2.3852E+00  1.9149E-01 -2.3755E-02  9.7717E-02  2.9008E-02  3.5347E-03
             9.1708E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1683.78113938697        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1250             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6715E-01  1.3976E+00  5.5723E-01  7.1957E-01  9.2190E-01  8.5668E-01  8.7263E-01  6.3997E-01  1.0492E+00  7.8495E-01
             9.9209E-01
 PARAMETER:  6.6601E-02  4.3475E-01 -4.8478E-01 -2.2910E-01  1.8686E-02 -5.4686E-02 -3.6242E-02 -3.4634E-01  1.4805E-01 -1.4214E-01
             9.2063E-02
 GRADIENT:   4.1631E+02  2.9815E+02  7.2392E+00  8.2224E+01  9.1649E+00  3.4149E+01  4.7380E+00  1.8441E-01  5.5554E+00  8.2395E-01
             8.0892E-01

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1683.78113938697        NO. OF FUNC. EVALS.:  55
 CUMULATIVE NO. OF FUNC. EVALS.:     1305
 NPARAMETR:  9.6715E-01  1.3976E+00  5.5723E-01  7.1957E-01  9.2190E-01  8.5668E-01  8.7263E-01  6.3997E-01  1.0492E+00  7.8495E-01
             9.9209E-01
 PARAMETER:  6.6601E-02  4.3475E-01 -4.8478E-01 -2.2910E-01  1.8686E-02 -5.4686E-02 -3.6242E-02 -3.4634E-01  1.4805E-01 -1.4214E-01
             9.2063E-02
 GRADIENT:   1.5657E+00 -1.1596E+00  3.3401E-01 -4.4863E-01  2.2035E-02  8.6349E-02 -1.0738E-03  6.4771E-03 -1.7182E-03 -3.2325E-03
             8.9277E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1305
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5203E-04 -1.6810E-02 -1.8526E-02  1.2816E-02 -2.7646E-02
 SE:             2.9815E-02  2.3626E-02  8.2772E-03  2.3234E-02  2.1049E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9593E-01  4.7676E-01  2.5207E-02  5.8124E-01  1.8906E-01

 ETASHRINKSD(%)  1.1575E-01  2.0851E+01  7.2270E+01  2.2162E+01  2.9482E+01
 ETASHRINKVR(%)  2.3137E-01  3.7354E+01  9.2311E+01  3.9412E+01  5.0272E+01
 EBVSHRINKSD(%)  5.5164E-01  2.0767E+01  7.4625E+01  2.2861E+01  2.8558E+01
 EBVSHRINKVR(%)  1.1002E+00  3.7222E+01  9.3561E+01  4.0496E+01  4.8960E+01
 RELATIVEINF(%)  9.8834E+01  2.9162E+00  4.0785E-01  2.7203E+00  6.1936E+00
 EPSSHRINKSD(%)  4.4723E+01
 EPSSHRINKVR(%)  6.9444E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1683.7811393869683     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -948.63031282323016     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.35
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.53
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1683.781       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.67E-01  1.40E+00  5.57E-01  7.20E-01  9.22E-01  8.57E-01  8.73E-01  6.40E-01  1.05E+00  7.85E-01  9.92E-01
 


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
+        1.60E+03
 
 TH 2
+       -8.93E+00  4.85E+02
 
 TH 3
+        1.21E+01  2.26E+02  6.08E+02
 
 TH 4
+       -2.25E+01  3.83E+02 -4.33E+02  1.17E+03
 
 TH 5
+       -5.85E+00 -3.46E+02 -6.80E+02  4.43E+02  1.09E+03
 
 TH 6
+       -3.45E-01 -1.46E+00  2.67E+00 -4.25E+00 -1.37E+00  2.66E+02
 
 TH 7
+        8.57E-01  1.93E+01 -2.28E+01 -9.88E+00 -4.63E+00 -5.44E-02  1.09E+02
 
 TH 8
+        3.14E-01 -1.17E+01 -3.06E+01  9.78E+00  6.86E+00  7.61E-02  5.43E+00  4.46E+00
 
 TH 9
+        2.36E+00 -2.01E+01 -3.19E+01  5.17E+01 -5.60E+00  1.97E-01  2.11E+01  5.09E+00  7.32E+01
 
 TH10
+       -6.11E-01 -1.27E+01 -3.46E+01 -1.74E+01 -8.60E+01  4.25E-01  1.86E+01  1.09E+01  9.43E+00  9.83E+01
 
 TH11
+       -9.71E+00 -1.30E+01 -1.94E+01 -2.35E+00 -7.99E+00  2.15E+00  7.82E+00  3.50E+00  8.57E+00  1.72E+01  2.14E+02
 
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
 #CPUT: Total CPU Time in Seconds,       21.948
Stop Time:
Wed Sep 29 18:51:39 CDT 2021
