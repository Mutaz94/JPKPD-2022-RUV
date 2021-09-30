Thu Sep 30 01:24:48 CDT 2021
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
$DATA ../../../../data/spa1/TD1/dat44.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m44.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2022.11802524946        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.7000E+02  6.6117E+01  1.0570E+01  1.2580E+02  2.8709E+01  2.4637E+01  5.4356E+00 -1.7869E+01  2.2781E+01 -3.7942E+00
            -3.9889E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2031.24262369852        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  9.3276E-01  9.8632E-01  9.5478E-01  9.8769E-01  9.6243E-01  1.0131E+00  9.9018E-01  1.0600E+00  9.4503E-01  1.0050E+00
             1.0742E+00
 PARAMETER:  3.0390E-02  8.6228E-02  5.3726E-02  8.7612E-02  6.1708E-02  1.1297E-01  9.0132E-02  1.5823E-01  4.3463E-02  1.0498E-01
             1.7153E-01
 GRADIENT:   1.5585E+00  1.2838E+01  3.9021E+00  1.8119E+01  1.4083E+00 -7.5079E+00  4.1742E-01 -9.3629E+00  7.6392E+00  5.9500E+00
             2.1570E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2033.24558226548        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      329
 NPARAMETR:  9.3458E-01  9.5133E-01  8.6540E-01  9.9677E-01  8.9425E-01  1.0226E+00  1.1264E+00  1.2000E+00  8.4552E-01  8.5772E-01
             1.0360E+00
 PARAMETER:  3.2337E-02  5.0109E-02 -4.4564E-02  9.6766E-02 -1.1769E-02  1.2231E-01  2.1906E-01  2.8230E-01 -6.7801E-02 -5.3475E-02
             1.3538E-01
 GRADIENT:   5.2580E+00  4.4880E+00 -4.2530E+00  9.4636E+00  4.0218E+00 -4.2103E+00  2.5745E-01  2.3040E+00 -1.1755E+00  1.2005E+00
            -4.9109E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2033.53253989647        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      511
 NPARAMETR:  9.2850E-01  7.9217E-01  1.0739E+00  1.1044E+00  9.2035E-01  1.0362E+00  1.2694E+00  1.2612E+00  7.9775E-01  9.0416E-01
             1.0501E+00
 PARAMETER:  2.5818E-02 -1.3298E-01  1.7130E-01  1.9927E-01  1.6998E-02  1.3551E-01  3.3855E-01  3.3206E-01 -1.2595E-01 -7.4347E-04
             1.4886E-01
 GRADIENT:  -3.3107E+00  6.5227E+00  4.9735E+00  2.9404E+00 -5.7927E+00  2.3255E+00  1.9037E-01 -1.6817E+00 -7.3357E-01 -8.0450E-01
             3.7641E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2034.21406844580        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      691
 NPARAMETR:  9.2715E-01  4.5324E-01  1.5184E+00  1.3454E+00  9.6403E-01  1.0200E+00  1.6295E+00  1.5776E+00  7.5544E-01  9.9839E-01
             1.0369E+00
 PARAMETER:  2.4364E-02 -6.9132E-01  5.1768E-01  3.9668E-01  6.3366E-02  1.1982E-01  5.8828E-01  5.5590E-01 -1.8046E-01  9.8385E-02
             1.3627E-01
 GRADIENT:   4.6513E+00  1.0963E+01  3.6201E+00  2.9381E+01 -9.0892E+00 -1.9451E+00  2.5634E+00  2.7011E-01  3.5516E+00  5.7122E-01
            -7.0760E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2034.43179393161        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      874
 NPARAMETR:  9.2426E-01  4.2478E-01  1.5132E+00  1.3383E+00  9.7193E-01  1.0312E+00  1.6577E+00  1.5789E+00  7.3325E-01  9.8905E-01
             1.0542E+00
 PARAMETER:  2.1240E-02 -7.5619E-01  5.1425E-01  3.9142E-01  7.1524E-02  1.3074E-01  6.0543E-01  5.5674E-01 -2.1027E-01  8.8993E-02
             1.5281E-01
 GRADIENT:  -1.2290E+00 -3.7776E+00 -4.1924E+00 -2.8617E+01  1.2585E+01  2.6954E+00  1.0527E+00  1.0644E+00 -2.6620E+00 -1.9250E+00
             6.1833E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2034.79151841515        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1051
 NPARAMETR:  9.2538E-01  4.1553E-01  1.5161E+00  1.3482E+00  9.6443E-01  1.0251E+00  1.4503E+00  1.5690E+00  7.4101E-01  9.9756E-01
             1.0459E+00
 PARAMETER:  2.2444E-02 -7.7819E-01  5.1613E-01  3.9876E-01  6.3780E-02  1.2481E-01  4.7178E-01  5.5044E-01 -1.9975E-01  9.7553E-02
             1.4488E-01
 GRADIENT:   1.6081E+00 -7.4800E-01 -1.6876E+00 -1.8142E+01  5.1953E+00  3.7846E-01  6.8225E-01  2.0555E-01 -6.0321E+00 -5.7244E-01
            -3.3014E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2035.88564852197        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1227
 NPARAMETR:  9.2422E-01  3.5556E-01  1.5473E+00  1.3871E+00  9.5374E-01  1.0241E+00  3.6767E-01  1.5911E+00  7.8314E-01  1.0074E+00
             1.0454E+00
 PARAMETER:  2.1198E-02 -9.3407E-01  5.3653E-01  4.2720E-01  5.2640E-02  1.2382E-01 -9.0057E-01  5.6444E-01 -1.4444E-01  1.0734E-01
             1.4443E-01
 GRADIENT:   7.3652E-01 -9.0731E-01 -9.4051E-01 -5.5177E+00 -3.3194E-01  1.3668E-01  1.8607E-01 -1.3290E-01  2.1842E+00  8.9343E-02
             3.1307E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2035.92795489294        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1407
 NPARAMETR:  9.2473E-01  3.3645E-01  1.5735E+00  1.4032E+00  9.5643E-01  1.0246E+00  7.3151E-02  1.6120E+00  7.7538E-01  1.0086E+00
             1.0454E+00
 PARAMETER:  2.1745E-02 -9.8931E-01  5.5333E-01  4.3878E-01  5.5447E-02  1.2429E-01 -2.5152E+00  5.7745E-01 -1.5441E-01  1.0855E-01
             1.4440E-01
 GRADIENT:   2.4230E+00  3.1346E-01 -1.1923E+00  2.6803E+00  1.6375E-01  3.7110E-01  7.8068E-03 -2.5216E-01  9.2335E-01 -2.0371E-01
             5.1696E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2035.93578106677        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     1576
 NPARAMETR:  9.2441E-01  3.3501E-01  1.5808E+00  1.4033E+00  9.5777E-01  1.0243E+00  5.8593E-02  1.6175E+00  7.7392E-01  1.0103E+00
             1.0454E+00
 PARAMETER:  2.1382E-02 -9.9363E-01  5.5784E-01  4.3879E-01  5.6836E-02  1.2398E-01 -2.7100E+00  5.8085E-01 -1.5623E-01  1.1026E-01
             1.4441E-01
 GRADIENT:  -8.1197E+04 -2.0170E-01 -7.7831E-01 -1.8508E+04 -8.1180E+04 -3.2526E-02  4.2807E-03 -2.5303E-01  3.7504E-01  7.3614E+04
             5.6213E+04
 NUMSIGDIG:         2.3         2.7         2.4         2.3         2.3         3.2         0.4         2.5         1.9         2.3
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1576
 NO. OF SIG. DIGITS IN FINAL EST.:  0.4

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.2834E-04 -8.3834E-04 -4.0073E-02 -8.0169E-03 -4.4136E-02
 SE:             2.9870E-02  3.8599E-04  1.9624E-02  2.9340E-02  2.0243E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9123E-01  2.9861E-02  4.1146E-02  7.8467E-01  2.9233E-02

 ETASHRINKSD(%)  1.0000E-10  9.8707E+01  3.4257E+01  1.7061E+00  3.2184E+01
 ETASHRINKVR(%)  1.0000E-10  9.9983E+01  5.6779E+01  3.3830E+00  5.4011E+01
 EBVSHRINKSD(%)  3.5503E-01  9.8739E+01  3.7901E+01  2.0840E+00  2.8886E+01
 EBVSHRINKVR(%)  7.0879E-01  9.9984E+01  6.1437E+01  4.1245E+00  4.9428E+01
 RELATIVEINF(%)  9.7973E+01  9.6721E-04  1.1510E+01  6.6029E+00  1.2424E+01
 EPSSHRINKSD(%)  3.4483E+01
 EPSSHRINKVR(%)  5.7076E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2035.9357810667727     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1116.9972478621000     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.53
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.90
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2035.936       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.24E-01  3.35E-01  1.58E+00  1.40E+00  9.58E-01  1.02E+00  6.02E-02  1.62E+00  7.74E-01  1.01E+00  1.05E+00
 


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
+        2.38E+07
 
 TH 2
+       -2.60E+01  4.62E+02
 
 TH 3
+       -1.88E+00  3.53E+01  2.60E+05
 
 TH 4
+        3.57E+06  9.89E+05 -2.25E+01  9.04E+02
 
 TH 5
+       -9.92E+03 -2.04E+02 -1.53E+02  3.44E+06  2.21E+07
 
 TH 6
+       -1.53E+01 -3.84E+00 -2.38E-01 -4.60E+00 -1.56E+01  1.26E+07
 
 TH 7
+        1.74E+02 -1.02E-01  6.00E-04  2.59E+01  1.68E+02  5.19E-03  1.20E+00
 
 TH 8
+       -4.76E-01 -3.02E+00  7.25E+02 -4.54E+00 -7.10E+00 -1.63E-01  1.42E-02  2.45E+01
 
 TH 9
+        4.57E+03 -1.09E+02  5.01E+00  6.81E+02  4.41E+03 -1.32E+07  5.47E-01 -2.46E-01  1.39E+07
 
 TH10
+        1.48E+04  7.59E+00 -2.76E+00 -2.96E+06 -1.90E+07  1.38E+01 -1.44E+02  1.10E+01 -3.79E+03  1.63E+07
 
 TH11
+        1.45E+07  4.05E+06 -4.13E+00 -2.18E+06  6.06E+03  1.16E+01 -1.06E+02  4.73E+00 -2.79E+03 -9.04E+03  3.66E+02
 
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
 #CPUT: Total CPU Time in Seconds,       31.505
Stop Time:
Thu Sep 30 01:25:21 CDT 2021
