Sat Sep 25 11:32:33 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat14.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m14.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1637.97921046582        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.1457E+01 -1.2272E+02 -7.3515E+01 -6.0358E+01  1.2970E+02  2.9576E+00 -3.8318E+00  1.1606E+01  3.3343E+01 -2.8191E+01
            -9.0305E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1652.66039674232        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0135E+00  1.1164E+00  1.1481E+00  9.3055E-01  1.1210E+00  9.9317E-01  1.0374E+00  9.1311E-01  7.8574E-01  1.2034E+00
             1.3366E+00
 PARAMETER:  1.1339E-01  2.1013E-01  2.3813E-01  2.8022E-02  2.1424E-01  9.3142E-02  1.3669E-01  9.1046E-03 -1.4113E-01  2.8517E-01
             3.9011E-01
 GRADIENT:   2.2395E+01 -1.1921E+02 -3.6638E+01 -1.0132E+02  7.4166E+01  3.1014E-01 -3.1419E+00  5.2131E+00  9.3293E-01 -1.9988E+00
             4.9428E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1659.71347958745        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0178E+00  1.4106E+00  1.1165E+00  8.1418E-01  1.2613E+00  1.0063E+00  8.6535E-01  6.5653E-01  7.7412E-01  1.4362E+00
             1.2404E+00
 PARAMETER:  1.1767E-01  4.4404E-01  2.1021E-01 -1.0558E-01  3.3211E-01  1.0626E-01 -4.4626E-02 -3.2078E-01 -1.5603E-01  4.6202E-01
             3.1547E-01
 GRADIENT:   3.4221E+01 -1.0822E+01 -2.5057E+01  1.2221E+01  5.9280E+01  4.5032E+00 -8.5572E+00  1.0025E+00 -1.0497E+01  7.3258E+00
             1.5327E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1664.62184687023        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.0031E+00  1.2640E+00  1.1105E+00  8.9265E-01  1.0823E+00  9.9651E-01  8.4668E-01  4.3909E-01  9.2508E-01  1.2094E+00
             1.1958E+00
 PARAMETER:  1.0313E-01  3.3431E-01  2.0485E-01 -1.3565E-02  1.7909E-01  9.6504E-02 -6.6436E-02 -7.2304E-01  2.2130E-02  2.9014E-01
             2.7886E-01
 GRADIENT:   2.9284E+00 -1.3868E+01 -3.0937E+00 -2.0384E+01  7.9819E+00  8.5915E-01 -4.2207E-01  5.6587E-01  9.5413E-01 -3.2607E+00
            -5.0679E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1665.13124336912        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  1.0027E+00  1.0646E+00  1.2372E+00  1.0346E+00  1.0321E+00  9.9719E-01  9.3239E-01  2.3056E-01  8.4940E-01  1.2331E+00
             1.2024E+00
 PARAMETER:  1.0273E-01  1.6260E-01  3.1281E-01  1.3399E-01  1.3164E-01  9.7181E-02  2.9998E-02 -1.3673E+00 -6.3222E-02  3.0950E-01
             2.8431E-01
 GRADIENT:   3.6189E+00 -6.1769E-01  2.5158E-01 -1.0853E+00 -3.9847E+00  1.5065E+00 -4.2947E-01  1.4795E-01  1.3129E+00  3.2029E-01
             3.3066E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1665.15829060258        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  1.0012E+00  1.0544E+00  1.2622E+00  1.0412E+00  1.0427E+00  9.9341E-01  9.5408E-01  1.8894E-01  8.3490E-01  1.2464E+00
             1.1947E+00
 PARAMETER:  1.0120E-01  1.5302E-01  3.3289E-01  1.4040E-01  1.4177E-01  9.3391E-02  5.2987E-02 -1.5663E+00 -8.0449E-02  3.2024E-01
             2.7789E-01
 GRADIENT:   4.3633E-01 -1.8905E+00 -3.9740E-01 -1.6927E+00  2.8322E-01  2.5830E-02 -4.6407E-02  8.7493E-02  6.6018E-02  2.7199E-01
             4.0735E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1665.21412193748        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      451
 NPARAMETR:  1.0014E+00  1.0695E+00  1.2654E+00  1.0331E+00  1.0487E+00  9.9364E-01  9.2895E-01  1.0000E-02  8.4498E-01  1.2549E+00
             1.1958E+00
 PARAMETER:  1.0138E-01  1.6722E-01  3.3537E-01  1.3256E-01  1.4753E-01  9.3624E-02  2.6297E-02 -4.5703E+00 -6.8439E-02  3.2703E-01
             2.7885E-01
 GRADIENT:   6.9845E-01  4.1391E-01  1.9399E-01 -2.3424E-02 -9.4745E-01  7.7056E-02 -4.3417E-01  0.0000E+00 -2.1594E-01  4.1383E-01
             9.5893E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1665.45122018277        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      631
 NPARAMETR:  1.0148E+00  1.0662E+00  1.3021E+00  1.0403E+00  1.0618E+00  1.0006E+00  9.2975E-01  1.0000E-02  8.4700E-01  1.2692E+00
             1.1983E+00
 PARAMETER:  1.1468E-01  1.6410E-01  3.6396E-01  1.3949E-01  1.5994E-01  1.0057E-01  2.7162E-02 -4.8913E+00 -6.6060E-02  3.3836E-01
             2.8088E-01
 GRADIENT:  -1.0448E-01 -8.1689E-02  2.6687E-02 -1.4781E-01 -6.4851E-02 -2.0349E-02 -9.5993E-03  0.0000E+00  4.6036E-02  2.0759E-02
             1.7598E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1665.45131825528        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      808
 NPARAMETR:  1.0148E+00  1.0580E+00  1.3055E+00  1.0456E+00  1.0598E+00  1.0006E+00  9.3593E-01  1.0000E-02  8.4323E-01  1.2685E+00
             1.1981E+00
 PARAMETER:  1.1468E-01  1.5641E-01  3.6661E-01  1.4463E-01  1.5806E-01  1.0055E-01  3.3780E-02 -4.6024E+00 -7.0518E-02  3.3786E-01
             2.8074E-01
 GRADIENT:  -3.0372E-03  1.5335E-02 -2.3230E-03  2.8491E-02 -5.5151E-05 -5.2027E-03 -3.6715E-03  0.0000E+00 -4.4430E-03  2.1599E-03
            -8.8464E-03

0ITERATION NO.:   42    OBJECTIVE VALUE:  -1665.45132014344        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:      865
 NPARAMETR:  1.0148E+00  1.0588E+00  1.3050E+00  1.0451E+00  1.0599E+00  1.0006E+00  9.3551E-01  1.0000E-02  8.4352E-01  1.2685E+00
             1.1981E+00
 PARAMETER:  1.1468E-01  1.5717E-01  3.6623E-01  1.4412E-01  1.5821E-01  1.0057E-01  3.3338E-02 -4.6277E+00 -7.0171E-02  3.3783E-01
             2.8076E-01
 GRADIENT:  -8.9601E-04 -2.9889E-04 -1.9914E-03  2.0010E-03  2.4720E-03 -1.5457E-03 -5.6997E-04  0.0000E+00 -4.2153E-04  1.3619E-03
            -2.4445E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      865
 NO. OF SIG. DIGITS IN FINAL EST.:  3.8
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -5.8954E-04 -1.1665E-02 -1.6878E-04 -2.9966E-03 -2.8853E-02
 SE:             2.9742E-02  1.7451E-02  9.0700E-05  2.3504E-02  2.4371E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8419E-01  5.0386E-01  6.2763E-02  8.9855E-01  2.3645E-01

 ETASHRINKSD(%)  3.5940E-01  4.1537E+01  9.9696E+01  2.1259E+01  1.8355E+01
 ETASHRINKVR(%)  7.1751E-01  6.5821E+01  9.9999E+01  3.7999E+01  3.3341E+01
 EBVSHRINKSD(%)  5.8481E-01  4.1247E+01  9.9692E+01  2.1509E+01  1.5063E+01
 EBVSHRINKVR(%)  1.1662E+00  6.5481E+01  9.9999E+01  3.8391E+01  2.7857E+01
 RELATIVEINF(%)  9.7484E+01  3.5789E-01  8.3639E-05  6.8196E-01  8.2411E+00
 EPSSHRINKSD(%)  4.0521E+01
 EPSSHRINKVR(%)  6.4622E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1665.4513201434415     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -930.30049357970336     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.74
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.99
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1665.451       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.06E+00  1.31E+00  1.05E+00  1.06E+00  1.00E+00  9.36E-01  1.00E-02  8.44E-01  1.27E+00  1.20E+00
 


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
+        1.07E+03
 
 TH 2
+       -1.54E+01  4.02E+02
 
 TH 3
+        6.25E+00  4.15E+01  6.67E+01
 
 TH 4
+       -1.66E+01  5.25E+02 -4.91E+01  8.80E+02
 
 TH 5
+       -3.68E+00 -1.45E+02 -1.44E+02  5.31E+01  4.45E+02
 
 TH 6
+       -3.46E+00 -1.98E+00  1.70E+00 -5.65E+00 -1.86E+00  1.94E+02
 
 TH 7
+       -1.26E+00  5.67E+00  7.82E+00 -9.81E+00 -9.13E+00  7.43E-03  3.27E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.64E+00 -1.63E+01 -4.71E+00  1.31E+01  5.06E+00 -1.49E+00  4.07E+01  0.00E+00  1.12E+02
 
 TH10
+       -7.92E-01 -3.32E-01 -1.02E+01 -8.57E+00 -4.06E+01  7.38E-01  4.27E+00  0.00E+00  5.32E-01  6.18E+01
 
 TH11
+       -9.37E+00 -2.00E+01 -1.79E+01 -6.21E+00  8.76E+00  1.61E-01  3.03E+00  0.00E+00  1.16E+01  1.78E+01  1.66E+02
 
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
 #CPUT: Total CPU Time in Seconds,       14.792
Stop Time:
Sat Sep 25 11:32:50 CDT 2021
