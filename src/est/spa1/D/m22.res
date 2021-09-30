Thu Sep 30 02:48:24 CDT 2021
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
$DATA ../../../../data/spa1/D/dat22.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   1458.26082655873        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.4040E+02  1.8987E+01 -8.1506E+01 -4.2145E+01  3.8806E+02 -1.2696E+03 -3.1174E+02 -4.2130E+01 -6.0423E+02 -5.8359E+02
            -4.2477E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1203.13886599957        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.0951E-01  6.0798E-01  2.0742E+00  1.5112E+00  1.4503E+00  2.7333E+00  1.5044E+00  1.0906E+00  4.5074E+00  1.4827E+00
             3.0234E+00
 PARAMETER:  5.1505E-03 -3.9762E-01  8.2958E-01  5.1289E-01  4.7179E-01  1.1055E+00  5.0840E-01  1.8671E-01  1.6057E+00  4.9385E-01
             1.2064E+00
 GRADIENT:  -1.4281E+01 -2.9771E+01 -6.3862E+01  7.2784E+01 -1.1607E+00  3.9835E+02  6.4450E+00 -5.1852E+00  2.3749E+02 -4.7996E+00
            -1.1585E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1248.85481976847        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  7.6403E-01  6.4863E-02  3.3238E+00  2.0817E+00  5.9338E+00  2.4971E+00  5.7458E+00  1.1273E+00  2.2636E+00  2.4627E+00
             3.3924E+00
 PARAMETER: -1.6915E-01 -2.6355E+00  1.3011E+00  8.3319E-01  1.8807E+00  1.0151E+00  1.8485E+00  2.1983E-01  9.1696E-01  1.0012E+00
             1.3215E+00
 GRADIENT:  -6.1028E+01  1.7351E+00  3.1485E+01  1.8260E+02 -4.4499E+00  2.9401E+02  1.2329E+00 -1.4100E+01  7.8763E+01 -1.7899E-01
            -1.4029E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1328.42331911907        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      294
 NPARAMETR:  9.1674E-01  4.1206E-02  3.2632E+00  1.6440E+00  7.9089E+00  2.6264E+00  1.5078E+01  5.1139E+00  1.5376E+00  1.9527E+00
             3.5134E+00
 PARAMETER:  1.3070E-02 -3.0892E+00  1.2827E+00  5.9714E-01  2.1680E+00  1.0656E+00  2.8132E+00  1.7320E+00  5.3023E-01  7.6922E-01
             1.3566E+00
 GRADIENT:  -3.2769E+01  7.0003E-01  2.2613E+01 -2.4319E+01 -5.9963E+00  8.0470E+01  2.5614E+00 -4.0063E+00  2.9594E+01  1.3720E+00
             1.0329E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1371.64619254919        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      472
 NPARAMETR:  8.6730E-01  1.5577E-02  7.3801E-01  1.2806E+00  1.0807E+01  1.5462E+00  1.0406E+01  4.7122E+00  8.6515E-01  8.9364E+00
             3.6395E+00
 PARAMETER: -4.2365E-02 -4.0620E+00 -2.0380E-01  3.4733E-01  2.4802E+00  5.3582E-01  2.4424E+00  1.6502E+00 -4.4850E-02  2.2901E+00
             1.3919E+00
 GRADIENT:  -2.0251E+01  1.4703E-01  2.1499E+01 -8.9671E+00  3.5941E-01 -4.4828E+01 -6.1325E-04  3.2319E+01  4.1977E+00  6.4583E+00
             3.6925E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1378.39251344083        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:      614
 NPARAMETR:  8.5409E-01  1.4298E-02  4.0801E-01  1.2011E+00  9.4813E+00  1.6185E+00  1.3877E+00  4.3742E+00  7.4241E-01  8.1302E+00
             3.6461E+00
 PARAMETER: -5.7721E-02 -4.1477E+00 -7.9646E-01  2.8323E-01  2.3493E+00  5.8149E-01  4.2765E-01  1.5757E+00 -1.9785E-01  2.1956E+00
             1.3937E+00
 GRADIENT:   3.3889E+01  3.1239E-01  1.8218E+01  6.9115E+01  3.9219E+00  4.3344E+01 -7.4229E-04  1.4031E+02  3.2898E+00  8.0982E+00
             5.4434E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1380.90795896281        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      690
 NPARAMETR:  8.4080E-01  1.7177E-02  3.3675E-01  1.1566E+00  9.5706E+00  1.6275E+00  3.3342E-02  4.3179E+00  7.0761E-01  7.5951E+00
             3.6498E+00
 PARAMETER: -7.3399E-02 -3.9642E+00 -9.8840E-01  2.4548E-01  2.3587E+00  5.8705E-01 -3.3010E+00  1.5628E+00 -2.4586E-01  2.1275E+00
             1.3947E+00
 GRADIENT:   3.8070E+01  4.2819E-01  1.8077E+01  7.9461E+01  1.3829E+00  4.6575E+01  2.7656E-06  1.3286E+02  2.6885E+00  1.0685E+01
             5.4027E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1386.07507282582        NO. OF FUNC. EVALS.: 145
 CUMULATIVE NO. OF FUNC. EVALS.:      835
 NPARAMETR:  8.2441E-01  1.3854E-02  3.2787E-01  1.1141E+00  9.6271E+00  1.6735E+00  1.0000E-02  4.2331E+00  6.4222E-01  7.0594E+00
             3.6434E+00
 PARAMETER: -9.3093E-02 -4.1792E+00 -1.0151E+00  2.0805E-01  2.3646E+00  6.1493E-01 -5.7525E+00  1.5429E+00 -3.4283E-01  2.0544E+00
             1.3929E+00
 GRADIENT:  -4.1950E-01  2.0487E-01 -5.0745E+00  1.6889E+01  1.7308E+00 -2.6500E+00  0.0000E+00  5.8627E+01 -1.1502E+00  5.9162E+00
             3.8318E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1428.22575327036        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1015
 NPARAMETR:  6.2736E-01  1.0000E-02  9.1311E-02  6.1290E-01  1.1978E+01  1.7334E+00  1.0000E-02  2.0790E+00  1.0983E-01  9.7535E-01
             3.4464E+00
 PARAMETER: -3.6623E-01 -1.0867E+01 -2.2935E+00 -3.8955E-01  2.5831E+00  6.5010E-01 -9.4537E+01  8.3191E-01 -2.1088E+00  7.5039E-02
             1.3373E+00
 GRADIENT:   1.7255E+01  0.0000E+00  4.7292E+01 -1.0814E+02 -2.0520E-01  1.7875E+01  0.0000E+00  3.0222E+01 -1.3849E-02 -1.8388E-03
            -6.9688E+00

0ITERATION NO.:   43    OBJECTIVE VALUE:  -1428.53882458938        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:     1130
 NPARAMETR:  6.2441E-01  1.0000E-02  9.0418E-02  6.1027E-01  1.1992E+01  1.7323E+00  1.0000E-02  2.0658E+00  1.1055E-01  9.5971E-01
             3.4439E+00
 PARAMETER: -3.6864E-01 -1.0929E+01 -2.3039E+00 -3.9375E-01  2.5849E+00  6.4963E-01 -9.5291E+01  8.2578E-01 -2.1235E+00  5.8847E-02
             1.3370E+00
 GRADIENT:   1.4978E+01  0.0000E+00 -6.2811E+02  3.9037E+03  6.1396E+02  1.2235E+03  0.0000E+00  9.7858E+02 -2.1931E-02 -7.8708E+03
             1.1584E+03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1130
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -9.9878E-03 -2.9924E-06 -4.9451E-02 -1.7530E-03 -2.9553E-04
 SE:             2.8385E-02  1.6276E-06  2.5281E-02  3.2343E-03  1.7320E-04
 N:                     100         100         100         100         100

 P VAL.:         7.2493E-01  6.5982E-02  5.0461E-02  5.8782E-01  8.7956E-02

 ETASHRINKSD(%)  4.9080E+00  9.9995E+01  1.5305E+01  8.9165E+01  9.9420E+01
 ETASHRINKVR(%)  9.5752E+00  1.0000E+02  2.8268E+01  9.8826E+01  9.9997E+01
 EBVSHRINKSD(%)  8.0190E-01  9.9996E+01  4.5503E+00  8.9945E+01  9.9603E+01
 EBVSHRINKVR(%)  1.5974E+00  1.0000E+02  8.8936E+00  9.8989E+01  9.9998E+01
 RELATIVEINF(%)  1.0965E+01  4.3498E-09  1.1338E+00  5.6113E-03  2.0048E-04
 EPSSHRINKSD(%)  1.8574E+01
 EPSSHRINKVR(%)  3.3697E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1428.5388245893805     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -509.60029138470782     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.82
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.41
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1428.539       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         6.26E-01  1.00E-02  9.04E-02  6.10E-01  1.20E+01  1.73E+00  1.00E-02  2.07E+00  1.08E-01  9.60E-01  3.45E+00
 


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
+        7.39E+05
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.36E+03  0.00E+00  9.26E+05
 
 TH 4
+        4.15E+02  0.00E+00 -1.24E+04  6.88E+05
 
 TH 5
+        9.78E+00  0.00E+00  7.64E+01  2.96E+01  4.15E+01
 
 TH 6
+        2.19E+02  0.00E+00 -1.66E+05  8.84E+02 -8.51E+00  3.09E+04
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        1.59E+02  0.00E+00  9.33E+02  5.76E+02 -6.31E+00  2.04E+04  0.00E+00  1.35E+04
 
 TH 9
+       -2.50E+00  0.00E+00  1.76E+02 -2.32E+02 -9.62E-01 -3.25E+01  0.00E+00 -1.45E+01  4.92E+00
 
 TH10
+        1.78E+06  0.00E+00  1.11E+03 -1.71E+06 -1.87E+00 -1.98E+02  0.00E+00 -1.15E+02  6.11E+02  4.28E+06
 
 TH11
+       -3.67E+04  0.00E+00  8.14E+02  3.53E+04 -5.31E+00 -1.42E+02  0.00E+00 -9.25E+01 -3.68E+04 -8.82E+04  1.87E+03
 
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
 #CPUT: Total CPU Time in Seconds,       25.294
Stop Time:
Thu Sep 30 02:48:50 CDT 2021
