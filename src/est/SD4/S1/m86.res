Sun Oct 24 02:51:37 CDT 2021
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
$DATA ../../../../data/SD4/S1/dat86.csv ignore=@
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

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       24 OCT 2021
Days until program expires : 175
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

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

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
 RAW OUTPUT FILE (FILE): m86.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1706.30010137939        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.7563E+02 -3.4758E+01 -4.1961E+01  2.5506E+01  6.0238E+01  5.3737E+01  6.6143E-01  1.0703E+01  1.6140E+01  1.2023E+01
             4.4407E-01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1712.40501922426        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      196
 NPARAMETR:  1.0362E+00  1.0747E+00  1.0740E+00  1.0067E+00  1.0125E+00  9.8674E-01  1.0018E+00  9.5438E-01  9.5961E-01  9.3760E-01
             1.0131E+00
 PARAMETER:  1.3551E-01  1.7200E-01  1.7135E-01  1.0663E-01  1.1242E-01  8.6650E-02  1.0180E-01  5.3305E-02  5.8767E-02  3.5563E-02
             1.1301E-01
 GRADIENT:   5.3320E+00  4.2978E+00  3.7381E-01  9.4593E+00  3.2113E+00 -9.0061E-01 -1.4085E+00  2.7503E+00 -1.2999E+00 -9.2281E-01
            -1.3412E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1713.63077140485        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  1.0423E+00  1.1445E+00  8.4183E-01  9.5443E-01  9.3851E-01  9.9915E-01  1.0235E+00  5.3346E-01  9.9595E-01  9.0031E-01
             1.0108E+00
 PARAMETER:  1.4140E-01  2.3495E-01 -7.2182E-02  5.3354E-02  3.6541E-02  9.9152E-02  1.2322E-01 -5.2838E-01  9.5941E-02 -5.0161E-03
             1.1079E-01
 GRADIENT:   1.4425E+01  6.6061E+00 -6.8716E+00  1.8073E+01 -1.0691E+00  3.1214E+00  1.5454E+00  2.1607E+00  4.5573E+00  4.0936E+00
             1.1335E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1714.74141900777        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      551
 NPARAMETR:  1.0350E+00  1.3286E+00  7.4469E-01  8.2248E-01  9.8766E-01  9.9175E-01  9.0244E-01  2.1596E-01  1.0845E+00  8.9816E-01
             1.0112E+00
 PARAMETER:  1.3438E-01  3.8415E-01 -1.9478E-01 -9.5426E-02  8.7585E-02  9.1714E-02 -2.6583E-03 -1.4326E+00  1.8113E-01 -7.4057E-03
             1.1111E-01
 GRADIENT:  -3.1085E+00 -4.2330E+00 -1.3311E+00 -3.2804E-01  3.4709E+00 -2.2182E-01 -8.8875E-01  2.6465E-01 -9.7339E-01 -7.6779E-01
            -1.7155E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1714.84824596077        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      727
 NPARAMETR:  1.0376E+00  1.4832E+00  6.7916E-01  7.2858E-01  1.0338E+00  9.9382E-01  8.3288E-01  9.0416E-02  1.2004E+00  9.1901E-01
             1.0123E+00
 PARAMETER:  1.3691E-01  4.9418E-01 -2.8690E-01 -2.1666E-01  1.3325E-01  9.3796E-02 -8.2869E-02 -2.3033E+00  2.8267E-01  1.5537E-02
             1.1219E-01
 GRADIENT:   1.7881E+00  5.9260E+00  2.3125E-01  4.9464E+00 -1.7544E+00  3.1288E-01  4.2568E-02  4.4661E-02  2.9144E-01  1.0435E-01
            -1.3551E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1714.88091789031        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      909
 NPARAMETR:  1.0371E+00  1.4879E+00  6.7538E-01  7.2050E-01  1.0374E+00  9.9332E-01  8.2954E-01  4.3197E-02  1.2060E+00  9.1948E-01
             1.0125E+00
 PARAMETER:  1.3646E-01  4.9739E-01 -2.9248E-01 -2.2782E-01  1.3674E-01  9.3294E-02 -8.6886E-02 -3.0420E+00  2.8734E-01  1.6048E-02
             1.1242E-01
 GRADIENT:   8.2615E-01 -1.6898E+00 -1.8802E-03 -2.4304E-01 -2.2758E-03  1.2867E-01 -2.0014E-02  1.0317E-02 -3.0778E-02  4.9330E-03
             1.9962E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1714.88623176452        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1086
 NPARAMETR:  1.0382E+00  1.4888E+00  6.7518E-01  7.2027E-01  1.0377E+00  9.9339E-01  8.2947E-01  1.3108E-02  1.2070E+00  9.1962E-01
             1.0125E+00
 PARAMETER:  1.3751E-01  4.9797E-01 -2.9277E-01 -2.2813E-01  1.3697E-01  9.3367E-02 -8.6973E-02 -4.2345E+00  2.8815E-01  1.6209E-02
             1.1239E-01
 GRADIENT:   3.2086E+00 -1.1809E+00  5.7522E-02  8.2305E-02 -1.0585E-01  1.5644E-01  2.3301E-02  9.6226E-04  3.8955E-02 -1.2971E-02
            -1.7808E-02

0ITERATION NO.:   33    OBJECTIVE VALUE:  -1714.88680949814        NO. OF FUNC. EVALS.: 102
 CUMULATIVE NO. OF FUNC. EVALS.:     1188
 NPARAMETR:  1.0382E+00  1.4886E+00  6.7675E-01  7.2145E-01  1.0363E+00  9.9331E-01  8.2987E-01  1.0000E-02  1.2060E+00  9.1900E-01
             1.0121E+00
 PARAMETER:  1.3755E-01  4.9718E-01 -2.9295E-01 -2.2844E-01  1.3703E-01  9.3324E-02 -8.7100E-02 -5.1923E+00  2.8816E-01  1.6301E-02
             1.1243E-01
 GRADIENT:   1.6433E-02 -3.5166E-01 -2.4981E-01 -5.2764E-01  6.1233E-01  3.9761E-03 -2.3523E-02  0.0000E+00  4.5758E-02  3.2772E-02
             5.4265E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1188
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -5.0443E-05 -2.3456E-02 -3.1858E-04  1.6124E-02 -2.9440E-02
 SE:             2.9834E-02  2.2145E-02  1.2384E-04  2.3714E-02  2.2770E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9865E-01  2.8951E-01  1.0096E-02  4.9654E-01  1.9604E-01

 ETASHRINKSD(%)  5.2265E-02  2.5813E+01  9.9585E+01  2.0557E+01  2.3717E+01
 ETASHRINKVR(%)  1.0450E-01  4.4963E+01  9.9998E+01  3.6887E+01  4.1809E+01
 EBVSHRINKSD(%)  4.4704E-01  2.5390E+01  9.9616E+01  2.1440E+01  2.2688E+01
 EBVSHRINKVR(%)  8.9208E-01  4.4334E+01  9.9999E+01  3.8283E+01  4.0228E+01
 RELATIVEINF(%)  9.8975E+01  2.8166E+00  1.5751E-04  3.5216E+00  9.2186E+00
 EPSSHRINKSD(%)  4.3814E+01
 EPSSHRINKVR(%)  6.8432E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1714.8868094981369     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -979.73598293439875     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     5.32
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1714.887       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.49E+00  6.75E-01  7.20E-01  1.04E+00  9.93E-01  8.29E-01  1.00E-02  1.21E+00  9.20E-01  1.01E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.00
 #CPUT: Total CPU Time in Seconds,       34.134
Stop Time:
Sun Oct 24 02:51:45 CDT 2021
