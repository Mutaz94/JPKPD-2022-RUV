Sun Oct 24 03:14:17 CDT 2021
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
$DATA ../../../../data/SD4/SL2/dat34.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m34.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1575.98013117633        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.9263E+02 -8.8317E+01 -3.5600E+01 -2.3119E+01  1.2641E+02 -2.5469E+00 -6.9070E+00 -1.2942E+01  3.7821E+01 -3.5029E+01
            -4.4488E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1594.75748326179        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      170
 NPARAMETR:  1.0007E+00  1.1025E+00  9.5558E-01  1.0095E+00  9.5523E-01  1.1356E+00  1.0230E+00  1.0895E+00  8.1905E-01  1.1081E+00
             1.1413E+00
 PARAMETER:  1.0072E-01  1.9757E-01  5.4560E-02  1.0942E-01  5.4195E-02  2.2720E-01  1.2275E-01  1.8574E-01 -9.9610E-02  2.0267E-01
             2.3215E-01
 GRADIENT:   9.1156E+01 -1.5058E+01 -1.2348E+01  1.2564E+01  1.8174E+01  1.7214E+01 -7.1483E+00 -4.0735E+00  1.0794E+01  1.9154E+00
             1.5975E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1597.78643034477        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      348
 NPARAMETR:  9.8842E-01  1.0422E+00  1.3019E+00  1.0703E+00  1.0301E+00  1.0489E+00  1.2556E+00  1.9315E+00  5.9757E-01  1.0896E+00
             1.1070E+00
 PARAMETER:  8.8353E-02  1.4131E-01  3.6385E-01  1.6793E-01  1.2966E-01  1.4775E-01  3.2762E-01  7.5830E-01 -4.1489E-01  1.8584E-01
             2.0161E-01
 GRADIENT:   8.5782E+01  1.8836E+01 -1.1167E+01  4.3616E+01 -9.1904E+00 -1.0993E+01 -3.2328E+00  1.4096E+01 -3.9602E+00 -6.2276E+00
             7.3298E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1598.82099275149        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:      537
 NPARAMETR:  9.8449E-01  1.0217E+00  1.4028E+00  1.0859E+00  1.0577E+00  1.0496E+00  1.2805E+00  1.9748E+00  5.9160E-01  1.1274E+00
             1.1073E+00
 PARAMETER:  8.4372E-02  1.2145E-01  4.3849E-01  1.8243E-01  1.5612E-01  1.4840E-01  3.4724E-01  7.8048E-01 -4.2493E-01  2.1991E-01
             2.0195E-01
 GRADIENT:   7.8560E+01  1.8859E+01 -1.0660E+01  4.3031E+01 -7.2093E+00 -1.0043E+01 -3.0280E+00  1.2646E+01 -3.9980E+00 -5.8383E+00
             7.6740E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1601.75673408325        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      719            RESET HESSIAN, TYPE II
 NPARAMETR:  9.4696E-01  1.0009E+00  1.4725E+00  1.0670E+00  1.0677E+00  1.0766E+00  1.3017E+00  1.7144E+00  6.2847E-01  1.1715E+00
             1.1050E+00
 PARAMETER:  4.5500E-02  1.0091E-01  4.8696E-01  1.6486E-01  1.6547E-01  1.7384E-01  3.6365E-01  6.3907E-01 -3.6447E-01  2.5827E-01
             1.9983E-01
 GRADIENT:   3.0944E+02  9.3524E+00  1.1260E+01  5.5219E+01 -7.4991E+00  5.9079E+01  1.3899E+01  1.0187E-01  7.7462E+00 -2.5896E+00
             3.1039E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1602.39158501057        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      876
 NPARAMETR:  9.4669E-01  1.0062E+00  1.4884E+00  1.0820E+00  1.1047E+00  1.0658E+00  1.2788E+00  1.7325E+00  6.4368E-01  1.2559E+00
             1.0972E+00
 PARAMETER:  4.5221E-02  1.0617E-01  4.9771E-01  1.7884E-01  1.9958E-01  1.6373E-01  3.4589E-01  6.4958E-01 -3.4055E-01  3.2782E-01
             1.9274E-01
 GRADIENT:   1.9182E+00 -7.4661E-01 -5.4203E+00  1.7673E+00  3.6833E+00 -8.6457E-01  3.0821E-01 -3.9185E-01  1.2132E+00  1.5480E+00
            -2.0069E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1602.96224344656        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1060             RESET HESSIAN, TYPE I
 NPARAMETR:  9.4579E-01  9.9616E-01  1.7903E+00  1.0930E+00  1.1672E+00  1.0684E+00  1.2870E+00  1.9814E+00  6.3427E-01  1.3045E+00
             1.1021E+00
 PARAMETER:  4.4267E-02  9.6153E-02  6.8241E-01  1.8893E-01  2.5464E-01  1.6620E-01  3.5234E-01  7.8380E-01 -3.5528E-01  3.6580E-01
             1.9718E-01
 GRADIENT:   3.1073E+02  2.0899E+01  3.9398E+00  1.0648E+02  1.2657E+01  5.4476E+01  1.0662E+01  1.0599E+00  8.5204E+00  3.8931E+00
             1.7436E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1602.98847036733        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1239
 NPARAMETR:  9.4509E-01  9.8846E-01  1.8021E+00  1.0977E+00  1.1679E+00  1.0667E+00  1.3073E+00  1.9814E+00  6.1746E-01  1.3097E+00
             1.1025E+00
 PARAMETER:  4.3529E-02  8.8390E-02  6.8897E-01  1.9325E-01  2.5520E-01  1.6454E-01  3.6794E-01  7.8380E-01 -3.8213E-01  3.6981E-01
             1.9761E-01
 GRADIENT:  -5.7666E-01  4.6455E-02 -1.9561E-01 -4.1740E+00  1.8297E-01 -2.1856E-01  1.7072E-01 -2.6234E+00 -6.6744E-02  9.6697E-02
             2.3736E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1603.01750837322        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1421
 NPARAMETR:  9.4401E-01  9.7435E-01  1.8189E+00  1.1072E+00  1.1668E+00  1.0657E+00  1.3235E+00  1.9815E+00  6.1422E-01  1.3107E+00
             1.1023E+00
 PARAMETER:  4.2385E-02  7.4012E-02  6.9822E-01  2.0186E-01  2.5427E-01  1.6362E-01  3.8025E-01  7.8386E-01 -3.8741E-01  3.7056E-01
             1.9738E-01
 GRADIENT:  -2.6136E+00  4.2192E-01 -1.6954E-01 -3.6596E+00  2.3868E-01 -5.5908E-01  2.9893E-01 -2.7553E+00 -5.4769E-02  2.1608E-01
             1.9040E-01

0ITERATION NO.:   41    OBJECTIVE VALUE:  -1603.01750837322        NO. OF FUNC. EVALS.:  26
 CUMULATIVE NO. OF FUNC. EVALS.:     1447
 NPARAMETR:  9.4401E-01  9.7435E-01  1.8189E+00  1.1072E+00  1.1668E+00  1.0657E+00  1.3235E+00  1.9815E+00  6.1422E-01  1.3107E+00
             1.1023E+00
 PARAMETER:  4.2385E-02  7.4012E-02  6.9822E-01  2.0186E-01  2.5427E-01  1.6362E-01  3.8025E-01  7.8386E-01 -3.8741E-01  3.7056E-01
             1.9738E-01
 GRADIENT:  -2.1868E+04 -4.0789E+03  6.2777E+03  2.1653E+04 -1.5961E+03  2.6726E+04  1.1497E+04  1.1516E+04 -1.1289E+04  1.0976E+03
            -2.2160E+04

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1447
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.3963E-03  3.4395E-03 -4.8632E-02 -1.5327E-02 -4.5882E-02
 SE:             2.9864E-02  2.2652E-02  1.6645E-02  1.8511E-02  2.1225E-02
 N:                     100         100         100         100         100

 P VAL.:         9.3605E-01  8.7931E-01  3.4812E-03  4.0766E-01  3.0643E-02

 ETASHRINKSD(%)  1.0000E-10  2.4111E+01  4.4237E+01  3.7986E+01  2.8893E+01
 ETASHRINKVR(%)  1.0000E-10  4.2409E+01  6.8905E+01  6.1543E+01  4.9438E+01
 EBVSHRINKSD(%)  4.9017E-01  2.4192E+01  5.3566E+01  3.8735E+01  2.2846E+01
 EBVSHRINKVR(%)  9.7793E-01  4.2531E+01  7.8439E+01  6.2466E+01  4.0472E+01
 RELATIVEINF(%)  9.8734E+01  2.6308E+00  5.1462E+00  1.6624E+00  2.1728E+01
 EPSSHRINKSD(%)  4.4485E+01
 EPSSHRINKVR(%)  6.9181E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1603.0175083732181     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -867.86668180947993     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.55
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1603.018       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.44E-01  9.74E-01  1.82E+00  1.11E+00  1.17E+00  1.07E+00  1.32E+00  1.98E+00  6.14E-01  1.31E+00  1.10E+00
 


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
 #CPUT: Total CPU Time in Seconds,       42.190
Stop Time:
Sun Oct 24 03:14:27 CDT 2021
