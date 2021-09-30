Wed Sep 29 14:04:01 CDT 2021
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
$DATA ../../../../data/spa/S1/dat15.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m15.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1698.41460764132        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.3004E+02 -7.4990E+01 -4.7201E+01 -3.2721E+01  1.1903E+02  6.9592E+01 -7.1686E+00 -2.3433E+00 -2.4851E+01  1.7368E+00
             3.1466E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1713.81157584627        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.0422E+00  1.0465E+00  1.0198E+00  1.0453E+00  9.2757E-01  8.4143E-01  1.0214E+00  1.0379E+00  1.1455E+00  9.3417E-01
             8.7841E-01
 PARAMETER:  1.4130E-01  1.4544E-01  1.1961E-01  1.4426E-01  2.4808E-02 -7.2650E-02  1.2113E-01  1.3717E-01  2.3581E-01  3.1901E-02
            -2.9646E-02
 GRADIENT:  -4.8615E+00 -3.8921E+00  5.8718E+00 -6.0697E-01 -7.0828E+00 -2.7346E+01  1.5805E+00 -5.8441E+00  7.7662E+00  4.3749E+00
            -2.1563E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1715.36456362296        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      368
 NPARAMETR:  1.0469E+00  1.3219E+00  9.3111E-01  8.9017E-01  9.7162E-01  8.6434E-01  9.5567E-01  1.4080E+00  1.3011E+00  8.5023E-01
             8.9724E-01
 PARAMETER:  1.4583E-01  3.7905E-01  2.8618E-02 -1.6339E-02  7.1213E-02 -4.5783E-02  5.4654E-02  4.4215E-01  3.6318E-01 -6.2252E-02
            -8.4296E-03
 GRADIENT:   5.4077E+00  2.9219E+01  1.1641E+01  1.7695E+01 -4.0245E+01 -1.5938E+01  7.9207E+00  2.1212E+00  1.0426E+01 -1.0904E+00
            -1.2230E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1717.72291815317        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      545
 NPARAMETR:  1.0460E+00  1.4601E+00  9.3089E-01  7.9043E-01  1.0628E+00  9.1207E-01  7.8122E-01  1.5029E+00  1.3602E+00  9.3924E-01
             9.2844E-01
 PARAMETER:  1.4498E-01  4.7851E-01  2.8391E-02 -1.3518E-01  1.6090E-01  7.9608E-03 -1.4690E-01  5.0741E-01  4.0761E-01  3.7316E-02
             2.5753E-02
 GRADIENT:   1.6861E+00  1.4048E+01  4.9814E+00  8.7513E+00 -4.1409E+00  5.8564E+00 -2.1603E+00 -2.0091E+00 -5.2749E+00 -1.6234E+00
             1.2225E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1718.57618212280        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      722
 NPARAMETR:  1.0462E+00  1.7777E+00  6.2525E-01  5.8399E-01  1.0923E+00  8.9466E-01  7.5731E-01  1.3047E+00  1.6818E+00  9.5782E-01
             9.2213E-01
 PARAMETER:  1.4520E-01  6.7530E-01 -3.6960E-01 -4.3788E-01  1.8826E-01 -1.1316E-02 -1.7798E-01  3.6601E-01  6.1987E-01  5.6900E-02
             1.8936E-02
 GRADIENT:  -7.3748E-02  3.4845E+01  7.0071E+00  1.3376E+01 -2.1080E+01 -2.1676E+00 -5.4315E-01 -1.6792E+00  9.8369E-01  1.9762E+00
            -1.5835E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1719.16259867054        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      902
 NPARAMETR:  1.0456E+00  1.9488E+00  4.7521E-01  4.5786E-01  1.1388E+00  8.8926E-01  7.3571E-01  1.2356E+00  1.9624E+00  9.7539E-01
             9.2225E-01
 PARAMETER:  1.4462E-01  7.6720E-01 -6.4399E-01 -6.8120E-01  2.2993E-01 -1.7370E-02 -2.0692E-01  3.1158E-01  7.7416E-01  7.5079E-02
             1.9059E-02
 GRADIENT:  -1.5738E+00  1.1659E+01  1.5662E+00  9.1449E+00 -9.9829E+00 -4.6916E+00  9.9933E-01 -4.9479E-01  4.1753E+00  2.3562E+00
            -1.8433E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1719.33009529155        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:     1043
 NPARAMETR:  1.0443E+00  1.9409E+00  4.7681E-01  4.5623E-01  1.1483E+00  9.0262E-01  7.2857E-01  1.2336E+00  1.9545E+00  9.5698E-01
             9.2928E-01
 PARAMETER:  1.4331E-01  7.6313E-01 -6.4063E-01 -6.8476E-01  2.3827E-01 -2.4550E-03 -2.1667E-01  3.0996E-01  7.7012E-01  5.6030E-02
             2.6657E-02
 GRADIENT:  -5.1824E+00 -4.5844E+00 -2.0928E-01  7.1762E+00  5.9720E+00  1.2431E+00 -1.1998E+00 -6.9079E-01  2.8882E+00 -1.7618E+00
             7.8630E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1719.41992806391        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:     1181
 NPARAMETR:  1.0434E+00  1.9308E+00  4.7682E-01  4.5361E-01  1.1434E+00  8.9845E-01  7.3164E-01  1.2409E+00  1.9441E+00  9.6448E-01
             9.2752E-01
 PARAMETER:  1.4252E-01  7.5795E-01 -6.4062E-01 -6.9052E-01  2.3403E-01 -7.0881E-03 -2.1247E-01  3.1587E-01  7.6480E-01  6.3830E-02
             2.4755E-02
 GRADIENT:   7.1906E+02  1.0905E+03  4.1358E+00  1.1161E+02  1.9668E+01  3.6425E+01  1.6963E+01 -1.3608E-01  4.2619E+01  1.2063E+00
             1.3629E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1719.45570551373        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1367            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0474E+00  1.9283E+00  4.7763E-01  4.5259E-01  1.1428E+00  8.9968E-01  7.3208E-01  1.2511E+00  1.9272E+00  9.6191E-01
             9.2650E-01
 PARAMETER:  1.4628E-01  7.5664E-01 -6.3893E-01 -6.9276E-01  2.3351E-01 -5.7191E-03 -2.1186E-01  3.2406E-01  7.5607E-01  6.1161E-02
             2.3655E-02
 GRADIENT:   7.5042E+02  1.0860E+03  4.5317E+00  1.1013E+02  1.9667E+01  3.6651E+01  1.6871E+01 -1.2651E-01  4.0435E+01  8.7738E-01
             8.2420E-01

0ITERATION NO.:   42    OBJECTIVE VALUE:  -1719.45570551373        NO. OF FUNC. EVALS.:  59
 CUMULATIVE NO. OF FUNC. EVALS.:     1426
 NPARAMETR:  1.0474E+00  1.9283E+00  4.7763E-01  4.5259E-01  1.1428E+00  8.9968E-01  7.3208E-01  1.2511E+00  1.9272E+00  9.6191E-01
             9.2650E-01
 PARAMETER:  1.4628E-01  7.5664E-01 -6.3893E-01 -6.9276E-01  2.3351E-01 -5.7191E-03 -2.1186E-01  3.2406E-01  7.5607E-01  6.1161E-02
             2.3655E-02
 GRADIENT:  -8.5315E-01  7.7674E+03 -9.1859E+03  8.4794E+03  8.1344E-01 -4.0670E-02 -1.8053E-01  1.8088E+04  7.7406E+03  6.5415E-02
            -5.8845E+04

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1426
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.2541E-04 -4.0934E-02 -2.3120E-02  3.8429E-02 -5.0129E-02
 SE:             2.9862E-02  2.3287E-02  9.0981E-03  2.3020E-02  2.1599E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7795E-01  7.8781E-02  1.1047E-02  9.5039E-02  2.0295E-02

 ETASHRINKSD(%)  1.0000E-10  2.1986E+01  6.9520E+01  2.2882E+01  2.7639E+01
 ETASHRINKVR(%)  1.0000E-10  3.9139E+01  9.0710E+01  4.0528E+01  4.7639E+01
 EBVSHRINKSD(%)  4.7687E-01  2.1429E+01  7.2363E+01  2.3333E+01  2.5296E+01
 EBVSHRINKVR(%)  9.5146E-01  3.8266E+01  9.2362E+01  4.1221E+01  4.4193E+01
 RELATIVEINF(%)  9.9015E+01  5.9141E+00  1.0089E+00  5.9050E+00  1.7267E+01
 EPSSHRINKSD(%)  4.6271E+01
 EPSSHRINKVR(%)  7.1131E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1719.4557055137320     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -984.30487894999385     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.04
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.20
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1719.456       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  1.93E+00  4.78E-01  4.53E-01  1.14E+00  9.00E-01  7.32E-01  1.25E+00  1.93E+00  9.62E-01  9.26E-01
 


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
+        1.24E+03
 
 TH 2
+        5.92E+01  6.94E+04
 
 TH 3
+       -3.08E+02  4.61E+02  1.57E+06
 
 TH 4
+        2.90E+02 -6.19E+02  4.07E+03  1.49E+06
 
 TH 5
+       -4.99E+00 -1.08E+02 -5.12E+02  4.76E+02  5.31E+02
 
 TH 6
+        1.14E+00  7.07E+01 -3.40E+02  3.31E+02 -9.90E-01  2.42E+02
 
 TH 7
+        7.08E-01  4.04E+01 -1.95E+02  1.75E+02 -1.62E+01 -6.99E-01  1.64E+02
 
 TH 8
+        2.33E+02  2.48E+05  8.57E+03 -3.22E+03  2.11E+02  2.56E+02  1.48E+02  8.88E+05
 
 TH 9
+       -6.55E+05 -7.42E+02  3.45E+03  3.20E+05  3.76E+05  7.27E+01 -6.47E+05 -2.59E+03  6.86E+04
 
 TH10
+        1.07E-02  2.04E+01 -1.80E+02  1.48E+02 -7.17E+01 -3.22E-02  1.59E+01  1.23E+02 -1.04E+06  8.24E+01
 
 TH11
+       -1.03E+03 -5.11E+02  2.37E+03 -5.06E+06 -9.25E+02 -1.12E+03 -6.08E+02 -1.78E+03 -1.08E+06  1.65E+07  1.71E+07
 
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
 #CPUT: Total CPU Time in Seconds,       25.298
Stop Time:
Wed Sep 29 14:04:28 CDT 2021
