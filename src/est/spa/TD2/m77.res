Wed Sep 29 19:18:32 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat77.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m77.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1687.34243771182        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9570E+02  3.7529E+01  1.2207E+01  7.2830E+01 -3.9397E+01  3.9898E+01  5.1505E+00 -1.4716E+00  4.5052E+01  9.8352E+00
            -1.4522E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1692.76783308668        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      194
 NPARAMETR:  1.0096E+00  1.1069E+00  1.0374E+00  9.6005E-01  1.0974E+00  1.0627E+00  1.0107E+00  1.0241E+00  7.3918E-01  9.4839E-01
             1.0830E+00
 PARAMETER:  1.0952E-01  2.0156E-01  1.3673E-01  5.9235E-02  1.9292E-01  1.6085E-01  1.1060E-01  1.2377E-01 -2.0221E-01  4.7011E-02
             1.7969E-01
 GRADIENT:  -2.9027E+01  3.4811E+01  6.8270E+00  4.0417E+01  3.4406E+00  1.1849E+01 -1.0962E+01 -8.4408E+00 -1.1390E+01 -4.2410E+00
             1.1380E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1694.17981204382        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  1.0120E+00  9.8647E-01  1.1434E+00  1.0348E+00  1.0786E+00  1.0436E+00  1.2497E+00  1.2654E+00  6.6303E-01  8.6729E-01
             1.0587E+00
 PARAMETER:  1.1194E-01  8.6382E-02  2.3405E-01  1.3421E-01  1.7564E-01  1.4264E-01  3.2291E-01  3.3542E-01 -3.1094E-01 -4.2386E-02
             1.5706E-01
 GRADIENT:  -2.1984E+01  3.0379E+01  8.0313E-01  3.5630E+01  5.6828E+00  5.2397E+00  1.5774E+00 -1.5183E+00 -7.6849E+00 -9.1995E+00
             1.5086E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1697.16348549244        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      555
 NPARAMETR:  1.0201E+00  7.9994E-01  2.2268E+00  1.1776E+00  1.3516E+00  1.0241E+00  1.2696E+00  2.1633E+00  7.2230E-01  1.2422E+00
             1.0442E+00
 PARAMETER:  1.1986E-01 -1.2322E-01  9.0055E-01  2.6345E-01  4.0126E-01  1.2386E-01  3.3871E-01  8.7164E-01 -2.2532E-01  3.1688E-01
             1.4324E-01
 GRADIENT:   4.0901E-01  1.1836E+01 -5.8127E+00  2.4320E+01  8.4556E+00 -1.0426E+00  2.3241E+00  2.3478E+00  2.7956E+00  1.1611E+00
            -3.8278E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1697.85164247186        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      731
 NPARAMETR:  1.0194E+00  5.8248E-01  3.4026E+00  1.3301E+00  1.4143E+00  1.0246E+00  1.3607E+00  2.7183E+00  6.9111E-01  1.2991E+00
             1.0504E+00
 PARAMETER:  1.1921E-01 -4.4046E-01  1.3245E+00  3.8527E-01  4.4663E-01  1.2427E-01  4.0798E-01  1.1000E+00 -2.6946E-01  3.6163E-01
             1.4915E-01
 GRADIENT:   2.9244E-01  7.1329E+00  6.4188E-02  1.5492E+01 -1.1627E+00 -3.5514E-01  4.6366E-04  6.7653E-02 -3.4571E-01  5.2541E-01
             5.6834E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1698.04471529604        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      907
 NPARAMETR:  1.0221E+00  3.9609E-01  5.6976E+00  1.4767E+00  1.5116E+00  1.0232E+00  1.4061E+00  3.6808E+00  6.8097E-01  1.3654E+00
             1.0530E+00
 PARAMETER:  1.2190E-01 -8.2612E-01  1.8400E+00  4.8980E-01  5.1314E-01  1.2297E-01  4.4084E-01  1.4031E+00 -2.8424E-01  4.1143E-01
             1.5163E-01
 GRADIENT:   4.1084E+00  8.5990E+00 -2.8655E-01  3.2168E+01  4.2782E+00 -1.0592E+00  8.3478E-01 -1.0128E+00  2.7001E+00  1.4064E+00
             2.2913E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1698.34517354329        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1090
 NPARAMETR:  1.0213E+00  3.4230E-01  6.6984E+00  1.5163E+00  1.5336E+00  1.0267E+00  1.0792E+00  4.0542E+00  6.8612E-01  1.3712E+00
             1.0576E+00
 PARAMETER:  1.2110E-01 -9.7208E-01  2.0019E+00  5.1631E-01  5.2765E-01  1.2635E-01  1.7625E-01  1.4998E+00 -2.7671E-01  4.1570E-01
             1.5601E-01
 GRADIENT:   1.4457E+00  7.4755E+00 -2.9908E+00  3.2874E+01  1.7127E+00  2.1610E-01  5.3612E-01  4.6609E+00  2.3981E+00 -2.5100E-01
             1.6341E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1698.63638172736        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1266
 NPARAMETR:  1.0207E+00  3.4203E-01  6.6997E+00  1.5161E+00  1.5335E+00  1.0267E+00  1.1610E-01  4.0534E+00  7.0649E-01  1.3711E+00
             1.0558E+00
 PARAMETER:  1.2049E-01 -9.7285E-01  2.0021E+00  5.1611E-01  5.2754E-01  1.2638E-01 -2.0533E+00  1.4996E+00 -2.4744E-01  4.1564E-01
             1.5431E-01
 GRADIENT:   1.7266E-01  7.9651E+00 -2.9594E+00  3.9927E+01  1.7021E+00  2.9926E-01  1.3156E-02  4.3653E+00  2.2772E+00 -3.0639E-01
             7.9470E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1698.86231008243        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1446
 NPARAMETR:  1.0216E+00  3.3820E-01  6.7284E+00  1.5113E+00  1.5314E+00  1.0240E+00  1.0000E-02  4.0374E+00  7.0642E-01  1.3703E+00
             1.0545E+00
 PARAMETER:  1.2139E-01 -9.8412E-01  2.0063E+00  5.1299E-01  5.2616E-01  1.2375E-01 -2.1724E+01  1.4956E+00 -2.4754E-01  4.1501E-01
             1.5304E-01
 GRADIENT:   2.4610E+00  4.8461E+00 -1.7316E+00  2.0473E+01  2.8294E+00 -7.1459E-01  0.0000E+00  2.4535E+00  3.0807E+00  3.0891E-01
             1.4532E+00

0ITERATION NO.:   42    OBJECTIVE VALUE:  -1698.86589138050        NO. OF FUNC. EVALS.:  69
 CUMULATIVE NO. OF FUNC. EVALS.:     1515
 NPARAMETR:  1.0216E+00  3.3865E-01  6.7080E+00  1.5124E+00  1.5301E+00  1.0253E+00  1.0000E-02  4.0464E+00  7.0605E-01  1.3694E+00
             1.0542E+00
 PARAMETER:  1.2138E-01 -9.8432E-01  2.0064E+00  5.1293E-01  5.2614E-01  1.2378E-01 -2.2086E+01  1.4955E+00 -2.4769E-01  4.1500E-01
             1.5299E-01
 GRADIENT:   3.1900E-02 -3.1901E+03  1.5676E+03 -6.0962E+03  5.9765E+03 -7.8372E-01  0.0000E+00 -2.1333E+03  6.3369E+03  7.5368E+03
             1.0254E+04

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1515
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.7300E-03 -1.8684E-04 -4.6073E-02 -2.2701E-02 -5.9903E-02
 SE:             2.9864E-02  6.7203E-05  1.6026E-02  2.8525E-02  2.0751E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5380E-01  5.4333E-03  4.0423E-03  4.2612E-01  3.8930E-03

 ETASHRINKSD(%)  1.0000E-10  9.9775E+01  4.6310E+01  4.4382E+00  3.0480E+01
 ETASHRINKVR(%)  1.0000E-10  9.9999E+01  7.1174E+01  8.6793E+00  5.1670E+01
 EBVSHRINKSD(%)  4.4443E-01  9.9784E+01  5.9004E+01  3.6667E+00  2.3554E+01
 EBVSHRINKVR(%)  8.8688E-01  1.0000E+02  8.3193E+01  7.1989E+00  4.1560E+01
 RELATIVEINF(%)  9.8426E+01  2.7944E-05  1.0719E+01  5.6734E+00  3.5086E+01
 EPSSHRINKSD(%)  4.3783E+01
 EPSSHRINKVR(%)  6.8397E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1698.8658913804961     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -963.71506481675794     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.98
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.50
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1698.866       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  3.38E-01  6.73E+00  1.51E+00  1.53E+00  1.02E+00  1.00E-02  4.04E+00  7.06E-01  1.37E+00  1.05E+00
 


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
+        1.01E+03
 
 TH 2
+       -9.39E+01  7.12E+05
 
 TH 3
+        2.10E+00 -7.53E+01  4.35E+02
 
 TH 4
+       -4.60E+01  3.05E+05 -2.13E+01  2.62E+05
 
 TH 5
+        3.64E+01 -1.31E+03  3.36E+01  1.25E+05  1.22E+05
 
 TH 6
+       -3.56E-02 -8.32E+01  2.11E+00 -3.79E+01  3.42E+01  1.88E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -6.85E+00  1.64E+02 -3.86E+01  3.55E+01 -9.57E+01 -5.17E+00  0.00E+00  2.25E+03
 
 TH 9
+        1.49E+02 -3.80E+02  5.34E+00 -1.30E+02  9.63E+01  1.56E+02  0.00E+00 -1.04E+01  2.56E+06
 
 TH10
+        4.87E+01 -1.83E+03 -5.89E+01  1.78E+05  1.71E+05  4.92E+01  0.00E+00  1.27E+02  1.40E+02  2.41E+05
 
 TH11
+        1.61E+02 -1.46E+06  4.65E+00  1.64E+02 -1.15E+01  1.76E+02  0.00E+00 -2.32E+01  5.14E+02 -1.09E+02  3.01E+06
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       29.588
Stop Time:
Wed Sep 29 19:19:30 CDT 2021
