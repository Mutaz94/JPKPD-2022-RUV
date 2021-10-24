Sun Oct 24 02:48:04 CDT 2021
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
$DATA ../../../../data/SD4/S1/dat62.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m62.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1688.08140241358        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3140E+02  4.4827E+01 -3.0332E+01  1.1422E+02  7.8897E+01  8.0220E+01 -2.5393E+00  3.7665E+00 -1.7267E+01  3.1776E+00
             2.1454E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1693.18759856779        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      192
 NPARAMETR:  1.0033E+00  9.3029E-01  1.0009E+00  1.0167E+00  9.1280E-01  8.6582E-01  1.0208E+00  9.7466E-01  1.2188E+00  9.1774E-01
             9.1437E-01
 PARAMETER:  1.0331E-01  2.7741E-02  1.0088E-01  1.1653E-01  8.7606E-03 -4.4075E-02  1.2058E-01  7.4329E-02  2.9787E-01  1.4159E-02
             1.0475E-02
 GRADIENT:  -3.9473E+01  2.6437E-01  6.4816E+00 -2.2930E+00 -9.1098E+00 -2.3301E+01  6.7100E+00  1.4593E+00  2.2380E+01  1.6342E+00
            -1.3705E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1694.49068568014        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      368
 NPARAMETR:  1.0073E+00  7.6207E-01  8.7593E-01  1.1353E+00  7.8444E-01  8.6855E-01  1.1094E+00  7.5419E-01  1.0766E+00  8.3663E-01
             9.1007E-01
 PARAMETER:  1.0725E-01 -1.7171E-01 -3.2474E-02  2.2693E-01 -1.4278E-01 -4.0933E-02  2.0384E-01 -1.8211E-01  1.7385E-01 -7.8377E-02
             5.7626E-03
 GRADIENT:  -2.8795E+01  1.4020E+01  1.5791E+00  2.9518E+01 -1.8407E+01 -2.1728E+01 -1.1332E-01  1.6650E+00  1.3483E+01  1.3439E+00
            -1.6053E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1698.31254748514        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      543
 NPARAMETR:  1.0200E+00  7.6444E-01  7.3445E-01  1.0925E+00  7.2850E-01  9.2004E-01  1.4651E+00  4.6557E-01  9.3658E-01  7.3834E-01
             9.5252E-01
 PARAMETER:  1.1978E-01 -1.6861E-01 -2.0864E-01  1.8849E-01 -2.1677E-01  1.6657E-02  4.8194E-01 -6.6449E-01  3.4477E-02 -2.0335E-01
             5.1356E-02
 GRADIENT:   3.1619E+00 -2.2210E+00 -5.9978E+00 -7.6706E-01  6.7974E+00  2.1391E+00  6.8427E-01  1.4389E+00 -1.8382E+00  2.5828E+00
             2.4729E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1698.84897853377        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      718
 NPARAMETR:  1.0203E+00  8.5993E-01  6.5271E-01  1.0270E+00  7.1710E-01  9.1738E-01  1.3294E+00  2.0979E-01  9.8486E-01  6.9651E-01
             9.4912E-01
 PARAMETER:  1.2005E-01 -5.0904E-02 -3.2663E-01  1.2664E-01 -2.3253E-01  1.3769E-02  3.8474E-01 -1.4617E+00  8.4739E-02 -2.6167E-01
             4.7782E-02
 GRADIENT:   4.6501E-01  8.6675E-01 -5.9953E-01  1.1854E+00 -2.1690E-01  3.7741E-01 -8.2553E-03  2.5919E-01 -5.5881E-01 -1.6223E-01
             4.1774E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1698.96462927991        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      893
 NPARAMETR:  1.0202E+00  8.7589E-01  6.4770E-01  1.0162E+00  7.2168E-01  9.1645E-01  1.3084E+00  4.9028E-02  9.9741E-01  7.0515E-01
             9.4776E-01
 PARAMETER:  1.1998E-01 -3.2516E-02 -3.3433E-01  1.1604E-01 -2.2617E-01  1.2747E-02  3.6881E-01 -2.9154E+00  9.7403E-02 -2.4934E-01
             4.6348E-02
 GRADIENT:   1.1445E-01  1.0803E-01  7.7051E-01  5.5263E-03 -4.3679E-01 -7.7501E-02 -2.0726E-01  1.2414E-02  2.8292E-01 -2.4551E-01
            -1.1013E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1698.97241887580        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1068
 NPARAMETR:  1.0202E+00  8.7450E-01  6.4474E-01  1.0164E+00  7.1953E-01  9.1669E-01  1.3122E+00  1.0000E-02  9.9483E-01  7.0329E-01
             9.4795E-01
 PARAMETER:  1.1996E-01 -3.4104E-02 -3.3890E-01  1.1625E-01 -2.2916E-01  1.3012E-02  3.7169E-01 -4.8309E+00  9.4815E-02 -2.5199E-01
             4.6547E-02
 GRADIENT:  -6.0385E-03 -5.0714E-02 -7.1601E-02 -7.0998E-02  9.2312E-02  1.7695E-02  2.2383E-02  0.0000E+00 -2.2793E-02  1.3904E-03
             1.7244E-02

0ITERATION NO.:   31    OBJECTIVE VALUE:  -1698.97241887580        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1090
 NPARAMETR:  1.0202E+00  8.7450E-01  6.4474E-01  1.0164E+00  7.1953E-01  9.1669E-01  1.3122E+00  1.0000E-02  9.9483E-01  7.0329E-01
             9.4795E-01
 PARAMETER:  1.1996E-01 -3.4104E-02 -3.3890E-01  1.1625E-01 -2.2916E-01  1.3012E-02  3.7169E-01 -4.8309E+00  9.4815E-02 -2.5199E-01
             4.6547E-02
 GRADIENT:  -6.0385E-03 -5.0714E-02 -7.1601E-02 -7.0998E-02  9.2312E-02  1.7695E-02  2.2383E-02  0.0000E+00 -2.2793E-02  1.3904E-03
             1.7244E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1090
 NO. OF SIG. DIGITS IN FINAL EST.:  2.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2101E-03  2.0569E-03 -5.0801E-04 -2.9737E-03 -6.6894E-03
 SE:             2.9856E-02  2.1819E-02  2.1247E-04  2.5999E-02  2.1817E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6767E-01  9.2489E-01  1.6805E-02  9.0894E-01  7.5914E-01

 ETASHRINKSD(%)  1.0000E-10  2.6904E+01  9.9288E+01  1.2899E+01  2.6910E+01
 ETASHRINKVR(%)  1.0000E-10  4.6570E+01  9.9995E+01  2.4134E+01  4.6579E+01
 EBVSHRINKSD(%)  4.5260E-01  2.6385E+01  9.9350E+01  1.2936E+01  2.7015E+01
 EBVSHRINKVR(%)  9.0316E-01  4.5808E+01  9.9996E+01  2.4199E+01  4.6732E+01
 RELATIVEINF(%)  9.8880E+01  4.9045E+00  5.1740E-04  1.0084E+01  4.1311E+00
 EPSSHRINKSD(%)  4.4856E+01
 EPSSHRINKVR(%)  6.9592E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1698.9724188757962     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -963.82159231205799     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     4.93
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1698.972       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  8.74E-01  6.45E-01  1.02E+00  7.20E-01  9.17E-01  1.31E+00  1.00E-02  9.95E-01  7.03E-01  9.48E-01
 


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
 #CPUT: Total CPU Time in Seconds,       31.461
Stop Time:
Sun Oct 24 02:48:12 CDT 2021
