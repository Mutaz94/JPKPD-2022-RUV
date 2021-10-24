Sat Oct 23 15:04:08 CDT 2021
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
$DATA ../../../../data/SD1/SL2/dat26.csv ignore=@
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
Current Date:       23 OCT 2021
Days until program expires : 176
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
 NO. OF DATA RECS IN DATA SET:      997
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      897
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
 RAW OUTPUT FILE (FILE): m26.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -567.636502207352        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3243E+02  5.0401E+01  3.7663E+02  1.0984E+02  1.5884E+02  5.9535E+01 -1.0466E+02 -7.9095E+02 -1.8685E+02 -5.7415E+01
            -5.4168E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2796.82735925719        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:       88
 NPARAMETR:  1.0049E+00  1.2157E+00  8.3369E-01  8.8793E-01  1.0399E+00  8.7875E-01  1.1034E+00  1.0860E+00  9.2748E-01  1.0312E+00
             2.5852E+00
 PARAMETER:  1.0489E-01  2.9534E-01 -8.1889E-02 -1.8858E-02  1.3908E-01 -2.9254E-02  1.9836E-01  1.8251E-01  2.4719E-02  1.3071E-01
             1.0498E+00
 GRADIENT:   4.5168E+01  3.4185E+01 -9.4707E+00  2.4321E+01  1.7381E+00 -1.5144E+01  1.0771E+01  1.0839E+01 -1.1635E+00 -7.6940E+00
             5.7707E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2803.23679201451        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      182
 NPARAMETR:  1.0080E+00  1.2615E+00  8.4279E-01  8.6296E-01  1.0880E+00  9.1033E-01  9.4878E-01  4.8469E-01  9.4534E-01  1.1544E+00
             2.5905E+00
 PARAMETER:  1.0801E-01  3.3229E-01 -7.1036E-02 -4.7388E-02  1.8434E-01  6.0510E-03  4.7420E-02 -6.2425E-01  4.3791E-02  2.4361E-01
             1.0518E+00
 GRADIENT:  -1.5140E+01 -2.5315E+01 -6.3081E+00  1.7516E+01  3.6506E-01 -8.7471E+00 -4.5227E+00  1.4563E+00 -2.7045E+00 -5.4408E-01
             3.7267E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2815.03303035993        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      360
 NPARAMETR:  1.0010E+00  1.8496E+00  9.2366E-01  5.4150E-01  1.5811E+00  9.6228E-01  7.4632E-01  5.6205E-01  1.2699E+00  1.4280E+00
             2.4738E+00
 PARAMETER:  1.0096E-01  7.1497E-01  2.0594E-02 -5.1341E-01  5.5812E-01  6.1551E-02 -1.9261E-01 -4.7616E-01  3.3890E-01  4.5627E-01
             1.0058E+00
 GRADIENT:  -2.8077E+01  7.2753E+01  7.5326E+00  4.4654E+01  8.1324E+00  1.0365E+01 -6.4617E+00 -5.0230E-01 -9.9663E+00 -7.5640E+00
            -4.6368E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2825.14378458592        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      535
 NPARAMETR:  1.0122E+00  2.1886E+00  4.2750E-01  2.7769E-01  1.7864E+00  9.3404E-01  6.6013E-01  9.4118E-02  1.9397E+00  1.6005E+00
             2.4924E+00
 PARAMETER:  1.1213E-01  8.8328E-01 -7.4979E-01 -1.1813E+00  6.8020E-01  3.1760E-02 -3.1532E-01 -2.2632E+00  7.6256E-01  5.7034E-01
             1.0132E+00
 GRADIENT:  -1.6634E+00 -3.0179E-01 -3.9477E+00  9.3149E+00  6.0985E+00  4.1341E-01 -6.3015E+00  1.3621E-02 -1.8989E+00 -1.9953E+00
            -5.3064E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2826.72647627644        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      710
 NPARAMETR:  1.0128E+00  2.3958E+00  2.9016E-01  1.3740E-01  1.9660E+00  9.3230E-01  6.5163E-01  1.0415E-02  2.7364E+00  1.7046E+00
             2.4950E+00
 PARAMETER:  1.1274E-01  9.7373E-01 -1.1373E+00 -1.8848E+00  7.7599E-01  2.9898E-02 -3.2827E-01 -4.4645E+00  1.1066E+00  6.3333E-01
             1.0143E+00
 GRADIENT:   2.5601E-01  4.9566E-01 -2.5578E+00  1.3037E+00  5.7003E+00 -8.9169E-02 -1.6682E+00  1.0469E-04 -1.5994E+00 -5.9295E-01
            -9.7066E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2827.23295690457        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      885
 NPARAMETR:  1.0124E+00  2.3942E+00  3.8170E-01  1.3297E-01  1.9659E+00  9.3280E-01  6.4524E-01  1.0000E-02  3.0547E+00  1.7045E+00
             2.4927E+00
 PARAMETER:  1.1230E-01  9.7307E-01 -8.6313E-01 -1.9176E+00  7.7596E-01  3.0432E-02 -3.3814E-01 -4.6446E+00  1.2167E+00  6.3326E-01
             1.0134E+00
 GRADIENT:  -2.7469E-01 -5.8437E+00 -2.8426E-01 -9.2837E-03 -5.5704E-01  1.8587E-01 -8.3675E-01  0.0000E+00  1.9531E-01  4.1318E-01
             7.3172E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2827.68289851537        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:     1002
 NPARAMETR:  1.0129E+00  2.3930E+00  4.4497E-01  1.3214E-01  1.9723E+00  9.3205E-01  6.5443E-01  1.0000E-02  3.0492E+00  1.7002E+00
             2.4916E+00
 PARAMETER:  1.1277E-01  9.7255E-01 -7.0976E-01 -1.9239E+00  7.7920E-01  2.9629E-02 -3.2399E-01 -4.6446E+00  1.2149E+00  6.3075E-01
             1.0129E+00
 GRADIENT:   1.1541E+00 -9.6198E+00  3.9233E-01  2.9744E-01 -2.0839E+00 -1.2303E-01  8.6740E-01  0.0000E+00 -8.4629E-01 -4.9256E-01
             4.8369E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2827.74601417263        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1178
 NPARAMETR:  1.0124E+00  2.4055E+00  4.3375E-01  1.2605E-01  1.9822E+00  9.3243E-01  6.5085E-01  1.0000E-02  3.1432E+00  1.7074E+00
             2.4907E+00
 PARAMETER:  1.1235E-01  9.7777E-01 -7.3528E-01 -1.9711E+00  7.8419E-01  3.0034E-02 -3.2947E-01 -4.6446E+00  1.2453E+00  6.3500E-01
             1.0126E+00
 GRADIENT:   2.2512E-02 -4.2273E+00  1.3282E-01  9.6353E-01 -1.1102E+00 -3.8903E-03 -2.8308E-02  0.0000E+00 -2.6147E-01 -1.0771E-01
            -8.3074E-01

0ITERATION NO.:   41    OBJECTIVE VALUE:  -2827.74601417263        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1200
 NPARAMETR:  1.0124E+00  2.4055E+00  4.3375E-01  1.2605E-01  1.9822E+00  9.3243E-01  6.5085E-01  1.0000E-02  3.1432E+00  1.7074E+00
             2.4907E+00
 PARAMETER:  1.1235E-01  9.7777E-01 -7.3528E-01 -1.9711E+00  7.8419E-01  3.0034E-02 -3.2947E-01 -4.6446E+00  1.2453E+00  6.3500E-01
             1.0126E+00
 GRADIENT:   2.2512E-02 -4.2273E+00  1.3282E-01  9.6353E-01 -1.1102E+00 -3.8903E-03 -2.8308E-02  0.0000E+00 -2.6147E-01 -1.0771E-01
            -8.3074E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1200
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3786E-03 -2.8994E-02 -2.4336E-05  3.5544E-02 -1.9979E-02
 SE:             2.9440E-02  2.5513E-02  1.3240E-05  1.6459E-02  2.6564E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6265E-01  2.5576E-01  6.6042E-02  3.0804E-02  4.5197E-01

 ETASHRINKSD(%)  1.3727E+00  1.4529E+01  9.9956E+01  4.4862E+01  1.1008E+01
 ETASHRINKVR(%)  2.7265E+00  2.6947E+01  1.0000E+02  6.9598E+01  2.0804E+01
 EBVSHRINKSD(%)  1.5629E+00  1.1985E+01  9.9916E+01  5.6543E+01  7.6690E+00
 EBVSHRINKVR(%)  3.1013E+00  2.2534E+01  1.0000E+02  8.1115E+01  1.4750E+01
 RELATIVEINF(%)  9.6834E+01  2.1251E+01  5.2561E-05  4.8997E+00  6.2373E+01
 EPSSHRINKSD(%)  1.6370E+01
 EPSSHRINKVR(%)  3.0061E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          897
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1648.5757285691827     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2827.7460141726287     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1179.1702856034460     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.50
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2827.746       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  2.41E+00  4.34E-01  1.26E-01  1.98E+00  9.32E-01  6.51E-01  1.00E-02  3.14E+00  1.71E+00  2.49E+00
 


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
 #CPUT: Total CPU Time in Seconds,       79.689
Stop Time:
Sat Oct 23 15:04:22 CDT 2021
