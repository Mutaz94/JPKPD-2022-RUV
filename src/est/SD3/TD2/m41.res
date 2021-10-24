Sun Oct 24 00:44:11 CDT 2021
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
$DATA ../../../../data/SD3/TD2/dat41.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m41.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2178.61452009286        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.5611E+02 -4.3709E+01 -1.9167E+00 -1.8356E+01  3.7390E+01  4.7178E+01  1.1719E+00  1.5473E-01  2.7017E+01 -1.2777E+01
             7.3091E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2189.66482875058        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      192
 NPARAMETR:  1.0368E+00  1.0853E+00  9.5074E-01  1.0200E+00  9.9473E-01  9.6479E-01  1.0030E+00  1.0147E+00  8.5734E-01  1.0899E+00
             8.9126E-01
 PARAMETER:  1.3611E-01  1.8189E-01  4.9480E-02  1.1984E-01  9.4713E-02  6.4155E-02  1.0296E-01  1.1456E-01 -5.3924E-02  1.8612E-01
            -1.5119E-02
 GRADIENT:   9.7249E+00 -4.1226E-01 -7.8660E+00  1.4010E+01  5.2743E+00 -5.6340E+00 -4.5601E+00  1.9725E+00 -5.1902E+00  8.3546E-01
            -1.3797E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2189.98625009842        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  1.0398E+00  1.2029E+00  9.1201E-01  9.5199E-01  1.0320E+00  9.7877E-01  9.8815E-01  9.4898E-01  8.9362E-01  1.1475E+00
             8.9260E-01
 PARAMETER:  1.3905E-01  2.8476E-01  7.8968E-03  5.0798E-02  1.3147E-01  7.8538E-02  8.8081E-02  4.7632E-02 -1.2480E-02  2.3755E-01
            -1.3622E-02
 GRADIENT:   1.4915E+01  1.4921E+01  7.4466E-01  1.7349E+01 -1.7928E+00 -2.3664E-01  1.2098E-01  1.2821E-01 -3.8688E+00  3.3850E+00
            -1.2002E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2190.81001971802        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      549
 NPARAMETR:  1.0367E+00  1.3776E+00  6.7417E-01  8.2315E-01  9.9782E-01  9.8460E-01  9.0271E-01  4.9175E-01  1.0034E+00  1.0995E+00
             9.0308E-01
 PARAMETER:  1.3603E-01  4.2031E-01 -2.9427E-01 -9.4616E-02  9.7815E-02  8.4482E-02 -2.3579E-03 -6.0978E-01  1.0340E-01  1.9489E-01
            -1.9465E-03
 GRADIENT:   2.6211E+00  5.1512E+00 -4.0410E+00  7.8573E+00 -1.3835E+00  1.3235E+00  2.0025E+00  6.8106E-01 -1.1952E+00  2.5788E+00
            -2.5306E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2190.95760442758        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      726
 NPARAMETR:  1.0363E+00  1.4631E+00  6.3383E-01  7.6437E-01  1.0297E+00  9.8155E-01  8.4700E-01  3.8779E-01  1.0713E+00  1.1136E+00
             9.0515E-01
 PARAMETER:  1.3564E-01  4.8052E-01 -3.5597E-01 -1.6870E-01  1.2926E-01  8.1382E-02 -6.6053E-02 -8.4729E-01  1.6890E-01  2.0757E-01
             3.5043E-04
 GRADIENT:   8.4609E-01 -5.6925E-01 -1.5173E+00  1.5520E+00  8.7051E-01 -4.7484E-02 -3.6423E-01  3.2040E-01 -6.1294E-02  3.8307E-01
            -6.4875E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2191.03558867606        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      905
 NPARAMETR:  1.0375E+00  1.4685E+00  6.2706E-01  7.6080E-01  1.0291E+00  9.8232E-01  8.5175E-01  2.2231E-01  1.0772E+00  1.1239E+00
             9.0629E-01
 PARAMETER:  1.3677E-01  4.8425E-01 -3.6671E-01 -1.7338E-01  1.2865E-01  8.2164E-02 -6.0459E-02 -1.4037E+00  1.7441E-01  2.1679E-01
             1.6059E-03
 GRADIENT:   3.4642E+00  2.0973E+00  9.9823E-01  1.2024E+00 -9.6082E-01  2.4522E-01  3.3068E-01  6.9151E-02  1.6220E-01  9.9436E-01
             1.3303E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2191.08222379729        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1085             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0386E+00  1.4757E+00  6.1035E-01  7.5329E-01  1.0232E+00  9.8215E-01  8.4745E-01  3.6339E-02  1.0826E+00  1.1088E+00
             9.0623E-01
 PARAMETER:  1.3790E-01  4.8915E-01 -3.9373E-01 -1.8330E-01  1.2293E-01  8.1985E-02 -6.5519E-02 -3.2149E+00  1.7938E-01  2.0329E-01
             1.5345E-03
 GRADIENT:   7.2445E+02  4.4842E+02  6.6395E+00  9.1483E+01  1.1362E+01  4.4851E+01  5.6137E+00  9.0517E-03  8.8520E+00  2.2597E+00
             9.1548E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2191.08521915685        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     1241
 NPARAMETR:  1.0377E+00  1.4750E+00  6.1050E-01  7.5380E-01  1.0234E+00  9.8212E-01  8.4804E-01  1.2787E-02  1.0828E+00  1.1101E+00
             9.0629E-01
 PARAMETER:  1.3697E-01  4.8862E-01 -3.9348E-01 -1.8263E-01  1.2314E-01  8.1963E-02 -6.4832E-02 -4.2593E+00  1.7958E-01  2.0446E-01
             1.6011E-03
 GRADIENT:   3.6248E+00 -1.6327E+00 -4.4150E-01 -3.9294E-01  2.0305E-01  1.0777E-01  2.3138E-02  3.8377E-04  1.1862E-01 -7.0508E-03
             2.5961E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2191.08613241279        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     1433             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0376E+00  1.4737E+00  6.1215E-01  7.5462E-01  1.0227E+00  9.8212E-01  8.4813E-01  1.0000E-02  1.0809E+00  1.1104E+00
             9.0623E-01
 PARAMETER:  1.3696E-01  4.8774E-01 -3.9078E-01 -1.8154E-01  1.2245E-01  8.1954E-02 -6.4721E-02 -5.5983E+00  1.7776E-01  2.0475E-01
             1.5430E-03
 GRADIENT:   7.1739E+02  4.4620E+02  7.1479E+00  9.0584E+01  1.0649E+01  4.4970E+01  5.5505E+00  0.0000E+00  8.6536E+00  2.4232E+00
             9.0986E-01

0ITERATION NO.:   44    OBJECTIVE VALUE:  -2191.08688850649        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:     1575
 NPARAMETR:  1.0376E+00  1.4729E+00  6.1328E-01  7.5514E-01  1.0222E+00  9.8211E-01  8.4835E-01  1.0000E-02  1.0801E+00  1.1104E+00
             9.0620E-01
 PARAMETER:  1.3695E-01  4.8732E-01 -3.9158E-01 -1.8100E-01  1.2316E-01  8.1945E-02 -6.3721E-02 -5.5983E+00  1.7861E-01  2.0499E-01
             1.6157E-03
 GRADIENT:   3.8217E-04  8.9071E-02 -5.7402E-01 -8.7815E-02  7.0461E-01 -1.4660E-03  4.7754E-02  0.0000E+00  1.3062E-01  2.5910E-02
             4.3292E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1575
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.9476E-04 -2.4235E-02 -3.7160E-04  1.8895E-02 -2.7319E-02
 SE:             2.9894E-02  2.2546E-02  1.5025E-04  2.4137E-02  2.3817E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9213E-01  2.8241E-01  1.3391E-02  4.3373E-01  2.5136E-01

 ETASHRINKSD(%)  1.0000E-10  2.4468E+01  9.9497E+01  1.9137E+01  2.0210E+01
 ETASHRINKVR(%)  1.0000E-10  4.2950E+01  9.9997E+01  3.4611E+01  3.6336E+01
 EBVSHRINKSD(%)  2.7847E-01  2.4212E+01  9.9582E+01  2.0585E+01  1.7893E+01
 EBVSHRINKVR(%)  5.5617E-01  4.2561E+01  9.9998E+01  3.6932E+01  3.2585E+01
 RELATIVEINF(%)  9.9301E+01  3.9993E+00  2.7904E-04  4.6086E+00  1.3491E+01
 EPSSHRINKSD(%)  3.4297E+01
 EPSSHRINKVR(%)  5.6831E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2191.0868885064870     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1272.1483553018143     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.17
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2191.087       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.47E+00  6.12E-01  7.55E-01  1.02E+00  9.82E-01  8.49E-01  1.00E-02  1.08E+00  1.11E+00  9.06E-01
 


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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       53.604
Stop Time:
Sun Oct 24 00:44:22 CDT 2021
