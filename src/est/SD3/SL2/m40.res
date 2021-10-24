Sat Oct 23 23:56:07 CDT 2021
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
$DATA ../../../../data/SD3/SL2/dat40.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m40.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2137.78294805807        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4454E+02  3.0535E+01 -2.9914E+01  1.1967E+02  5.7360E+01  8.9595E+01  1.5509E+00  1.0903E+00  1.0068E+01 -2.9872E-01
            -1.5730E-01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2140.93671137148        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      175
 NPARAMETR:  9.9867E-01  1.0516E+00  1.0335E+00  9.6000E-01  9.9231E-01  7.8540E-01  1.0076E+00  1.0040E+00  9.9242E-01  9.8927E-01
             1.0037E+00
 PARAMETER:  9.8673E-02  1.5031E-01  1.3292E-01  5.9173E-02  9.2281E-02 -1.4156E-01  1.0752E-01  1.0402E-01  9.2396E-02  8.9213E-02
             1.0368E-01
 GRADIENT:  -3.2895E+01  6.2420E+00  7.9435E+00  6.8663E+00 -1.1058E+01 -4.2790E+01 -1.2673E+00 -2.9016E+00 -1.9232E+00 -1.6671E+00
             7.8568E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2141.33075572670        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  1.0010E+00  1.3395E+00  9.1903E-01  8.1136E-01  1.0538E+00  7.9945E-01  9.2820E-01  1.2262E+00  1.0810E+00  9.9085E-01
             9.8873E-01
 PARAMETER:  1.0104E-01  3.9231E-01  1.5561E-02 -1.0904E-01  1.5238E-01 -1.2383E-01  2.5489E-02  3.0390E-01  1.7787E-01  9.0804E-02
             8.8667E-02
 GRADIENT:  -3.0366E+01  6.4555E+01  1.4471E+01  4.0782E+01 -2.6092E+01 -3.5605E+01  1.9725E-01 -2.2875E+00 -6.2339E+00 -3.1928E+00
            -5.3040E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2144.92276885667        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      527
 NPARAMETR:  1.0109E+00  1.3207E+00  8.0530E-01  7.8458E-01  1.0243E+00  8.7600E-01  8.7927E-01  9.9451E-01  1.1393E+00  9.8964E-01
             9.9352E-01
 PARAMETER:  1.1086E-01  3.7814E-01 -1.1654E-01 -1.4261E-01  1.2405E-01 -3.2386E-02 -2.8660E-02  9.4499E-02  2.3038E-01  8.9587E-02
             9.3494E-02
 GRADIENT:   2.8544E+00  6.2178E+00  2.5470E+00  9.2611E+00 -7.2246E+00  3.2668E+00 -7.2849E-01 -2.2458E-01 -8.0343E-01  5.7611E-01
            -7.5163E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2146.28069686198        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      703
 NPARAMETR:  1.0115E+00  1.7107E+00  5.2193E-01  5.3615E-01  1.1219E+00  8.6498E-01  7.6650E-01  8.1420E-01  1.4590E+00  1.0336E+00
             9.9154E-01
 PARAMETER:  1.1140E-01  6.3688E-01 -5.5023E-01 -5.2335E-01  2.1502E-01 -4.5048E-02 -1.6592E-01 -1.0555E-01  4.7772E-01  1.3303E-01
             9.1507E-02
 GRADIENT:  -1.1006E+00  1.7936E+01 -2.5075E+00  1.5496E+01  9.8636E-01 -2.5813E+00 -1.4362E+00  4.7577E-01 -8.2346E-01 -7.3214E-01
             1.0107E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2146.38364857833        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      879
 NPARAMETR:  1.0112E+00  1.8210E+00  4.6012E-01  4.6188E-01  1.1701E+00  8.6526E-01  7.4680E-01  7.0279E-01  1.6014E+00  1.0742E+00
             9.9015E-01
 PARAMETER:  1.1112E-01  6.9937E-01 -6.7628E-01 -6.7245E-01  2.5713E-01 -4.4722E-02 -1.9196E-01 -2.5269E-01  5.7085E-01  1.7159E-01
             9.0100E-02
 GRADIENT:  -2.4200E+00  1.5088E+01 -1.3006E+00  1.1169E+01  1.7703E+00 -2.6076E+00 -1.3142E+00  3.1265E-01 -6.2189E-01 -4.1509E-01
             7.5237E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2146.65495945202        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1064             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0132E+00  1.8134E+00  4.4762E-01  4.4820E-01  1.1736E+00  8.7146E-01  7.4913E-01  3.5392E-01  1.6365E+00  1.0841E+00
             9.9024E-01
 PARAMETER:  1.1315E-01  6.9519E-01 -7.0380E-01 -7.0253E-01  2.6012E-01 -3.7584E-02 -1.8884E-01 -9.3868E-01  5.9255E-01  1.8074E-01
             9.0192E-02
 GRADIENT:   5.4545E+02  9.8564E+02  3.9149E+00  1.2831E+02  2.2786E+01  4.2183E+01  1.2888E+01  8.7393E-02  2.5008E+01  1.7356E+00
             5.5301E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2146.66883519011        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1240
 NPARAMETR:  1.0123E+00  1.8225E+00  4.4294E-01  4.4707E-01  1.1710E+00  8.7114E-01  7.5061E-01  3.3202E-01  1.6337E+00  1.0833E+00
             9.9071E-01
 PARAMETER:  1.1226E-01  7.0020E-01 -7.1433E-01 -7.0504E-01  2.5782E-01 -3.7953E-02 -1.8687E-01 -1.0026E+00  5.9087E-01  1.8000E-01
             9.0671E-02
 GRADIENT:   1.1605E+00 -7.7311E+00  6.4379E-01  1.3947E+00  1.1886E+00  7.2940E-02 -1.5359E-01  2.6674E-02  5.5207E-02  8.7442E-02
            -1.7398E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2146.68863094989        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1420
 NPARAMETR:  1.0124E+00  1.8230E+00  4.3830E-01  4.4532E-01  1.1692E+00  8.7115E-01  7.5173E-01  2.0990E-01  1.6350E+00  1.0825E+00
             9.9102E-01
 PARAMETER:  1.1235E-01  7.0050E-01 -7.2485E-01 -7.0896E-01  2.5636E-01 -3.7938E-02 -1.8538E-01 -1.4611E+00  5.9162E-01  1.7929E-01
             9.0977E-02
 GRADIENT:   1.3806E+00 -9.6485E+00  7.1950E-01  5.3529E-01  8.0521E-01  7.4434E-02 -1.7326E-02  1.1091E-02 -2.5340E-02  6.2011E-02
            -7.8651E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2146.70548172096        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1604
 NPARAMETR:  1.0129E+00  1.8203E+00  4.3389E-01  4.4485E-01  1.1658E+00  8.7132E-01  7.5227E-01  2.1752E-02  1.6364E+00  1.0806E+00
             9.9118E-01
 PARAMETER:  1.1280E-01  6.9899E-01 -7.3497E-01 -7.1002E-01  2.5338E-01 -3.7752E-02 -1.8467E-01 -3.7281E+00  5.9247E-01  1.7754E-01
             9.1137E-02
 GRADIENT:   2.6949E+00 -1.4212E+01  9.3642E-02 -4.0678E-01  6.2671E-01  1.4872E-01  6.3772E-02  2.2440E-04  3.6885E-01  3.0052E-01
             9.3914E-02

0ITERATION NO.:   48    OBJECTIVE VALUE:  -2146.70653434569        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:     1703
 NPARAMETR:  1.0128E+00  1.8187E+00  4.3165E-01  4.4546E-01  1.1634E+00  8.7130E-01  7.5353E-01  1.0000E-02  1.6339E+00  1.0798E+00
             9.9143E-01
 PARAMETER:  1.1271E-01  6.9920E-01 -7.3555E-01 -7.0936E-01  2.5281E-01 -3.7792E-02 -1.8483E-01 -4.6308E+00  5.9110E-01  1.7504E-01
             9.1020E-02
 GRADIENT:  -2.3374E-02  8.9499E-01  2.5164E-01 -1.0047E-01  5.4047E-01 -3.4048E-03 -1.0103E-01  0.0000E+00  4.7540E-03 -9.0563E-02
            -9.0027E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1703
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.3321E-04 -3.5214E-02 -2.2061E-04  3.4838E-02 -3.8774E-02
 SE:             2.9855E-02  2.4210E-02  9.7046E-05  2.2908E-02  2.2931E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8842E-01  1.4579E-01  2.3012E-02  1.2832E-01  9.0861E-02

 ETASHRINKSD(%)  1.0000E-10  1.8894E+01  9.9675E+01  2.3255E+01  2.3178E+01
 ETASHRINKVR(%)  1.0000E-10  3.4219E+01  9.9999E+01  4.1101E+01  4.0984E+01
 EBVSHRINKSD(%)  4.1355E-01  1.8350E+01  9.9724E+01  2.5498E+01  2.0514E+01
 EBVSHRINKVR(%)  8.2539E-01  3.3332E+01  9.9999E+01  4.4495E+01  3.6820E+01
 RELATIVEINF(%)  9.9137E+01  8.0018E+00  1.9997E-04  6.4282E+00  1.8435E+01
 EPSSHRINKSD(%)  3.4198E+01
 EPSSHRINKVR(%)  5.6701E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2146.7065343456920     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1227.7680011410193     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.40
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2146.707       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.82E+00  4.34E-01  4.45E-01  1.17E+00  8.71E-01  7.52E-01  1.00E-02  1.63E+00  1.08E+00  9.91E-01
 


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
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,      134.961
Stop Time:
Sat Oct 23 23:56:27 CDT 2021
