Wed Sep 29 16:52:59 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat79.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m79.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1592.67420153708        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9575E+02 -2.9187E+01 -1.9320E+01  2.3861E+01  4.3600E+01  5.6638E+01 -1.3510E+01  1.0411E-02  1.6136E+01 -6.1585E+00
            -1.7622E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1620.32192522438        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:       88
 NPARAMETR:  1.0320E+00  1.0816E+00  1.0352E+00  9.5460E-01  1.0494E+00  8.8175E-01  1.1251E+00  9.8840E-01  8.8257E-01  9.5943E-01
             1.4192E+00
 PARAMETER:  1.3147E-01  1.7844E-01  1.3460E-01  5.3541E-02  1.4819E-01 -2.5847E-02  2.1789E-01  8.8330E-02 -2.4915E-02  5.8581E-02
             4.5010E-01
 GRADIENT:   3.1731E+02 -1.7575E+01 -6.6686E+00 -1.5486E+01  2.9559E+01 -8.6242E+00 -5.4256E+00 -7.9060E-01  2.0177E+00 -2.3803E+00
             2.1875E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1621.57642400594        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:      223
 NPARAMETR:  1.0270E+00  1.0920E+00  1.0968E+00  1.0037E+00  1.0386E+00  9.4101E-01  1.3795E+00  1.0158E+00  7.9893E-01  9.2798E-01
             1.4276E+00
 PARAMETER:  1.2662E-01  1.8798E-01  1.9239E-01  1.0366E-01  1.3785E-01  3.9197E-02  4.2171E-01  1.1572E-01 -1.2448E-01  2.5257E-02
             4.5602E-01
 GRADIENT:   1.6709E+01  2.5204E+01  8.2838E+00  7.4531E+00 -4.1800E+00  1.4548E+00  1.4133E+01 -3.5147E+00 -3.8744E+00 -4.0958E+00
             1.8199E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1623.58324745232        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      399
 NPARAMETR:  1.0164E+00  1.0578E+00  1.3225E+00  1.0193E+00  1.1178E+00  9.3460E-01  1.2354E+00  1.4430E+00  8.4253E-01  1.0379E+00
             1.3601E+00
 PARAMETER:  1.1627E-01  1.5621E-01  3.7955E-01  1.1916E-01  2.1139E-01  3.2358E-02  3.1137E-01  4.6674E-01 -7.1341E-02  1.3721E-01
             4.0757E-01
 GRADIENT:  -5.4038E+00  2.1896E-02 -3.4140E-01 -8.9646E-01  3.0784E-01 -1.0750E+00 -1.9617E-01  3.0977E-01 -3.5126E-01  1.7885E-01
             9.4310E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1623.75540196455        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      579
 NPARAMETR:  1.0153E+00  8.7343E-01  1.4106E+00  1.1417E+00  1.0744E+00  9.3333E-01  1.4450E+00  1.3968E+00  7.7600E-01  1.0196E+00
             1.3605E+00
 PARAMETER:  1.1521E-01 -3.5327E-02  4.4399E-01  2.3254E-01  1.7173E-01  3.0999E-02  4.6814E-01  4.3418E-01 -1.5360E-01  1.1944E-01
             4.0789E-01
 GRADIENT:  -4.4093E+00  4.5675E+00 -1.4150E+00  8.8413E+00 -6.3447E-02 -9.5681E-01 -7.5095E-01  3.2489E-01 -4.8853E-02  2.6225E-01
             1.3781E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1624.12665355258        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      757
 NPARAMETR:  1.0142E+00  6.3432E-01  2.0303E+00  1.3047E+00  1.1442E+00  9.3159E-01  1.8342E+00  1.8588E+00  7.0018E-01  1.0918E+00
             1.3666E+00
 PARAMETER:  1.1414E-01 -3.5520E-01  8.0818E-01  3.6600E-01  2.3468E-01  2.9134E-02  7.0660E-01  7.1995E-01 -2.5641E-01  1.8785E-01
             4.1232E-01
 GRADIENT:  -5.5787E-01  4.3465E+00  3.1540E+00  2.8711E+00 -5.0322E+00 -1.2304E-01 -1.5036E-02 -8.7989E-01 -5.4044E-01  8.8409E-01
             1.6594E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1624.53815073984        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      936
 NPARAMETR:  1.0147E+00  3.5824E-01  3.0413E+00  1.5085E+00  1.2145E+00  9.3121E-01  2.6903E+00  2.5051E+00  6.4352E-01  1.1176E+00
             1.3530E+00
 PARAMETER:  1.1462E-01 -9.2656E-01  1.2123E+00  5.1113E-01  2.9431E-01  2.8732E-02  1.0897E+00  1.0183E+00 -3.4081E-01  2.1121E-01
             4.0229E-01
 GRADIENT:   5.6348E+00  9.5231E+00  2.8371E+00  3.1432E+01 -3.5982E+00  1.5257E+00  1.7426E+00 -1.7603E+00 -3.0529E+00 -2.6307E+00
            -6.0381E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1625.94902507754        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1115
 NPARAMETR:  1.0098E+00  1.6172E-01  3.5497E+00  1.6202E+00  1.2282E+00  9.2270E-01  3.9034E+00  2.8081E+00  6.3730E-01  1.1530E+00
             1.3730E+00
 PARAMETER:  1.0975E-01 -1.7219E+00  1.3669E+00  5.8253E-01  3.0553E-01  1.9544E-02  1.4618E+00  1.1325E+00 -3.5051E-01  2.4236E-01
             4.1701E-01
 GRADIENT:  -2.3167E+00  3.3767E-01  9.1370E-01 -1.2687E+01 -2.1582E+00 -3.7701E-01  2.0491E-01 -4.8808E-01  7.3333E-01  5.9858E-01
             2.5916E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1626.10099736543        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1291
 NPARAMETR:  1.0101E+00  8.8643E-02  3.9838E+00  1.6867E+00  1.2614E+00  9.2220E-01  5.0214E+00  3.0718E+00  6.2706E-01  1.1765E+00
             1.3677E+00
 PARAMETER:  1.1000E-01 -2.3231E+00  1.4822E+00  6.2280E-01  3.3223E-01  1.9008E-02  1.7137E+00  1.2223E+00 -3.6671E-01  2.6251E-01
             4.1310E-01
 GRADIENT:  -2.1952E+00  4.8164E-01 -2.0378E+00  3.5446E+01  1.0688E+00 -3.6335E-01 -2.6049E+00  4.6861E-01  9.9371E-01  9.8497E-02
            -2.9387E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1626.47625550485        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1468
 NPARAMETR:  1.0106E+00  6.8453E-02  4.6059E+00  1.6909E+00  1.2954E+00  9.2068E-01  5.5632E+00  3.3594E+00  6.1739E-01  1.2084E+00
             1.3751E+00
 PARAMETER:  1.1053E-01 -2.5816E+00  1.6273E+00  6.2529E-01  3.5882E-01  1.7356E-02  1.8162E+00  1.3118E+00 -3.8226E-01  2.8933E-01
             4.1852E-01
 GRADIENT:  -1.9954E+00  1.1250E-01 -9.7651E-01 -7.3268E+00  1.7623E-01 -3.4433E-01  9.7973E-01  9.0261E-01 -7.9789E-01  6.2756E-01
             2.0120E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1626.52661834133        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:     1631
 NPARAMETR:  1.0125E+00  6.2545E-02  5.4562E+00  1.7027E+00  1.3331E+00  9.2101E-01  5.7696E+00  3.6745E+00  6.1317E-01  1.2227E+00
             1.3716E+00
 PARAMETER:  1.1243E-01 -2.6719E+00  1.7967E+00  6.3223E-01  3.8747E-01  1.7720E-02  1.8526E+00  1.4014E+00 -3.8911E-01  3.0105E-01
             4.1597E-01
 GRADIENT:   4.6151E-01  4.4059E-01 -1.2573E-02 -2.1785E-01  1.7482E-01  1.1802E-01  9.8464E-01 -5.7581E-02 -7.8969E-01 -4.4114E-01
            -6.2444E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1631
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.0695E-03  3.2378E-02 -4.6828E-02 -3.8479E-02 -5.6524E-02
 SE:             2.9636E-02  1.4139E-02  1.7399E-02  2.4954E-02  1.9192E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7121E-01  2.2021E-02  7.1139E-03  1.2307E-01  3.2275E-03

 ETASHRINKSD(%)  7.1637E-01  5.2632E+01  4.1712E+01  1.6401E+01  3.5705E+01
 ETASHRINKVR(%)  1.4276E+00  7.7563E+01  6.6026E+01  3.0113E+01  5.8662E+01
 EBVSHRINKSD(%)  9.4973E-01  6.7241E+01  5.2300E+01  1.1120E+01  3.0197E+01
 EBVSHRINKVR(%)  1.8904E+00  8.9268E+01  7.7247E+01  2.1003E+01  5.1275E+01
 RELATIVEINF(%)  9.7754E+01  5.6575E+00  1.2654E+01  3.9954E+01  2.7331E+01
 EPSSHRINKSD(%)  4.2329E+01
 EPSSHRINKVR(%)  6.6740E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1626.5266183413264     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -891.37579177758823     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.36
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.89
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1626.527       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  6.25E-02  5.46E+00  1.70E+00  1.33E+00  9.21E-01  5.77E+00  3.67E+00  6.13E-01  1.22E+00  1.37E+00
 


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
+        1.26E+03
 
 TH 2
+        5.38E+01  6.39E+04
 
 TH 3
+       -6.94E-01 -5.36E+00  1.12E+00
 
 TH 4
+       -3.91E+01 -1.00E+04 -1.64E+00  2.70E+03
 
 TH 5
+       -2.70E+00  1.54E+02 -3.11E+00 -8.01E+01  2.55E+02
 
 TH 6
+       -1.04E+00  1.35E+00 -1.60E-01 -1.36E+01  5.65E-02  2.27E+02
 
 TH 7
+        1.06E+00  1.14E+03 -2.15E-01 -2.07E+02  6.17E+00  8.88E-01  2.32E+01
 
 TH 8
+       -3.00E-01 -1.40E+01 -2.40E+00  7.58E-01 -1.62E+01  8.38E-01  2.83E-01  8.05E+00
 
 TH 9
+       -1.27E+00 -2.27E+03  2.27E+00  3.09E+02 -4.94E+01 -9.96E+00 -6.82E+02  6.11E+00  8.14E+02
 
 TH10
+        1.56E+00 -1.85E+02  2.00E+00  3.92E+01 -3.70E+01 -1.17E+00 -4.64E+00 -3.01E+00  5.79E+01  5.93E+01
 
 TH11
+       -1.36E+01 -4.91E+01  1.33E+00 -1.59E+01  1.17E+01  1.42E+00 -9.68E-01 -5.96E+00  1.57E+01  1.78E+01  1.35E+02
 
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
 #CPUT: Total CPU Time in Seconds,       33.328
Stop Time:
Wed Sep 29 16:53:36 CDT 2021
