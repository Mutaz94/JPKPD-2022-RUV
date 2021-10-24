Sun Oct 24 00:21:48 CDT 2021
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
$DATA ../../../../data/SD3/SL3/dat97.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m97.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2002.84090446902        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.4466E+02 -5.9499E+01 -5.2468E+01  3.1919E+01  9.4757E+01  3.3785E+01 -1.1444E+01  6.7982E+00  2.8907E+01 -1.4169E+01
            -1.1424E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2019.47291478848        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      171
 NPARAMETR:  1.0037E+00  1.1361E+00  1.0736E+00  9.7689E-01  1.0447E+00  1.0861E+00  1.0654E+00  9.6264E-01  8.9064E-01  1.0427E+00
             1.1549E+00
 PARAMETER:  1.0368E-01  2.2763E-01  1.7102E-01  7.6617E-02  1.4370E-01  1.8261E-01  1.6334E-01  6.1925E-02 -1.5818E-02  1.4179E-01
             2.4404E-01
 GRADIENT:   1.3043E+02 -4.4442E+00 -1.5664E+01  2.8260E+01  2.5363E+01  8.9797E+00 -6.4934E+00  3.3998E+00  4.3311E+00 -1.0820E+01
             1.3641E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2022.10633637519        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      349
 NPARAMETR:  9.9218E-01  1.1351E+00  1.2637E+00  9.9451E-01  1.1160E+00  1.0629E+00  1.1450E+00  8.6048E-01  8.3038E-01  1.2256E+00
             1.1423E+00
 PARAMETER:  9.2148E-02  2.2676E-01  3.3404E-01  9.4493E-02  2.0975E-01  1.6097E-01  2.3540E-01 -5.0262E-02 -8.5873E-02  3.0346E-01
             2.3305E-01
 GRADIENT:   1.1530E+02  2.1958E+01  8.5337E-01  3.4029E+01  4.6756E+00  3.1888E+00 -8.4851E-01 -8.4124E-01 -2.3784E+00  5.6529E-01
             1.7373E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2026.33593527111        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      526
 NPARAMETR:  9.3468E-01  1.2247E+00  1.0218E+00  9.1388E-01  1.0583E+00  1.0436E+00  1.0596E+00  5.8479E-01  9.0001E-01  1.1622E+00
             1.1332E+00
 PARAMETER:  3.2453E-02  3.0273E-01  1.2161E-01  9.9390E-03  1.5663E-01  1.4272E-01  1.5791E-01 -4.3649E-01 -5.3474E-03  2.5032E-01
             2.2506E-01
 GRADIENT:  -3.4862E+00  5.1592E+00 -2.5186E-01  1.0703E+01 -3.2388E+00  2.5928E+00 -1.4369E+00  2.0741E-01  5.5558E-01  8.6268E-01
            -1.1380E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2026.97717923621        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      702
 NPARAMETR:  9.3886E-01  1.5731E+00  7.5894E-01  6.8908E-01  1.1313E+00  1.0352E+00  8.8823E-01  2.6509E-01  1.0698E+00  1.1813E+00
             1.1339E+00
 PARAMETER:  3.6908E-02  5.5302E-01 -1.7583E-01 -2.7239E-01  2.2339E-01  1.3459E-01 -1.8526E-02 -1.2277E+00  1.6746E-01  2.6659E-01
             2.2570E-01
 GRADIENT:   1.4397E+00  1.1502E+01 -1.3142E+00  1.0563E+01  1.5640E+00 -1.3166E+00 -3.0440E-01  9.5096E-02 -1.5242E+00 -8.2781E-01
             1.3295E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2027.07247741317        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      879
 NPARAMETR:  9.3889E-01  1.6869E+00  6.9925E-01  6.1409E-01  1.1697E+00  1.0363E+00  8.3068E-01  1.7809E-01  1.1849E+00  1.2105E+00
             1.1309E+00
 PARAMETER:  3.6943E-02  6.2292E-01 -2.5774E-01 -3.8761E-01  2.5675E-01  1.3563E-01 -8.5516E-02 -1.6255E+00  2.6963E-01  2.9102E-01
             2.2301E-01
 GRADIENT:   8.9747E-01  1.1016E+01  8.4520E-01  7.3440E+00 -2.4365E+00 -1.0263E+00 -1.2781E+00  4.0956E-02 -1.6090E-01 -2.6463E-01
            -7.8322E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2027.15184296905        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1060
 NPARAMETR:  9.3865E-01  1.7172E+00  6.7278E-01  5.8759E-01  1.1824E+00  1.0388E+00  8.2463E-01  6.8975E-02  1.2139E+00  1.2173E+00
             1.1312E+00
 PARAMETER:  3.6684E-02  6.4067E-01 -2.9634E-01 -4.3173E-01  2.6754E-01  1.3808E-01 -9.2819E-02 -2.5740E+00  2.9388E-01  2.9663E-01
             2.2329E-01
 GRADIENT:   2.2235E-01 -1.0687E+00  4.0657E-01  1.1425E+00 -4.3032E-01 -5.8375E-02 -7.3678E-03  6.7351E-03 -1.3861E-02 -9.3151E-02
            -6.0306E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2027.15832847896        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1242
 NPARAMETR:  9.3893E-01  1.7184E+00  6.6880E-01  5.8503E-01  1.1830E+00  1.0399E+00  8.2422E-01  1.9680E-02  1.2155E+00  1.2179E+00
             1.1312E+00
 PARAMETER:  3.6983E-02  6.4142E-01 -3.0227E-01 -4.3610E-01  2.6807E-01  1.3917E-01 -9.3312E-02 -3.8282E+00  2.9515E-01  2.9716E-01
             2.2324E-01
 GRADIENT:   8.2101E-01 -4.3626E+00  7.1205E-02 -1.1391E-01  2.2244E-01  3.7878E-01  5.9945E-02  5.8262E-04 -1.1516E-02  8.2889E-02
             5.9274E-02

0ITERATION NO.:   38    OBJECTIVE VALUE:  -2027.15879673192        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:     1338
 NPARAMETR:  9.3912E-01  1.7181E+00  6.6851E-01  5.8500E-01  1.1827E+00  1.0403E+00  8.2397E-01  1.0000E-02  1.2159E+00  1.2175E+00
             1.1311E+00
 PARAMETER:  3.7190E-02  6.4119E-01 -3.0270E-01 -4.3614E-01  2.6780E-01  1.3949E-01 -9.3623E-02 -4.9489E+00  2.9547E-01  2.9681E-01
             2.2317E-01
 GRADIENT:   1.2453E+00 -4.8726E+00  5.5447E-02 -3.2803E-01  2.1555E-01  5.0629E-01  1.4226E-02  0.0000E+00  1.8893E-02  7.5482E-02
             1.4821E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1338
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.4646E-04 -2.2618E-02 -2.5509E-04  2.0814E-02 -3.1663E-02
 SE:             2.9831E-02  2.3752E-02  9.1453E-05  2.1249E-02  2.3573E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8806E-01  3.4096E-01  5.2813E-03  3.2730E-01  1.7921E-01

 ETASHRINKSD(%)  6.2886E-02  2.0429E+01  9.9694E+01  2.8815E+01  2.1026E+01
 ETASHRINKVR(%)  1.2573E-01  3.6685E+01  9.9999E+01  4.9326E+01  3.7631E+01
 EBVSHRINKSD(%)  4.0439E-01  1.9743E+01  9.9740E+01  3.2012E+01  1.7833E+01
 EBVSHRINKVR(%)  8.0714E-01  3.5588E+01  9.9999E+01  5.3776E+01  3.2485E+01
 RELATIVEINF(%)  9.9035E+01  2.9868E+00  6.9928E-05  1.9266E+00  1.6418E+01
 EPSSHRINKSD(%)  3.2722E+01
 EPSSHRINKVR(%)  5.4736E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2027.1587967319178     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1108.2202635272452     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.94
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2027.159       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.39E-01  1.72E+00  6.69E-01  5.85E-01  1.18E+00  1.04E+00  8.24E-01  1.00E-02  1.22E+00  1.22E+00  1.13E+00
 


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
 #CPUT: Total CPU Time in Seconds,       45.394
Stop Time:
Sun Oct 24 00:21:57 CDT 2021
