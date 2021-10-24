Sat Oct 23 18:06:40 CDT 2021
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
$DATA ../../../../data/SD2/A3/dat78.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      800
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

 TOT. NO. OF OBS RECS:      700
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
 RAW OUTPUT FILE (FILE): m78.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -366.848594147945        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6136E+02  1.5890E+02  2.3789E+02 -1.2254E+00  2.9055E+02  4.3691E+01 -1.8545E+02 -2.1077E+02 -1.6776E+01 -2.1754E+02
            -4.5058E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2120.69467353091        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.8221E-01  1.0362E+00  9.7277E-01  1.0912E+00  9.3498E-01  8.3638E-01  1.0809E+00  9.2628E-01  8.6401E-01  8.5597E-01
             3.0989E+00
 PARAMETER:  8.2049E-02  1.3560E-01  7.2393E-02  1.8724E-01  3.2767E-02 -7.8668E-02  1.7783E-01  2.3424E-02 -4.6175E-02 -5.5522E-02
             1.2311E+00
 GRADIENT:  -5.9281E+00  2.1014E+01 -3.4996E+00  4.3370E+01 -6.7375E+00 -3.7031E+01 -9.6150E+00  7.0732E+00  3.6141E+00  1.2824E+01
             1.1659E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2127.56313469771        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.7787E-01  8.2596E-01  7.0958E-01  1.2052E+00  7.1571E-01  9.0331E-01  1.4512E+00  5.6523E-01  7.4817E-01  6.2524E-01
             3.0259E+00
 PARAMETER:  7.7622E-02 -9.1206E-02 -2.4309E-01  2.8662E-01 -2.3448E-01 -1.6885E-03  4.7239E-01 -4.7053E-01 -1.9012E-01 -3.6963E-01
             1.2072E+00
 GRADIENT:  -9.8828E+00  2.4532E+01 -3.2748E+01  1.1711E+02  4.7555E+01 -9.8281E+00  6.8087E+00  3.1252E+00 -6.0248E+00  2.1080E+00
             9.1047E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2133.98428023612        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      294
 NPARAMETR:  9.9104E-01  7.5831E-01  6.8796E-01  1.2074E+00  6.7598E-01  9.3591E-01  1.3935E+00  2.6484E-01  8.5420E-01  6.6268E-01
             2.8500E+00
 PARAMETER:  9.0995E-02 -1.7667E-01 -2.7402E-01  2.8848E-01 -2.9160E-01  3.3760E-02  4.3183E-01 -1.2286E+00 -5.7585E-02 -3.1147E-01
             1.1473E+00
 GRADIENT:  -3.6613E+00  2.3057E+00 -1.7293E+01  3.5171E+01  4.3173E+01  1.2485E+00 -1.5620E+00  4.0944E-01  7.5648E+00 -8.3861E+00
             2.2837E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2179.32415379279        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      475
 NPARAMETR:  9.8681E-01  3.4548E-01  2.9869E-01  1.2430E+00  2.9619E-01  9.3516E-01  1.1729E+00  1.0000E-02  1.1155E+00  8.6297E-01
             2.3542E+00
 PARAMETER:  8.6717E-02 -9.6282E-01 -1.1084E+00  3.1749E-01 -1.1168E+00  3.2962E-02  2.5946E-01 -6.1133E+00  2.0930E-01 -4.7370E-02
             9.5619E-01
 GRADIENT:  -5.8191E-01  1.3095E+01  3.5441E+01  7.6856E+01 -2.3156E+01  1.8047E+00  7.6578E+00  0.0000E+00  5.2945E+00  7.8749E+00
            -1.2191E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2194.78215015993        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      650
 NPARAMETR:  9.9218E-01  2.7981E-01  2.0025E-01  1.0712E+00  2.2763E-01  9.2495E-01  7.1549E-01  1.0000E-02  1.1980E+00  8.5387E-01
             2.5149E+00
 PARAMETER:  9.2150E-02 -1.1736E+00 -1.5082E+00  1.6882E-01 -1.3800E+00  2.1987E-02 -2.3479E-01 -7.5996E+00  2.8063E-01 -5.7978E-02
             1.0222E+00
 GRADIENT:   3.1022E+00 -6.6544E-01 -1.2993E+00 -8.7249E+00  3.1712E+00  8.6157E-01  2.0421E+00  0.0000E+00  1.9888E+00 -7.9964E-01
             2.7547E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2195.12936387381        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      826
 NPARAMETR:  9.9104E-01  2.8211E-01  2.0459E-01  1.0823E+00  2.3033E-01  9.2157E-01  4.8034E-01  1.0000E-02  1.1816E+00  8.7294E-01
             2.5252E+00
 PARAMETER:  9.1001E-02 -1.1655E+00 -1.4867E+00  1.7913E-01 -1.3683E+00  1.8321E-02 -6.3325E-01 -7.3527E+00  2.6691E-01 -3.5885E-02
             1.0263E+00
 GRADIENT:  -2.4012E-01 -4.1752E-02  1.3759E-02 -2.3953E-02 -3.0283E-02 -4.8078E-02  1.8912E-02  0.0000E+00 -7.1166E-02  2.0459E-02
             1.1752E-01

0ITERATION NO.:   32    OBJECTIVE VALUE:  -2195.12941525148        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:      883
 NPARAMETR:  9.9115E-01  2.8215E-01  2.0461E-01  1.0824E+00  2.3035E-01  9.2167E-01  4.7638E-01  1.0000E-02  1.1819E+00  8.7315E-01
             2.5252E+00
 PARAMETER:  9.1106E-02 -1.1653E+00 -1.4867E+00  1.7916E-01 -1.3682E+00  1.8435E-02 -6.4154E-01 -7.3482E+00  2.6714E-01 -3.5644E-02
             1.0263E+00
 GRADIENT:   3.5748E-02  1.4708E-02  3.9191E-03 -5.4681E-03 -1.5448E-02 -1.6397E-03  3.4696E-03  0.0000E+00  1.9070E-03  6.1327E-03
             4.8363E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      883
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5445E-03 -2.0826E-04  1.7813E-04 -6.9638E-03  3.4760E-04
 SE:             2.9207E-02  8.0981E-03  2.1329E-04  2.7427E-02  2.7704E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5783E-01  9.7948E-01  4.0365E-01  7.9957E-01  9.8999E-01

 ETASHRINKSD(%)  2.1519E+00  7.2870E+01  9.9285E+01  8.1161E+00  7.1872E+00
 ETASHRINKVR(%)  4.2576E+00  9.2640E+01  9.9995E+01  1.5573E+01  1.3858E+01
 EBVSHRINKSD(%)  2.0706E+00  7.3060E+01  9.9333E+01  6.4120E+00  7.6471E+00
 EBVSHRINKVR(%)  4.0984E+00  9.2743E+01  9.9996E+01  1.2413E+01  1.4709E+01
 RELATIVEINF(%)  9.5865E+01  1.5951E+00  4.6253E-04  3.7227E+01  6.2702E+00
 EPSSHRINKSD(%)  2.1005E+01
 EPSSHRINKVR(%)  3.7597E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          700
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1286.5139464865417     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2195.1294152514847     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -908.61546876494299     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.30
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2195.129       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.91E-01  2.82E-01  2.05E-01  1.08E+00  2.30E-01  9.22E-01  4.76E-01  1.00E-02  1.18E+00  8.73E-01  2.53E+00
 


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
 #CPUT: Total CPU Time in Seconds,       53.105
Stop Time:
Sat Oct 23 18:06:51 CDT 2021
