Thu Sep 30 23:09:14 CDT 2021
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
$DATA ../../../../data/SD1/TD2/dat4.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m4.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3292.70614825717        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0273E+02 -1.7532E+01 -2.8420E+01  2.3589E+02  1.8019E+02  2.7199E+01 -3.2305E+01 -6.0029E+01 -9.7060E+01 -3.3301E+01
            -9.1609E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3448.35320236262        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      128
 NPARAMETR:  1.0230E+00  1.2360E+00  9.6174E-01  8.4234E-01  9.8775E-01  1.1380E+00  1.0944E+00  1.1353E+00  1.4300E+00  1.1691E+00
             1.4947E+00
 PARAMETER:  1.2272E-01  3.1190E-01  6.0988E-02 -7.1568E-02  8.7679E-02  2.2924E-01  1.9016E-01  2.2687E-01  4.5765E-01  2.5625E-01
             5.0189E-01
 GRADIENT:   7.1335E+01  2.7469E+01 -1.9261E+01 -4.9860E+00 -7.4829E+01  1.8603E+01  2.1435E+01  9.4347E-01  1.9256E+01 -1.2405E+00
             1.6740E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3464.75794454301        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      307
 NPARAMETR:  1.0152E+00  1.5605E+00  1.4501E+00  6.7479E-01  1.2859E+00  1.1092E+00  6.9482E-01  2.1174E+00  1.4661E+00  1.3720E+00
             1.4184E+00
 PARAMETER:  1.1510E-01  5.4503E-01  4.7160E-01 -2.9335E-01  3.5146E-01  2.0367E-01 -2.6410E-01  8.5019E-01  4.8260E-01  4.1630E-01
             4.4951E-01
 GRADIENT:   6.1338E+01  9.8023E+01  2.2875E+01 -3.5008E+00 -5.6174E+01  9.7236E+00 -8.1156E+00 -1.0855E+01  2.5559E+00 -1.1099E+01
             5.4842E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3472.52100604339        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      491
 NPARAMETR:  9.6823E-01  1.2604E+00  1.3509E+00  8.1878E-01  1.1482E+00  1.1265E+00  9.9394E-01  2.2200E+00  1.2232E+00  1.4394E+00
             1.3706E+00
 PARAMETER:  6.7710E-02  3.3145E-01  4.0080E-01 -9.9938E-02  2.3819E-01  2.1910E-01  9.3922E-02  8.9752E-01  3.0148E-01  4.6426E-01
             4.1526E-01
 GRADIENT:  -2.1389E+01 -3.4012E+01 -4.0673E+00 -3.0663E+01 -2.5545E+01  1.7243E+01  1.7232E+01 -5.2731E+00  7.9767E+00  2.6156E+01
             4.5275E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3475.24507196372        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:      681             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7002E-01  1.2739E+00  1.3869E+00  8.2259E-01  1.1690E+00  1.1163E+00  9.5721E-01  2.3491E+00  1.2022E+00  1.4111E+00
             1.3640E+00
 PARAMETER:  6.9559E-02  3.4208E-01  4.2706E-01 -9.5302E-02  2.5618E-01  2.1001E-01  5.6266E-02  9.5402E-01  2.8411E-01  4.4435E-01
             4.1042E-01
 GRADIENT:   2.2876E+02  1.8288E+02  4.9756E+00  1.8356E+01  4.5673E+01  9.0596E+01  1.7609E+01  5.6426E+00  2.0773E+01  3.7822E+01
             4.2528E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3476.43775943174        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:      868            RESET HESSIAN, TYPE II
 NPARAMETR:  9.8088E-01  1.2728E+00  1.3869E+00  8.2218E-01  1.1690E+00  1.0770E+00  8.3708E-01  2.3380E+00  1.2022E+00  1.4142E+00
             1.3640E+00
 PARAMETER:  8.0690E-02  3.4123E-01  4.2706E-01 -9.5802E-02  2.5618E-01  1.7415E-01 -7.7837E-02  9.4929E-01  2.8411E-01  4.4656E-01
             4.1042E-01
 GRADIENT:   2.4839E+02  1.7440E+02  6.2131E+00  1.6475E+01  4.7993E+01  6.2039E+01  5.0192E+00  4.3151E+00  1.6717E+01  3.6351E+01
             3.7386E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3476.56802544593        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1054
 NPARAMETR:  9.8068E-01  1.2734E+00  1.3876E+00  8.2228E-01  1.1687E+00  1.0754E+00  8.3132E-01  2.4918E+00  1.2026E+00  1.4150E+00
             1.3633E+00
 PARAMETER:  8.0488E-02  3.4166E-01  4.2759E-01 -9.5678E-02  2.5586E-01  1.7272E-01 -8.4741E-02  1.0130E+00  2.8446E-01  4.4711E-01
             4.0990E-01
 GRADIENT:   5.6282E-01 -3.3624E+01 -7.7391E+00 -2.4093E+01 -1.9040E+01 -1.9065E-01  1.3597E+00  1.8853E+00  4.4326E+00  1.7483E+01
             3.2001E+01

0ITERATION NO.:   32    OBJECTIVE VALUE:  -3476.58264397071        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1111
 NPARAMETR:  9.8068E-01  1.2734E+00  1.3876E+00  8.2228E-01  1.1687E+00  1.0759E+00  8.3132E-01  2.4418E+00  1.2026E+00  1.4150E+00
             1.3633E+00
 PARAMETER:  8.0488E-02  3.4167E-01  4.2759E-01 -9.5678E-02  2.5586E-01  1.7317E-01 -8.4740E-02  9.9273E-01  2.8447E-01  4.4711E-01
             4.0989E-01
 GRADIENT:   5.5312E-01 -3.3348E+01 -6.3884E+00 -2.3582E+01 -1.8127E+01 -9.5299E-03  1.1841E+00 -3.1850E-01  3.7995E+00  1.7697E+01
             3.0993E+01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1111
 NO. OF SIG. DIGITS IN FINAL EST.:  2.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.6050E-04 -3.3213E-02 -2.1195E-02  3.7513E-02 -2.3874E-02
 SE:             2.9874E-02  1.9529E-02  2.0402E-02  2.5803E-02  2.4019E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7702E-01  8.9000E-02  2.9888E-01  1.4600E-01  3.2024E-01

 ETASHRINKSD(%)  1.0000E-10  3.4576E+01  3.1650E+01  1.3556E+01  1.9532E+01
 ETASHRINKVR(%)  1.0000E-10  5.7197E+01  5.3282E+01  2.5274E+01  3.5249E+01
 EBVSHRINKSD(%)  4.2503E-01  3.4634E+01  3.4297E+01  1.4076E+01  1.2282E+01
 EBVSHRINKVR(%)  8.4825E-01  5.7272E+01  5.6831E+01  2.6171E+01  2.3055E+01
 RELATIVEINF(%)  9.9148E+01  1.7499E+01  3.6386E+01  3.4536E+01  4.6718E+01
 EPSSHRINKSD(%)  2.2007E+01
 EPSSHRINKVR(%)  3.9171E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3476.5826439707116     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1822.4932842023009     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.05
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3476.583       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.81E-01  1.27E+00  1.39E+00  8.22E-01  1.17E+00  1.08E+00  8.31E-01  2.44E+00  1.20E+00  1.41E+00  1.36E+00
 


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
 #CPUT: Total CPU Time in Seconds,       21.059
Stop Time:
Thu Sep 30 23:09:36 CDT 2021
