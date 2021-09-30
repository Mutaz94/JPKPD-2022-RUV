Wed Sep 29 12:52:16 CDT 2021
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
$DATA ../../../../data/spa/A2/dat51.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m51.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1234.11899243194        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.2360E+02  4.8813E+01  3.4866E+01  4.0161E+01  6.0196E+01  5.4307E+01 -1.0932E+00 -9.6429E+00 -1.7856E+01 -5.3174E+01
            -6.9362E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1352.07325239084        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0869E+00  9.5197E-01  9.0936E-01  1.1098E+00  8.9491E-01  8.6221E-01  9.2362E-01  9.7820E-01  9.8861E-01  1.0676E+00
             3.6654E+00
 PARAMETER:  1.8331E-01  5.0783E-02  4.9874E-03  2.0419E-01 -1.1037E-02 -4.8254E-02  2.0543E-02  7.7958E-02  8.8547E-02  1.6539E-01
             1.3989E+00
 GRADIENT:   3.3602E+02  1.4655E+00 -1.1550E+01  2.2892E+01 -2.1183E+01 -4.5544E+01  1.0561E+01  8.6259E+00  2.1680E+01  3.0982E+01
             2.0791E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1389.23848717530        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0467E+00  6.3677E-01  3.2002E-01  1.2254E+00  3.9652E-01  9.1560E-01  2.6623E-01  2.5422E-01  1.1963E+00  6.5993E-01
             2.7457E+00
 PARAMETER:  1.4565E-01 -3.5135E-01 -1.0394E+00  3.0324E-01 -8.2503E-01  1.1821E-02 -1.2234E+00 -1.2696E+00  2.7923E-01 -3.1563E-01
             1.1100E+00
 GRADIENT:   2.5057E+02  5.4668E+01 -6.6505E+00  1.6697E+02 -1.3137E+01 -2.4544E+01 -3.8016E-02  1.3405E+00  4.7194E+01  8.3381E+00
             1.1965E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1438.26760566864        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      339
 NPARAMETR:  9.7237E-01  8.0222E-01  7.8680E-01  1.0765E+00  7.5668E-01  9.2579E-01  3.3400E-01  3.7614E-01  1.0864E+00  9.5944E-01
             2.1699E+00
 PARAMETER:  7.1979E-02 -1.2037E-01 -1.3978E-01  1.7375E-01 -1.7882E-01  2.2890E-02 -9.9660E-01 -8.7779E-01  1.8284E-01  5.8593E-02
             8.7470E-01
 GRADIENT:   4.0667E+01 -1.6427E+01  1.2889E+01 -6.2310E+01 -2.4394E+01 -7.1430E+00 -2.6944E-01  7.3894E-01  7.4820E+00  5.8616E+00
             2.8122E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1444.02161311725        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      519
 NPARAMETR:  9.5011E-01  6.1909E-01  6.7019E-01  1.2395E+00  6.3401E-01  9.5760E-01  1.3276E+00  6.0538E-02  8.9019E-01  8.3414E-01
             2.0119E+00
 PARAMETER:  4.8822E-02 -3.7951E-01 -3.0019E-01  3.1469E-01 -3.5569E-01  5.6673E-02  3.8337E-01 -2.7045E+00 -1.6317E-02 -8.1359E-02
             7.9908E-01
 GRADIENT:  -1.8607E+01  1.5367E+01 -3.8022E+00  3.6604E+01  1.3614E+00  4.2304E+00  1.5872E+00  3.1658E-02  2.4154E+00  3.2262E+00
             1.7312E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1445.71647924319        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      694
 NPARAMETR:  9.5932E-01  4.5875E-01  5.6245E-01  1.2828E+00  5.1946E-01  9.4459E-01  1.6399E+00  3.8435E-02  8.3432E-01  7.1692E-01
             1.9657E+00
 PARAMETER:  5.8468E-02 -6.7925E-01 -4.7546E-01  3.4907E-01 -5.5496E-01  4.3000E-02  5.9466E-01 -3.1588E+00 -8.1139E-02 -2.3279E-01
             7.7583E-01
 GRADIENT:   5.7641E+00  3.1546E+00 -4.4767E+00  1.3409E+01  3.6912E+00 -1.0055E+00 -1.3494E+00  1.5713E-02 -1.7389E+00 -2.5018E+00
            -6.6045E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1446.40324707691        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      869
 NPARAMETR:  9.5320E-01  3.0446E-01  6.2538E-01  1.3759E+00  5.2082E-01  9.3854E-01  2.2148E+00  1.1676E-02  8.0522E-01  8.0480E-01
             1.9824E+00
 PARAMETER:  5.2071E-02 -1.0892E+00 -3.6940E-01  4.1909E-01 -5.5235E-01  3.6570E-02  8.9515E-01 -4.3503E+00 -1.1664E-01 -1.1717E-01
             7.8430E-01
 GRADIENT:   1.1803E+00  3.2167E+00  7.7741E+00  9.6031E+00 -1.1959E+01 -1.3061E+00  3.3595E-02  1.3613E-03  1.9968E-01  3.0330E-01
             1.8178E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1446.51397661580        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1050             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5086E-01  2.4153E-01  6.2968E-01  1.4009E+00  5.1509E-01  9.4039E-01  2.5957E+00  1.0000E-02  7.9234E-01  8.1316E-01
             1.9789E+00
 PARAMETER:  4.9612E-02 -1.3208E+00 -3.6254E-01  4.3711E-01 -5.6341E-01  3.8538E-02  1.0539E+00 -4.9136E+00 -1.3277E-01 -1.0682E-01
             7.8252E-01
 GRADIENT:   9.3376E+01  8.3511E+00  7.2720E+00  1.4747E+02  2.5819E+01  1.0652E+01  3.8291E+00  0.0000E+00  4.5664E+00  3.7748E-01
             6.4302E+00

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1446.51473804035        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:     1143
 NPARAMETR:  9.5076E-01  2.4204E-01  6.2965E-01  1.4016E+00  5.1501E-01  9.4029E-01  2.5986E+00  1.0000E-02  7.9128E-01  8.1307E-01
             1.9785E+00
 PARAMETER:  4.9502E-02 -1.3187E+00 -3.6259E-01  4.3760E-01 -5.6357E-01  3.8430E-02  1.0550E+00 -4.9136E+00 -1.3410E-01 -1.0694E-01
             7.8235E-01
 GRADIENT:   1.6505E-01  7.7046E-02  2.0032E-01 -9.7303E-01  1.0557E-01  1.6487E-02  4.0881E-02  0.0000E+00 -1.1484E-01  6.3240E-03
            -1.6374E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1143
 NO. OF SIG. DIGITS IN FINAL EST.:  2.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.2333E-04  2.1215E-02 -1.6209E-04 -1.8820E-02 -6.4586E-03
 SE:             2.9386E-02  1.3843E-02  2.0029E-04  2.5968E-02  2.1719E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8308E-01  1.2540E-01  4.1834E-01  4.6862E-01  7.6619E-01

 ETASHRINKSD(%)  1.5535E+00  5.3624E+01  9.9329E+01  1.3004E+01  2.7237E+01
 ETASHRINKVR(%)  3.0829E+00  7.8492E+01  9.9995E+01  2.4317E+01  4.7056E+01
 EBVSHRINKSD(%)  1.6907E+00  6.1808E+01  9.9275E+01  1.1758E+01  2.3514E+01
 EBVSHRINKVR(%)  3.3527E+00  8.5414E+01  9.9995E+01  2.2133E+01  4.1499E+01
 RELATIVEINF(%)  9.4834E+01  2.5409E+00  2.6473E-04  1.9491E+01  2.6059E+00
 EPSSHRINKSD(%)  3.6777E+01
 EPSSHRINKVR(%)  6.0029E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1446.5147380403509     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -711.36391147661277     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.32
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.63
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1446.515       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.51E-01  2.42E-01  6.30E-01  1.40E+00  5.15E-01  9.40E-01  2.60E+00  1.00E-02  7.91E-01  8.13E-01  1.98E+00
 


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
+        1.34E+03
 
 TH 2
+       -5.13E+01  6.94E+02
 
 TH 3
+        1.94E+01  2.34E+02  1.55E+03
 
 TH 4
+       -3.26E+01  2.89E+02 -1.78E+02  7.21E+02
 
 TH 5
+        2.51E+01 -6.04E+02 -2.30E+03 -1.15E+01  3.76E+03
 
 TH 6
+       -2.02E-01 -6.88E+00  7.48E+00 -8.48E+00 -1.15E+00  2.09E+02
 
 TH 7
+        1.39E+00  5.54E+01 -1.03E+01 -7.34E+00  8.33E+00  2.46E-01  8.02E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        5.18E+00 -9.50E+01  5.33E+01  6.55E+00 -4.15E+01  1.68E+00 -8.37E+00  0.00E+00  2.24E+02
 
 TH10
+       -3.32E+00  2.16E+01 -5.10E+01 -1.63E+01 -1.02E+01  1.03E-01  2.73E+00  0.00E+00  7.48E+00  8.66E+01
 
 TH11
+       -1.35E+01 -8.81E+00 -2.67E+01 -1.05E+01  9.18E+00  2.68E+00  1.36E-01  0.00E+00  1.28E+01  2.26E+01  6.53E+01
 
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
 #CPUT: Total CPU Time in Seconds,       21.012
Stop Time:
Wed Sep 29 12:52:47 CDT 2021
