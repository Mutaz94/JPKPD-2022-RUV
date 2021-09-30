Wed Sep 29 14:10:05 CDT 2021
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
$DATA ../../../../data/spa/S1/dat29.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m29.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1579.78196749711        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.8160E+02 -1.2269E+01 -9.0868E+00  1.4600E+01  2.5526E+01  1.1862E+01 -2.6808E+01  3.0188E-01 -2.2967E+01 -6.2881E-01
            -7.0764E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1593.40435164851        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.7201E-01  1.0832E+00  1.0572E+00  9.9001E-01  1.0447E+00  1.0968E+00  1.2157E+00  9.8383E-01  1.1292E+00  9.7217E-01
             1.1537E+00
 PARAMETER:  7.1610E-02  1.7993E-01  1.5561E-01  8.9959E-02  1.4368E-01  1.9240E-01  2.9529E-01  8.3695E-02  2.2151E-01  7.1778E-02
             2.4300E-01
 GRADIENT:  -1.4976E+01 -3.8627E+00  4.2939E+00 -4.7331E-01 -5.9824E-01  3.0804E+00 -4.8968E+00 -2.1143E-01  5.5466E+00 -2.6103E+00
            -2.9838E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1593.90695069113        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      354
 NPARAMETR:  9.8445E-01  1.1175E+00  9.6308E-01  9.6666E-01  1.0196E+00  1.0891E+00  1.3204E+00  9.3217E-01  1.1051E+00  9.4382E-01
             1.1500E+00
 PARAMETER:  8.4332E-02  2.1108E-01  6.2380E-02  6.6094E-02  1.1936E-01  1.8539E-01  3.7795E-01  2.9759E-02  1.9992E-01  4.2184E-02
             2.3979E-01
 GRADIENT:   7.1506E+00  6.2264E-01 -5.8164E-01  2.2593E+00  1.5431E+00  1.9190E-01  4.6769E+00  1.4639E+00  5.1673E+00  5.8673E-01
            -1.7531E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1594.57095595313        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      536
 NPARAMETR:  9.8191E-01  1.2750E+00  7.2596E-01  8.5872E-01  9.6964E-01  1.0885E+00  1.1976E+00  4.5567E-01  1.1181E+00  8.9750E-01
             1.1540E+00
 PARAMETER:  8.1747E-02  3.4294E-01 -2.2026E-01 -5.2307E-02  6.9169E-02  1.8481E-01  2.8034E-01 -6.8599E-01  2.1165E-01 -8.1371E-03
             2.4326E-01
 GRADIENT:  -1.5892E+00  8.0445E+00  1.8771E+00  5.5598E+00 -8.5801E+00 -6.7183E-01 -5.6361E-01  2.5120E-01 -1.1383E+00  1.9439E+00
            -3.3999E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1594.81128976831        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      717             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8420E-01  1.4303E+00  6.3222E-01  7.4707E-01  1.0113E+00  1.0921E+00  1.1039E+00  1.1386E-01  1.2185E+00  8.9294E-01
             1.1562E+00
 PARAMETER:  8.4075E-02  4.5788E-01 -3.5851E-01 -1.9160E-01  1.1119E-01  1.8812E-01  1.9882E-01 -2.0728E+00  2.9766E-01 -1.3233E-02
             2.4516E-01
 GRADIENT:   3.2803E+02  2.5334E+02  5.4905E+00  5.6943E+01  4.1216E+00  7.4378E+01  1.1815E+01  3.8593E-02  1.2658E+01 -3.6821E-02
             1.4970E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1594.82350586605        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      895
 NPARAMETR:  9.8351E-01  1.4319E+00  6.2986E-01  7.4795E-01  1.0114E+00  1.0911E+00  1.1035E+00  7.1858E-02  1.2193E+00  8.9694E-01
             1.1567E+00
 PARAMETER:  8.3369E-02  4.5903E-01 -3.6226E-01 -1.9042E-01  1.1136E-01  1.8720E-01  1.9848E-01 -2.5331E+00  2.9825E-01 -8.7629E-03
             2.4558E-01
 GRADIENT:   4.4187E-01 -4.3988E-01  4.2238E-01 -6.8321E-01 -7.0140E-01  1.0736E-01 -6.1021E-02  7.5083E-03 -1.6791E-01  1.6340E-02
            -5.8943E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1594.82425485248        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1071
 NPARAMETR:  9.8318E-01  1.4343E+00  6.2812E-01  7.4707E-01  1.0119E+00  1.0907E+00  1.1029E+00  4.0867E-02  1.2211E+00  8.9688E-01
             1.1568E+00
 PARAMETER:  8.3041E-02  4.6070E-01 -3.6503E-01 -1.9159E-01  1.1182E-01  1.8680E-01  1.9795E-01 -3.0974E+00  2.9975E-01 -8.8349E-03
             2.4563E-01
 GRADIENT:  -2.0039E-01 -5.6722E-04  4.9305E-02  2.2388E-01 -1.1757E-01 -5.8451E-02  4.6809E-02  2.5396E-03 -9.6233E-03  1.1505E-02
            -1.0598E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1594.82724256065        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1252
 NPARAMETR:  9.8397E-01  1.4326E+00  6.2812E-01  7.4704E-01  1.0118E+00  1.0920E+00  1.1030E+00  1.0000E-02  1.2215E+00  8.9671E-01
             1.1568E+00
 PARAMETER:  8.3843E-02  4.5952E-01 -3.6502E-01 -1.9164E-01  1.1172E-01  1.8805E-01  1.9806E-01 -4.9710E+00  3.0010E-01 -9.0238E-03
             2.4565E-01
 GRADIENT:   1.3051E+00 -1.4424E+00 -1.7464E-01 -4.4499E-01  4.2124E-01  4.3591E-01 -3.1569E-03  0.0000E+00  9.9685E-02  2.7746E-02
             3.0939E-02

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1594.82724256065        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1274
 NPARAMETR:  9.8397E-01  1.4326E+00  6.2812E-01  7.4704E-01  1.0118E+00  1.0920E+00  1.1030E+00  1.0000E-02  1.2215E+00  8.9671E-01
             1.1568E+00
 PARAMETER:  8.3843E-02  4.5952E-01 -3.6502E-01 -1.9164E-01  1.1172E-01  1.8805E-01  1.9806E-01 -4.9710E+00  3.0010E-01 -9.0238E-03
             2.4565E-01
 GRADIENT:   1.3051E+00 -1.4424E+00 -1.7464E-01 -4.4499E-01  4.2124E-01  4.3591E-01 -3.1569E-03  0.0000E+00  9.9685E-02  2.7746E-02
             3.0939E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1274
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0963E-04 -1.5227E-02 -3.4207E-04  1.1624E-02 -2.7657E-02
 SE:             2.9774E-02  2.4252E-02  1.2802E-04  2.3069E-02  2.0923E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9706E-01  5.3009E-01  7.5411E-03  6.1434E-01  1.8622E-01

 ETASHRINKSD(%)  2.5236E-01  1.8754E+01  9.9571E+01  2.2716E+01  2.9906E+01
 ETASHRINKVR(%)  5.0408E-01  3.3991E+01  9.9998E+01  4.0272E+01  5.0868E+01
 EBVSHRINKSD(%)  5.2258E-01  1.8256E+01  9.9622E+01  2.3743E+01  2.9204E+01
 EBVSHRINKVR(%)  1.0424E+00  3.3179E+01  9.9999E+01  4.1849E+01  4.9879E+01
 RELATIVEINF(%)  9.8859E+01  5.0797E+00  1.5117E-04  4.3352E+00  7.8104E+00
 EPSSHRINKSD(%)  4.3172E+01
 EPSSHRINKVR(%)  6.7706E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1594.8272425606458     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -859.67641599690762     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.41
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.23
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1594.827       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.84E-01  1.43E+00  6.28E-01  7.47E-01  1.01E+00  1.09E+00  1.10E+00  1.00E-02  1.22E+00  8.97E-01  1.16E+00
 


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
+        9.53E+02
 
 TH 2
+       -4.23E+00  3.09E+02
 
 TH 3
+        1.01E+01  1.54E+02  4.19E+02
 
 TH 4
+       -1.18E+01  2.28E+02 -3.07E+02  7.96E+02
 
 TH 5
+       -3.39E+00 -2.07E+02 -4.34E+02  3.11E+02  6.79E+02
 
 TH 6
+        1.35E+00 -1.30E+00  1.37E+00 -3.89E+00 -2.63E+00  1.61E+02
 
 TH 7
+        6.81E-01  1.90E+01 -1.14E+01 -1.59E+01 -5.39E-02 -4.43E-01  7.61E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.23E+00 -1.60E+01 -3.54E+01  4.51E+01 -1.15E+00 -2.81E-02  1.29E+01  0.00E+00  5.50E+01
 
 TH10
+       -9.33E-01 -1.15E+01 -3.43E+01 -1.33E+01 -6.48E+01  7.37E-01  1.12E+01  0.00E+00  1.01E+01  7.12E+01
 
 TH11
+       -7.61E+00 -1.35E+01 -3.22E+01  1.80E+00 -2.57E+00  3.64E+00  8.27E+00  0.00E+00  5.73E+00  1.83E+01  1.61E+02
 
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
 #CPUT: Total CPU Time in Seconds,       23.645
Stop Time:
Wed Sep 29 14:10:31 CDT 2021
