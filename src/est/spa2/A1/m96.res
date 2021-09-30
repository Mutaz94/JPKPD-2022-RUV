Thu Sep 30 05:23:18 CDT 2021
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
$DATA ../../../../data/spa2/A1/dat96.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      700
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

 TOT. NO. OF OBS RECS:      600
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
 RAW OUTPUT FILE (FILE): m96.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2156.16285983211        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.9231E+02  3.2592E+01  4.2979E+01  8.9486E+01  7.2138E+01  2.9088E+01  1.0567E+01 -1.2693E+02  2.3730E+01  7.3919E+00
            -4.5579E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2234.85403349386        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0575E+00  1.0124E+00  8.1489E-01  9.8833E-01  8.9893E-01  1.0803E+00  9.1479E-01  1.8031E+00  8.7101E-01  8.9732E-01
             1.5106E+00
 PARAMETER:  1.5586E-01  1.1230E-01 -1.0471E-01  8.8263E-02 -6.5484E-03  1.7719E-01  1.0943E-02  6.8949E-01 -3.8098E-02 -8.3413E-03
             5.1250E-01
 GRADIENT:   4.4755E+02  2.0490E+01 -2.9942E+01  2.6290E+01 -1.8687E+01  3.7392E+01  1.1803E+00  2.0408E+01 -6.0006E+00  1.4308E+01
             1.1305E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2246.95954486125        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0204E+00  1.0293E+00  1.1622E+00  9.9495E-01  1.0323E+00  1.0537E+00  1.1454E+00  2.2088E+00  7.5969E-01  9.2104E-01
             1.4771E+00
 PARAMETER:  1.2020E-01  1.2891E-01  2.5034E-01  9.4939E-02  1.3182E-01  1.5233E-01  2.3579E-01  8.9246E-01 -1.7485E-01  1.7746E-02
             4.9006E-01
 GRADIENT:   3.3512E+02  3.0481E+01  1.8513E+00  2.4661E+01 -2.0781E+01  3.6734E+01  1.4466E+01  3.1524E+01 -5.9832E+00  1.2482E+01
             9.2447E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2254.28477351727        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:      259
 NPARAMETR:  9.7841E-01  1.0156E+00  1.0873E+00  9.9197E-01  9.9904E-01  1.0010E+00  1.0201E+00  2.0211E+00  8.1052E-01  8.7861E-01
             1.4237E+00
 PARAMETER:  7.8175E-02  1.1552E-01  1.8370E-01  9.1939E-02  9.9041E-02  1.0099E-01  1.1986E-01  8.0364E-01 -1.1008E-01 -2.9418E-02
             4.5328E-01
 GRADIENT:   2.7352E+01 -7.8965E+00 -2.1864E+00 -2.3491E+01 -2.9258E+01 -1.6748E+00 -2.2590E-01  9.7699E+00 -1.2483E+01  7.0372E+00
             5.4418E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2255.11596354999        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:      426
 NPARAMETR:  9.7706E-01  1.0176E+00  1.0923E+00  9.9159E-01  1.0026E+00  9.9569E-01  1.0137E+00  2.0172E+00  8.5119E-01  8.7708E-01
             1.4176E+00
 PARAMETER:  7.6792E-02  1.1747E-01  1.8828E-01  9.1549E-02  1.0263E-01  9.5683E-02  1.1365E-01  8.0172E-01 -6.1124E-02 -3.1154E-02
             4.4898E-01
 GRADIENT:   2.2969E+02  2.0737E+01  2.3762E+00  1.8326E+01 -1.3777E+01  1.5494E+01  4.9892E+00  2.1366E+01 -2.3477E+00  7.5099E+00
             5.6695E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2255.24090253388        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      586
 NPARAMETR:  9.7709E-01  1.0176E+00  1.0922E+00  9.9159E-01  1.0027E+00  9.9983E-01  1.0039E+00  2.0172E+00  8.7588E-01  8.7710E-01
             1.4176E+00
 PARAMETER:  7.6824E-02  1.1747E-01  1.8821E-01  9.1550E-02  1.0266E-01  9.9834E-02  1.0390E-01  8.0172E-01 -3.2521E-02 -3.1129E-02
             4.4896E-01
 GRADIENT:   2.2969E+02  1.9846E+01  2.4637E+00  1.9920E+01 -1.2848E+01  1.7124E+01  5.5728E+00  2.1565E+01  2.0021E+00  7.5270E+00
             5.7531E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2255.45880256195        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      705             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7709E-01  1.0170E+00  1.1072E+00  9.9159E-01  1.0414E+00  1.0150E+00  9.9852E-01  2.0174E+00  8.8638E-01  8.7711E-01
             1.4176E+00
 PARAMETER:  7.6824E-02  1.1689E-01  2.0188E-01  9.1550E-02  1.4059E-01  1.1492E-01  9.8517E-02  8.0179E-01 -2.0608E-02 -3.1128E-02
             4.4894E-01
 GRADIENT:   2.2916E+02 -4.4566E-01 -1.0148E+01  2.6887E+01  3.4487E+01  2.5820E+01  4.7471E+00  1.8847E+01  4.5491E+00  4.3510E+00
             5.6438E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2256.03633778793        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      822
 NPARAMETR:  9.7734E-01  1.0368E+00  1.1135E+00  9.9309E-01  1.0373E+00  9.9681E-01  9.6823E-01  2.0215E+00  8.7746E-01  8.4439E-01
             1.4159E+00
 PARAMETER:  7.7080E-02  1.3612E-01  2.0747E-01  9.3062E-02  1.3657E-01  9.6809E-02  6.7719E-02  8.0385E-01 -3.0725E-02 -6.9143E-02
             4.4776E-01
 GRADIENT:   2.3046E+02  2.9178E+01 -4.5823E+00  4.4921E+01  2.0092E+01  1.5833E+01  1.7345E+00  1.8236E+01  1.3600E-01 -4.2607E-02
             5.2553E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2256.12759238936        NO. OF FUNC. EVALS.: 171
 CUMULATIVE NO. OF FUNC. EVALS.:      993
 NPARAMETR:  9.7734E-01  1.0414E+00  1.1134E+00  9.8858E-01  1.0373E+00  1.0026E+00  9.6565E-01  2.0222E+00  8.8833E-01  8.4979E-01
             1.4159E+00
 PARAMETER:  7.7079E-02  1.4055E-01  2.0747E-01  8.8512E-02  1.3657E-01  1.0198E-01  6.5779E-02  8.0384E-01 -1.9408E-02 -6.3321E-02
             4.4774E-01
 GRADIENT:   2.5227E+01 -2.0491E+01  1.1702E+04 -1.8145E+01  4.0724E+00 -1.3900E+00  3.9910E-01 -3.0765E+03 -1.0700E+00 -3.5862E-01
             8.2550E+01
 NUMSIGDIG:         5.6         5.5         2.3         5.7         6.2         1.2         1.0         2.3         0.9         1.2
                    4.4

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      993
 NO. OF SIG. DIGITS IN FINAL EST.:  0.9

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -9.8860E-03 -1.5238E-02 -2.2957E-02  7.7534E-03 -3.8836E-02
 SE:             2.9883E-02  1.9934E-02  2.0234E-02  2.4703E-02  1.9553E-02
 N:                     100         100         100         100         100

 P VAL.:         7.4078E-01  4.4462E-01  2.5655E-01  7.5362E-01  4.7013E-02

 ETASHRINKSD(%)  1.0000E-10  3.3218E+01  3.2213E+01  1.7243E+01  3.4494E+01
 ETASHRINKVR(%)  1.0000E-10  5.5401E+01  5.4049E+01  3.1513E+01  5.7090E+01
 EBVSHRINKSD(%)  6.1904E-01  3.3239E+01  3.1390E+01  1.8764E+01  3.3773E+01
 EBVSHRINKVR(%)  1.2342E+00  5.5430E+01  5.2927E+01  3.4007E+01  5.6140E+01
 RELATIVEINF(%)  9.8732E+01  8.8457E+00  2.0800E+01  1.6509E+01  1.2001E+01
 EPSSHRINKSD(%)  3.2175E+01
 EPSSHRINKVR(%)  5.3998E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2256.1275923893609     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1153.4013525437538     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.05
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.29
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2256.128       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.77E-01  1.04E+00  1.11E+00  9.89E-01  1.04E+00  1.00E+00  9.66E-01  2.02E+00  8.87E-01  8.49E-01  1.42E+00
 


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
+        6.35E+06
 
 TH 2
+       -4.24E+06  2.83E+06
 
 TH 3
+       -6.29E+03 -1.47E+03  2.27E+06
 
 TH 4
+        6.28E+06  4.07E+03 -6.30E+03  6.21E+06
 
 TH 5
+       -8.95E+02  2.93E+06 -4.54E+03  1.18E+02  3.02E+06
 
 TH 6
+       -6.08E+06 -2.21E+00  1.22E+02 -2.26E+00 -1.17E+00  1.94E+02
 
 TH 7
+       -4.62E+01  4.29E+06 -2.72E+06 -5.18E+01 -4.43E+06 -3.92E-01  4.73E+01
 
 TH 8
+        1.59E+04  2.16E+02 -3.93E+02  1.57E+04  2.69E+05 -1.72E+01  2.07E-01  2.40E+04
 
 TH 9
+        7.00E+06 -1.73E+01 -2.59E+02  2.58E+01  8.68E+00 -4.66E-02  3.54E+01  4.02E+01  1.10E+02
 
 TH10
+       -3.23E+02 -4.88E+06  3.10E+06 -3.12E+02  5.04E+06  4.69E-01 -7.39E+06  6.15E+00  2.44E+00  8.41E+06
 
 TH11
+       -9.82E+05 -4.21E+03  8.28E+05 -6.22E+03  1.20E+02  2.49E+00 -9.87E+05 -1.18E+05  1.22E+01  1.12E+06  3.02E+05
 
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
 #CPUT: Total CPU Time in Seconds,       27.429
Stop Time:
Thu Sep 30 05:23:47 CDT 2021
