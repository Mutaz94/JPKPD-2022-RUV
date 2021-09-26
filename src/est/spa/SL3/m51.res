Sat Sep 25 11:45:29 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat51.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
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

$OMEGA  0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       25 SEP 2021
Days until program expires : 204
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 NO. OF SIG. FIGURES REQUIRED:            3
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1549.44109722550        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.3174E+02 -1.7636E+01 -2.4598E+01  3.6164E+01  1.6940E+02  1.1653E+01  5.0283E+00 -2.4989E+01 -2.5613E+00 -6.6286E+01
            -8.8627E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1573.83585206607        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.6548E-01  9.4393E-01  8.9260E-01  1.0337E+00  8.5787E-01  9.6565E-01  9.3414E-01  1.1412E+00  9.8999E-01  1.1262E+00
             1.1867E+00
 PARAMETER:  6.4875E-02  4.2296E-02 -1.3621E-02  1.3316E-01 -5.3308E-02  6.5045E-02  3.1873E-02  2.3209E-01  8.9939E-02  2.1882E-01
             2.7114E-01
 GRADIENT:   3.7600E+01  2.2210E+00 -2.1193E+01  3.2015E+01  3.6018E+01  1.1789E+00  6.1898E+00 -1.7078E+00  2.1811E+00  4.1942E+00
             6.1824E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1576.95678055906        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.5308E-01  8.4368E-01  1.2260E+00  1.1023E+00  9.4714E-01  9.5515E-01  4.8121E-01  1.4446E+00  1.0469E+00  1.2283E+00
             1.2044E+00
 PARAMETER:  5.1942E-02 -6.9983E-02  3.0374E-01  1.9744E-01  4.5693E-02  5.4115E-02 -6.3146E-01  4.6783E-01  1.4588E-01  3.0562E-01
             2.8602E-01
 GRADIENT:   1.2330E+01  1.2924E+00 -6.4818E+00  2.9087E+01  2.3521E+01 -2.4409E+00  2.8789E+00 -5.3602E+00  1.0339E+01  2.8222E+00
             1.2530E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1578.84063923652        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.4834E-01  9.9355E-01  1.1484E+00  9.8803E-01  9.5266E-01  9.5789E-01  2.7869E-01  1.5766E+00  1.1388E+00  1.2081E+00
             1.1782E+00
 PARAMETER:  4.6958E-02  9.3531E-02  2.3839E-01  8.7956E-02  5.1498E-02  5.6974E-02 -1.1777E+00  5.5527E-01  2.2996E-01  2.8907E-01
             2.6398E-01
 GRADIENT:   2.3900E-01  2.9295E+00 -7.4863E-01  5.5513E+00  4.4384E+00 -1.7579E+00  1.4646E+00 -2.4870E+00  1.6579E+00  1.5718E+00
             3.0064E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1579.64367092816        NO. OF FUNC. EVALS.: 130
 CUMULATIVE NO. OF FUNC. EVALS.:      360
 NPARAMETR:  9.6228E-01  9.9081E-01  1.1494E+00  9.8768E-01  9.5293E-01  9.7813E-01  8.3775E-02  1.5789E+00  1.1445E+00  1.2082E+00
             1.1782E+00
 PARAMETER:  6.1550E-02  9.0763E-02  2.3922E-01  8.7600E-02  5.1790E-02  7.7889E-02 -2.3796E+00  5.5671E-01  2.3493E-01  2.8915E-01
             2.6399E-01
 GRADIENT:   3.4967E+01  3.5100E+00 -9.0620E-01  5.1691E+00  5.3275E+00  6.5632E+00  1.6486E-01 -2.8500E+00 -1.0617E+00  1.1282E+00
             2.8294E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1579.67265309012        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      499
 NPARAMETR:  9.5965E-01  9.9157E-01  1.1494E+00  9.8768E-01  9.5293E-01  9.7225E-01  7.4700E-02  1.5789E+00  1.1512E+00  1.2082E+00
             1.1782E+00
 PARAMETER:  5.8815E-02  9.1531E-02  2.3922E-01  8.7600E-02  5.1790E-02  7.1853E-02 -2.4943E+00  5.5671E-01  2.4078E-01  2.8915E-01
             2.6399E-01
 GRADIENT:  -2.6905E+00 -7.9005E-02 -1.0352E+00  1.7703E+00  4.4245E+00 -2.1853E-01  9.4299E-02 -3.1239E+00 -8.8542E-01  8.1538E-01
             2.7611E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1579.72083940833        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:      664
 NPARAMETR:  9.5967E-01  9.9154E-01  1.1491E+00  9.8758E-01  9.5284E-01  9.7275E-01  1.0000E-02  1.5797E+00  1.1514E+00  1.2079E+00
             1.1785E+00
 PARAMETER:  5.8830E-02  9.1507E-02  2.3898E-01  8.7500E-02  5.1690E-02  7.2373E-02 -4.6731E+00  5.5727E-01  2.4097E-01  2.8886E-01
             2.6425E-01
 GRADIENT:  -2.6731E+00  1.9321E-01 -1.0850E+00  1.9209E+00  4.4408E+00 -2.7733E-02  0.0000E+00 -3.1184E+00 -1.2523E+00  7.6443E-01
             2.8669E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1579.76217189727        NO. OF FUNC. EVALS.: 151
 CUMULATIVE NO. OF FUNC. EVALS.:      815
 NPARAMETR:  9.6063E-01  9.9155E-01  1.1499E+00  9.8637E-01  9.4863E-01  9.7271E-01  1.0000E-02  1.5877E+00  1.1555E+00  1.2072E+00
             1.1705E+00
 PARAMETER:  5.9833E-02  9.1512E-02  2.3967E-01  8.6272E-02  4.7261E-02  7.2329E-02 -4.6731E+00  5.6230E-01  2.4451E-01  2.8834E-01
             2.5746E-01
 GRADIENT:   3.1653E+01  3.4182E+00  7.2382E-01  5.0525E+00  6.6478E-01  4.4776E+00  0.0000E+00 -2.6075E+00  1.3824E+00  1.4708E+00
             3.3135E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1579.76267672503        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      993
 NPARAMETR:  9.6069E-01  9.9226E-01  1.1499E+00  9.8564E-01  9.4864E-01  9.7289E-01  1.0000E-02  1.5877E+00  1.1569E+00  1.2072E+00
             1.1700E+00
 PARAMETER:  5.9897E-02  9.2233E-02  2.3967E-01  8.5537E-02  4.7270E-02  7.2513E-02 -4.6731E+00  5.6230E-01  2.4572E-01  2.8834E-01
             2.5697E-01
 GRADIENT:   3.6733E-02 -2.8863E-01  6.1053E-01  1.5215E-01 -1.6685E-01  2.0796E-02  0.0000E+00 -2.8661E+00  1.5349E-02  1.1564E+00
            -6.7871E-02

0ITERATION NO.:   41    OBJECTIVE VALUE:  -1579.76267672503        NO. OF FUNC. EVALS.:  30
 CUMULATIVE NO. OF FUNC. EVALS.:     1023
 NPARAMETR:  9.6079E-01  9.9236E-01  1.1502E+00  9.8574E-01  9.4854E-01  9.7279E-01  1.0000E-02  1.5886E+00  1.1566E+00  1.2069E+00
             1.1703E+00
 PARAMETER:  5.9897E-02  9.2233E-02  2.3967E-01  8.5537E-02  4.7270E-02  7.2513E-02 -4.6731E+00  5.6230E-01  2.4572E-01  2.8834E-01
             2.5697E-01
 GRADIENT:  -1.6763E+06 -8.3814E+05 -6.9939E+05 -1.6763E+06  1.6763E+06  8.3814E+05  0.0000E+00 -1.4908E+05  6.8218E+05  2.9068E+05
            -3.2630E+05

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1023
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2975E-03 -8.0292E-04 -3.7986E-02 -2.7121E-03 -3.3094E-02
 SE:             2.9818E-02  2.1781E-04  1.6195E-02  2.8871E-02  2.2785E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6529E-01  2.2763E-04  1.8998E-02  9.2516E-01  1.4637E-01

 ETASHRINKSD(%)  1.0752E-01  9.9270E+01  4.5745E+01  3.2797E+00  2.3669E+01
 ETASHRINKVR(%)  2.1492E-01  9.9995E+01  7.0564E+01  6.4517E+00  4.1735E+01
 EBVSHRINKSD(%)  6.0206E-01  9.9333E+01  5.0181E+01  3.5936E+00  2.0502E+01
 EBVSHRINKVR(%)  1.2005E+00  9.9996E+01  7.5180E+01  7.0581E+00  3.6801E+01
 RELATIVEINF(%)  9.8560E+01  4.3726E-04  6.6951E+00  1.0845E+01  1.8711E+01
 EPSSHRINKSD(%)  4.4658E+01
 EPSSHRINKVR(%)  6.9373E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1579.7626767250254     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -844.61185016128718     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.26
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.01
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1579.763       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.61E-01  9.92E-01  1.15E+00  9.86E-01  9.49E-01  9.73E-01  1.00E-02  1.59E+00  1.16E+00  1.21E+00  1.17E+00
 


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
+        4.54E+09
 
 TH 2
+       -1.13E+04  4.26E+09
 
 TH 3
+       -4.05E+03 -6.98E+04  5.52E+08
 
 TH 4
+       -1.14E+04  6.63E+04 -7.03E+04  4.31E+09
 
 TH 5
+        1.18E+04 -6.86E+04  7.29E+04 -4.48E+09  9.31E+09
 
 TH 6
+       -1.46E+04 -1.41E+04 -5.08E+03 -1.42E+04  1.48E+04  4.43E+09
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -1.25E+03  1.71E+05  6.16E+04  1.72E+05 -1.79E+05 -1.57E+03  0.00E+00  5.26E+07
 
 TH 9
+       -1.10E+04 -1.08E+04 -3.84E+03 -1.07E+04  1.12E+04  4.93E+03  0.00E+00 -1.18E+03  5.19E+08
 
 TH10
+        3.21E+03  5.66E+03  2.03E+03  5.70E+03 -5.98E+03  4.02E+03  0.00E+00  6.37E+02  3.04E+03  3.46E+08
 
 TH11
+        1.23E+06  1.19E+06  4.29E+05  1.20E+06 -1.25E+06 -4.66E+03  0.00E+00  1.32E+05 -3.51E+03 -3.40E+05  4.64E+08
 
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
 #CPUT: Total CPU Time in Seconds,       18.344
Stop Time:
Sat Sep 25 11:45:49 CDT 2021
