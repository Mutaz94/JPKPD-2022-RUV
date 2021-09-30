Wed Sep 29 13:13:47 CDT 2021
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
$DATA ../../../../data/spa/A2/dat100.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m100.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -954.280993232289        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0750E+02  1.2240E+01  3.7301E+01  1.2185E+01  7.8921E+01  2.9239E+01 -1.8416E+01 -2.5564E+01 -2.5032E+01 -5.4204E+01
            -1.2107E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1322.31626736027        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.2150E+00  1.0370E+00  9.6000E-01  1.0655E+00  9.3295E-01  1.2928E+00  1.0254E+00  1.0237E+00  1.0230E+00  1.0542E+00
             2.6580E+00
 PARAMETER:  2.9476E-01  1.3635E-01  5.9174E-02  1.6345E-01  3.0592E-02  3.5678E-01  1.2505E-01  1.2340E-01  1.2272E-01  1.5277E-01
             1.0776E+00
 GRADIENT:   4.6192E+02  2.0318E+01  7.9472E+00  2.8763E+01 -2.2496E+01  3.1386E+01  4.6937E+00  2.5630E+00  1.1392E+01  7.4145E+00
             3.0340E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1360.12603159006        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.0156E+00  1.1725E+00  7.9965E-01  9.4153E-01  9.4020E-01  9.5275E-01  1.1438E+00  1.2696E+00  5.2772E-01  5.8003E-01
             2.4493E+00
 PARAMETER:  1.1549E-01  2.5917E-01 -1.2359E-01  3.9754E-02  3.8341E-02  5.1597E-02  2.3432E-01  3.3871E-01 -5.3919E-01 -4.4468E-01
             9.9581E-01
 GRADIENT:   1.8878E+02  1.9111E+01 -1.2974E-01  1.9428E+01 -5.6207E+00 -1.5556E+01 -2.5434E+00 -1.2863E-02 -5.8075E+00  3.1298E+00
            -9.8324E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1372.87008935519        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  9.4301E-01  1.2093E+00  1.1994E+00  9.1365E-01  1.1478E+00  9.9263E-01  1.1488E+00  1.6404E+00  2.6323E-01  3.8676E-01
             2.8650E+00
 PARAMETER:  4.1323E-02  2.9004E-01  2.8180E-01  9.6871E-03  2.3789E-01  9.2606E-02  2.3871E-01  5.9494E-01 -1.2347E+00 -8.4996E-01
             1.1526E+00
 GRADIENT:  -6.0387E+00  2.0014E+00  2.7269E-01  2.8110E+00 -7.8085E+00  8.1257E+00 -4.3287E+00 -9.6910E-01 -4.7984E-01  1.3143E+00
            -1.5963E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1376.50011907689        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  9.4622E-01  1.1363E+00  2.8665E+00  9.9036E-01  1.4396E+00  9.7055E-01  1.2252E+00  3.3970E+00  2.8642E-01  3.2793E-01
             2.8903E+00
 PARAMETER:  4.4724E-02  2.2782E-01  1.1531E+00  9.0318E-02  4.6435E-01  7.0103E-02  3.0312E-01  1.3229E+00 -1.1503E+00 -1.0150E+00
             1.1614E+00
 GRADIENT:  -1.2183E+00 -4.1673E+00 -4.0737E+00  1.0518E+00  8.2519E+00 -3.6153E-01  2.9496E-01  2.1690E+00  4.7741E-02  3.1483E-01
            -2.9337E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1378.23376572126        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:      443             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8733E-01  1.1406E+00  4.7694E+00  9.9697E-01  1.4832E+00  1.0010E+00  1.1740E+00  4.7327E+00  3.8713E-01  1.1429E-01
             2.9178E+00
 PARAMETER:  8.7245E-02  2.3158E-01  1.6622E+00  9.6968E-02  4.9424E-01  1.0097E-01  2.6045E-01  1.6545E+00 -8.4898E-01 -2.0690E+00
             1.1708E+00
 GRADIENT:   8.8502E+01 -1.1341E+01  9.0683E-01 -1.6228E+01  6.4913E+00  1.0406E+01  2.1905E+00 -1.1258E+00  9.7090E-01  5.4247E-02
             7.9217E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1378.77046117835        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      579
 NPARAMETR:  9.7124E-01  1.1438E+00  4.7409E+00  9.9755E-01  1.4786E+00  9.8284E-01  1.1737E+00  4.7592E+00  3.8714E-01  1.8242E-02
             2.9187E+00
 PARAMETER:  7.0818E-02  2.3435E-01  1.6562E+00  9.7548E-02  4.9108E-01  8.2687E-02  2.6014E-01  1.6601E+00 -8.4896E-01 -3.9041E+00
             1.1711E+00
 GRADIENT:   6.9155E+00 -1.3603E+01 -2.0385E+00 -1.1538E+01 -5.7981E-01 -2.9689E-01  9.4181E-01  6.1460E-01 -2.3531E-02  7.5292E-04
            -4.4649E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1378.88776820583        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      754
 NPARAMETR:  9.6949E-01  1.1536E+00  4.8645E+00  9.9917E-01  1.4715E+00  9.8456E-01  1.1722E+00  4.8104E+00  3.8739E-01  1.0000E-02
             2.9405E+00
 PARAMETER:  6.9011E-02  2.4291E-01  1.6820E+00  9.9170E-02  4.8630E-01  8.4436E-02  2.5884E-01  1.6708E+00 -8.4833E-01 -7.8418E+00
             1.1786E+00
 GRADIENT:   8.5409E-01 -2.8473E+00 -3.2278E-01 -6.5089E-01 -4.0183E+00  1.2320E-01  1.5223E+00 -1.7635E+00 -2.0200E-01  0.0000E+00
             3.5811E+00

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1378.92721236246        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      890
 NPARAMETR:  9.7043E-01  1.1568E+00  4.9877E+00  9.9949E-01  1.4784E+00  9.8647E-01  1.1673E+00  4.8364E+00  3.8755E-01  1.0000E-02
             2.9249E+00
 PARAMETER:  7.0489E-02  2.4583E-01  1.7056E+00  9.9572E-02  4.9051E-01  8.7039E-02  2.5725E-01  1.6773E+00 -8.4719E-01 -8.6061E+00
             1.1724E+00
 GRADIENT:   3.1851E+00  1.0083E+02 -1.6136E+01  1.2338E+02 -3.0077E+01  6.7455E-01  1.0565E+00  1.7009E+01  2.6744E+01  0.0000E+00
            -2.8203E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      890
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2645E-03  6.3622E-03 -1.8861E-02 -1.6663E-02 -7.1724E-05
 SE:             2.9051E-02  2.2886E-02  9.6722E-03  9.8718E-03  1.2912E-04
 N:                     100         100         100         100         100

 P VAL.:         9.6528E-01  7.8101E-01  5.1177E-02  9.1433E-02  5.7857E-01

 ETASHRINKSD(%)  2.6755E+00  2.3330E+01  6.7597E+01  6.6928E+01  9.9567E+01
 ETASHRINKVR(%)  5.2794E+00  4.1218E+01  8.9500E+01  8.9062E+01  9.9998E+01
 EBVSHRINKSD(%)  2.7569E+00  2.2858E+01  7.7020E+01  6.7787E+01  9.9538E+01
 EBVSHRINKVR(%)  5.4378E+00  4.0490E+01  9.4719E+01  8.9623E+01  9.9998E+01
 RELATIVEINF(%)  9.3175E+01  2.9938E+00  2.5202E+00  5.1981E-01  9.8230E-04
 EPSSHRINKSD(%)  2.3169E+01
 EPSSHRINKVR(%)  4.0971E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1378.9272123624582     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -643.77638579872007     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.36
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.23
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1378.927       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.71E-01  1.16E+00  4.98E+00  1.00E+00  1.48E+00  9.87E-01  1.17E+00  4.84E+00  3.88E-01  1.00E-02  2.92E+00
 


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
+        1.15E+03
 
 TH 2
+        2.08E+02  8.56E+03
 
 TH 3
+       -4.26E+00 -1.50E+01  2.88E+01
 
 TH 4
+        7.88E+02  2.36E+04 -4.10E+01  6.56E+04
 
 TH 5
+       -5.71E+01 -2.97E+02  2.60E+01 -5.85E+02  1.58E+03
 
 TH 6
+       -8.43E-01 -1.55E+02  1.23E+00 -4.86E+02  1.86E+01  1.81E+02
 
 TH 7
+        1.93E+04  2.73E+02 -3.34E+00  8.14E+02 -4.51E+01  2.72E+00  6.53E+03
 
 TH 8
+        4.59E+00  1.09E+01 -3.14E+01  4.23E+01 -2.36E+02  5.67E-02  2.95E+00  4.96E+01
 
 TH 9
+        1.83E+04  3.37E+02 -4.56E+00  1.14E+03 -5.23E+01 -4.30E+02  6.06E+03  3.99E+00  5.63E+03
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.30E+01 -3.63E+01  7.27E+01 -1.02E+02  5.11E+01  3.98E+00  1.42E+00 -8.79E+01 -3.95E+00  0.00E+00  3.11E+02
 
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
 #CPUT: Total CPU Time in Seconds,       18.650
Stop Time:
Wed Sep 29 13:14:08 CDT 2021
