Sat Sep 18 12:03:56 CDT 2021
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
$DATA ../../../../data/spa/SL2/dat7.csv ignore=@
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
Current Date:       18 SEP 2021
Days until program expires : 211
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
 RAW OUTPUT FILE (FILE): m7.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1155.04340610915        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.6508E+01 -9.7017E+00  5.1172E+01  3.3410E+00  7.5021E+01 -3.1344E+01  1.4853E+01 -2.0783E+02 -5.3408E+00 -1.0885E+01
            -7.6693E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1649.01525971233        NO. OF FUNC. EVALS.:  84
 CUMULATIVE NO. OF FUNC. EVALS.:       97
 NPARAMETR:  9.6558E-01  1.0080E+00  9.8710E-01  9.8638E-01  9.7427E-01  1.0644E+00  9.5631E-01  1.0784E+00  9.5258E-01  9.7803E-01
             1.1306E+00
 PARAMETER:  6.4975E-02  1.0793E-01  8.7012E-02  8.6291E-02  7.3928E-02  1.6240E-01  5.5325E-02  1.7549E-01  5.1421E-02  7.7787E-02
             2.2279E-01
 GRADIENT:  -1.8201E+01 -2.3506E+00 -1.1068E+01  1.5892E+01  2.6450E+01 -6.0500E+00  9.1236E+00  5.7284E-01  9.0744E+00  1.9159E+00
             2.4786E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1650.13325294425        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      281
 NPARAMETR:  9.6415E-01  1.0535E+00  9.1309E-01  9.4968E-01  9.5098E-01  1.0970E+00  7.3839E-01  1.0708E+00  1.0364E+00  9.5116E-01
             1.1346E+00
 PARAMETER:  6.3492E-02  1.5215E-01  9.0796E-03  4.8373E-02  4.9737E-02  1.9254E-01 -2.0328E-01  1.6838E-01  1.3571E-01  4.9926E-02
             2.2632E-01
 GRADIENT:  -5.5620E+01 -1.5017E+01 -1.1873E+01  5.6965E+00  1.6143E+01 -4.6649E+00  2.7535E+00  1.1435E+00  1.1437E+01 -3.8203E-02
             2.7003E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1650.44799005435        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:      473             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6502E-01  1.0770E+00  9.0555E-01  9.3597E-01  9.5743E-01  1.1034E+00  6.8377E-01  1.0678E+00  1.0531E+00  9.5837E-01
             1.1328E+00
 PARAMETER:  6.4398E-02  1.7416E-01  7.9105E-04  3.3827E-02  5.6493E-02  1.9837E-01 -2.8014E-01  1.6564E-01  1.5175E-01  5.7476E-02
             2.2473E-01
 GRADIENT:  -1.7248E+01 -5.2946E+00 -1.0772E+01  1.3635E+01  1.5582E+01  1.0267E+01  1.9770E+00  3.1390E-01  1.0374E+01 -7.1062E-01
             2.6139E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1650.45887075843        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      612
 NPARAMETR:  9.6502E-01  1.0770E+00  9.0555E-01  9.3597E-01  9.5743E-01  1.0948E+00  6.6055E-01  1.0678E+00  1.0531E+00  9.5837E-01
             1.1328E+00
 PARAMETER:  6.4398E-02  1.7416E-01  7.9104E-04  3.3827E-02  5.6493E-02  1.9061E-01 -3.1468E-01  1.6564E-01  1.5175E-01  5.7476E-02
             2.2473E-01
 GRADIENT:  -5.4511E+01 -1.1540E+01 -1.1163E+01  8.4650E+00  1.5273E+01 -5.5114E+00  8.8784E-01  2.0340E-02  8.1652E+00 -1.4138E+00
             2.5733E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1650.49707715025        NO. OF FUNC. EVALS.: 121
 CUMULATIVE NO. OF FUNC. EVALS.:      733
 NPARAMETR:  9.6507E-01  1.0769E+00  9.0551E-01  9.3592E-01  9.5747E-01  1.1101E+00  6.6066E-01  1.0665E+00  1.0532E+00  9.5842E-01
             1.1328E+00
 PARAMETER:  6.4448E-02  1.7408E-01  7.4102E-04  3.3777E-02  5.6543E-02  2.0449E-01 -3.1452E-01  1.6436E-01  1.5183E-01  5.7526E-02
             2.2473E-01
 GRADIENT:  -5.2947E+01 -1.1808E+01 -1.1104E+01  8.2027E+00  1.5418E+01  1.0694E-02  8.9633E-01 -2.6261E-02  8.1967E+00 -1.4334E+00
             2.5758E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1650.52690546973        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:      931
 NPARAMETR:  9.6531E-01  1.0770E+00  9.0572E-01  9.3588E-01  9.5739E-01  1.1101E+00  6.6038E-01  1.0673E+00  1.0529E+00  9.5849E-01
             1.1325E+00
 PARAMETER:  6.4694E-02  1.7416E-01  9.7324E-04  3.3733E-02  5.6452E-02  2.0448E-01 -3.1493E-01  1.6515E-01  1.5153E-01  5.7602E-02
             2.2440E-01
 GRADIENT:  -5.2504E+01 -1.1606E+01 -1.0959E+01  8.1807E+00  1.5098E+01  3.1086E-02  8.8419E-01 -2.2319E-02  8.0982E+00 -1.4058E+00
             2.5628E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1650.53217091490        NO. OF FUNC. EVALS.: 148
 CUMULATIVE NO. OF FUNC. EVALS.:     1079
 NPARAMETR:  9.6531E-01  1.0770E+00  9.0571E-01  9.3587E-01  9.5739E-01  1.1094E+00  6.5593E-01  1.0673E+00  1.0529E+00  9.5849E-01
             1.1325E+00
 PARAMETER:  6.4692E-02  1.7416E-01  9.7168E-04  3.3731E-02  5.6454E-02  2.0365E-01 -3.2172E-01  1.6514E-01  1.5153E-01  5.7600E-02
             2.2440E-01
 GRADIENT:  -4.9015E+01  1.1188E+06  1.9485E+06  9.7423E+05 -1.9484E+06 -2.8956E-01 -6.0564E+05  1.1798E+06 -9.4540E+00 -6.2611E+00
            -5.4931E+02
 NUMSIGDIG:         8.2         3.3         3.3         3.3         3.3         2.1         3.3         3.3         8.7         9.1
                    6.8

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1079
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.6742E-02 -2.4092E-02 -2.0919E-02  4.2322E-03 -3.9708E-02
 SE:             2.9769E-02  1.3991E-02  1.3850E-02  2.5714E-02  2.2393E-02
 N:                     100         100         100         100         100

 P VAL.:         3.6902E-01  8.5069E-02  1.3095E-01  8.6927E-01  7.6195E-02

 ETASHRINKSD(%)  2.7005E-01  5.3130E+01  5.3600E+01  1.3854E+01  2.4980E+01
 ETASHRINKVR(%)  5.3937E-01  7.8032E+01  7.8470E+01  2.5789E+01  4.3720E+01
 EBVSHRINKSD(%)  4.3825E-01  5.3240E+01  5.6259E+01  1.1733E+01  2.3797E+01
 EBVSHRINKVR(%)  8.7458E-01  7.8135E+01  8.0868E+01  2.2089E+01  4.1932E+01
 RELATIVEINF(%)  9.8724E+01  7.3461E-01  2.2539E+00  3.4388E+00  8.8680E+00
 EPSSHRINKSD(%)  4.7170E+01
 EPSSHRINKVR(%)  7.2090E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1650.5321709149018     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -915.38134435116365     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.00
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.91
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1650.532       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.65E-01  1.08E+00  9.06E-01  9.36E-01  9.57E-01  1.11E+00  6.56E-01  1.07E+00  1.05E+00  9.58E-01  1.13E+00
 


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
+        1.05E+10
 
 TH 2
+        4.96E+03  1.38E+09
 
 TH 3
+        1.05E+01  8.69E+04  5.94E+09
 
 TH 4
+       -5.39E+09 -4.12E+04  1.74E+05  5.56E+09
 
 TH 5
+       -1.32E+00  2.37E+04 -1.70E+05  4.82E+04  5.31E+09
 
 TH 6
+        4.40E+01  5.83E+03  1.21E+04  1.17E+04 -1.14E+04  1.60E+02
 
 TH 7
+        4.23E-01  4.21E+03  8.76E+03  8.46E+03 -8.29E+03 -5.19E+03  1.09E+09
 
 TH 8
+        2.86E+09 -1.47E+09 -3.05E+09 -2.95E+09  2.89E+09  6.21E+03  4.50E+03  1.57E+09
 
 TH 9
+        3.14E+00 -2.19E+04 -4.53E+04 -3.26E+09  3.19E+09  3.03E+00  1.45E+09  1.73E+09  3.83E+09
 
 TH10
+       -3.04E-01 -6.56E+03 -1.36E+04 -5.43E+09  5.31E+09 -9.36E-01  2.41E+09 -6.98E+03  6.37E+09  1.06E+10
 
 TH11
+       -5.70E+00  1.02E+09 -1.70E+01  2.05E+09 -6.83E+00  9.76E-01  2.40E+00  3.35E+01  1.39E+00  1.12E+01  1.51E+09
 
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
 #CPUT: Total CPU Time in Seconds,       19.960
Stop Time:
Sat Sep 18 12:04:18 CDT 2021
