Wed Sep 29 12:37:50 CDT 2021
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
$DATA ../../../../data/spa/A2/dat18.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m18.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1032.20317722845        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2865E+02 -4.2260E+01  1.7767E+01 -5.4624E+01  1.1194E+02  4.5312E+01 -5.7957E+01 -1.4271E+01 -6.6865E+01 -6.7337E+01
            -1.0469E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1359.71947267160        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.1924E+00  1.0830E+00  1.0117E+00  1.0627E+00  9.2548E-01  1.2044E+00  1.2157E+00  9.8386E-01  1.1762E+00  1.0081E+00
             2.7390E+00
 PARAMETER:  2.7601E-01  1.7971E-01  1.1158E-01  1.6083E-01  2.2559E-02  2.8601E-01  2.9534E-01  8.3730E-02  2.6225E-01  1.0803E-01
             1.1076E+00
 GRADIENT:   3.9917E+02  1.2575E+01  1.4269E+01  5.8378E+00 -4.0644E+01  4.4469E+01  4.5687E-01  3.2421E+00  1.7217E+01  1.2759E+01
             6.1731E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1395.13527161155        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0586E+00  9.5259E-01  5.6193E-01  1.0629E+00  6.6809E-01  9.7372E-01  1.3376E+00  2.6958E-01  8.9080E-01  2.0238E-01
             2.6694E+00
 PARAMETER:  1.5691E-01  5.1433E-02 -4.7638E-01  1.6100E-01 -3.0334E-01  7.3369E-02  3.9089E-01 -1.2109E+00 -1.5632E-02 -1.4976E+00
             1.0818E+00
 GRADIENT:   1.8925E+02 -7.5766E+00  1.2622E+01  1.0143E+01 -3.6692E+00  1.2271E+01 -2.4393E+01 -1.8017E-01 -2.7121E+00 -3.5611E-01
             2.5585E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1401.13312074820        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:      293
 NPARAMETR:  1.0112E+00  8.5994E-01  5.1425E-01  1.0939E+00  6.0659E-01  9.2771E-01  1.5830E+00  1.9389E-01  9.0786E-01  3.1278E-01
             2.4526E+00
 PARAMETER:  1.1111E-01 -5.0894E-02 -5.6505E-01  1.8977E-01 -3.9991E-01  2.4968E-02  5.5935E-01 -1.5405E+00  3.3361E-03 -1.0623E+00
             9.9714E-01
 GRADIENT:   9.2313E+00 -7.9548E+00 -4.8292E+00 -3.8772E+00  1.0615E+01 -6.3441E+00 -7.5814E+00  1.2410E-02  5.2857E-02 -9.8405E-01
            -4.2385E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1403.12836345176        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      468
 NPARAMETR:  9.9735E-01  6.5852E-01  6.7616E-01  1.2480E+00  6.3206E-01  9.2962E-01  2.0904E+00  1.4056E-01  8.3855E-01  6.7505E-01
             2.4167E+00
 PARAMETER:  9.7347E-02 -3.1777E-01 -2.9133E-01  3.2155E-01 -3.5877E-01  2.7017E-02  8.3738E-01 -1.8621E+00 -7.6086E-02 -2.9296E-01
             9.8240E-01
 GRADIENT:  -1.8240E+01  5.2144E+00  3.5498E+00  8.3857E-01 -5.3439E+00 -5.1088E+00  1.1886E-01  6.8921E-02 -3.7201E-01  8.7647E-01
             5.8591E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1403.72842788804        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      643
 NPARAMETR:  1.0024E+00  5.0791E-01  6.2558E-01  1.3103E+00  5.6304E-01  9.4204E-01  2.4661E+00  1.3473E-01  8.2198E-01  6.6603E-01
             2.3614E+00
 PARAMETER:  1.0243E-01 -5.7746E-01 -3.6908E-01  3.7028E-01 -4.7441E-01  4.0292E-02  1.0026E+00 -1.9045E+00 -9.6034E-02 -3.0641E-01
             9.5927E-01
 GRADIENT:  -2.0200E+00 -4.7590E-01 -1.5111E+00 -3.7001E-01  2.2052E+00 -6.0899E-01  3.9965E-01  9.5039E-02  4.6156E-02  1.2224E-01
             1.1231E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1403.77163171163        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      819
 NPARAMETR:  1.0041E+00  5.1390E-01  6.0961E-01  1.3042E+00  5.5379E-01  9.4389E-01  2.4247E+00  6.3855E-02  8.2467E-01  6.5703E-01
             2.3538E+00
 PARAMETER:  1.0413E-01 -5.6572E-01 -3.9493E-01  3.6562E-01 -4.9096E-01  4.2254E-02  9.8571E-01 -2.6511E+00 -9.2777E-02 -3.2003E-01
             9.5604E-01
 GRADIENT:   1.1957E+00 -1.4510E-02  3.6186E-01  3.2978E-01 -6.6797E-01 -1.3199E-01 -4.0349E-01  2.1076E-02 -2.6095E-01 -1.0008E-01
            -1.9344E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1403.78281937042        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      994
 NPARAMETR:  1.0035E+00  5.1249E-01  6.1531E-01  1.3062E+00  5.5725E-01  9.4410E-01  2.4389E+00  1.4605E-02  8.2501E-01  6.6469E-01
             2.3539E+00
 PARAMETER:  1.0347E-01 -5.6847E-01 -3.8563E-01  3.6712E-01 -4.8474E-01  4.2475E-02  9.9155E-01 -4.1264E+00 -9.2359E-02 -3.0844E-01
             9.5607E-01
 GRADIENT:  -1.9644E-02 -2.1586E-02 -5.3520E-02  3.4068E-02  1.0200E-01  1.9100E-03 -5.2055E-04  1.1154E-03  4.2972E-03  1.3714E-03
            -3.8331E-02

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1403.78314468665        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1086
 NPARAMETR:  1.0035E+00  5.1253E-01  6.1524E-01  1.3062E+00  5.5720E-01  9.4390E-01  2.4394E+00  1.0000E-02  8.2502E-01  6.6467E-01
             2.3539E+00
 PARAMETER:  1.0352E-01 -5.6841E-01 -3.8575E-01  3.6713E-01 -4.8483E-01  4.2269E-02  9.9175E-01 -4.9083E+00 -9.2346E-02 -3.0846E-01
             9.5605E-01
 GRADIENT:   9.4804E-02  1.9729E-02 -3.4990E-02  7.4478E-02  5.9044E-02 -7.3582E-02  3.5295E-02  0.0000E+00  8.7671E-03  9.4311E-04
            -4.6882E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1086
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.8484E-03  3.5182E-02 -3.3587E-04 -3.4401E-02  7.8392E-03
 SE:             2.9175E-02  2.0853E-02  1.7224E-04  2.2898E-02  1.6587E-02
 N:                     100         100         100         100         100

 P VAL.:         9.4948E-01  9.1579E-02  5.1171E-02  1.3301E-01  6.3649E-01

 ETASHRINKSD(%)  2.2615E+00  3.0140E+01  9.9423E+01  2.3289E+01  4.4431E+01
 ETASHRINKVR(%)  4.4719E+00  5.1195E+01  9.9997E+01  4.1154E+01  6.9121E+01
 EBVSHRINKSD(%)  2.4305E+00  3.1740E+01  9.9387E+01  2.1921E+01  4.2428E+01
 EBVSHRINKVR(%)  4.8020E+00  5.3405E+01  9.9996E+01  3.9037E+01  6.6855E+01
 RELATIVEINF(%)  9.4218E+01  9.8986E+00  1.9934E-04  1.8028E+01  1.5590E+00
 EPSSHRINKSD(%)  3.4112E+01
 EPSSHRINKVR(%)  5.6587E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1403.7831446866467     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -668.63231812290849     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.23
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.42
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1403.783       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  5.13E-01  6.15E-01  1.31E+00  5.57E-01  9.44E-01  2.44E+00  1.00E-02  8.25E-01  6.65E-01  2.35E+00
 


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
+        1.18E+03
 
 TH 2
+       -3.02E+01  3.59E+02
 
 TH 3
+        3.07E+01  2.62E+02  1.23E+03
 
 TH 4
+       -3.70E+01  1.99E+02 -2.86E+02  6.02E+02
 
 TH 5
+        7.17E+00 -5.41E+02 -1.74E+03  2.54E+02  2.69E+03
 
 TH 6
+        3.94E-01 -6.10E+00  1.13E+01 -1.22E+01 -5.25E+00  2.01E+02
 
 TH 7
+        2.32E+00  3.19E+01 -7.21E+00 -8.00E+00  1.50E+00  1.64E-01  1.18E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.96E+00 -1.55E+01 -1.81E+01 -2.05E+01  4.29E+01 -1.87E-01  7.13E+00  0.00E+00  1.18E+02
 
 TH10
+       -4.46E+00  6.22E-01 -5.76E+01 -2.23E+01  2.20E+01 -8.65E-01  1.57E+00  0.00E+00  9.80E-01  5.57E+01
 
 TH11
+       -1.51E+01 -5.46E+00 -3.09E+01 -1.10E+01  1.07E+01  4.38E+00  1.90E+00  0.00E+00  1.22E+01  2.15E+01  4.82E+01
 
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
 #CPUT: Total CPU Time in Seconds,       19.700
Stop Time:
Wed Sep 29 12:38:14 CDT 2021
