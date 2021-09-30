Thu Sep 30 08:50:55 CDT 2021
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
$DATA ../../../../data/spa2/D/dat26.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m26.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   11982.0001308543        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.1296E+02  1.6689E+02 -1.2720E+01 -5.2204E+01  2.1319E+02 -9.4117E+02 -4.9513E+02 -2.6077E+02 -8.6237E+02 -3.4017E+02
            -2.5805E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -806.555802589814        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2690E+00  1.3456E+00  9.8367E-01  1.9673E+00  1.1238E+00  2.4307E+00  2.4200E+00  9.7606E-01  2.3408E+00  1.2953E+00
             1.2976E+01
 PARAMETER:  3.3819E-01  3.9688E-01  8.3539E-02  7.7667E-01  2.1667E-01  9.8817E-01  9.8376E-01  7.5773E-02  9.5050E-01  3.5875E-01
             2.6631E+00
 GRADIENT:  -4.5730E+01 -2.9131E+01 -5.5528E+01  7.9474E+01  4.8353E+01  6.9946E+01 -3.9975E+00  2.9116E+00  3.0125E+01  1.8026E+01
             4.5483E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -860.740219240126        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.4338E+00  1.6614E+00  1.2275E+01  2.5213E+00  4.9614E+00  1.7153E+00  1.0150E+01  7.4627E-01  3.3601E+00  3.4911E+00
             1.1690E+01
 PARAMETER:  4.6034E-01  6.0769E-01  2.6075E+00  1.0248E+00  1.7017E+00  6.3959E-01  2.4174E+00 -1.9267E-01  1.3120E+00  1.3502E+00
             2.5587E+00
 GRADIENT:   6.8286E+01  3.1213E+01  1.5459E+00  8.0800E+00 -3.1984E+01  3.8676E+00  1.0067E+02  1.7600E-02  7.9839E+01  2.5008E+01
             3.6887E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1039.24981283294        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  9.8366E-01  5.6282E-01  9.9979E+00  1.7556E+00  2.0903E+00  2.4840E+00  3.5622E+00  4.0191E-01  1.4624E+00  1.0696E+00
             7.0703E+00
 PARAMETER:  8.3529E-02 -4.7480E-01  2.4024E+00  6.6280E-01  8.3732E-01  1.0099E+00  1.3704E+00 -8.1152E-01  4.8008E-01  1.6730E-01
             2.0559E+00
 GRADIENT:  -6.8316E+01 -1.7424E+01 -2.5207E+00  3.5833E+01  1.4914E+01  9.5784E+01 -4.9065E+01  5.2339E-02 -3.0819E+01  1.3764E+01
             1.1444E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1069.77878715183        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  1.1252E+00  7.4832E-01  7.2779E+00  1.3892E+00  1.9255E+00  1.7814E+00  3.8385E+00  1.4192E+00  1.5860E+00  7.7579E-01
             6.6154E+00
 PARAMETER:  2.1799E-01 -1.8992E-01  2.0848E+00  4.2875E-01  7.5518E-01  6.7743E-01  1.4451E+00  4.5012E-01  5.6119E-01 -1.5387E-01
             1.9894E+00
 GRADIENT:  -5.2773E+00 -2.5595E+01 -3.1148E+00 -1.9470E+01  1.1673E+01  4.0906E+00  8.0185E+00  4.5523E-01  9.7575E+00  7.9872E+00
             1.8542E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1072.37091092187        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      377
 NPARAMETR:  1.1695E+00  1.1006E+00  5.3369E+00  1.2142E+00  1.9059E+00  1.7535E+00  3.4122E+00  1.3479E+00  1.6029E+00  6.5040E-01
             6.6573E+00
 PARAMETER:  2.5653E-01  1.9581E-01  1.7746E+00  2.9412E-01  7.4493E-01  6.6162E-01  1.3274E+00  3.9852E-01  5.7182E-01 -3.3016E-01
             1.9957E+00
 GRADIENT:   2.1231E+01 -1.8398E+01 -3.6793E+00 -7.8906E+00  8.6464E+00 -4.1997E+00  7.1583E+00  4.4431E-01  1.0201E+01  5.6530E+00
             1.9267E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1079.49861832070        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      558
 NPARAMETR:  1.1978E+00  1.0701E+00  6.6760E+00  1.3160E+00  1.9354E+00  1.9384E+00  4.2295E+00  3.4046E-01  1.1785E+00  4.4046E-01
             6.6644E+00
 PARAMETER:  2.8047E-01  1.6776E-01  1.9985E+00  3.7463E-01  7.6031E-01  7.6188E-01  1.5421E+00 -9.7746E-01  2.6426E-01 -7.1994E-01
             1.9968E+00
 GRADIENT:   4.7516E+00  2.5778E+00 -1.3173E+00  8.9320E+00 -2.5076E+00  3.8831E+00 -2.5984E-01  2.1516E-02 -5.8946E+00  2.1983E+00
            -1.1207E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1081.37467983635        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      734
 NPARAMETR:  1.1654E+00  1.0307E+00  1.6167E+01  1.3767E+00  2.1068E+00  1.9407E+00  4.2709E+00  1.6265E-02  1.3664E+00  1.4189E-01
             6.6651E+00
 PARAMETER:  2.5306E-01  1.3028E-01  2.8829E+00  4.1967E-01  8.4519E-01  7.6304E-01  1.5518E+00 -4.0188E+00  4.1217E-01 -1.8527E+00
             1.9969E+00
 GRADIENT:  -1.0643E+01  2.5778E+00 -2.2947E-01  1.2093E+00  1.7811E+00  4.0173E+00  1.6343E+00  6.7237E-06 -1.0002E+00  2.2104E-01
            -5.6224E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1082.57631689049        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      910
 NPARAMETR:  1.1895E+00  4.8445E-01  6.4726E+01  1.8352E+00  2.1713E+00  1.9072E+00  5.3550E+00  1.0000E-02  1.6831E+00  6.3134E-02
             6.7195E+00
 PARAMETER:  2.7355E-01 -6.2474E-01  4.2702E+00  7.0717E-01  8.7533E-01  7.4566E-01  1.7780E+00 -7.2528E+00  6.2066E-01 -2.6625E+00
             2.0050E+00
 GRADIENT:  -4.1622E-01  3.1009E-01  1.1262E-01  2.8279E+00  1.4377E+00 -1.7731E-01 -3.0677E-01  0.0000E+00  5.6956E-01  4.3867E-02
            -3.4547E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1082.72114240666        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1094
 NPARAMETR:  1.1920E+00  4.8736E-01  2.3932E+01  1.8124E+00  2.0995E+00  1.9133E+00  5.3081E+00  1.0000E-02  1.6615E+00  1.5012E-02
             6.7215E+00
 PARAMETER:  2.7559E-01 -6.1876E-01  3.2752E+00  6.9465E-01  8.4169E-01  7.4882E-01  1.7692E+00 -7.0619E+00  6.0770E-01 -4.0989E+00
             2.0053E+00
 GRADIENT:   1.1950E+00 -2.8125E-01 -1.2314E-01  2.7297E-01  2.8304E+00  8.7948E-01 -7.8619E-01  0.0000E+00 -8.8049E-01  2.4425E-03
             5.5547E-01

0ITERATION NO.:   49    OBJECTIVE VALUE:  -1082.75899238229        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1221
 NPARAMETR:  1.1886E+00  4.8292E-01  2.1523E+01  1.8134E+00  2.0591E+00  1.9090E+00  5.3277E+00  1.0000E-02  1.6781E+00  1.0000E-02
             6.7163E+00
 PARAMETER:  2.7277E-01 -6.2790E-01  3.1691E+00  6.9519E-01  8.2226E-01  7.4657E-01  1.7729E+00 -7.0619E+00  6.1769E-01 -4.6959E+00
             2.0045E+00
 GRADIENT:  -4.0543E-02 -2.3235E-02 -4.6626E-03 -2.1873E-01 -1.6056E-02 -1.9259E-03 -8.2400E-02  0.0000E+00 -6.0611E-02  0.0000E+00
             5.2730E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1221
 NO. OF SIG. DIGITS IN FINAL EST.:  2.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.7299E-03  4.8830E-02  3.7444E-06 -6.5323E-02 -2.5515E-05
 SE:             2.9272E-02  2.0927E-02  4.5514E-06  1.8826E-02  8.1616E-05
 N:                     100         100         100         100         100

 P VAL.:         8.1816E-01  1.9628E-02  4.1068E-01  5.2091E-04  7.5457E-01

 ETASHRINKSD(%)  1.9349E+00  2.9892E+01  9.9985E+01  3.6930E+01  9.9727E+01
 ETASHRINKVR(%)  3.8324E+00  5.0849E+01  1.0000E+02  6.0222E+01  9.9999E+01
 EBVSHRINKSD(%)  3.1344E+00  3.9785E+01  9.9977E+01  2.5811E+01  9.9639E+01
 EBVSHRINKVR(%)  6.1705E+00  6.3741E+01  1.0000E+02  4.4960E+01  9.9999E+01
 RELATIVEINF(%)  9.2609E+01  1.2673E+01  6.7809E-07  1.8452E+01  1.6087E-04
 EPSSHRINKSD(%)  1.3117E+01
 EPSSHRINKVR(%)  2.4514E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1082.7589923822884     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       19.967247463318699     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.19
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    11.70
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1082.759       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.19E+00  4.83E-01  2.15E+01  1.81E+00  2.06E+00  1.91E+00  5.33E+00  1.00E-02  1.68E+00  1.00E-02  6.72E+00
 


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
+        1.98E+02
 
 TH 2
+       -1.80E+00  7.07E+01
 
 TH 3
+        1.05E-02  1.85E-02  2.09E-03
 
 TH 4
+       -5.08E+00  4.31E+01  2.14E-02  7.06E+01
 
 TH 5
+       -3.16E+00 -1.03E+01 -1.98E-01 -8.44E+00  3.10E+01
 
 TH 6
+       -9.36E-01 -6.88E-01  8.10E-03 -1.18E-01 -5.17E-01  4.67E+01
 
 TH 7
+        4.36E-01  9.07E+00 -4.80E-03 -2.60E+00 -1.07E-01 -1.04E-01  3.12E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -4.35E-01 -1.45E+00 -6.98E-04 -1.77E+01  3.46E+00 -1.83E+00  1.87E+00  0.00E+00  2.08E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -6.57E+00 -2.80E+00  2.40E-03 -5.43E+00  1.18E-01  1.86E+00  1.99E-01  0.00E+00  1.53E+00  0.00E+00  1.53E+01
 
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
 #CPUT: Total CPU Time in Seconds,       35.978
Stop Time:
Thu Sep 30 08:51:32 CDT 2021
