Sat Sep 25 12:44:28 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat19.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m19.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1739.55908947512        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -2.0705E+01 -4.9083E+01 -7.2268E+00 -4.0879E+01  3.1189E+00  1.9567E+01  6.4973E+00  5.1973E+00  4.4391E+01 -1.4996E+00
             2.7480E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1748.85394499899        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0193E+00  1.0913E+00  9.9035E-01  9.6968E-01  1.0070E+00  9.6563E-01  9.9078E-01  9.7103E-01  7.5402E-01  1.0675E+00
             9.7677E-01
 PARAMETER:  1.1912E-01  1.8738E-01  9.0305E-02  6.9211E-02  1.0696E-01  6.5023E-02  9.0739E-02  7.0603E-02 -1.8234E-01  1.6534E-01
             7.6494E-02
 GRADIENT:   2.9324E+01  9.2459E-01  1.9019E+01 -2.5650E+01 -4.1235E+01  7.9693E+00 -4.4905E+00  1.7623E+00  1.5842E+00  4.9408E+00
             1.5323E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1749.79930853691        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0135E+00  9.4233E-01  9.2603E-01  1.0617E+00  9.2828E-01  9.9291E-01  1.2155E+00  7.3691E-01  6.5343E-01  1.0241E+00
             9.5810E-01
 PARAMETER:  1.1339E-01  4.0596E-02  2.3147E-02  1.5990E-01  2.5577E-02  9.2889E-02  2.9517E-01 -2.0529E-01 -3.2552E-01  1.2385E-01
             5.7197E-02
 GRADIENT:   1.6997E+01  1.5659E+00  4.3410E+00  3.9201E+00 -1.7059E+01  1.8583E+01 -2.7205E-01  1.6103E+00 -1.3339E+00  6.0397E+00
             8.9458E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1750.21584850492        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      253
 NPARAMETR:  1.0103E+00  1.0348E+00  8.2526E-01  9.9978E-01  9.2560E-01  9.5646E-01  1.1201E+00  5.7895E-01  6.9728E-01  9.8433E-01
             9.4581E-01
 PARAMETER:  1.1028E-01  1.3420E-01 -9.2057E-02  9.9781E-02  2.2684E-02  5.5482E-02  2.1342E-01 -4.4655E-01 -2.6057E-01  8.4203E-02
             4.4286E-02
 GRADIENT:  -5.1705E+01 -5.9147E+00 -2.5995E+00 -6.4417E+00  4.1382E-01 -1.8301E+00 -6.2684E-01  4.5310E-01 -1.1766E+00  1.1682E+00
             3.1467E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1750.87437057987        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      431
 NPARAMETR:  1.0325E+00  1.1453E+00  7.6770E-01  9.3439E-01  9.4541E-01  9.6051E-01  1.0241E+00  5.0282E-01  7.4594E-01  9.8744E-01
             9.3948E-01
 PARAMETER:  1.3202E-01  2.3568E-01 -1.6436E-01  3.2136E-02  4.3859E-02  5.9710E-02  1.2381E-01 -5.8753E-01 -1.9312E-01  8.7357E-02
             3.7573E-02
 GRADIENT:  -1.3629E+00  1.9729E-01 -3.8993E-02  6.6954E-01 -3.5834E-01  3.8577E-01 -8.5412E-02  4.0240E-02 -1.3811E-01  1.2823E-01
            -3.0699E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1750.90706511318        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      607
 NPARAMETR:  1.0340E+00  1.2845E+00  6.7701E-01  8.4489E-01  9.6536E-01  9.6133E-01  9.2408E-01  3.6626E-01  8.0695E-01  9.8701E-01
             9.4037E-01
 PARAMETER:  1.3344E-01  3.5040E-01 -2.9006E-01 -6.8545E-02  6.4741E-02  6.0565E-02  2.1043E-02 -9.0441E-01 -1.1450E-01  8.6924E-02
             3.8518E-02
 GRADIENT:  -2.9366E-01  1.5848E+00  2.2238E-01  1.1586E+00 -6.4562E-01  1.3029E-01 -9.2468E-02  4.3201E-02 -5.3915E-02  5.6875E-03
            -3.6103E-03

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1750.91061034832        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      782
 NPARAMETR:  1.0343E+00  1.3192E+00  6.5248E-01  8.2082E-01  9.7134E-01  9.6128E-01  9.0290E-01  3.0884E-01  8.2433E-01  9.8775E-01
             9.4047E-01
 PARAMETER:  1.3376E-01  3.7703E-01 -3.2697E-01 -9.7454E-02  7.0917E-02  6.0513E-02 -2.1463E-03 -1.0749E+00 -9.3179E-02  8.7674E-02
             3.8626E-02
 GRADIENT:   4.4943E-02 -1.0413E+00 -2.6508E-01 -7.5049E-01  5.6134E-01 -1.3286E-02  4.2003E-02  4.1005E-02  5.2699E-02  3.2946E-02
             8.9730E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1750.91860314148        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      961
 NPARAMETR:  1.0345E+00  1.3286E+00  6.3157E-01  8.1393E-01  9.6223E-01  9.6157E-01  8.9704E-01  1.7489E-01  8.2939E-01  9.7910E-01
             9.4037E-01
 PARAMETER:  1.3394E-01  3.8413E-01 -3.5955E-01 -1.0588E-01  6.1494E-02  6.0813E-02 -8.6599E-03 -1.6436E+00 -8.7071E-02  7.8875E-02
             3.8523E-02
 GRADIENT:   5.5727E-02  3.0660E-01 -1.8472E-02  3.0429E-01 -6.2388E-02  7.3166E-03 -1.5893E-02  4.3547E-03 -2.2718E-02  4.6393E-03
             1.0424E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1750.92012942341        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1136
 NPARAMETR:  1.0346E+00  1.3436E+00  6.1993E-01  8.0371E-01  9.6375E-01  9.6170E-01  8.8822E-01  8.1042E-02  8.3712E-01  9.7854E-01
             9.4039E-01
 PARAMETER:  1.3402E-01  3.9534E-01 -3.7814E-01 -1.1852E-01  6.3079E-02  6.0945E-02 -1.8540E-02 -2.4128E+00 -7.7785E-02  7.8307E-02
             3.8541E-02
 GRADIENT:   6.7071E-02 -1.4310E-01 -2.6989E-02 -4.7343E-02  1.4385E-01  3.8819E-03 -2.1469E-02  6.4197E-04 -1.9833E-02 -2.7416E-02
             9.8135E-03

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1750.92046930098        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1311
 NPARAMETR:  1.0346E+00  1.3477E+00  6.1636E-01  8.0097E-01  9.6373E-01  9.6173E-01  8.8597E-01  1.8185E-02  8.3926E-01  9.7813E-01
             9.4034E-01
 PARAMETER:  1.3402E-01  3.9837E-01 -3.8392E-01 -1.2193E-01  6.3058E-02  6.0982E-02 -2.1067E-02 -3.9072E+00 -7.5241E-02  7.7884E-02
             3.8490E-02
 GRADIENT:  -5.2027E-03 -9.6078E-02 -7.1890E-02 -1.0368E-02  5.5321E-02  3.3063E-04  1.3260E-02  5.0254E-05  9.7237E-03  2.7147E-02
             1.1073E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1750.92049474215        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     1479
 NPARAMETR:  1.0347E+00  1.3473E+00  6.1667E-01  8.0121E-01  9.6369E-01  9.6175E-01  8.8605E-01  1.0000E-02  8.3907E-01  9.7811E-01
             9.4035E-01
 PARAMETER:  1.3402E-01  3.9813E-01 -3.8344E-01 -1.2162E-01  6.3019E-02  6.0977E-02 -2.0877E-02 -5.0543E+00 -7.5419E-02  7.7878E-02
             3.8496E-02
 GRADIENT:  -3.7110E-02  9.7067E-03 -1.2380E-03  1.5508E-03  2.7967E-07 -1.1054E-03  1.8361E-03  0.0000E+00  6.5734E-04  2.8230E-04
            -6.5742E-05

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1479
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.4503E-04 -1.4130E-02 -3.9616E-04  8.8902E-03 -2.1926E-02
 SE:             2.9843E-02  2.3327E-02  1.5797E-04  2.2337E-02  2.3441E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9612E-01  5.4470E-01  1.2148E-02  6.9063E-01  3.4960E-01

 ETASHRINKSD(%)  2.1735E-02  2.1852E+01  9.9471E+01  2.5168E+01  2.1470E+01
 ETASHRINKVR(%)  4.3466E-02  3.8928E+01  9.9997E+01  4.4002E+01  3.8330E+01
 EBVSHRINKSD(%)  4.0945E-01  2.1670E+01  9.9551E+01  2.6684E+01  1.9478E+01
 EBVSHRINKVR(%)  8.1723E-01  3.8644E+01  9.9998E+01  4.6248E+01  3.5162E+01
 RELATIVEINF(%)  9.8999E+01  2.5933E+00  1.6986E-04  2.1107E+00  7.7320E+00
 EPSSHRINKSD(%)  4.4171E+01
 EPSSHRINKVR(%)  6.8831E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1750.9204947421483     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1015.7696681784101     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.04
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.73
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1750.920       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.35E+00  6.17E-01  8.01E-01  9.64E-01  9.62E-01  8.86E-01  1.00E-02  8.39E-01  9.78E-01  9.40E-01
 


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
+        1.11E+03
 
 TH 2
+       -7.01E+00  4.95E+02
 
 TH 3
+        1.65E+01  1.90E+02  6.68E+02
 
 TH 4
+       -2.10E+01  4.83E+02 -4.67E+02  1.33E+03
 
 TH 5
+       -2.57E+00 -2.59E+02 -5.98E+02  4.03E+02  8.12E+02
 
 TH 6
+       -6.35E-01 -1.65E+00  5.64E+00 -6.37E+00 -3.04E+00  2.11E+02
 
 TH 7
+       -1.65E-01  2.39E+01 -2.20E+01 -1.37E+01  5.04E-01 -2.12E+00  1.03E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -3.57E-01 -2.02E+01 -4.06E+01  4.28E+01  1.22E+01 -1.68E+00  3.47E+01  0.00E+00  9.27E+01
 
 TH10
+       -7.54E-01 -1.78E+01 -5.82E+01  3.20E-01 -6.11E+01  8.57E-01  8.58E+00  0.00E+00  1.78E+01  8.75E+01
 
 TH11
+       -7.87E+00 -1.75E+01 -3.98E+01  5.37E+00  6.15E+00  3.36E+00  9.30E+00  0.00E+00  1.13E+01  1.54E+01  2.40E+02
 
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
 #CPUT: Total CPU Time in Seconds,       22.827
Stop Time:
Sat Sep 25 12:44:52 CDT 2021
