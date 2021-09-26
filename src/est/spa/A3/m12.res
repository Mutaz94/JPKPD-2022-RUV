Sat Sep 25 09:04:30 CDT 2021
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
$DATA ../../../../data/spa/A3/dat12.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m12.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   214.380110694884        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.1854E+02 -1.4210E+00  8.7798E+01 -1.3481E+02  5.6066E+01  2.5681E+01 -5.5495E+01 -1.1313E+01 -1.1266E+02 -9.9805E+01
            -3.4333E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1204.81314243628        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0707E+00  1.0142E+00  9.9283E-01  1.1659E+00  1.1001E+00  7.3831E-01  9.7870E-01  9.3395E-01  1.0226E+00  9.3548E-01
             5.4268E+00
 PARAMETER:  1.6829E-01  1.1406E-01  9.2809E-02  2.5347E-01  1.9542E-01 -2.0339E-01  7.8466E-02  3.1670E-02  1.2240E-01  3.3301E-02
             1.7913E+00
 GRADIENT:   1.3978E+02 -2.9048E+01 -2.5170E+01 -1.9376E+01  1.2105E+01 -2.6801E+01  9.9823E+00  5.9836E+00  1.9580E+01  1.4256E+01
             1.9132E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1241.32143944413        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0055E+00  8.6622E-01  6.9114E-01  1.2250E+00  7.4831E-01  8.3656E-01  3.8719E-01  7.4828E-02  1.0429E+00  6.0957E-01
             4.6543E+00
 PARAMETER:  1.0553E-01 -4.3612E-02 -2.6941E-01  3.0293E-01 -1.8994E-01 -7.8462E-02 -8.4885E-01 -2.4926E+00  1.4199E-01 -3.9501E-01
             1.6378E+00
 GRADIENT:  -2.9641E+01  1.1313E+01 -1.3878E+01  3.4062E+01 -6.7045E+00  1.9605E-01  9.1816E-01  1.0140E-01  8.7107E+00  1.3933E+01
             1.2440E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1260.47558646117        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  9.9151E-01  7.1258E-01  9.3469E-01  1.2821E+00  8.3767E-01  8.4686E-01  5.6193E-01  1.0551E-01  9.8467E-01  2.8599E-01
             3.8858E+00
 PARAMETER:  9.1478E-02 -2.3887E-01  3.2455E-02  3.4854E-01 -7.7129E-02 -6.6221E-02 -4.7637E-01 -2.1490E+00  8.4555E-02 -1.1518E+00
             1.4573E+00
 GRADIENT:  -4.7352E+00  4.6377E+00  5.3524E+00  2.7739E-01 -1.0836E+01  2.4896E+00  9.3950E-01  1.6393E-01  7.1088E-01  2.1377E+00
            -2.7658E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1262.17639060830        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  9.9045E-01  5.0318E-01  9.2500E-01  1.4018E+00  7.7231E-01  8.4157E-01  2.1649E-01  9.3767E-02  9.1284E-01  5.7952E-02
             3.9053E+00
 PARAMETER:  9.0402E-02 -5.8681E-01  2.2040E-02  4.3777E-01 -1.5837E-01 -7.2480E-02 -1.4302E+00 -2.2669E+00  8.8017E-03 -2.7481E+00
             1.4623E+00
 GRADIENT:   3.5028E-01  2.5575E-01  7.0012E-01 -1.1982E-01 -1.0157E+00  9.9033E-01  6.1267E-02  1.5545E-01  1.9757E-01  8.8133E-02
             1.3103E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1262.36812730351        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  9.8863E-01  3.5193E-01  8.7255E-01  1.4793E+00  6.9956E-01  8.3754E-01  5.7087E-02  1.1282E-02  8.6138E-01  1.0000E-02
             3.8981E+00
 PARAMETER:  8.8561E-02 -9.4434E-01 -3.6337E-02  4.9154E-01 -2.5731E-01 -7.7284E-02 -2.7632E+00 -4.3846E+00 -4.9218E-02 -4.7404E+00
             1.4605E+00
 GRADIENT:   2.2995E+00 -6.5758E-02  4.8464E-01 -5.3327E-01 -5.7074E-01 -6.0707E-01  1.9774E-03  2.5430E-03 -1.7599E-01  0.0000E+00
             8.2324E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1262.37991331387        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      553
 NPARAMETR:  9.8798E-01  3.1995E-01  8.7764E-01  1.5015E+00  6.9524E-01  8.3933E-01  3.8153E-02  1.0000E-02  8.4955E-01  1.0000E-02
             3.8987E+00
 PARAMETER:  8.7909E-02 -1.0396E+00 -3.0519E-02  5.0649E-01 -2.6349E-01 -7.5153E-02 -3.1662E+00 -5.0260E+00 -6.3044E-02 -5.3065E+00
             1.4606E+00
 GRADIENT:  -6.0387E-01  1.2458E-01  8.1134E-02  7.0597E-01 -1.4433E-01 -1.3092E-01  6.7022E-04  0.0000E+00 -1.3745E-01  0.0000E+00
            -1.9265E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1262.38064058112        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:      738             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8796E-01  3.0686E-01  8.7026E-01  1.5071E+00  6.8770E-01  8.3976E-01  1.7795E-02  1.0000E-02  8.4662E-01  1.0000E-02
             3.8984E+00
 PARAMETER:  8.7890E-02 -1.0814E+00 -3.8960E-02  5.1020E-01 -2.7441E-01 -7.4644E-02 -3.9288E+00 -5.3669E+00 -6.6499E-02 -5.5851E+00
             1.4606E+00
 GRADIENT:   2.4410E+00  1.8076E-01  8.7055E-02  3.4847E+00  1.0943E-01  1.6663E-01  1.5640E-04  0.0000E+00  5.1154E-02  0.0000E+00
             1.1850E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1262.38068351289        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      814
 NPARAMETR:  9.8795E-01  3.0716E-01  8.7022E-01  1.5069E+00  6.8773E-01  8.3974E-01  1.0000E-02  1.0000E-02  8.4672E-01  1.0000E-02
             3.8983E+00
 PARAMETER:  8.7873E-02 -1.0804E+00 -3.9012E-02  5.1003E-01 -2.7436E-01 -7.4662E-02 -5.3565E+00 -5.3669E+00 -6.6389E-02 -5.5851E+00
             1.4605E+00
 GRADIENT:   2.3971E+00  1.7449E-01  1.0201E-01  3.3906E+00  9.3027E-02  1.6252E-01  0.0000E+00  0.0000E+00  4.3906E-02  0.0000E+00
             1.1581E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1262.38069426567        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:      978
 NPARAMETR:  9.8792E-01  3.0705E-01  8.7018E-01  1.5070E+00  6.8775E-01  8.3973E-01  1.0000E-02  1.0000E-02  8.4671E-01  1.0000E-02
             3.8983E+00
 PARAMETER:  8.7849E-02 -1.0807E+00 -3.9058E-02  5.1015E-01 -2.7433E-01 -7.4677E-02 -4.6879E+00 -5.3669E+00 -6.6395E-02 -5.5851E+00
             1.4605E+00
 GRADIENT:  -7.7678E-03 -1.1229E-03 -2.2556E-03 -4.2996E-03 -3.7810E-03 -6.3378E-04  0.0000E+00  0.0000E+00 -6.6735E-04  0.0000E+00
            -2.6641E-03

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1262.38070541036        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1157
 NPARAMETR:  9.8798E-01  3.0956E-01  8.7111E-01  1.5058E+00  6.8892E-01  8.3973E-01  1.0000E-02  1.0000E-02  8.4749E-01  1.0000E-02
             3.8984E+00
 PARAMETER:  8.7903E-02 -1.0726E+00 -3.7984E-02  5.0933E-01 -2.7263E-01 -7.4672E-02 -7.1016E+00 -5.3669E+00 -6.5475E-02 -5.5851E+00
             1.4606E+00
 GRADIENT:   8.1666E-04 -1.0754E-03 -5.9259E-03 -6.5211E-03  4.0690E-03  2.3737E-03  0.0000E+00  0.0000E+00  5.0683E-03  0.0000E+00
            -5.7719E-03

0ITERATION NO.:   51    OBJECTIVE VALUE:  -1262.38070541036        NO. OF FUNC. EVALS.:  37
 CUMULATIVE NO. OF FUNC. EVALS.:     1194
 NPARAMETR:  9.8798E-01  3.0958E-01  8.7114E-01  1.5058E+00  6.8891E-01  8.3971E-01  1.0000E-02  1.0000E-02  8.4741E-01  1.0000E-02
             3.8985E+00
 PARAMETER:  8.7903E-02 -1.0726E+00 -3.7984E-02  5.0933E-01 -2.7263E-01 -7.4672E-02 -7.1016E+00 -5.3669E+00 -6.5475E-02 -5.5851E+00
             1.4606E+00
 GRADIENT:   6.9041E-04 -1.1829E-03 -5.6437E-03 -2.3619E-03  3.9553E-03  2.4716E-03  0.0000E+00  0.0000E+00  5.4861E-03  0.0000E+00
            -5.1939E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1194
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.4744E-04 -1.0019E-04  1.0317E-04 -1.3963E-02 -1.5587E-05
 SE:             2.7959E-02  4.5244E-05  1.1873E-04  2.4107E-02  1.7191E-04
 N:                     100         100         100         100         100

 P VAL.:         9.8153E-01  2.6806E-02  3.8488E-01  5.6243E-01  9.2776E-01

 ETASHRINKSD(%)  6.3336E+00  9.9848E+01  9.9602E+01  1.9239E+01  9.9424E+01
 ETASHRINKVR(%)  1.2266E+01  1.0000E+02  9.9998E+01  3.4777E+01  9.9997E+01
 EBVSHRINKSD(%)  6.1036E+00  9.9851E+01  9.9491E+01  1.8888E+01  9.9315E+01
 EBVSHRINKVR(%)  1.1835E+01  1.0000E+02  9.9997E+01  3.4208E+01  9.9995E+01
 RELATIVEINF(%)  7.8880E+01  6.8693E-06  1.0838E-04  4.0267E+00  1.0423E-04
 EPSSHRINKSD(%)  1.9332E+01
 EPSSHRINKVR(%)  3.4927E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1262.3807054103629     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -527.22987884662473     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.44
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.08
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1262.381       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.88E-01  3.10E-01  8.71E-01  1.51E+00  6.89E-01  8.40E-01  1.00E-02  1.00E-02  8.47E-01  1.00E-02  3.90E+00
 


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
+        1.43E+03
 
 TH 2
+       -1.11E+02  2.75E+02
 
 TH 3
+        8.34E+00  1.53E+02  3.52E+02
 
 TH 4
+       -1.16E+02  2.89E+02  1.97E+01  4.44E+02
 
 TH 5
+        4.86E+01 -3.65E+02 -5.86E+02 -1.57E+02  1.08E+03
 
 TH 6
+       -3.19E+00 -1.54E+01  1.33E+01 -2.63E+01 -6.59E+00  2.23E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.02E+01 -5.31E+01  6.95E-01 -1.10E+01  2.81E+01  1.63E+00  0.00E+00  0.00E+00  1.13E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.25E+01 -1.26E+01 -4.49E+00 -9.69E+00  1.01E+01  5.86E+00  0.00E+00  0.00E+00  1.28E+01  0.00E+00  2.92E+01
 
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
 #CPUT: Total CPU Time in Seconds,       19.580
Stop Time:
Sat Sep 25 09:04:51 CDT 2021
