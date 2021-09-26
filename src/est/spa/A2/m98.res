Sat Sep 25 08:59:27 CDT 2021
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
$DATA ../../../../data/spa/A2/dat98.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m98.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -715.909063519614        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.2944E+02 -1.8962E+01  2.0010E+00 -3.4062E+01  4.2530E+01  1.9885E+01 -1.7660E+01 -6.4266E+00 -4.6757E+01 -2.5934E+01
            -1.7430E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1279.83383095042        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.7167E-01  9.7451E-01  1.2226E+00  1.0594E+00  1.0706E+00  8.8486E-01  1.0672E+00  8.4923E-01  1.1848E+00  1.0896E+00
             1.8336E+00
 PARAMETER:  7.1256E-02  7.4182E-02  3.0097E-01  1.5772E-01  1.6824E-01 -2.2329E-02  1.6506E-01 -6.3428E-02  2.6960E-01  1.8585E-01
             7.0626E-01
 GRADIENT:   3.6565E+01  2.9764E+00 -1.4709E+00  4.7267E+00  9.6999E+00 -1.8008E+01  3.7716E+00  2.1725E+00  9.4536E+00 -1.1721E+01
            -3.6713E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1314.83752850109        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.7417E-01  7.3823E-01  1.2388E+00  1.2227E+00  9.3921E-01  8.4781E-01  3.6073E-01  1.0687E-01  1.2244E+00  1.0203E+00
             2.0638E+00
 PARAMETER:  7.3835E-02 -2.0351E-01  3.1415E-01  3.0107E-01  3.7284E-02 -6.5098E-02 -9.1963E-01 -2.1362E+00  3.0248E-01  1.2007E-01
             8.2455E-01
 GRADIENT:   4.2756E+01  2.2198E+01  1.7883E+01  3.7396E+01 -3.4114E+01 -3.2773E+01 -1.6772E-01  5.2227E-02  2.1394E+01 -9.0786E+00
            -2.4641E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1356.52994203963        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.7223E-01  6.5010E-01  1.6229E+00  1.2806E+00  1.0943E+00  8.9921E-01  5.2420E-01  1.4816E-02  8.3182E-01  1.0028E+00
             2.9565E+00
 PARAMETER:  7.1835E-02 -3.3063E-01  5.8419E-01  3.4732E-01  1.9014E-01 -6.2375E-03 -5.4587E-01 -4.1121E+00 -8.4143E-02  1.0279E-01
             1.1840E+00
 GRADIENT:   7.4915E+00 -9.5390E+00 -3.5211E+00 -4.1125E+01 -7.0068E+00  2.6158E+00 -7.0963E-01  1.0595E-03 -1.2750E+01  2.1706E+01
            -8.3182E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1365.98027964490        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  9.6587E-01  4.7689E-01  2.1187E+00  1.4189E+00  1.1397E+00  8.8582E-01  1.2541E+00  3.3428E-02  7.7159E-01  3.3620E-01
             3.0321E+00
 PARAMETER:  6.5275E-02 -6.4047E-01  8.5079E-01  4.4989E-01  2.3074E-01 -2.1238E-02  3.2639E-01 -3.2983E+00 -1.5930E-01 -9.9006E-01
             1.2093E+00
 GRADIENT:  -3.3169E+00 -1.5425E+00  3.1827E+00 -8.4730E+00 -3.4457E+00  8.7837E-02  6.1795E-01  3.0572E-03  2.0125E+00  1.4113E+00
             1.2093E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1366.76633349887        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  9.6575E-01  4.1453E-01  1.7157E+00  1.4506E+00  1.0369E+00  8.8669E-01  1.3843E+00  1.7954E-02  7.6320E-01  1.1417E-01
             3.0391E+00
 PARAMETER:  6.5154E-02 -7.8060E-01  6.3982E-01  4.7201E-01  1.3622E-01 -2.0261E-02  4.2522E-01 -3.9199E+00 -1.7024E-01 -2.0700E+00
             1.2116E+00
 GRADIENT:  -1.1925E+00 -3.0665E-01 -1.4702E-01  6.2622E-02 -3.2911E-02  8.5715E-02  1.0229E-01  1.6055E-03  1.1274E-01  1.3812E-01
             8.2195E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1366.83853233426        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      446
 NPARAMETR:  9.6514E-01  3.8346E-01  1.7434E+00  1.4729E+00  1.0343E+00  8.8392E-01  1.4662E+00  1.6091E-02  7.5414E-01  5.7946E-02
             3.0420E+00
 PARAMETER:  6.4518E-02 -8.5851E-01  6.5586E-01  4.8722E-01  1.3374E-01 -2.3388E-02  4.8267E-01 -4.0295E+00 -1.8218E-01 -2.7482E+00
             1.2125E+00
 GRADIENT:  -2.2834E+00  4.9246E-01  3.0174E-01  4.2326E+00 -1.1005E+00 -8.9713E-01  1.3829E-01  1.2550E-03  1.1187E-01  3.4480E-02
             1.1305E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1366.86622735103        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      561
 NPARAMETR:  9.6730E-01  3.7547E-01  1.7711E+00  1.4806E+00  1.0407E+00  8.8590E-01  1.4679E+00  1.5775E-02  7.5282E-01  4.7889E-02
             3.0430E+00
 PARAMETER:  6.6752E-02 -8.7957E-01  6.7163E-01  4.9244E-01  1.3985E-01 -2.1153E-02  4.8381E-01 -4.0493E+00 -1.8392E-01 -2.9389E+00
             1.2129E+00
 GRADIENT:  -6.2605E-01  2.4853E-01  1.0675E-01 -3.1041E-01 -6.9865E-01 -4.9782E-01  9.2551E-02  1.1503E-03  1.7821E-01  2.2807E-02
             1.5099E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1366.92878900801        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      736
 NPARAMETR:  9.6663E-01  2.6684E-01  1.8931E+00  1.5549E+00  1.0417E+00  8.8697E-01  1.6258E+00  1.0000E-02  7.2719E-01  1.0000E-02
             3.0430E+00
 PARAMETER:  6.6058E-02 -1.2211E+00  7.3822E-01  5.4141E-01  1.4089E-01 -1.9940E-02  5.8601E-01 -4.7433E+00 -2.1856E-01 -6.4691E+00
             1.2128E+00
 GRADIENT:   1.5254E-01  6.1128E-01  2.4513E-01  4.0898E+00 -8.2770E-01  1.4252E-01 -1.2064E-01  0.0000E+00 -7.0896E-01  0.0000E+00
            -4.8220E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1366.97337566555        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      911
 NPARAMETR:  9.6487E-01  1.7070E-01  1.8808E+00  1.6115E+00  1.0130E+00  8.8628E-01  2.1939E+00  1.0000E-02  7.0956E-01  1.0000E-02
             3.0426E+00
 PARAMETER:  6.4242E-02 -1.6679E+00  7.3170E-01  5.7719E-01  1.1293E-01 -2.0727E-02  8.8568E-01 -5.8049E+00 -2.4311E-01 -1.2010E+01
             1.2127E+00
 GRADIENT:  -8.5262E-01  5.3957E-02 -8.5324E-02 -3.0932E-01  6.7279E-02  1.3094E-01  3.0158E-02  0.0000E+00 -1.2437E-02  0.0000E+00
             7.2414E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1366.99386605380        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1087
 NPARAMETR:  9.6450E-01  1.0429E-01  1.9378E+00  1.6564E+00  1.0090E+00  8.8535E-01  2.8158E+00  1.0000E-02  6.9568E-01  1.0000E-02
             3.0429E+00
 PARAMETER:  6.3854E-02 -2.1606E+00  7.6158E-01  6.0466E-01  1.0901E-01 -2.1772E-02  1.1352E+00 -7.1144E+00 -2.6286E-01 -1.8765E+01
             1.2128E+00
 GRADIENT:  -1.6025E-01  1.8129E-01  1.5459E-01  3.0565E+00 -6.2159E-01 -8.7091E-02  1.1527E-02  0.0000E+00 -2.9122E-01  0.0000E+00
            -1.0010E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1367.00606108298        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1263
 NPARAMETR:  9.6405E-01  6.1832E-02  1.9600E+00  1.6825E+00  1.0039E+00  8.8512E-01  3.2939E+00  1.0000E-02  6.8933E-01  1.0000E-02
             3.0422E+00
 PARAMETER:  6.3393E-02 -2.6833E+00  7.7296E-01  6.2025E-01  1.0384E-01 -2.2031E-02  1.2921E+00 -8.7498E+00 -2.7204E-01 -2.6511E+01
             1.2126E+00
 GRADIENT:   2.1261E-01  3.3697E-02  1.8551E-02  1.0831E+00 -1.2889E-01 -5.7503E-02 -3.0386E-02  0.0000E+00 -2.7195E-01  0.0000E+00
            -2.5118E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1367.01430693624        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1438
 NPARAMETR:  9.6357E-01  3.3611E-02  1.9609E+00  1.6999E+00  9.9664E-01  8.8495E-01  4.7029E+00  1.0000E-02  6.8459E-01  1.0000E-02
             3.0426E+00
 PARAMETER:  6.2889E-02 -3.2929E+00  7.7342E-01  6.3056E-01  9.6639E-02 -2.2222E-02  1.6482E+00 -1.0548E+01 -2.7893E-01 -3.5801E+01
             1.2127E+00
 GRADIENT:  -8.6841E-02  2.4324E-02  6.2413E-02  1.3291E+00 -3.2488E-01 -5.1721E-02 -8.1685E-03  0.0000E+00 -1.4649E-01  0.0000E+00
            -3.1357E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1367.01774465238        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1614
 NPARAMETR:  9.6359E-01  1.6669E-02  1.9720E+00  1.7105E+00  9.9537E-01  8.8445E-01  7.6182E+00  1.0000E-02  6.8160E-01  1.0000E-02
             3.0426E+00
 PARAMETER:  6.2908E-02 -3.9942E+00  7.7905E-01  6.3681E-01  9.5362E-02 -2.2793E-02  2.1305E+00 -1.2577E+01 -2.8331E-01 -4.6737E+01
             1.2127E+00
 GRADIENT:   5.6470E-01  1.5470E-02 -3.6746E-02  5.9137E-01  4.2302E-02 -1.9575E-01  1.5276E-02  0.0000E+00  1.9231E-02  0.0000E+00
            -6.1472E-03

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1367.02001987249        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1792
 NPARAMETR:  9.6326E-01  1.0000E-02  1.9719E+00  1.7143E+00  9.9357E-01  8.8484E-01  9.2733E+00  1.0000E-02  6.8070E-01  1.0000E-02
             3.0425E+00
 PARAMETER:  6.2572E-02 -4.5193E+00  7.7902E-01  6.3903E-01  9.3554E-02 -2.2349E-02  2.3271E+00 -1.4300E+01 -2.8464E-01 -5.5064E+01
             1.2127E+00
 GRADIENT:  -2.2171E-02  0.0000E+00 -1.5364E-02  1.5107E-01  8.1746E-03 -2.0976E-02  2.9298E-03  0.0000E+00 -1.1966E-02  0.0000E+00
             3.9068E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1367.02008205911        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1967
 NPARAMETR:  9.6327E-01  1.0000E-02  1.9728E+00  1.7143E+00  9.9374E-01  8.8490E-01  8.8765E+00  1.0000E-02  6.8079E-01  1.0000E-02
             3.0424E+00
 PARAMETER:  6.2576E-02 -4.5097E+00  7.7948E-01  6.3900E-01  9.3720E-02 -2.2282E-02  2.2834E+00 -1.4312E+01 -2.8450E-01 -5.4924E+01
             1.2126E+00
 GRADIENT:   5.3888E-03  4.7491E-06  1.7708E-03 -2.8459E-02 -1.7264E-04  3.6556E-03 -2.9383E-04  0.0000E+00  2.4250E-04  0.0000E+00
            -3.3885E-03

0ITERATION NO.:   77    OBJECTIVE VALUE:  -1367.02008223914        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     2024
 NPARAMETR:  9.6327E-01  1.0000E-02  1.9728E+00  1.7143E+00  9.9373E-01  8.8489E-01  8.9146E+00  1.0000E-02  6.8078E-01  1.0000E-02
             3.0424E+00
 PARAMETER:  6.2575E-02 -4.5100E+00  7.7944E-01  6.3900E-01  9.3706E-02 -2.2290E-02  2.2877E+00 -1.4308E+01 -2.8451E-01 -5.4927E+01
             1.2126E+00
 GRADIENT:  -2.2313E-03  0.0000E+00  7.9046E-05  2.8155E-03 -6.0435E-04 -2.1618E-04 -2.3963E-05  0.0000E+00  1.3157E-04  0.0000E+00
            -1.7419E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2024
 NO. OF SIG. DIGITS IN FINAL EST.:  3.8
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.7248E-04 -1.0315E-03  7.2086E-05 -1.0621E-02 -1.8626E-05
 SE:             2.8764E-02  1.3563E-03  6.6594E-05  2.5612E-02  1.5758E-04
 N:                     100         100         100         100         100

 P VAL.:         9.9244E-01  4.4696E-01  2.7904E-01  6.7836E-01  9.0591E-01

 ETASHRINKSD(%)  3.6365E+00  9.5456E+01  9.9777E+01  1.4196E+01  9.9472E+01
 ETASHRINKVR(%)  7.1408E+00  9.9794E+01  1.0000E+02  2.6377E+01  9.9997E+01
 EBVSHRINKSD(%)  3.6394E+00  9.5437E+01  9.9735E+01  1.3982E+01  9.9432E+01
 EBVSHRINKVR(%)  7.1463E+00  9.9792E+01  9.9999E+01  2.6010E+01  9.9997E+01
 RELATIVEINF(%)  7.7319E+01  9.6979E-04  3.8020E-05  4.1664E-01  1.2636E-04
 EPSSHRINKSD(%)  2.2111E+01
 EPSSHRINKVR(%)  3.9333E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1367.0200822391396     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -631.86925567540140     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.72
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.87
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1367.020       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.63E-01  1.00E-02  1.97E+00  1.71E+00  9.94E-01  8.85E-01  8.91E+00  1.00E-02  6.81E-01  1.00E-02  3.04E+00
 


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
+        0.00E+00  5.41E+02
 
 TH 3
+       -2.17E+00  0.00E+00  1.85E+01
 
 TH 4
+       -8.50E+01  0.00E+00  3.33E+00  5.98E+02
 
 TH 5
+       -2.56E+01  0.00E+00 -7.64E+01 -1.33E+02  3.44E+02
 
 TH 6
+        7.94E+00  0.00E+00  2.94E+00 -1.67E+01 -5.03E+00  1.90E+02
 
 TH 7
+        3.01E-02  0.00E+00 -3.28E-03 -8.62E-03 -7.16E-02  3.86E-04 -8.24E-05
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.01E+01  0.00E+00  2.97E+00 -1.83E+01  9.98E+00 -5.89E+00  5.38E-02  0.00E+00  2.28E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.85E+01  0.00E+00  1.99E-01 -1.03E+01  2.87E+00  5.60E+00  5.23E-03  0.00E+00  1.71E+01  0.00E+00  4.59E+01
 
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
 #CPUT: Total CPU Time in Seconds,       28.656
Stop Time:
Sat Sep 25 08:59:57 CDT 2021
