Thu Sep 30 03:20:45 CDT 2021
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
$DATA ../../../../data/spa1/D/dat60.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      600
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E22.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m60.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   16145.4355582307        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.5667E+02  3.3275E+02 -2.9974E+01 -3.0016E+01  3.3534E+02 -1.9996E+03 -8.1419E+02 -3.0217E+02 -1.7691E+03 -7.1836E+02
            -3.0586E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -650.163064697154        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.8146E+00  1.0183E+00  9.6660E-01  2.4743E+00  1.1656E+00  3.4360E+00  1.3385E+00  9.6695E-01  2.9251E+00  1.2370E+00
             1.1870E+01
 PARAMETER:  6.9589E-01  1.1815E-01  6.6034E-02  1.0060E+00  2.5327E-01  1.3343E+00  3.9158E-01  6.6393E-02  1.1733E+00  3.1266E-01
             2.5740E+00
 GRADIENT:   7.7110E+01  3.9606E+01 -4.6373E+01  9.6551E+01 -5.4369E+00  1.0018E+02  1.8742E+00  7.9021E+00 -7.6780E+00  6.5403E+00
             1.2975E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -697.494773580661        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.6452E+00  1.3713E+00  5.5217E+00  2.5779E+00  9.7216E+00  3.1927E+00  3.7567E+00  6.1050E-01  4.0818E+00  1.2663E+01
             9.7403E+00
 PARAMETER:  5.9786E-01  4.1575E-01  1.8087E+00  1.0470E+00  2.3744E+00  1.2609E+00  1.4236E+00 -3.9347E-01  1.5065E+00  2.6387E+00
             2.3763E+00
 GRADIENT:   8.1862E+01  2.3928E+01  3.3979E+00  7.4817E+01 -2.7414E+00  7.6037E+01  1.5367E+01 -6.5703E-02  5.0197E+01  1.6124E+01
             9.0253E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -744.737451502313        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.1198E+00  7.6326E-01  2.2201E+00  1.3887E+00  3.9423E+00  2.4084E+00  1.1437E+00  3.5224E+00  3.1157E+00  5.2620E+00
             9.5810E+00
 PARAMETER:  2.1313E-01 -1.7016E-01  8.9755E-01  4.2837E-01  1.4718E+00  9.7898E-01  2.3430E-01  1.3591E+00  1.2365E+00  1.7605E+00
             2.3598E+00
 GRADIENT:  -2.9055E+01  1.8723E+00  5.1149E+00 -2.3531E+01 -1.8665E+01  5.5663E+00  3.8319E+00  1.7334E+01  3.3851E+01  1.6173E+00
             1.4332E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -827.802494547524        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  8.6931E-01  9.1371E-02  1.3871E-01  9.6821E-01  1.3662E+01  1.6294E+00  2.3124E+00  2.2453E+00  9.6618E-01  1.7590E+00
             6.9350E+00
 PARAMETER: -4.0057E-02 -2.2928E+00 -1.8754E+00  6.7692E-02  2.7146E+00  5.8819E-01  9.3828E-01  9.0884E-01  6.5597E-02  6.6477E-01
             2.0366E+00
 GRADIENT:   8.2648E+01  7.6947E+01 -4.8134E+01  1.3053E+02 -5.3063E+00 -5.0192E+01  3.8890E+00 -1.0459E+00 -4.5590E+00  1.8467E-02
            -2.6008E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -876.856924386772        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      442
 NPARAMETR:  6.6057E-01  3.5423E-02  7.1824E-02  5.6054E-01  6.6404E+01  1.6659E+00  3.4585E-01  1.5234E+00  5.5712E-01  1.3080E+00
             6.8466E+00
 PARAMETER: -3.1465E-01 -3.2404E+00 -2.5335E+00 -4.7885E-01  4.2958E+00  6.1035E-01 -9.6175E-01  5.2091E-01 -4.8497E-01  3.6851E-01
             2.0238E+00
 GRADIENT:   4.7017E+01  7.7829E+00  3.0194E+01 -3.7692E+01 -1.4089E-02 -3.4807E+01  9.9340E-02 -3.6096E+01 -4.8334E+00  9.7569E-05
            -2.4231E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -907.356146916338        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      617
 NPARAMETR:  5.6897E-01  2.1128E-02  4.8231E-02  4.5531E-01  2.9109E+02  1.7589E+00  6.5874E-02  1.6650E+00  3.6755E-01  1.3393E+00
             8.2239E+00
 PARAMETER: -4.6392E-01 -3.7572E+00 -2.9317E+00 -6.8678E-01  5.7736E+00  6.6469E-01 -2.6200E+00  6.0982E-01 -9.0091E-01  3.9214E-01
             2.2070E+00
 GRADIENT:   1.0114E-01 -7.4839E-01  1.0373E+00 -1.6354E+00  1.0389E-02 -6.1099E-01  4.8207E-04  7.3212E-01  1.1550E+00 -1.9074E-06
             6.7049E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -907.667004896903        NO. OF FUNC. EVALS.: 123
 CUMULATIVE NO. OF FUNC. EVALS.:      740
 NPARAMETR:  5.6182E-01  2.3178E-02  4.6075E-02  4.4156E-01  2.2005E+02  1.7613E+00  7.6913E-02  1.6983E+00  2.2545E-01  1.3264E+00
             8.2262E+00
 PARAMETER: -4.7657E-01 -3.6646E+00 -2.9775E+00 -7.1744E-01  5.4939E+00  6.6607E-01 -2.4651E+00  6.2962E-01 -1.3896E+00  3.8245E-01
             2.2073E+00
 GRADIENT:   5.1886E+01  2.0638E+00  6.9395E+01  2.4795E+01  7.4278E-03  1.9967E+01  3.7282E-03  4.9604E+00  4.4777E-01  7.9490E-06
             2.1596E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -907.746263851605        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      812
 NPARAMETR:  5.5813E-01  2.1968E-02  4.5361E-02  4.3626E-01  1.5892E+01  1.7568E+00  1.0000E-02  1.7416E+00  4.0270E-02  1.2820E+00
             8.2130E+00
 PARAMETER: -4.8316E-01 -3.7182E+00 -2.9931E+00 -7.2951E-01  2.8658E+00  6.6347E-01 -2.1789E+01  6.5483E-01 -3.1121E+00  3.4842E-01
             2.2057E+00
 GRADIENT:   5.1201E+01  3.7003E+00  7.1877E+01  2.0530E+01 -9.5735E-02  1.9193E+01  0.0000E+00  8.4941E+00  2.4516E-02  9.8532E-03
             2.0837E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -907.754171960381        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      886
 NPARAMETR:  5.5640E-01  2.0844E-02  4.4801E-02  4.3315E-01  1.0996E+01  1.7573E+00  1.0000E-02  1.7344E+00  1.3350E-02  9.2327E-01
             8.2086E+00
 PARAMETER: -4.8626E-01 -3.7707E+00 -3.0055E+00 -7.3667E-01  2.4976E+00  6.6376E-01 -3.2028E+01  6.5065E-01 -4.2163E+00  2.0163E-02
             2.2052E+00
 GRADIENT:   5.1540E+01  3.9339E+00  7.2534E+01  2.0543E+01 -6.5841E-01  1.9185E+01  0.0000E+00  8.6293E+00  3.5815E-03  3.1341E-02
             2.0369E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -907.754735882440        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:      984
 NPARAMETR:  5.5531E-01  2.0552E-02  4.4509E-02  4.3140E-01  1.0643E+01  1.7575E+00  1.0000E-02  1.7303E+00  1.0000E-02  8.7886E-01
             8.2074E+00
 PARAMETER: -4.8823E-01 -3.7848E+00 -3.0121E+00 -7.4072E-01  2.4649E+00  6.6390E-01 -3.7880E+01  6.4829E-01 -4.7727E+00 -2.9127E-02
             2.2050E+00
 GRADIENT:  -9.9349E-01  2.6468E+00 -1.2037E+00 -6.7883E+00 -9.2330E-01 -5.7403E-01  0.0000E+00  4.0068E+00  0.0000E+00  3.4802E-02
            -1.7930E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -907.891180768366        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1164
 NPARAMETR:  5.6111E-01  2.0488E-02  4.6061E-02  4.4194E-01  1.2997E+01  1.7617E+00  2.4421E-02  1.7127E+00  1.0751E-01  6.3544E-01
             8.2202E+00
 PARAMETER: -4.7784E-01 -3.7879E+00 -2.9778E+00 -7.1659E-01  2.6647E+00  6.6627E-01 -3.6123E+00  6.3809E-01 -2.1302E+00 -3.5344E-01
             2.2066E+00
 GRADIENT:  -1.5162E+00  5.1477E-01 -1.4701E+00 -8.9588E-01  5.3201E-02 -3.8275E-01  2.7472E-04  1.3018E-01  5.0788E-03  3.3035E-03
            -6.3021E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -907.904594688807        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1346
 NPARAMETR:  5.6221E-01  1.9965E-02  4.6118E-02  4.4251E-01  1.2366E+01  1.7633E+00  2.2218E-02  1.7101E+00  1.0519E-01  3.2176E-01
             8.2238E+00
 PARAMETER: -4.7588E-01 -3.8138E+00 -2.9766E+00 -7.1528E-01  2.6149E+00  6.6719E-01 -3.7069E+00  6.3655E-01 -2.1520E+00 -1.0339E+00
             2.2070E+00
 GRADIENT:  -5.6817E-01  1.2739E-01 -1.7000E+00 -4.1318E-01  9.2316E-02 -1.9803E-01  1.8405E-04 -2.0875E-01  3.2269E-03  9.3921E-04
            -3.3513E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -907.912760494210        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:     1488
 NPARAMETR:  5.6327E-01  1.9441E-02  4.6037E-02  4.4258E-01  1.1253E+01  1.7666E+00  1.0000E-02  1.7143E+00  8.4310E-02  1.0000E-02
             8.2247E+00
 PARAMETER: -4.7399E-01 -3.8404E+00 -2.9783E+00 -7.1514E-01  2.5207E+00  6.6906E-01 -5.0252E+00  6.3902E-01 -2.3733E+00 -1.1551E+01
             2.2071E+00
 GRADIENT:   5.1425E+01  4.1239E-01  6.9146E+01  2.6921E+01  5.4616E-02  2.0113E+01  0.0000E+00  4.5871E+00  7.4340E-02  0.0000E+00
             2.1648E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -907.913251664596        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     1658
 NPARAMETR:  5.6293E-01  1.9482E-02  4.6033E-02  4.4235E-01  1.1199E+01  1.7659E+00  1.0000E-02  1.7135E+00  7.2619E-02  1.0000E-02
             8.2258E+00
 PARAMETER: -4.7460E-01 -3.8382E+00 -2.9784E+00 -7.1565E-01  2.5158E+00  6.6866E-01 -4.8927E+00  6.3855E-01 -2.5225E+00 -9.7988E+00
             2.2073E+00
 GRADIENT:   3.7107E-01  5.5724E-02 -2.3354E+00 -6.7722E-02 -1.2399E-02  1.9690E-01  0.0000E+00 -4.5274E-02  6.9626E-04  0.0000E+00
            -2.7514E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -907.913451083133        NO. OF FUNC. EVALS.: 172
 CUMULATIVE NO. OF FUNC. EVALS.:     1830
 NPARAMETR:  5.6285E-01  1.9459E-02  4.5989E-02  4.4219E-01  1.1238E+01  1.7660E+00  1.0000E-02  1.7131E+00  7.1151E-02  1.0000E-02
             8.2257E+00
 PARAMETER: -4.7474E-01 -3.8395E+00 -2.9794E+00 -7.1602E-01  2.5193E+00  6.6870E-01 -4.8927E+00  6.3832E-01 -2.5429E+00 -9.7988E+00
             2.2073E+00
 GRADIENT:   4.4077E-01  3.4364E-02 -2.5782E+00  2.7947E-01  7.7354E-04  2.0262E-01  0.0000E+00 -7.2498E-02  3.7143E-04  0.0000E+00
            -3.1763E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -907.913706426207        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     2021
 NPARAMETR:  5.6286E-01  1.9453E-02  4.5973E-02  4.4202E-01  1.1241E+01  1.7663E+00  1.0000E-02  1.7130E+00  7.4086E-02  1.0000E-02
             8.2261E+00
 PARAMETER: -4.7473E-01 -3.8398E+00 -2.9797E+00 -7.1639E-01  2.5196E+00  6.6888E-01 -4.8927E+00  6.3825E-01 -2.5025E+00 -9.7988E+00
             2.2073E+00
 GRADIENT:   5.5408E-01  3.2090E-02 -2.4571E+00  8.4503E-05  6.3475E-04  2.7854E-01  0.0000E+00 -2.7297E-04  9.0833E-04  0.0000E+00
            -2.3952E-01

0ITERATION NO.:   84    OBJECTIVE VALUE:  -907.913949953483        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:     2136
 NPARAMETR:  5.6260E-01  1.9444E-02  4.5940E-02  4.4182E-01  1.1235E+01  1.7658E+00  1.0000E-02  1.7132E+00  7.6344E-02  1.0000E-02
             8.2258E+00
 PARAMETER: -4.7519E-01 -3.8402E+00 -2.9804E+00 -7.1685E-01  2.5191E+00  6.6862E-01 -4.8927E+00  6.3838E-01 -2.4725E+00 -9.7988E+00
             2.2073E+00
 GRADIENT:   4.1269E-01  3.4196E-02 -2.4745E+00  3.1243E-02 -2.1480E-03  2.0691E-01  0.0000E+00  1.0293E-01  1.5104E-03  0.0000E+00
            -2.2179E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2136
 NO. OF SIG. DIGITS IN FINAL EST.:  2.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.7228E-03 -2.8953E-05  7.6531E-03 -3.0902E-03  4.4377E-06
 SE:             2.9079E-02  6.1725E-06  2.5489E-02  1.9108E-03  2.4260E-06
 N:                     100         100         100         100         100

 P VAL.:         7.9056E-01  2.7257E-06  7.6398E-01  1.0582E-01  6.7366E-02

 ETASHRINKSD(%)  2.5828E+00  9.9979E+01  1.4610E+01  9.3599E+01  9.9992E+01
 ETASHRINKVR(%)  5.0989E+00  1.0000E+02  2.7085E+01  9.9590E+01  1.0000E+02
 EBVSHRINKSD(%)  2.7226E+00  9.9961E+01  1.3042E+01  9.3867E+01  9.9990E+01
 EBVSHRINKVR(%)  5.3710E+00  1.0000E+02  2.4383E+01  9.9624E+01  1.0000E+02
 RELATIVEINF(%)  1.7889E+01  1.4645E-06  2.9674E+00  1.1564E-02  8.1989E-08
 EPSSHRINKSD(%)  1.1993E+01
 EPSSHRINKVR(%)  2.2548E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -907.91394995348253     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       11.024583251190165     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    39.28
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     9.67
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -907.914       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         5.63E-01  1.94E-02  4.59E-02  4.42E-01  1.12E+01  1.77E+00  1.00E-02  1.71E+00  7.63E-02  1.00E-02  8.23E+00
 


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
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        5.11E+00
 
 TH 2
+        5.42E+00  5.76E+00
 
 TH 3
+       -8.12E+02 -8.61E+02  1.29E+05
 
 TH 4
+        1.21E+02  1.29E+02 -1.93E+04  2.88E+03
 
 TH 5
+        1.62E-02  1.72E-02 -2.58E+00  3.85E-01  5.14E-05
 
 TH 6
+       -4.18E-01 -4.44E-01  6.65E+01 -9.93E+00 -1.33E-03  3.43E-02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        2.36E-01  2.51E-01 -3.75E+01  5.61E+00  7.50E-04 -1.93E-02  0.00E+00  1.09E-02
 
 TH 9
+       -2.41E-01 -2.56E-01  3.83E+01 -5.73E+00 -7.66E-04  1.98E-02  0.00E+00 -1.12E-02  1.14E-02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -7.49E-01 -7.95E-01  1.19E+02 -1.78E+01 -2.38E-03  6.13E-02  0.00E+00 -3.46E-02  3.54E-02  0.00E+00  1.10E-01
 
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
+        1.05E+03
 
 TH 2
+       -3.66E+02  1.97E+04
 
 TH 3
+       -9.14E+02 -1.04E+03  1.32E+05
 
 TH 4
+       -5.17E+02 -5.44E+02 -1.97E+04  3.59E+03
 
 TH 5
+        2.23E-01 -7.53E+00 -2.57E+00  7.37E-01  1.09E-02
 
 TH 6
+        4.63E+00  3.72E+01  6.59E+01 -2.61E+01  1.71E-03  5.33E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        7.87E+00 -5.44E+01 -4.94E+01 -6.63E+01 -1.81E-02  3.17E+00  0.00E+00  3.84E+01
 
 TH 9
+        9.23E-01 -1.50E+01  3.84E+01 -1.23E+01 -1.87E-03  5.79E-01  0.00E+00  4.32E+00  9.30E-01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.30E+01 -1.96E+01  1.22E+02 -1.60E+01  2.72E-03  1.08E+00  0.00E+00  1.76E+00  2.84E-01  0.00E+00  7.54E+00
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.07E+03
 
 TH 2
+       -3.11E+02  5.10E+02
 
 TH 3
+       -1.97E+03  4.19E+02  1.36E+05
 
 TH 4
+       -5.10E+02  2.39E+02 -2.02E+04  3.79E+03
 
 TH 5
+        1.22E-01 -1.43E-01 -5.14E+00  7.98E-01  3.33E-04
 
 TH 6
+        6.31E+01 -8.59E+00 -2.23E+02 -4.42E+01  1.10E-02  5.88E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        5.92E+01 -8.66E+01  1.16E+03 -2.99E+02 -7.88E-02  7.92E+00  0.00E+00  6.88E+01
 
 TH 9
+        9.96E+00 -1.47E+01  2.30E+02 -5.29E+01 -1.51E-02  1.31E+00  0.00E+00  1.14E+01  2.20E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -7.40E+01  2.16E+01  4.49E+02 -6.62E+01 -3.47E-02  1.87E+01  0.00E+00 -1.27E+00 -8.05E-01  0.00E+00  1.54E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.02
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       49.029
Stop Time:
Thu Sep 30 03:21:35 CDT 2021
