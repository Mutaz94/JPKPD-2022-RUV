Sat Sep 25 13:38:58 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat60.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1664.09715572865        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.7186E+01 -8.0362E+01 -3.8417E+01 -5.4156E+01  3.6838E+01 -2.7005E+00 -1.7877E+01  7.1272E+00  4.7960E+00  1.5994E+01
            -2.2264E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1670.20802718225        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:      111
 NPARAMETR:  1.0124E+00  1.2109E+00  1.1527E+00  9.1627E-01  1.0895E+00  1.0048E+00  1.2316E+00  9.5526E-01  9.2731E-01  8.0793E-01
             1.1172E+00
 PARAMETER:  1.1232E-01  2.9136E-01  2.4207E-01  1.2555E-02  1.8575E-01  1.0479E-01  3.0828E-01  5.4226E-02  2.4534E-02 -1.1328E-01
             2.1086E-01
 GRADIENT:   1.8274E+01  1.5685E+01  3.7711E+01 -3.1186E+01 -2.0495E+01 -9.4381E+00  7.2919E+00 -1.0060E+01 -6.6998E+00 -1.1442E+01
             7.9599E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1671.12264315593        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:      285
 NPARAMETR:  1.0241E+00  1.2132E+00  1.1341E+00  9.2301E-01  1.0856E+00  1.0282E+00  1.2373E+00  9.6161E-01  9.6436E-01  8.9406E-01
             1.1276E+00
 PARAMETER:  1.2382E-01  2.9323E-01  2.2582E-01  1.9886E-02  1.8216E-01  1.2783E-01  3.1293E-01  6.0856E-02  6.3705E-02 -1.1985E-02
             2.2010E-01
 GRADIENT:   4.0887E+01  1.8713E+01  2.8016E+01 -1.7115E+01 -2.6789E+01 -9.0242E-01  1.0302E+01 -5.8629E+00  8.3416E-01  7.9733E-01
             1.8245E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1671.43086146537        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      465            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0233E+00  1.2121E+00  1.1300E+00  9.2355E-01  1.0864E+00  1.0306E+00  1.2339E+00  9.7265E-01  9.5880E-01  8.8779E-01
             1.1255E+00
 PARAMETER:  1.2302E-01  2.9235E-01  2.2219E-01  2.0470E-02  1.8290E-01  1.3013E-01  3.1021E-01  7.2266E-02  5.7926E-02 -1.9015E-02
             2.1820E-01
 GRADIENT:   8.3019E+01  2.9220E+01  2.6200E+01 -1.0179E+01 -2.1424E+01  9.3443E+00  1.1269E+01 -5.6485E+00  6.9382E-01  2.0579E-01
             1.7738E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1671.71869445153        NO. OF FUNC. EVALS.:  88
 CUMULATIVE NO. OF FUNC. EVALS.:      553
 NPARAMETR:  1.0215E+00  1.2106E+00  1.1268E+00  9.2385E-01  1.0866E+00  1.0294E+00  1.2301E+00  9.8236E-01  9.5831E-01  8.8729E-01
             1.1229E+00
 PARAMETER:  1.2130E-01  2.9108E-01  2.1936E-01  2.0791E-02  1.8307E-01  1.2900E-01  3.0712E-01  8.2203E-02  5.7416E-02 -1.9581E-02
             2.1590E-01
 GRADIENT:   7.9088E+01  2.7513E+01  2.4446E+01 -9.5892E+00 -1.8789E+01  8.9975E+00  1.0701E+01 -5.3224E+00  7.4394E-01  3.0074E-01
             1.7197E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1671.73066932379        NO. OF FUNC. EVALS.: 199
 CUMULATIVE NO. OF FUNC. EVALS.:      752             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0214E+00  1.2102E+00  1.1265E+00  9.2376E-01  1.0868E+00  1.0327E+00  1.2301E+00  9.8236E-01  9.5151E-01  8.8734E-01
             1.1231E+00
 PARAMETER:  1.2116E-01  2.9077E-01  2.1911E-01  2.0699E-02  1.8324E-01  1.3215E-01  3.0710E-01  8.2202E-02  5.0299E-02 -1.9524E-02
             2.1609E-01
 GRADIENT:   7.8503E+01  2.7231E+01  2.4429E+01 -1.0286E+01 -1.8388E+01  1.0469E+01  1.0488E+01 -5.3813E+00 -1.8779E-01  2.4224E-01
             1.7161E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1671.75563332232        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:      915            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0213E+00  1.2101E+00  1.1262E+00  9.2381E-01  1.0869E+00  1.0302E+00  1.2298E+00  9.8331E-01  9.5652E-01  8.8734E-01
             1.1229E+00
 PARAMETER:  1.2110E-01  2.9071E-01  2.1882E-01  2.0749E-02  1.8329E-01  1.2974E-01  3.0687E-01  8.3166E-02  5.5549E-02 -1.9532E-02
             2.1592E-01
 GRADIENT:   7.8527E+01  2.6974E+01  2.4093E+01 -9.8034E+00 -1.8058E+01  9.3560E+00  1.0600E+01 -5.2796E+00  5.2455E-01  3.0197E-01
             1.7226E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1671.87083138292        NO. OF FUNC. EVALS.: 148
 CUMULATIVE NO. OF FUNC. EVALS.:     1063             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0191E+00  1.2043E+00  1.1219E+00  9.2237E-01  1.0901E+00  1.0379E+00  1.2354E+00  9.8293E-01  9.5801E-01  8.8871E-01
             1.1265E+00
 PARAMETER:  1.1895E-01  2.8591E-01  2.1499E-01  1.9188E-02  1.8624E-01  1.3724E-01  3.1141E-01  8.2779E-02  5.7100E-02 -1.7990E-02
             2.1909E-01
 GRADIENT:   7.2577E+01  1.9706E+01  2.0667E+01 -1.3480E+01 -1.0767E+01  1.2946E+01  1.1249E+01 -4.7556E+00  1.0754E+00  4.7515E-01
             1.8950E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1671.88290367496        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:     1197
 NPARAMETR:  1.0191E+00  1.2043E+00  1.1219E+00  9.2237E-01  1.0901E+00  1.0298E+00  1.2354E+00  9.8293E-01  9.5801E-01  8.8871E-01
             1.1265E+00
 PARAMETER:  1.1895E-01  2.8591E-01  2.1499E-01  1.9188E-02  1.8624E-01  1.2933E-01  3.1141E-01  8.2779E-02  5.7100E-02 -1.7990E-02
             2.1909E-01
 GRADIENT:   3.0549E+01  8.2226E+00  2.0430E+01 -1.9142E+01 -1.1653E+01 -4.0268E-05  9.5941E+00 -4.7602E+00  5.8100E-01  4.3245E-01
             1.8733E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1673.03668524537        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     1353             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0161E+00  1.2021E+00  1.0618E+00  9.2364E-01  1.0907E+00  1.0295E+00  1.1621E+00  1.0125E+00  9.5134E-01  8.8205E-01
             1.1199E+00
 PARAMETER:  1.1594E-01  2.8407E-01  1.6000E-01  2.0569E-02  1.8683E-01  1.2912E-01  2.5019E-01  1.1244E-01  5.0116E-02 -2.5501E-02
             2.1328E-01
 GRADIENT:   6.5025E+01  4.8898E+00  1.2254E-01  2.9215E+00  2.1923E+01  9.2765E+00  1.0534E+00 -1.5201E+00 -1.1688E-01  1.4961E-01
             1.8398E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1673.03834106493        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1531
 NPARAMETR:  1.0161E+00  1.2021E+00  1.0623E+00  9.2364E-01  1.0907E+00  1.0300E+00  1.1625E+00  1.0125E+00  9.5581E-01  8.8087E-01
             1.1199E+00
 PARAMETER:  1.1594E-01  2.8407E-01  1.6044E-01  2.0569E-02  1.8683E-01  1.2958E-01  2.5058E-01  1.1244E-01  5.4808E-02 -2.6849E-02
             2.1328E-01
 GRADIENT:   2.3240E+01 -6.8834E+00  4.1755E-03 -2.5607E+00  2.0992E+01  2.0945E-03  9.4552E-03 -1.5276E+00 -2.9765E-03 -3.3706E-03
             1.8264E+01

0ITERATION NO.:   54    OBJECTIVE VALUE:  -1673.21636837097        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1658
 NPARAMETR:  1.0138E+00  1.2023E+00  1.0622E+00  9.2363E-01  1.0889E+00  1.0288E+00  1.1620E+00  1.0183E+00  9.5565E-01  8.8073E-01
             1.1142E+00
 PARAMETER:  1.1371E-01  2.8422E-01  1.6035E-01  2.0563E-02  1.8520E-01  1.2837E-01  2.5016E-01  1.1812E-01  5.4631E-02 -2.6904E-02
             2.0813E-01
 GRADIENT:   1.0346E+06  8.2784E+05  7.3370E+05  2.3529E+06 -1.2704E+06 -4.0580E-01 -4.7028E+05  1.9918E+06 -1.1764E+06  1.6113E-01
            -1.1308E+06
 NUMSIGDIG:         3.3         3.3         3.3         3.3         3.3         1.8         3.3         3.3         3.3         1.5
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1658
 NO. OF SIG. DIGITS IN FINAL EST.:  1.5
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -8.4750E-03 -2.8101E-03 -2.8177E-02  1.3340E-03 -3.8987E-02
 SE:             2.9851E-02  2.2369E-02  1.2605E-02  2.1800E-02  2.0370E-02
 N:                     100         100         100         100         100

 P VAL.:         7.7648E-01  9.0003E-01  2.5394E-02  9.5121E-01  5.5629E-02

 ETASHRINKSD(%)  1.0000E-10  2.5062E+01  5.7771E+01  2.6967E+01  3.1758E+01
 ETASHRINKVR(%)  1.0000E-10  4.3843E+01  8.2167E+01  4.6661E+01  5.3430E+01
 EBVSHRINKSD(%)  4.9465E-01  2.5251E+01  6.2293E+01  2.7962E+01  2.9778E+01
 EBVSHRINKVR(%)  9.8685E-01  4.4126E+01  8.5781E+01  4.8105E+01  5.0689E+01
 RELATIVEINF(%)  9.8557E+01  2.5785E+00  1.6105E+00  2.4657E+00  9.5385E+00
 EPSSHRINKSD(%)  4.5038E+01
 EPSSHRINKVR(%)  6.9792E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1673.2163683709744     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -938.06554180723617     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.59
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.06
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1673.216       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.20E+00  1.06E+00  9.24E-01  1.09E+00  1.03E+00  1.16E+00  1.02E+00  9.56E-01  8.81E-01  1.11E+00
 


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
+        4.43E+09
 
 TH 2
+        2.01E+03  5.04E+08
 
 TH 3
+        3.00E+09  4.99E+04  2.03E+09
 
 TH 4
+        7.44E+03 -3.27E+04  1.84E+05  6.89E+09
 
 TH 5
+       -3.41E+03  1.29E+04 -8.46E+04  4.86E+04  1.45E+09
 
 TH 6
+        1.08E+04  3.66E+03  7.34E+03  1.35E+04 -6.20E+03  1.85E+02
 
 TH 7
+        6.93E+03  2.37E+03  4.69E+03  8.64E+03 -3.97E+03 -4.30E+03  6.96E+08
 
 TH 8
+        5.72E+03 -1.69E+05  1.41E+05 -6.25E+05 -2.42E+09  1.04E+04  6.65E+03  4.07E+09
 
 TH 9
+       -9.11E+04 -3.08E+04 -6.17E+04 -1.14E+05  5.21E+04 -4.66E+09  2.12E+09 -8.73E+04  6.44E+09
 
 TH10
+       -2.36E+04 -7.97E+03 -1.60E+04 -2.95E+04  1.34E+04 -5.06E+09  2.30E+09 -2.26E+04  6.99E+09  7.58E+09
 
 TH11
+       -2.97E+03 -4.23E+05 -7.34E+04 -1.56E+06  1.26E+09 -5.39E+03 -3.44E+03 -2.11E+09  4.53E+04  1.18E+04  1.09E+09
 
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
 #CPUT: Total CPU Time in Seconds,       28.722
Stop Time:
Sat Sep 25 13:39:29 CDT 2021
