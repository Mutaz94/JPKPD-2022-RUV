Sat Sep 18 15:45:49 CDT 2021
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
$DATA ../../../../data/spa/D/dat97.csv ignore=@
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
 (E4.0,E3.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m97.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   25813.1103822308        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.3609E+02  6.5467E+02 -9.8891E+00  6.4420E+02  1.5125E+02 -3.1695E+03 -1.1957E+03 -4.9828E+01 -1.7935E+03 -8.0416E+02
            -4.7536E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -389.502824160237        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0316E+00  9.0314E-01  8.0931E-01  1.2777E+00  1.8626E+00  2.1849E+00  9.0981E-01  9.4015E-01  6.7370E-01  8.6137E-01
             1.5065E+01
 PARAMETER:  1.3113E-01 -1.8816E-03 -1.1157E-01  3.4505E-01  7.2199E-01  8.8155E-01  5.4827E-03  3.8287E-02 -2.9497E-01 -4.9235E-02
             2.8124E+00
 GRADIENT:  -2.2278E+01 -2.6891E+00 -5.9607E+00 -2.1424E+01 -4.3694E+00  6.1746E+01  1.5012E+00  3.3037E+00  5.7075E+00  3.6440E-02
             9.2186E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -410.975095794846        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.1090E+00  6.7748E-01  6.8124E-01  1.4880E+00  3.9485E+00  1.6815E+00  6.1182E-01  1.5115E-01  2.0111E-01  1.2027E+00
             1.7010E+01
 PARAMETER:  2.0348E-01 -2.8938E-01 -2.8384E-01  4.9740E-01  1.4733E+00  6.1971E-01 -3.9132E-01 -1.7895E+00 -1.5039E+00  2.8459E-01
             2.9338E+00
 GRADIENT:  -1.9165E+01  1.9829E+01 -2.8318E+00  2.9307E+01 -1.0509E+01 -9.3624E+00  2.2831E-01  1.1766E-01  7.8477E-01  1.3343E+00
             2.5139E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -441.968191049498        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  1.0080E+00  2.6083E-01  2.6497E-01  1.2870E+00  2.2061E+01  1.8756E+00  8.5149E-02  2.0862E-02  1.5573E-01  6.0905E+00
             1.6127E+01
 PARAMETER:  1.0801E-01 -1.2439E+00 -1.2282E+00  3.5233E-01  3.1938E+00  7.2895E-01 -2.3634E+00 -3.7698E+00 -1.7596E+00  1.9067E+00
             2.8805E+00
 GRADIENT:   1.9651E+01  1.8782E+01  7.4098E+00 -8.9073E+00 -4.5693E+00 -6.0724E-01  5.3117E-02  7.4992E-03  9.3915E-01  5.3664E+00
             8.2740E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -463.861538363950        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  7.5728E-01  1.0067E-01  1.1574E-01  9.5331E-01  6.0958E+01  2.2287E+00  1.2724E-01  1.0000E-02  5.4464E-02  1.0656E+01
             1.4626E+01
 PARAMETER: -1.7802E-01 -2.1959E+00 -2.0564E+00  5.2190E-02  4.2102E+00  9.0142E-01 -1.9616E+00 -4.7685E+00 -2.8102E+00  2.4661E+00
             2.7828E+00
 GRADIENT:  -1.0161E+01 -4.0718E+00  1.9765E+01  1.9511E+01 -9.8534E-01  1.8212E+01  1.9117E-01  0.0000E+00  1.4035E-01  2.7028E+00
            -1.4215E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -479.200300340497        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      379
 NPARAMETR:  4.2681E-01  1.5124E-02  2.1525E-02  2.7462E-01  1.0936E+03  1.8390E+00  9.0958E-02  1.4073E-02  1.0000E-02  7.8470E+01
             1.3493E+01
 PARAMETER: -7.5141E-01 -4.0915E+00 -3.7385E+00 -1.1924E+00  7.0973E+00  7.0921E-01 -2.2974E+00 -4.1635E+00 -6.4786E+00  4.4627E+00
             2.7022E+00
 GRADIENT:   8.3592E+00  2.4738E+00  1.8572E+01 -1.9266E+01 -4.1421E-03 -4.0161E+01  3.0039E-03  4.4647E-03  0.0000E+00  3.1980E-03
            -2.6806E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -480.011050863396        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      453
 NPARAMETR:  3.5239E-01  1.0000E-02  1.3518E-02  1.9224E-01  2.4966E+03  1.8943E+00  6.1552E-02  1.6928E-02  1.0000E-02  1.4805E+02
             1.3238E+01
 PARAMETER: -9.4301E-01 -4.6318E+00 -4.2037E+00 -1.5490E+00  7.9227E+00  7.3884E-01 -2.6879E+00 -3.9788E+00 -7.5617E+00  5.0975E+00
             2.6831E+00
 GRADIENT:   8.4371E+00  0.0000E+00  1.4106E+00  2.0693E+00  1.4121E-04 -3.1537E+01  7.3016E-04  8.2619E-03  0.0000E+00  2.1455E-04
            -3.2671E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -481.055996162489        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:      619
 NPARAMETR:  3.5889E-01  1.0000E-02  1.4813E-02  2.0993E-01  2.0706E+03  1.9070E+00  3.4353E-02  1.8530E-02  1.0000E-02  1.2692E+02
             1.3965E+01
 PARAMETER: -9.2473E-01 -4.5294E+00 -4.1122E+00 -1.4610E+00  7.7356E+00  7.4551E-01 -3.2711E+00 -3.8883E+00 -7.5889E+00  4.9436E+00
             2.7366E+00
 GRADIENT:  -1.7448E+01  0.0000E+00 -9.8378E+00  1.7655E+01 -9.0891E-04 -3.0909E+01  8.2220E-05  1.1305E-02  0.0000E+00  4.4759E-04
             3.6685E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -483.023222873211        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      798
 NPARAMETR:  3.9471E-01  1.0862E-02  1.7451E-02  2.3818E-01  1.4912E+03  2.0946E+00  3.4504E-02  1.0000E-02  1.0000E-02  1.1043E+02
             1.3896E+01
 PARAMETER: -8.2960E-01 -4.4225E+00 -3.9484E+00 -1.3347E+00  7.4073E+00  8.3938E-01 -3.2667E+00 -4.5149E+00 -7.1944E+00  4.8044E+00
             2.7316E+00
 GRADIENT:   1.1998E-01  1.0646E+00 -1.2872E+00 -4.9664E-02 -3.6035E-04 -1.7146E-01  5.2682E-05  0.0000E+00  0.0000E+00  4.4659E-04
            -2.7716E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -483.130571854175        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      975
 NPARAMETR:  4.0702E-01  1.0000E-02  1.8771E-02  2.5225E-01  1.1546E+03  2.0931E+00  3.3925E-02  1.0000E-02  1.0000E-02  8.3590E+01
             1.3981E+01
 PARAMETER: -7.9889E-01 -4.5213E+00 -3.8754E+00 -1.2774E+00  7.1515E+00  8.3867E-01 -3.2836E+00 -4.5930E+00 -7.0990E+00  4.5259E+00
             2.7377E+00
 GRADIENT:   7.0151E-01  0.0000E+00  5.1867E-01 -1.1349E+00  5.0613E-04 -5.9570E-01  7.7449E-06  0.0000E+00  0.0000E+00  3.2648E-05
            -2.3894E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -483.132543251099        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1151
 NPARAMETR:  4.0835E-01  1.0000E-02  1.8993E-02  2.5483E-01  1.1170E+03  2.0977E+00  3.3973E-02  1.0000E-02  1.0000E-02  7.6650E+01
             1.3986E+01
 PARAMETER: -7.9563E-01 -4.5183E+00 -3.8637E+00 -1.2671E+00  7.1184E+00  8.4086E-01 -3.2822E+00 -4.6082E+00 -7.0727E+00  4.4393E+00
             2.7380E+00
 GRADIENT:   6.9458E-03  0.0000E+00 -2.6141E-02  4.2008E-02  5.5589E-04 -4.9529E-03  6.8921E-06  0.0000E+00  0.0000E+00  1.2818E-05
            -2.6008E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -483.144186400940        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1338
 NPARAMETR:  4.0381E-01  1.0000E-02  1.8449E-02  2.4894E-01  2.5473E+01  2.0987E+00  3.2437E-02  1.0000E-02  1.0000E-02  1.0000E-02
             1.3981E+01
 PARAMETER: -8.0681E-01 -4.5183E+00 -3.8927E+00 -1.2905E+00  3.3376E+00  8.4133E-01 -3.3285E+00 -4.6082E+00 -7.0727E+00 -5.6588E+01
             2.7377E+00
 GRADIENT:  -6.7957E-01  0.0000E+00  9.3914E-01 -1.4095E+00  8.6888E-03  2.3940E-01  1.3545E-05  0.0000E+00  0.0000E+00  0.0000E+00
             5.7966E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -483.159127372196        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1518
 NPARAMETR:  4.1330E-01  1.0000E-02  1.9333E-02  2.5873E-01  1.2852E+01  2.0958E+00  3.2382E-02  1.0000E-02  1.0000E-02  1.0000E-02
             1.3975E+01
 PARAMETER: -7.8357E-01 -4.5183E+00 -3.8459E+00 -1.2520E+00  2.6535E+00  8.3995E-01 -3.3301E+00 -4.6082E+00 -7.0727E+00 -5.4723E+01
             2.7373E+00
 GRADIENT:   1.0036E+00  0.0000E+00 -4.1898E-01 -3.9334E-02  2.6291E-03 -5.9126E-01  1.7595E-05  0.0000E+00  0.0000E+00  0.0000E+00
            -1.0138E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -483.161289696588        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1706
 NPARAMETR:  4.1473E-01  1.0000E-02  1.9567E-02  2.6133E-01  1.1319E+01  2.1000E+00  3.1851E-02  1.0000E-02  1.0000E-02  1.0000E-02
             1.3994E+01
 PARAMETER: -7.8014E-01 -4.5183E+00 -3.8339E+00 -1.2420E+00  2.5265E+00  8.4192E-01 -3.3467E+00 -4.6082E+00 -7.0727E+00 -5.3519E+01
             2.7386E+00
 GRADIENT:  -1.5702E-02  0.0000E+00  2.7061E-02 -4.2825E-02 -3.6812E-04  1.3561E-03  2.0156E-05  0.0000E+00  0.0000E+00  0.0000E+00
             8.6188E-03

0ITERATION NO.:   70    OBJECTIVE VALUE:  -483.161293070327        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1884
 NPARAMETR:  4.1478E-01  1.0000E-02  1.9570E-02  2.6138E-01  1.1333E+01  2.1000E+00  2.2243E-02  1.0000E-02  1.0000E-02  1.0000E-02
             1.3994E+01
 PARAMETER: -7.8000E-01 -4.5183E+00 -3.8337E+00 -1.2418E+00  2.5277E+00  8.4192E-01 -3.7057E+00 -4.6082E+00 -7.0727E+00 -5.3519E+01
             2.7386E+00
 GRADIENT:   1.3865E-02  0.0000E+00 -4.4558E-02  3.7788E-02  2.7287E-04 -4.4482E-03  9.7582E-06  0.0000E+00  0.0000E+00  0.0000E+00
            -6.9319E-04

0ITERATION NO.:   75    OBJECTIVE VALUE:  -483.161293514825        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     2072
 NPARAMETR:  4.1481E-01  1.0000E-02  1.9574E-02  2.6141E-01  1.1280E+01  2.1000E+00  1.7033E-02  1.0000E-02  1.0000E-02  1.0000E-02
             1.3994E+01
 PARAMETER: -7.7994E-01 -4.5183E+00 -3.8336E+00 -1.2416E+00  2.5230E+00  8.4194E-01 -3.9726E+00 -4.6082E+00 -7.0727E+00 -5.3519E+01
             2.7386E+00
 GRADIENT:   2.1209E-04  0.0000E+00  8.7207E-03 -2.4417E-02 -1.8016E-04  4.2598E-03  5.7906E-06  0.0000E+00  0.0000E+00  0.0000E+00
             8.3510E-03

0ITERATION NO.:   79    OBJECTIVE VALUE:  -483.161293674503        NO. OF FUNC. EVALS.: 128
 CUMULATIVE NO. OF FUNC. EVALS.:     2200
 NPARAMETR:  4.1481E-01  1.0000E-02  1.9574E-02  2.6141E-01  1.1288E+01  2.1000E+00  1.6925E-02  1.0000E-02  1.0000E-02  1.0000E-02
             1.3994E+01
 PARAMETER: -7.7994E-01 -4.5183E+00 -3.8336E+00 -1.2416E+00  2.5237E+00  8.4194E-01 -3.9790E+00 -4.6082E+00 -7.0727E+00 -5.3519E+01
             2.7386E+00
 GRADIENT:   7.7091E-04  0.0000E+00  5.3805E-03 -1.9775E-02 -1.1293E-04  4.9250E-03  5.7101E-06  0.0000E+00  0.0000E+00  0.0000E+00
             5.1451E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2200
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.0729E-03  6.3666E-06  8.5439E-05 -1.9306E-04 -5.9372E-07
 SE:             2.9114E-02  2.9908E-06  2.3490E-04  2.6952E-04  1.3007E-06
 N:                     100         100         100         100         100

 P VAL.:         9.4324E-01  3.3274E-02  7.1607E-01  4.7381E-01  6.4806E-01

 ETASHRINKSD(%)  2.4650E+00  9.9990E+01  9.9213E+01  9.9097E+01  9.9996E+01
 ETASHRINKVR(%)  4.8693E+00  1.0000E+02  9.9994E+01  9.9992E+01  1.0000E+02
 EBVSHRINKSD(%)  2.6398E+00  9.9981E+01  9.9156E+01  9.8978E+01  9.9995E+01
 EBVSHRINKVR(%)  5.2100E+00  1.0000E+02  9.9993E+01  9.9990E+01  1.0000E+02
 RELATIVEINF(%)  9.5692E+00  6.5313E-07  4.9280E-05  6.8428E-05  2.7259E-08
 EPSSHRINKSD(%)  5.3642E+00
 EPSSHRINKVR(%)  1.0441E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -483.16129367450321     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       251.98953288923497     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    28.07
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.41
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -483.161       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.15E-01  1.00E-02  1.96E-02  2.61E-01  1.13E+01  2.10E+00  1.69E-02  1.00E-02  1.00E-02  1.00E-02  1.40E+01
 


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
+        1.37E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -9.45E+03  0.00E+00  1.31E+06
 
 TH 4
+       -1.61E+02  0.00E+00 -1.16E+05  1.12E+04
 
 TH 5
+        2.24E-01  0.00E+00 -1.03E+01  9.51E-01  3.62E-04
 
 TH 6
+       -9.12E-01  0.00E+00  9.24E+02 -9.16E+01 -2.76E-03  3.67E+01
 
 TH 7
+        3.33E-02  0.00E+00  6.74E-01  2.91E-01  1.17E-03  1.02E-02 -2.36E-01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.36E+01  0.00E+00  2.53E+02 -1.63E+01 -4.85E-03  5.57E-01  2.21E-03  0.00E+00  0.00E+00  0.00E+00  1.70E+00
 
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
 #CPUT: Total CPU Time in Seconds,       34.550
Stop Time:
Sat Sep 18 15:46:25 CDT 2021
