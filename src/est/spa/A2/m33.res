Sat Sep 18 09:47:52 CDT 2021
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
$DATA ../../../../data/spa/A2/dat33.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m33.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -343.172744859204        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.7164E+02  3.6144E+01  4.7586E+01  2.5800E+01  6.9533E+01  1.8799E+01 -2.7571E+01 -2.3799E+01 -1.8434E+01 -6.1307E+01
            -2.4314E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1224.17957574226        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.5755E-01  1.0396E+00  9.7854E-01  1.0566E+00  9.8351E-01  8.2904E-01  9.9067E-01  1.0015E+00  9.1736E-01  9.8031E-01
             5.3416E+00
 PARAMETER:  5.6628E-02  1.3887E-01  7.8308E-02  1.5505E-01  8.3368E-02 -8.7486E-02  9.0623E-02  1.0152E-01  1.3750E-02  8.0110E-02
             1.7755E+00
 GRADIENT:  -1.2275E+02  1.0402E+01 -1.1235E+01  2.3068E+01 -1.7681E+01 -2.0513E+01  1.2563E+01  4.8351E+00  2.0019E+01  1.9362E+01
             2.7933E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1298.74372664063        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.5474E-01  8.0008E-01  1.9238E+00  1.1689E+00  1.2180E+00  8.6333E-01  4.6419E-01  5.3453E-01  9.4902E-01  4.9435E-01
             3.9770E+00
 PARAMETER:  5.3685E-02 -1.2305E-01  7.5428E-01  2.5604E-01  2.9721E-01 -4.6963E-02 -6.6747E-01 -5.2637E-01  4.7676E-02 -6.0450E-01
             1.4805E+00
 GRADIENT:  -3.5022E+01 -1.5525E+01 -8.7093E-01 -9.7629E+00  3.1878E+00 -9.5778E+00  9.4835E-01  4.4635E-01  2.5830E+01  3.1490E+00
             1.2506E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1312.43449751492        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      235
 NPARAMETR:  9.6082E-01  7.0727E-01  2.3455E+00  1.2340E+00  1.2528E+00  9.0298E-01  1.0024E+00  5.1213E-01  7.0039E-01  2.8856E-01
             3.4289E+00
 PARAMETER:  6.0034E-02 -2.4635E-01  9.5252E-01  3.1026E-01  3.2536E-01 -2.0577E-03  1.0241E-01 -5.6917E-01 -2.5612E-01 -1.1428E+00
             1.3322E+00
 GRADIENT:   5.4974E+00 -4.3041E+00 -9.1681E-01 -7.3366E+00  1.3533E+00  1.5612E+00  9.3880E-02  1.7752E-01  1.2439E+00  6.4859E-01
             2.9062E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1313.07229405600        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      307
 NPARAMETR:  9.6118E-01  8.6254E-01  4.1004E+00  1.1612E+00  1.4378E+00  9.0107E-01  9.2097E-01  5.0851E-02  6.8864E-01  4.9561E-02
             3.4363E+00
 PARAMETER:  6.0409E-02 -4.7876E-02  1.5111E+00  2.4944E-01  4.6311E-01 -4.1747E-03  1.7669E-02 -2.8789E+00 -2.7304E-01 -2.9046E+00
             1.3344E+00
 GRADIENT:  -3.5643E+00  4.0623E+00 -5.8156E-01  1.0299E+01  2.1039E+00 -4.4593E-01 -1.8920E-01  2.5454E-04 -3.0789E-01  1.5553E-02
            -3.4215E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1313.21478082499        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      380
 NPARAMETR:  9.6367E-01  1.0893E+00  6.1752E+00  1.0108E+00  1.5019E+00  9.0180E-01  8.3288E-01  1.1721E-02  7.2067E-01  1.0000E-02
             3.4376E+00
 PARAMETER:  6.2991E-02  1.8557E-01  1.9205E+00  1.1076E-01  5.0672E-01 -3.3622E-03 -8.2871E-02 -4.3464E+00 -2.2757E-01 -4.8907E+00
             1.3348E+00
 GRADIENT:  -9.5247E-01  3.4690E+00 -1.1368E-01  4.5836E+00 -3.1524E-01 -3.0975E-01 -7.8798E-02  2.5579E-06 -3.2838E-01  0.0000E+00
            -3.9455E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1313.24955352381        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      497
 NPARAMETR:  9.6555E-01  1.1790E+00  7.7401E+00  9.5227E-01  1.5278E+00  9.0266E-01  7.9003E-01  1.0000E-02  7.5120E-01  1.0000E-02
             3.4419E+00
 PARAMETER:  6.4942E-02  2.6469E-01  2.1464E+00  5.1097E-02  5.2382E-01 -2.4112E-03 -1.3569E-01 -4.9483E+00 -1.8608E-01 -5.8439E+00
             1.3360E+00
 GRADIENT:  -1.3233E+00  2.5193E+00 -9.7801E-02  3.1467E+00  1.5207E-01 -3.7570E-01  1.9363E-02  0.0000E+00  6.8707E-02  0.0000E+00
            -4.6128E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1313.28695485017        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      672
 NPARAMETR:  9.6624E-01  1.3251E+00  1.2321E+01  8.5401E-01  1.5522E+00  9.0353E-01  7.4655E-01  1.0000E-02  7.8816E-01  1.0000E-02
             3.4421E+00
 PARAMETER:  6.5659E-02  3.8148E-01  2.6113E+00 -5.7813E-02  5.3965E-01 -1.4447E-03 -1.9229E-01 -6.1538E+00 -1.3806E-01 -7.8013E+00
             1.3361E+00
 GRADIENT:  -1.4385E+00  2.2921E+00 -4.3494E-02  2.2208E+00 -1.6100E-01 -1.8996E-01  1.0447E-02  0.0000E+00  4.3084E-02  0.0000E+00
            -5.7831E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1313.30883268005        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      847
 NPARAMETR:  9.6746E-01  1.4673E+00  2.2744E+01  7.5646E-01  1.5660E+00  9.0457E-01  7.1658E-01  1.0000E-02  8.1714E-01  1.0000E-02
             3.4457E+00
 PARAMETER:  6.6919E-02  4.8345E-01  3.2243E+00 -1.7910E-01  5.4851E-01 -2.9160E-04 -2.3327E-01 -7.6998E+00 -1.0195E-01 -1.0290E+01
             1.3371E+00
 GRADIENT:   4.7696E-01  1.4581E+00 -1.8498E-02  1.0269E+00 -2.7063E-01  1.2369E-01 -9.0708E-02  0.0000E+00 -2.4465E-01  0.0000E+00
            -9.4412E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1313.31791537275        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1022
 NPARAMETR:  9.6723E-01  1.5461E+00  3.4957E+01  7.0170E-01  1.5712E+00  9.0413E-01  6.8830E-01  1.0000E-02  8.7070E-01  1.0000E-02
             3.4445E+00
 PARAMETER:  6.6680E-02  5.3572E-01  3.6541E+00 -2.5425E-01  5.5184E-01 -7.8527E-04 -2.7354E-01 -8.7665E+00 -3.8456E-02 -1.1991E+01
             1.3368E+00
 GRADIENT:  -3.6904E-01  5.2000E-01 -9.8995E-03  4.1348E-01  8.7582E-02 -4.1346E-02 -2.7514E-02  0.0000E+00  9.0601E-02  0.0000E+00
            -9.7401E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1313.32251112548        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1199
 NPARAMETR:  9.6745E-01  1.6268E+00  5.8004E+01  6.4618E-01  1.5734E+00  9.0434E-01  6.7752E-01  1.0000E-02  8.9147E-01  1.0000E-02
             3.4453E+00
 PARAMETER:  6.6906E-02  5.8662E-01  4.1605E+00 -3.3668E-01  5.5327E-01 -5.4530E-04 -2.8932E-01 -1.0060E+01 -1.4882E-02 -1.3980E+01
             1.3370E+00
 GRADIENT:  -2.0168E-01  1.3734E+00 -5.1729E-03  7.6048E-01 -2.1554E-01 -1.1647E-04  1.6836E-02  0.0000E+00 -5.2301E-02  0.0000E+00
            -3.4143E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1313.32519937155        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1375
 NPARAMETR:  9.6751E-01  1.6882E+00  9.1828E+01  6.0301E-01  1.5744E+00  9.0418E-01  6.6118E-01  1.0000E-02  9.3135E-01  1.0000E-02
             3.4450E+00
 PARAMETER:  6.6970E-02  6.2363E-01  4.6199E+00 -4.0582E-01  5.5387E-01 -7.2279E-04 -3.1373E-01 -1.1234E+01  2.8875E-02 -1.5758E+01
             1.3369E+00
 GRADIENT:  -3.4909E-02  6.8008E-01 -2.7058E-03  3.6950E-01 -1.5206E-01 -5.0447E-02 -7.3046E-02  0.0000E+00 -6.5146E-03  0.0000E+00
            -1.1245E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1313.32674887721        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1550
 NPARAMETR:  9.6757E-01  1.7615E+00  1.7221E+02  5.5170E-01  1.5753E+00  9.0428E-01  6.4796E-01  1.0000E-02  9.6975E-01  1.0000E-02
             3.4452E+00
 PARAMETER:  6.7029E-02  6.6617E-01  5.2487E+00 -4.9476E-01  5.5441E-01 -6.2014E-04 -3.3393E-01 -1.2863E+01  6.9283E-02 -1.8168E+01
             1.3370E+00
 GRADIENT:   4.7659E-02  4.7558E-01 -1.2148E-03  2.1909E-01 -1.2944E-01 -1.9661E-02 -1.0598E-01  0.0000E+00 -4.3854E-02  0.0000E+00
            -9.7445E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1313.32739818865        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1726
 NPARAMETR:  9.6758E-01  1.8096E+00  2.7615E+02  5.1788E-01  1.5754E+00  9.0420E-01  6.3985E-01  1.0000E-02  1.0028E+00  1.0000E-02
             3.4450E+00
 PARAMETER:  6.7040E-02  6.9309E-01  5.7210E+00 -5.5800E-01  5.5450E-01 -7.0371E-04 -3.4652E-01 -1.4097E+01  1.0280E-01 -1.9964E+01
             1.3369E+00
 GRADIENT:   1.2010E-01  6.7774E-02 -6.4577E-04  3.5941E-02 -7.2833E-02 -2.9123E-02  1.1369E-02  0.0000E+00  1.3859E-02  0.0000E+00
            -2.5710E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1313.32762545954        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1902
 NPARAMETR:  9.6754E-01  1.8492E+00  4.1979E+02  4.9048E-01  1.5760E+00  9.0432E-01  6.3342E-01  1.0000E-02  1.0269E+00  1.0000E-02
             3.4452E+00
 PARAMETER:  6.6997E-02  7.1473E-01  6.1397E+00 -6.1236E-01  5.5487E-01 -5.7344E-04 -3.5662E-01 -1.5197E+01  1.2658E-01 -2.1546E+01
             1.3370E+00
 GRADIENT:  -8.9711E-02  5.0365E-01 -3.9955E-04  1.8935E-01 -7.9264E-02 -5.4020E-03 -2.5579E-02  0.0000E+00 -1.9417E-02  0.0000E+00
            -3.6244E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1313.32787613971        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2080
 NPARAMETR:  9.6754E-01  1.8927E+00  7.0371E+02  4.5979E-01  1.5761E+00  9.0430E-01  6.2579E-01  1.0000E-02  1.0630E+00  1.0000E-02
             3.4452E+00
 PARAMETER:  6.7005E-02  7.3801E-01  6.6564E+00 -6.7699E-01  5.5493E-01 -5.9162E-04 -3.6874E-01 -1.6563E+01  1.6110E-01 -2.3491E+01
             1.3370E+00
 GRADIENT:  -2.4084E-02  5.5751E-02 -2.0593E-04  2.3689E-02 -1.5112E-04 -1.8419E-04  2.6003E-03  0.0000E+00  3.4608E-03  0.0000E+00
             3.3381E-03

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1313.32792258823        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2257
 NPARAMETR:  9.6753E-01  1.9169E+00  9.5323E+02  4.4298E-01  1.5761E+00  9.0426E-01  6.2178E-01  1.0000E-02  1.0836E+00  1.0000E-02
             3.4452E+00
 PARAMETER:  6.6987E-02  7.5071E-01  6.9599E+00 -7.1423E-01  5.5493E-01 -6.3806E-04 -3.7518E-01 -1.7368E+01  1.8027E-01 -2.4630E+01
             1.3370E+00
 GRADIENT:  -9.9713E-02  2.9010E-01 -1.4147E-04  1.0155E-01 -4.3453E-02 -2.0510E-02  1.3753E-03  0.0000E+00  5.8348E-03  0.0000E+00
            -4.5296E-03

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1313.32798776193        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2434
 NPARAMETR:  9.6755E-01  1.9459E+00  1.4085E+03  4.2268E-01  1.5760E+00  9.0429E-01  6.1726E-01  1.0000E-02  1.1081E+00  1.0000E-02
             3.4452E+00
 PARAMETER:  6.7012E-02  7.6570E-01  7.3502E+00 -7.6113E-01  5.5492E-01 -6.0092E-04 -3.8247E-01 -1.8409E+01  2.0260E-01 -2.6091E+01
             1.3370E+00
 GRADIENT:  -3.8255E-02  2.3906E-01 -8.7086E-05  7.7436E-02 -4.6371E-02 -9.0943E-03 -6.7377E-03  0.0000E+00 -2.0073E-03  0.0000E+00
            -1.5651E-02

0ITERATION NO.:   88    OBJECTIVE VALUE:  -1313.32800003913        NO. OF FUNC. EVALS.: 103
 CUMULATIVE NO. OF FUNC. EVALS.:     2537
 NPARAMETR:  9.6757E-01  1.9581E+00  1.6743E+03  4.1353E-01  1.5765E+00  9.0429E-01  6.1538E-01  1.0000E-02  1.1193E+00  1.0000E-02
             3.4453E+00
 PARAMETER:  6.7028E-02  7.7219E-01  7.5283E+00 -7.8225E-01  5.5494E-01 -5.9632E-04 -3.8558E-01 -1.8884E+01  2.1266E-01 -2.6756E+01
             1.3370E+00
 GRADIENT:  -5.9002E-03  3.1636E-01  2.3168E-04  9.5071E-02 -5.0735E-02  1.9841E-03 -4.9433E-03  0.0000E+00 -8.0322E-03  0.0000E+00
            -2.1617E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2537
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1877E-03 -9.4654E-03  3.6810E-09 -5.2159E-03 -4.3304E-05
 SE:             2.8633E-02  2.0514E-02  2.8712E-09  1.1264E-02  1.0983E-04
 N:                     100         100         100         100         100

 P VAL.:         9.6691E-01  6.4450E-01  1.9983E-01  6.4331E-01  6.9338E-01

 ETASHRINKSD(%)  4.0768E+00  3.1277E+01  1.0000E+02  6.2265E+01  9.9632E+01
 ETASHRINKVR(%)  7.9874E+00  5.2772E+01  1.0000E+02  8.5761E+01  9.9999E+01
 EBVSHRINKSD(%)  4.0086E+00  3.1241E+01  1.0000E+02  6.2247E+01  9.9589E+01
 EBVSHRINKVR(%)  7.8566E+00  5.2722E+01  1.0000E+02  8.5747E+01  9.9998E+01
 RELATIVEINF(%)  6.1748E+01  1.1051E-06  0.0000E+00  3.3305E-07  9.0100E-05
 EPSSHRINKSD(%)  1.9735E+01
 EPSSHRINKVR(%)  3.5576E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1313.3280000391289     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -578.17717347539076     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    30.64
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     6.17
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1313.328       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.68E-01  1.96E+00  1.68E+03  4.14E-01  1.58E+00  9.04E-01  6.15E-01  1.00E-02  1.12E+00  1.00E-02  3.45E+00
 


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
+        1.77E+03
 
 TH 2
+       -7.86E+01  3.88E+02
 
 TH 3
+       -5.44E-05  4.47E-05  8.82E-10
 
 TH 4
+       -9.64E+01  5.61E+02  9.25E-05  8.13E+02
 
 TH 5
+       -2.80E+01 -8.75E+01 -4.07E-05 -1.28E+02  2.23E+01
 
 TH 6
+        8.78E+01  1.15E+00  4.98E-04  2.46E+01 -2.46E+01  3.75E+02
 
 TH 7
+       -7.70E+01  1.78E+01 -4.85E-05  1.80E+01  3.26E+00 -1.11E+02  7.71E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.13E+01  1.87E+00  5.18E-05  4.01E+00 -1.85E+00  2.53E+01  6.00E-01  0.00E+00  3.22E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -4.67E+01  5.23E+00  4.78E-05  7.34E+00 -8.52E-01  6.04E+00  1.58E+01  0.00E+00  3.74E+00  0.00E+00  7.52E+00
 
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
+        1.50E+03
 
 TH 2
+       -8.94E+01  3.56E+02
 
 TH 3
+        8.09E-05  4.57E-05  4.03E-09
 
 TH 4
+       -1.16E+02  5.07E+02 -1.47E-04  7.38E+02
 
 TH 5
+       -1.75E+01 -7.34E+01  3.64E-05 -1.06E+02  7.28E+01
 
 TH 6
+        5.94E+01 -2.40E+01  7.21E-04 -1.65E+01 -5.85E+00  3.25E+02
 
 TH 7
+        3.62E+00 -1.42E+01  3.98E-04 -2.11E+01  4.00E+00  4.51E+00  1.26E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.05E+00 -3.35E+00 -2.40E-04 -9.87E+00  6.96E-02  4.10E+01  1.87E+01  0.00E+00  3.94E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.94E+01 -9.95E+00 -9.58E-06 -1.46E+01  9.38E-01  2.52E+00  1.52E+01  0.00E+00  2.34E+00  0.00E+00  3.79E+01
 
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
+        1.34E+03
 
 TH 2
+       -1.44E+02  3.63E+02
 
 TH 3
+        4.92E-08 -8.84E-07  3.94E-14
 
 TH 4
+       -2.06E+02  5.18E+02 -1.26E-06  7.40E+02
 
 TH 5
+       -1.11E+00 -6.53E+01 -1.08E-06 -9.33E+01  6.08E+01
 
 TH 6
+       -2.09E+01 -4.06E+01  5.41E-07 -5.80E+01 -1.04E+01  1.33E+02
 
 TH 7
+        2.64E+01 -5.64E+01  6.71E-07 -8.06E+01  2.99E+00  2.53E+00  1.40E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.38E+00 -9.35E+00  1.11E-07 -1.34E+01  4.92E-01  4.25E-01  2.32E+01  0.00E+00  3.85E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        1.45E+02 -9.39E+01  3.33E-07 -1.34E+02  9.43E+00  4.63E+01  5.12E+01  0.00E+00  8.49E+00  0.00E+00  2.23E+02
 
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
 #CPUT: Total CPU Time in Seconds,       36.877
Stop Time:
Sat Sep 18 09:48:30 CDT 2021
