Sat Sep 25 10:44:58 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat79.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m79.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1657.86480549968        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.6398E+00 -4.5566E+01 -4.0076E+00 -4.2782E+01 -1.5649E+01  2.6715E+01 -1.3025E+01  3.6170E-01  1.4502E+01  8.5679E+00
            -7.6946E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1670.35686990768        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0122E+00  1.1667E+00  1.0607E+00  9.4039E-01  1.1111E+00  9.0577E-01  1.1481E+00  1.0057E+00  8.7603E-01  9.2209E-01
             1.1776E+00
 PARAMETER:  1.1212E-01  2.5422E-01  1.5893E-01  3.8538E-02  2.0537E-01  1.0348E-03  2.3812E-01  1.0570E-01 -3.2357E-02  1.8892E-02
             2.6351E-01
 GRADIENT:   2.0453E+01  1.5845E+01  1.0098E+01  8.8471E+00  5.0562E+00 -1.0556E+01 -1.9933E+00 -6.6623E+00 -2.2113E+00 -4.6101E+00
            -3.0298E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1670.85976999937        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0125E+00  1.1264E+00  1.0579E+00  9.5903E-01  1.0799E+00  9.1190E-01  1.2147E+00  1.1165E+00  8.3625E-01  8.6186E-01
             1.1780E+00
 PARAMETER:  1.1243E-01  2.1902E-01  1.5628E-01  5.8167E-02  1.7686E-01  7.7760E-03  2.9449E-01  2.1020E-01 -7.8829E-02 -4.8664E-02
             2.6386E-01
 GRADIENT:   2.1910E+01  1.0476E+01  8.4847E+00 -4.2622E-01 -1.6415E+00 -7.6789E+00  1.4581E+00 -3.8382E+00 -2.7850E+00 -5.5141E+00
            -1.9533E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1671.34963481771        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      317
 NPARAMETR:  1.0163E+00  1.0908E+00  1.1922E+00  9.8938E-01  1.1246E+00  9.3403E-01  1.2212E+00  1.3504E+00  8.4044E-01  9.3611E-01
             1.1769E+00
 PARAMETER:  1.1617E-01  1.8687E-01  2.7577E-01  8.9327E-02  2.1745E-01  3.1758E-02  2.9981E-01  4.0044E-01 -7.3829E-02  3.3980E-02
             2.6291E-01
 GRADIENT:  -1.7609E+00  5.6714E-01  2.2637E-01  1.0610E+00 -2.0413E-01 -2.6755E-01 -2.7709E-01 -1.6745E-02 -1.2856E-01 -4.8585E-01
            -3.7483E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1671.47847442857        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      496
 NPARAMETR:  1.0197E+00  1.3162E+00  9.6596E-01  8.4173E-01  1.1401E+00  9.3750E-01  1.0567E+00  1.2231E+00  9.4219E-01  9.3851E-01
             1.1780E+00
 PARAMETER:  1.1955E-01  3.7477E-01  6.5368E-02 -7.2300E-02  2.3109E-01  3.5466E-02  1.5513E-01  3.0138E-01  4.0451E-02  3.6543E-02
             2.6385E-01
 GRADIENT:   1.8126E+00  3.3526E+00  2.0598E+00  1.0141E+00 -2.0631E+00  3.3615E-01  7.1309E-01 -8.5417E-01 -7.8534E-01  3.3851E-01
             2.5582E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1671.86918019187        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      674
 NPARAMETR:  1.0200E+00  1.6290E+00  6.3232E-01  6.3968E-01  1.1446E+00  9.3766E-01  8.8796E-01  9.8829E-01  1.1424E+00  8.9703E-01
             1.1751E+00
 PARAMETER:  1.1979E-01  5.8794E-01 -3.5836E-01 -3.4679E-01  2.3504E-01  3.5630E-02 -1.8834E-02  8.8218E-02  2.3313E-01 -8.6619E-03
             2.6134E-01
 GRADIENT:  -3.0669E+00  2.0658E+01  3.2320E+00  9.3754E+00 -8.1540E+00 -6.0098E-01 -1.7128E+00 -2.5395E-01 -3.3056E-01 -1.3971E+00
            -1.3030E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1672.53957663246        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      851
 NPARAMETR:  1.0205E+00  1.8833E+00  4.2632E-01  4.6092E-01  1.2154E+00  9.3644E-01  7.9638E-01  8.8074E-01  1.3899E+00  9.2989E-01
             1.1777E+00
 PARAMETER:  1.2027E-01  7.3304E-01 -7.5256E-01 -6.7452E-01  2.9507E-01  3.4328E-02 -1.2767E-01 -2.6991E-02  4.2923E-01  2.7315E-02
             2.6357E-01
 GRADIENT:  -3.1313E+00  5.9319E+00 -2.9500E-01  3.4948E+00 -7.2710E-01 -1.1231E+00 -1.2873E+00  4.1031E-01 -1.6318E-01 -1.1866E+00
            -8.2904E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1672.74931485684        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1026
 NPARAMETR:  1.0225E+00  2.0327E+00  3.3396E-01  3.5713E-01  1.2850E+00  9.3883E-01  7.5973E-01  7.7599E-01  1.6272E+00  9.8913E-01
             1.1830E+00
 PARAMETER:  1.2230E-01  8.0935E-01 -9.9673E-01 -9.2965E-01  3.5079E-01  3.6876E-02 -1.7479E-01 -1.5361E-01  5.8685E-01  8.9070E-02
             2.6805E-01
 GRADIENT:   2.1430E+00 -1.4944E+00  8.5384E-02 -9.5381E-01 -2.1735E+00  1.3676E-01  6.1215E-01  3.9793E-01  2.3524E-01  1.7149E-01
             2.2508E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1672.84441580844        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1202
 NPARAMETR:  1.0225E+00  2.0358E+00  3.1664E-01  3.5660E-01  1.2789E+00  9.3925E-01  7.5761E-01  4.7167E-01  1.6158E+00  9.8385E-01
             1.1808E+00
 PARAMETER:  1.2223E-01  8.1089E-01 -1.0500E+00 -9.3115E-01  3.4600E-01  3.7328E-02 -1.7758E-01 -6.5148E-01  5.7985E-01  8.3719E-02
             2.6621E-01
 GRADIENT:   1.6089E+00  1.2578E+00 -1.4076E+00  2.3369E+00  3.9557E+00  1.8551E-01 -4.6315E-01  1.4680E-01 -4.6017E-01 -9.1168E-02
            -4.7849E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1672.92279238270        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1377
 NPARAMETR:  1.0218E+00  2.0450E+00  3.0929E-01  3.4838E-01  1.2798E+00  9.3880E-01  7.5752E-01  9.0882E-02  1.6529E+00  9.8474E-01
             1.1824E+00
 PARAMETER:  1.2160E-01  8.1538E-01 -1.0735E+00 -9.5447E-01  3.4673E-01  3.6845E-02 -1.7771E-01 -2.2982E+00  6.0251E-01  8.4625E-02
             2.6759E-01
 GRADIENT:   1.7887E-01  3.3650E-01 -1.2336E-01  2.4157E-01  2.9417E-01 -7.5354E-03  2.8787E-01  4.9135E-03  8.7026E-02 -5.5976E-02
             7.0452E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1672.92524588445        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1560
 NPARAMETR:  1.0218E+00  2.0442E+00  3.0984E-01  3.4859E-01  1.2798E+00  9.3882E-01  7.5682E-01  3.1602E-02  1.6526E+00  9.8571E-01
             1.1823E+00
 PARAMETER:  1.2154E-01  8.1502E-01 -1.0717E+00 -9.5387E-01  3.4668E-01  3.6873E-02 -1.7863E-01 -3.3545E+00  6.0234E-01  8.5610E-02
             2.6745E-01
 GRADIENT:   2.4271E-02 -2.5572E-01 -2.2390E-02 -1.8671E-02  6.7539E-02  3.6882E-03  1.3782E-02  7.1215E-04  2.0727E-02  3.3850E-02
             2.3554E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1672.92542574908        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1738
 NPARAMETR:  1.0218E+00  2.0443E+00  3.0974E-01  3.4861E-01  1.2796E+00  9.3884E-01  7.5679E-01  1.8114E-02  1.6522E+00  9.8494E-01
             1.1822E+00
 PARAMETER:  1.2156E-01  8.1507E-01 -1.0720E+00 -9.5379E-01  3.4657E-01  3.6892E-02 -1.7866E-01 -3.9111E+00  6.0208E-01  8.4827E-02
             2.6742E-01
 GRADIENT:   7.3883E-02  1.2841E-01 -3.1135E-03  4.7577E-02  1.1261E-01  8.8212E-03 -1.2588E-02  1.5179E-04 -4.1362E-03 -5.9374E-02
            -2.4987E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1672.92550774670        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1913
 NPARAMETR:  1.0218E+00  2.0443E+00  3.0960E-01  3.4859E-01  1.2794E+00  9.3883E-01  7.5681E-01  1.1191E-02  1.6519E+00  9.8513E-01
             1.1822E+00
 PARAMETER:  1.2155E-01  8.1504E-01 -1.0725E+00 -9.5387E-01  3.4641E-01  3.6880E-02 -1.7864E-01 -4.3927E+00  6.0191E-01  8.5020E-02
             2.6742E-01
 GRADIENT:   7.4471E-02 -6.9607E-03 -4.2209E-03 -3.5794E-03  4.5868E-03  8.7986E-03 -5.5173E-04  1.0480E-04 -5.1643E-03 -7.4054E-03
            -5.6489E-03

0ITERATION NO.:   64    OBJECTIVE VALUE:  -1672.92551599784        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     2040
 NPARAMETR:  1.0218E+00  2.0443E+00  3.0964E-01  3.4859E-01  1.2794E+00  9.3882E-01  7.5682E-01  1.0000E-02  1.6520E+00  9.8515E-01
             1.1823E+00
 PARAMETER:  1.2154E-01  8.1504E-01 -1.0724E+00 -9.5385E-01  3.4640E-01  3.6871E-02 -1.7863E-01 -4.5247E+00  6.0201E-01  8.5035E-02
             2.6742E-01
 GRADIENT:   1.5967E-02  2.2589E-02  3.7186E-03 -4.5550E-03 -3.3691E-02  1.7248E-03  7.8338E-04  0.0000E+00  4.3454E-03 -1.8601E-03
            -7.8947E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2040
 NO. OF SIG. DIGITS IN FINAL EST.:  3.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3752E-04 -2.7978E-02 -2.0156E-04  3.4972E-02 -4.2498E-02
 SE:             2.9765E-02  2.6199E-02  6.9181E-05  2.0256E-02  2.0645E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9631E-01  2.8557E-01  3.5741E-03  8.4260E-02  3.9543E-02

 ETASHRINKSD(%)  2.8239E-01  1.2230E+01  9.9768E+01  3.2139E+01  3.0836E+01
 ETASHRINKVR(%)  5.6399E-01  2.2964E+01  9.9999E+01  5.3949E+01  5.2164E+01
 EBVSHRINKSD(%)  6.7956E-01  1.2420E+01  9.9801E+01  3.5525E+01  2.8882E+01
 EBVSHRINKVR(%)  1.3545E+00  2.3297E+01  1.0000E+02  5.8429E+01  4.9423E+01
 RELATIVEINF(%)  9.8587E+01  6.9955E+00  3.4201E-05  2.6435E+00  1.3775E+01
 EPSSHRINKSD(%)  4.2929E+01
 EPSSHRINKVR(%)  6.7429E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1672.9255159978422     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -937.77468943410406     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.32
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.18
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1672.926       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  2.04E+00  3.10E-01  3.49E-01  1.28E+00  9.39E-01  7.57E-01  1.00E-02  1.65E+00  9.85E-01  1.18E+00
 


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
+        1.19E+03
 
 TH 2
+       -9.52E+00  3.69E+02
 
 TH 3
+        7.90E+00  1.31E+02  5.36E+02
 
 TH 4
+       -2.75E+01  3.21E+02 -5.92E+02  1.48E+03
 
 TH 5
+       -6.37E+00 -1.12E+02 -2.98E+02  3.25E+02  3.62E+02
 
 TH 6
+       -3.16E+00 -2.05E+00  1.50E+00 -7.67E+00 -1.67E+00  2.26E+02
 
 TH 7
+        1.70E+00  7.98E+00 -2.31E+01 -2.22E+01 -1.59E+00 -6.99E-01  2.21E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.45E+00 -1.16E+01 -4.33E+01  7.93E+01 -2.81E+00 -4.81E-01  9.03E+00  0.00E+00  2.42E+01
 
 TH10
+        4.23E-01 -1.18E+01 -2.75E+01  6.22E+00 -5.72E+01 -8.51E-01  7.14E+00  0.00E+00  4.64E+00  6.75E+01
 
 TH11
+       -9.41E+00 -1.66E+01 -2.55E+01  1.32E+01 -7.73E+00  3.01E+00  1.19E+01  0.00E+00  4.29E+00  1.61E+01  1.56E+02
 
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
 #CPUT: Total CPU Time in Seconds,       31.565
Stop Time:
Sat Sep 25 10:45:31 CDT 2021
