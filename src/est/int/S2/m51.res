Sat Sep 25 03:19:19 CDT 2021
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
$DATA ../../../../data/int/S2/dat51.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m51.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3360.74499473139        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.5333E+02  5.5474E+01  9.5973E+01  1.4537E+01  9.4687E+01  1.9860E+01 -1.7567E+01 -5.1415E+02 -1.5356E+02 -6.7503E+01
            -1.7424E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3637.61240877387        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.3380E-01  1.0444E+00  1.0601E+00  9.8056E-01  9.8999E-01  9.3812E-01  1.0863E+00  2.4498E+00  9.8374E-01  1.1644E+00
             1.0840E+00
 PARAMETER:  3.1508E-02  1.4344E-01  1.5833E-01  8.0370E-02  8.9937E-02  3.6124E-02  1.8277E-01  9.9600E-01  8.3609E-02  2.5220E-01
             1.8069E-01
 GRADIENT:  -9.6611E+00  4.5637E+01  1.4441E+00  2.3855E+01 -7.1815E+01  9.2781E-01  3.4074E+01 -5.1243E+01 -1.8113E+01  1.4204E+00
             6.1742E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3646.36770403617        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      163
 NPARAMETR:  9.3359E-01  1.0369E+00  1.3539E+00  1.0191E+00  1.0957E+00  9.4969E-01  1.0283E+00  2.5022E+00  9.2742E-01  1.2233E+00
             1.0953E+00
 PARAMETER:  3.1283E-02  1.3622E-01  4.0300E-01  1.1893E-01  1.9138E-01  4.8378E-02  1.2790E-01  1.0172E+00  2.4652E-02  3.0152E-01
             1.9107E-01
 GRADIENT:  -1.0193E+01  3.4912E+01  1.4636E+01  8.5649E+01 -6.2626E+00  5.6499E+00  1.7272E+01 -8.4311E+01 -1.3571E+01  1.5522E+01
             8.2535E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3654.14658796889        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:      297             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5091E-01  1.0141E+00  1.3823E+00  9.8344E-01  1.1094E+00  9.4088E-01  8.9542E-01  2.6655E+00  9.7825E-01  1.1310E+00
             1.0935E+00
 PARAMETER:  4.9668E-02  1.1396E-01  4.2374E-01  8.3303E-02  2.0384E-01  3.9063E-02 -1.0468E-02  1.0804E+00  7.8014E-02  2.2306E-01
             1.8942E-01
 GRADIENT:   3.5262E+01 -2.9793E+01  1.7973E+01  2.9039E+00  1.4935E+01  2.8378E+00  3.9003E+00 -5.6439E+01 -6.0480E+00 -3.1334E+00
             7.3826E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3656.45014177213        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  9.4131E-01  1.0731E+00  1.3823E+00  9.4775E-01  1.1211E+00  9.3497E-01  7.5358E-01  2.6656E+00  1.0585E+00  1.1933E+00
             1.0935E+00
 PARAMETER:  3.9513E-02  1.7056E-01  4.2375E-01  4.6333E-02  2.1435E-01  3.2756E-02 -1.8292E-01  1.0804E+00  1.5688E-01  2.7673E-01
             1.8943E-01
 GRADIENT:   9.4620E+00  2.7372E+00  2.6282E+01 -5.2329E+00 -3.5994E+00  6.4582E-02 -1.7293E-01 -5.5644E+01 -5.6246E-01 -4.6754E+00
             7.1590E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3656.92607706859        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      525
 NPARAMETR:  9.5247E-01  1.0954E+00  1.3823E+00  9.4087E-01  1.1428E+00  9.4652E-01  7.4476E-01  2.6658E+00  1.0770E+00  1.2356E+00
             1.0935E+00
 PARAMETER:  5.1300E-02  1.9110E-01  4.2375E-01  3.9048E-02  2.3348E-01  4.5035E-02 -1.9470E-01  1.0805E+00  1.7417E-01  3.1153E-01
             1.8942E-01
 GRADIENT:   1.8109E+00 -1.0993E+00  2.5978E+01 -1.7877E+00 -1.1195E+00  4.4931E-01  1.5746E-01 -6.2140E+01  2.1408E-01  5.3619E-01
             6.8623E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3657.42321948578        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      706
 NPARAMETR:  9.4334E-01  1.0920E+00  1.3814E+00  9.4857E-01  1.1303E+00  9.7331E-01  7.5638E-01  2.7042E+00  1.0455E+00  1.2262E+00
             1.0931E+00
 PARAMETER:  4.1666E-02  1.8797E-01  4.2310E-01  4.7200E-02  2.2246E-01  7.2952E-02 -1.7922E-01  1.0948E+00  1.4449E-01  3.0394E-01
             1.8906E-01
 GRADIENT:  -2.0701E+01  1.0082E+01  2.4033E+01  5.0240E+00 -1.2318E+01  1.0889E+01 -7.5390E-01 -5.7191E+01 -5.3380E+00 -6.1150E-01
             6.9499E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3659.54478611111        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:      891
 NPARAMETR:  9.3796E-01  1.0929E+00  1.3780E+00  9.5054E-01  1.1282E+00  9.9099E-01  7.7641E-01  2.8520E+00  1.0211E+00  1.2290E+00
             1.0916E+00
 PARAMETER:  3.5957E-02  1.8879E-01  4.2064E-01  4.9278E-02  2.2061E-01  9.0948E-02 -1.5308E-01  1.1480E+00  1.2089E-01  3.0619E-01
             1.8766E-01
 GRADIENT:  -3.2767E+01  1.5343E+01  1.8238E+01  7.3987E+00 -1.8680E+01  1.7121E+01 -8.6901E-01 -3.9508E+01 -8.1776E+00 -6.8647E-01
             6.8804E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3660.36448935071        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     1060
 NPARAMETR:  9.5107E-01  1.0929E+00  1.3781E+00  9.4554E-01  1.1283E+00  9.4691E-01  7.7640E-01  2.8524E+00  1.0537E+00  1.2290E+00
             1.0917E+00
 PARAMETER:  4.9829E-02  1.8888E-01  4.2073E-01  4.3999E-02  2.2075E-01  4.5454E-02 -1.5309E-01  1.1482E+00  1.5234E-01  3.0616E-01
             1.8770E-01
 GRADIENT:  -1.7362E+00  8.4408E+00  1.8731E+01  3.6885E-01 -1.9241E+01  5.7231E-01  1.1135E+00 -3.8557E+01 -6.9446E-01 -1.0490E+00
             6.8613E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3660.40924011058        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1239
 NPARAMETR:  9.5391E-01  1.0930E+00  1.3780E+00  9.4472E-01  1.1284E+00  9.4135E-01  7.7639E-01  2.8566E+00  1.0657E+00  1.2289E+00
             1.0916E+00
 PARAMETER:  5.2811E-02  1.8889E-01  4.2067E-01  4.3131E-02  2.2083E-01  3.9563E-02 -1.5310E-01  1.1496E+00  1.6359E-01  3.0614E-01
             1.8764E-01
 GRADIENT:   5.6839E+00  6.7897E+00  1.8692E+01 -3.9730E-01 -1.9376E+01 -1.7517E+00  1.7306E+00 -3.7735E+01  2.0262E+00 -1.2094E+00
             6.8506E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3660.43413520123        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1424
 NPARAMETR:  9.5158E-01  1.0930E+00  1.3782E+00  9.4437E-01  1.1285E+00  9.4610E-01  7.5218E-01  2.8561E+00  1.0625E+00  1.2288E+00
             1.0916E+00
 PARAMETER:  5.0367E-02  1.8894E-01  4.2077E-01  4.2758E-02  2.2090E-01  4.4594E-02 -1.8477E-01  1.1495E+00  1.6064E-01  3.0606E-01
             1.8768E-01
 GRADIENT:  -4.1777E-01  6.3307E+00  1.9257E+01 -2.2790E-01 -1.8416E+01  2.4469E-01 -2.9351E-01 -3.7654E+01 -4.2129E-01 -2.1536E+00
             6.7724E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -3661.12815649187        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1605
 NPARAMETR:  9.4935E-01  1.0935E+00  1.3775E+00  9.4141E-01  1.1305E+00  9.4480E-01  7.0945E-01  2.9248E+00  1.0497E+00  1.2277E+00
             1.0909E+00
 PARAMETER:  4.8021E-02  1.8943E-01  4.2026E-01  3.9623E-02  2.2263E-01  4.3220E-02 -2.4327E-01  1.1732E+00  1.4848E-01  3.0516E-01
             1.8697E-01
 GRADIENT:  -6.2577E+00  3.7840E+00  1.7486E+01 -3.8274E+00 -1.7779E+01 -3.1636E-01 -4.2237E+00 -2.9855E+01 -6.3081E+00 -5.1058E+00
             6.5023E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -3661.20205462693        NO. OF FUNC. EVALS.: 205
 CUMULATIVE NO. OF FUNC. EVALS.:     1810
 NPARAMETR:  9.5217E-01  1.0936E+00  1.3776E+00  9.4234E-01  1.1305E+00  9.4575E-01  7.0960E-01  2.9240E+00  1.0726E+00  1.2277E+00
             1.0909E+00
 PARAMETER:  5.0985E-02  1.8947E-01  4.2036E-01  4.0606E-02  2.2269E-01  4.4228E-02 -2.4305E-01  1.1729E+00  1.7006E-01  3.0514E-01
             1.8701E-01
 GRADIENT:   1.1187E+00  2.7167E+00  1.7467E+01 -1.0536E+00 -1.7834E+01  1.0064E-01 -3.0281E+00 -2.9313E+01 -6.7483E-01 -5.3318E+00
             6.5181E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -3661.22987441760        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1990
 NPARAMETR:  9.5085E-01  1.0936E+00  1.3776E+00  9.4405E-01  1.1306E+00  9.4443E-01  7.0959E-01  2.9271E+00  1.0833E+00  1.2277E+00
             1.0909E+00
 PARAMETER:  4.9605E-02  1.8949E-01  4.2032E-01  4.2426E-02  2.2277E-01  4.2824E-02 -2.4307E-01  1.1740E+00  1.7997E-01  3.0512E-01
             1.8697E-01
 GRADIENT:  -2.2816E+00  3.3519E+00  1.7235E+01  2.2916E+00 -1.7779E+01 -4.6141E-01 -2.5541E+00 -2.8662E+01  1.9602E+00 -5.4626E+00
             6.5190E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -3661.23588378517        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2175
 NPARAMETR:  9.5180E-01  1.0936E+00  1.3776E+00  9.4291E-01  1.1306E+00  9.4557E-01  7.0957E-01  2.9269E+00  1.0753E+00  1.2276E+00
             1.0909E+00
 PARAMETER:  5.0603E-02  1.8950E-01  4.2035E-01  4.1219E-02  2.2279E-01  4.4031E-02 -2.4310E-01  1.1739E+00  1.7265E-01  3.0509E-01
             1.8698E-01
 GRADIENT:   1.7847E-01  3.0141E+00  1.7319E+01  7.5901E-02 -1.7806E+01  2.0564E-02 -2.9155E+00 -2.8917E+01  3.9402E-02 -5.4042E+00
             6.5156E+01

0ITERATION NO.:   71    OBJECTIVE VALUE:  -3661.23588378517        NO. OF FUNC. EVALS.:  26
 CUMULATIVE NO. OF FUNC. EVALS.:     2201
 NPARAMETR:  9.5174E-01  1.0937E+00  1.3779E+00  9.4287E-01  1.1308E+00  9.4553E-01  7.0949E-01  2.9254E+00  1.0752E+00  1.2275E+00
             1.0910E+00
 PARAMETER:  5.0603E-02  1.8950E-01  4.2035E-01  4.1219E-02  2.2279E-01  4.4031E-02 -2.4310E-01  1.1739E+00  1.7265E-01  3.0509E-01
             1.8698E-01
 GRADIENT:   1.0162E-01 -6.3645E+03 -2.8562E+03  7.6815E-02 -2.7278E+03  1.9554E-02  4.9548E+03  1.0089E+03  2.7217E-02  3.9466E+03
            -6.4043E+03
 NUMSIGDIG:         3.1         3.3         3.3         3.3         3.3         3.3         3.3         3.3         2.9         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2201
 NO. OF SIG. DIGITS IN FINAL EST.:  2.9

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.4572E-04 -4.7999E-02 -5.1027E-02  1.8389E-02 -3.1915E-02
 SE:             2.9879E-02  1.8247E-02  2.6871E-02  2.7211E-02  2.5488E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8810E-01  8.5256E-03  5.7573E-02  4.9916E-01  2.1051E-01

 ETASHRINKSD(%)  1.0000E-10  3.8871E+01  9.9775E+00  8.8400E+00  1.4613E+01
 ETASHRINKVR(%)  1.0000E-10  6.2632E+01  1.8960E+01  1.6899E+01  2.7091E+01
 EBVSHRINKSD(%)  3.4007E-01  4.0546E+01  1.8814E+01  9.2439E+00  1.4754E+01
 EBVSHRINKVR(%)  6.7898E-01  6.4652E+01  3.4088E+01  1.7633E+01  2.7331E+01
 RELATIVEINF(%)  9.9318E+01  1.7474E+01  5.9099E+01  4.9206E+01  4.7797E+01
 EPSSHRINKSD(%)  2.4672E+01
 EPSSHRINKVR(%)  4.3257E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3661.2358837851666     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2007.1465240167558     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    65.98
 Elapsed covariance  time in seconds:    15.50
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3661.236       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.52E-01  1.09E+00  1.38E+00  9.43E-01  1.13E+00  9.46E-01  7.10E-01  2.93E+00  1.08E+00  1.23E+00  1.09E+00
 


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
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         2.71E-02  4.08E-04  1.14E-03  3.38E-02  4.97E-04  7.62E-02  3.43E-04  6.82E-03  1.45E-06  7.43E-04  4.01E-04
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        7.33E-04
 
 TH 2
+        1.45E-06  1.67E-07
 
 TH 3
+        3.90E-06  4.66E-07  1.30E-06
 
 TH 4
+        2.25E-05  2.12E-06  5.68E-06  1.15E-03
 
 TH 5
+        1.73E-06  2.03E-07  5.67E-07  2.53E-06  2.47E-07
 
 TH 6
+       -1.15E-04  2.44E-07  3.55E-07 -2.48E-05  3.07E-07  5.81E-03
 
 TH 7
+       -1.25E-06 -1.40E-07 -3.92E-07 -1.82E-06 -1.70E-07 -3.25E-07  1.18E-07
 
 TH 8
+       -1.80E-05 -2.76E-06 -7.75E-06 -3.02E-05 -3.37E-06  1.45E-05  2.32E-06  4.65E-05
 
 TH 9
+       -1.12E-09 -4.28E-10 -1.20E-09 -1.80E-09 -5.20E-10 -1.87E-09  3.59E-10  7.26E-09  2.11E-12
 
 TH10
+       -2.65E-06 -3.03E-07 -8.48E-07 -3.89E-06 -3.69E-07 -5.91E-07  2.55E-07  5.03E-06  7.78E-10  5.52E-07
 
 TH11
+        1.38E-06  1.64E-07  4.58E-07  2.03E-06  1.99E-07  2.31E-07 -1.38E-07 -2.72E-06 -4.21E-10 -2.98E-07  1.61E-07
 
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
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        2.71E-02
 
 TH 2
+        1.31E-01  4.08E-04
 
 TH 3
+        1.26E-01  1.00E+00  1.14E-03
 
 TH 4
+        2.45E-02  1.53E-01  1.47E-01  3.38E-02
 
 TH 5
+        1.29E-01  1.00E+00  1.00E+00  1.51E-01  4.97E-04
 
 TH 6
+       -5.57E-02  7.85E-03  4.07E-03 -9.61E-03  8.12E-03  7.62E-02
 
 TH 7
+       -1.35E-01 -1.00E+00 -1.00E+00 -1.57E-01 -1.00E+00 -1.24E-02  3.43E-04
 
 TH 8
+       -9.75E-02 -9.93E-01 -9.95E-01 -1.31E-01 -9.94E-01  2.79E-02  9.93E-01  6.82E-03
 
 TH 9
+       -2.86E-02 -7.22E-01 -7.22E-01 -3.65E-02 -7.21E-01 -1.69E-02  7.21E-01  7.32E-01  1.45E-06
 
 TH10
+       -1.32E-01 -1.00E+00 -1.00E+00 -1.55E-01 -1.00E+00 -1.04E-02  1.00E+00  9.94E-01  7.21E-01  7.43E-04
 
 TH11
+        1.27E-01  1.00E+00  1.00E+00  1.49E-01  1.00E+00  7.54E-03 -1.00E+00 -9.93E-01 -7.22E-01 -1.00E+00  4.01E-04
 
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
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.05E+04
 
 TH 2
+        1.32E+07  2.05E+11
 
 TH 3
+        1.30E+07 -2.38E+10  5.39E+10
 
 TH 4
+        7.30E+03  8.67E+06  1.30E+07  6.97E+03
 
 TH 5
+        1.47E+07  5.36E+10  3.10E+09  1.11E+07  6.97E+10
 
 TH 6
+        2.93E+03  5.92E+06  4.14E+06  2.34E+03  4.77E+06  1.12E+03
 
 TH 7
+        9.38E+07  1.29E+11  1.31E+11  7.37E+07  1.59E+11  2.95E+07  1.11E+12
 
 TH 8
+        5.02E+03 -7.08E+08  5.46E+08  4.52E+04 -1.76E+08 -1.41E+03  3.00E+08  1.01E+07
 
 TH 9
+       -1.09E+07  1.27E+11 -9.83E+10 -1.59E+07  9.82E+09 -1.56E+06 -9.32E+10 -1.60E+09  1.28E+12
 
 TH10
+        1.31E+06  8.09E+10 -4.87E+09  2.27E+06  1.97E+10  1.36E+06 -8.09E+10 -2.26E+08  2.85E+10  1.01E+11
 
 TH11
+        1.40E+07  4.09E+10 -2.17E+10  8.29E+06  1.93E+10  3.92E+06  1.05E+11 -6.13E+08  8.84E+10  2.13E+10  1.15E+11
 
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
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,       81.592
Stop Time:
Sat Sep 25 03:20:42 CDT 2021
