Wed Sep 29 07:06:03 CDT 2021
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
$DATA ../../../../data/int/TD2/dat18.csv ignore=@
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
Current Date:       29 SEP 2021
Days until program expires : 200
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
 RAW OUTPUT FILE (FILE): m18.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3347.04940041073        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4351E+02 -8.2170E+00 -1.3830E+01  9.8165E+01  1.4777E+02  4.0938E+01 -6.0072E+01 -4.4335E+01 -7.7694E+01 -4.5922E+01
            -8.7883E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3489.72787769252        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:      106
 NPARAMETR:  9.9315E-01  1.0651E+00  9.9652E-01  9.9358E-01  9.0096E-01  1.0061E+00  1.2093E+00  1.0272E+00  1.2271E+00  1.1429E+00
             1.3552E+00
 PARAMETER:  9.3124E-02  1.6305E-01  9.6509E-02  9.3560E-02 -4.2913E-03  1.0605E-01  2.9004E-01  1.2681E-01  3.0463E-01  2.3354E-01
             4.0396E-01
 GRADIENT:  -1.2912E+01  2.0869E+01  1.0617E+01 -1.9844E+00 -9.2096E+01  4.5771E+00 -5.2395E+00 -1.2448E+00  7.1278E-01 -3.3017E-01
            -1.0912E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3499.30129971212        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      283
 NPARAMETR:  1.0079E+00  1.4716E+00  1.3698E+00  7.7685E-01  1.3314E+00  9.6829E-01  1.0714E+00  1.4244E+00  1.1346E+00  1.5398E+00
             1.4132E+00
 PARAMETER:  1.0788E-01  4.8633E-01  4.1466E-01 -1.5251E-01  3.8626E-01  6.7777E-02  1.6893E-01  4.5376E-01  2.2630E-01  5.3168E-01
             4.4584E-01
 GRADIENT:   1.7289E+01  1.8608E+01  3.7684E+00  3.3219E+00 -1.7805E+01 -1.0530E+01 -1.4298E+01 -9.8903E+00 -7.7195E+00  1.7888E+01
             5.1919E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3504.10592786914        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      464
 NPARAMETR:  1.0024E+00  1.4704E+00  1.6479E+00  7.8260E-01  1.4302E+00  1.0336E+00  1.1384E+00  2.3302E+00  1.0411E+00  1.5806E+00
             1.3740E+00
 PARAMETER:  1.0238E-01  4.8552E-01  5.9952E-01 -1.4513E-01  4.5785E-01  1.3301E-01  2.2963E-01  9.4594E-01  1.4027E-01  5.5782E-01
             4.1775E-01
 GRADIENT:   5.0068E+00  6.7482E+00 -8.2587E+00  1.0030E+01 -1.8164E+00  1.4489E+01 -7.7526E+00  1.5125E+00  3.6855E+00  1.1061E+01
             1.3382E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3504.80955248587        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      646
 NPARAMETR:  1.0016E+00  1.4626E+00  1.6688E+00  7.8672E-01  1.4347E+00  9.8771E-01  1.2074E+00  2.3506E+00  9.7053E-01  1.5742E+00
             1.3721E+00
 PARAMETER:  1.0162E-01  4.8019E-01  6.1212E-01 -1.3988E-01  4.6095E-01  8.7631E-02  2.8850E-01  9.5469E-01  7.0091E-02  5.5372E-01
             4.1635E-01
 GRADIENT:   3.5813E+00  9.5697E+00 -9.3814E+00  5.6656E+00  7.8111E-01 -2.6882E+00  2.5296E+00 -8.6204E-01 -8.0658E-01  8.8645E+00
             1.0957E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3504.89523467875        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      824
 NPARAMETR:  1.0015E+00  1.4610E+00  1.6790E+00  7.8663E-01  1.4353E+00  9.8570E-01  1.1769E+00  2.3438E+00  9.9562E-01  1.5705E+00
             1.3698E+00
 PARAMETER:  1.0153E-01  4.7909E-01  6.1821E-01 -1.3999E-01  4.6137E-01  8.5598E-02  2.6291E-01  9.5178E-01  9.5606E-02  5.5138E-01
             4.1466E-01
 GRADIENT:   3.3721E+00  6.7007E+00 -8.1604E+00  6.1709E+00  8.6412E-01 -3.4941E+00 -2.6791E+00 -4.0584E-01  7.9559E-01  8.7615E+00
             7.1512E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3505.40695156075        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1008
 NPARAMETR:  1.0012E+00  1.4598E+00  1.6934E+00  7.8678E-01  1.4372E+00  9.9530E-01  1.1925E+00  2.3378E+00  9.8726E-01  1.5624E+00
             1.3677E+00
 PARAMETER:  1.0116E-01  4.7828E-01  6.2673E-01 -1.3981E-01  4.6268E-01  9.5290E-02  2.7609E-01  9.4921E-01  8.7176E-02  5.4622E-01
             4.1312E-01
 GRADIENT:   2.7089E+00  6.5654E+00 -1.3013E+01  5.1507E+00 -1.3247E+00  2.8433E-01 -7.7636E-02 -3.5501E-01 -4.7356E-01  7.5883E+00
             2.3000E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3505.95640436290        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1188
 NPARAMETR:  1.0002E+00  1.4479E+00  1.7772E+00  7.8668E-01  1.4546E+00  9.9690E-01  1.1992E+00  2.3990E+00  9.7317E-01  1.4833E+00
             1.3654E+00
 PARAMETER:  1.0018E-01  4.7010E-01  6.7502E-01 -1.3993E-01  4.7470E-01  9.6891E-02  2.8163E-01  9.7507E-01  7.2808E-02  4.9428E-01
             4.1143E-01
 GRADIENT:   5.8426E-01 -6.8549E-01 -1.0371E+01 -7.4280E+00  4.9849E+00  9.1667E-01  1.5025E-01 -2.0986E-03 -7.7499E-01 -7.6731E+00
            -1.9006E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3506.64560255498        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:     1322
 NPARAMETR:  9.9150E-01  1.4145E+00  2.0451E+00  8.0548E-01  1.4317E+00  9.8914E-01  1.1955E+00  2.4310E+00  9.7152E-01  1.5401E+00
             1.3659E+00
 PARAMETER:  9.1460E-02  4.4675E-01  8.1544E-01 -1.1632E-01  4.5888E-01  8.9078E-02  2.7857E-01  9.8830E-01  7.1102E-02  5.3184E-01
             4.1181E-01
 GRADIENT:   2.1883E+02  2.6810E+02  1.3858E+01  2.4380E+01  7.6760E+01  2.0351E+01  1.6723E+01  4.5480E+00  4.6149E+00  3.2865E+01
             1.0549E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3509.51936041675        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1506
 NPARAMETR:  9.9947E-01  1.3111E+00  2.6055E+00  9.0976E-01  1.5009E+00  9.9276E-01  1.3575E+00  2.9676E+00  7.5433E-01  1.4536E+00
             1.3526E+00
 PARAMETER:  9.9471E-02  3.7087E-01  1.0576E+00  5.4216E-03  5.0604E-01  9.2731E-02  4.0564E-01  1.1878E+00 -1.8193E-01  4.7402E-01
             4.0204E-01
 GRADIENT:  -7.1477E-01  2.5690E+01 -8.3679E+00  2.5315E+01  2.1970E+01 -7.3030E-01  9.3188E-02 -2.9450E+00 -2.3296E+00 -1.1855E+01
            -3.2864E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3512.39269389567        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1685
 NPARAMETR:  9.9893E-01  1.1332E+00  3.6021E+00  1.0130E+00  1.4493E+00  9.9472E-01  1.5249E+00  3.5144E+00  6.4486E-01  1.5406E+00
             1.3409E+00
 PARAMETER:  9.8930E-02  2.2503E-01  1.3815E+00  1.1290E-01  4.7111E-01  9.4707E-02  5.2193E-01  1.3569E+00 -3.3873E-01  5.3214E-01
             3.9332E-01
 GRADIENT:  -7.5769E-01  1.1332E+01 -2.9673E+00  9.2694E+00 -1.9660E+01  1.7548E-01 -3.2416E+00  9.3144E+00  3.5820E-01  2.8850E+00
            -6.4748E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -3512.57978940296        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     1846
 NPARAMETR:  9.9951E-01  1.1293E+00  3.5930E+00  1.0081E+00  1.4540E+00  9.9353E-01  1.5457E+00  3.4858E+00  6.2810E-01  1.5451E+00
             1.3436E+00
 PARAMETER:  9.9505E-02  2.2157E-01  1.3790E+00  1.0807E-01  4.7434E-01  9.3506E-02  5.3549E-01  1.3487E+00 -3.6505E-01  5.3507E-01
             3.9535E-01
 GRADIENT:   5.3759E-01  6.1448E+00 -2.6763E+00 -1.9256E+00 -1.7952E+01 -2.6031E-01 -7.2106E-01  7.4977E+00 -2.6121E-01  3.2561E+00
             3.2660E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -3513.07118374223        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2022
 NPARAMETR:  9.9836E-01  1.1104E+00  3.7336E+00  1.0175E+00  1.4870E+00  9.9392E-01  1.5565E+00  3.3503E+00  6.3304E-01  1.5216E+00
             1.3429E+00
 PARAMETER:  9.8360E-02  2.0475E-01  1.4174E+00  1.1732E-01  4.9676E-01  9.3899E-02  5.4243E-01  1.3090E+00 -3.5723E-01  5.1979E-01
             3.9486E-01
 GRADIENT:  -2.0449E+00  1.7143E+00 -3.3015E-01 -4.4018E+00  2.4247E-01 -1.6109E-01 -1.0069E+00 -5.2445E-01  2.0199E-01 -2.1668E-01
             8.4101E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -3513.19572769701        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2206            RESET HESSIAN, TYPE II
 NPARAMETR:  9.9914E-01  1.0807E+00  3.8525E+00  1.0398E+00  1.4819E+00  9.9441E-01  1.6073E+00  3.3244E+00  6.1217E-01  1.5190E+00
             1.3414E+00
 PARAMETER:  9.9139E-02  1.7762E-01  1.4487E+00  1.3901E-01  4.9335E-01  9.4389E-02  5.7456E-01  1.3013E+00 -3.9075E-01  5.1804E-01
             3.9372E-01
 GRADIENT:   2.4624E+02  6.9805E+01  2.7453E+01  9.1897E+01  1.0419E+02  2.2915E+01  5.0964E+01  2.3395E+01  7.8731E+00  3.0608E+01
             6.9931E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -3513.25759125514        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     2362
 NPARAMETR:  9.9920E-01  1.0485E+00  3.8489E+00  1.0560E+00  1.4694E+00  9.9483E-01  1.6411E+00  3.3243E+00  6.0513E-01  1.5022E+00
             1.3410E+00
 PARAMETER:  9.9200E-02  1.4733E-01  1.4478E+00  1.5444E-01  4.8482E-01  9.4818E-02  5.9540E-01  1.3013E+00 -4.0231E-01  5.0694E-01
             3.9343E-01
 GRADIENT:   4.2887E-02  7.0242E-03 -1.8615E+00  1.1666E+00  1.9911E-01  2.6454E-01  7.8920E-02 -2.2601E+00  3.7378E-01 -5.5094E-02
             4.6702E-01

0ITERATION NO.:   74    OBJECTIVE VALUE:  -3513.30552158670        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:     2496
 NPARAMETR:  9.9977E-01  1.0504E+00  3.8933E+00  1.0536E+00  1.4703E+00  9.9430E-01  1.6427E+00  3.3722E+00  5.9928E-01  1.5066E+00
             1.3406E+00
 PARAMETER:  9.9602E-02  1.4945E-01  1.4617E+00  1.5242E-01  4.8465E-01  9.4290E-02  5.9617E-01  1.3177E+00 -4.0794E-01  5.0900E-01
             3.9327E-01
 GRADIENT:  -4.7081E+03  3.1502E+03  3.2544E+02  3.0880E+03 -9.7147E+02  6.2506E-03 -1.1644E-01  3.5234E+02  3.3133E-01 -9.3025E+02
             9.6522E-01
 NUMSIGDIG:         2.3         2.3         2.3         2.3         2.3         3.8         3.1         2.3         1.5         2.3
                    2.9

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2496
 NO. OF SIG. DIGITS IN FINAL EST.:  1.5

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.0310E-04 -1.1305E-03 -4.5692E-02 -6.4546E-03 -3.1330E-02
 SE:             2.9792E-02  2.5834E-02  2.1539E-02  1.7693E-02  2.5019E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8653E-01  9.6510E-01  3.3890E-02  7.1525E-01  2.1049E-01

 ETASHRINKSD(%)  1.9197E-01  1.3454E+01  2.7842E+01  4.0727E+01  1.6182E+01
 ETASHRINKVR(%)  3.8357E-01  2.5097E+01  4.7933E+01  6.4867E+01  2.9746E+01
 EBVSHRINKSD(%)  5.0037E-01  1.3299E+01  2.8699E+01  4.4834E+01  1.3184E+01
 EBVSHRINKVR(%)  9.9824E-01  2.4830E+01  4.9161E+01  6.9567E+01  2.4630E+01
 RELATIVEINF(%)  9.8997E+01  2.2527E+01  3.1565E+01  8.2557E+00  6.1651E+01
 EPSSHRINKSD(%)  2.0341E+01
 EPSSHRINKVR(%)  3.6544E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3513.3055215867021     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1859.2161618182913     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    77.70
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    15.01
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3513.306       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.05E+00  3.90E+00  1.05E+00  1.47E+00  9.94E-01  1.64E+00  3.38E+00  6.02E-01  1.51E+00  1.34E+00
 


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
+        2.36E+06
 
 TH 2
+       -7.50E+05  4.78E+05
 
 TH 3
+       -4.13E+04  2.08E+00  7.27E+02
 
 TH 4
+       -7.33E+05  4.67E+05 -6.75E+00  4.57E+05
 
 TH 5
+       -6.33E+02 -1.05E+05  3.50E+00  4.54E+02  2.33E+04
 
 TH 6
+        7.73E+00 -4.58E+00 -1.36E-01 -4.42E+00 -3.29E-01  1.97E+02
 
 TH 7
+        1.10E+01  1.50E+01 -4.97E-01 -4.83E+01  2.85E+00 -2.36E-01  4.42E+01
 
 TH 8
+        6.31E+02 -1.71E+00 -1.33E+01  1.63E+01  1.09E+01 -7.86E-02 -4.60E-01  5.87E+02
 
 TH 9
+        1.60E+01 -3.01E+01 -6.99E-01 -3.46E+01  6.72E+04 -3.16E-01  2.46E+01  4.50E+00  7.39E+01
 
 TH10
+        1.54E+05 -9.83E+04 -2.85E+01 -9.61E+04  2.16E+04  5.97E-01  6.70E-01 -3.77E+01  5.30E+00  2.03E+04
 
 TH11
+       -2.23E+05  1.61E+02  5.70E+00  1.75E+02 -3.13E+04  2.30E+00  3.78E+00  1.13E+01 -9.10E+04 -2.33E+01  4.30E+04
 
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
 #CPUT: Total CPU Time in Seconds,       92.791
Stop Time:
Wed Sep 29 07:07:37 CDT 2021
