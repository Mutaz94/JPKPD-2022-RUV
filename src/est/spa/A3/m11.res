Sat Sep 25 09:04:08 CDT 2021
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
$DATA ../../../../data/spa/A3/dat11.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m11.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -13.0195381075263        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.3390E+02  1.1772E+00  5.9448E+01 -9.5831E+01  1.0354E+02  1.5786E+01 -4.0968E+01 -2.4409E+01 -1.0203E+02 -7.5509E+01
            -2.9270E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1064.93531252662        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0135E+00  1.0420E+00  9.4371E-01  1.2376E+00  9.0891E-01  8.0970E-01  1.0117E+00  1.0043E+00  1.1316E+00  9.9976E-01
             7.9424E+00
 PARAMETER:  1.1342E-01  1.4116E-01  4.2069E-02  3.1317E-01  4.4915E-03 -1.1109E-01  1.1161E-01  1.0425E-01  2.2365E-01  9.9761E-02
             2.1722E+00
 GRADIENT:  -5.3674E+01 -4.5732E+00 -8.9354E+00  6.8075E+00 -1.2651E+01  1.3048E+01  1.1676E+01  5.5773E+00  2.9778E+01  1.6205E+01
             3.7941E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1218.47101420722        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  8.7222E-01  6.3956E-01  2.1762E-01  1.1531E+00  3.4781E-01  8.2611E-01  6.0498E-01  1.0616E-01  8.6009E-01  3.3497E-02
             4.0556E+00
 PARAMETER: -3.6716E-02 -3.4697E-01 -1.4250E+00  2.4244E-01 -9.5611E-01 -9.1023E-02 -4.0256E-01 -2.1428E+00 -5.0712E-02 -3.2963E+00
             1.5001E+00
 GRADIENT:  -2.3584E+02 -1.7319E+01 -7.4969E+01  1.4352E+02  7.0106E+01 -4.4695E+01  9.3867E-01  1.4017E-01 -2.2107E+01  4.5258E-02
             7.0613E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1255.43397265304        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.4195E-01  1.2174E+00  8.5343E-01  9.3850E-01  9.8653E-01  8.3596E-01  6.6719E-01  1.6245E-01  8.7695E-01  3.8689E-02
             3.6954E+00
 PARAMETER:  4.0195E-02  2.9675E-01 -5.8488E-02  3.6531E-02  8.6437E-02 -7.9173E-02 -3.0468E-01 -1.7174E+00 -3.1303E-02 -3.1522E+00
             1.4071E+00
 GRADIENT:   4.4810E+00  7.8027E+00  4.0774E+00 -1.2228E+01 -1.6663E+01 -7.5966E-01 -3.8651E+00  1.4993E-01 -1.4882E+01  3.0351E-02
            -1.5880E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1259.10013914892        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  9.3722E-01  1.0496E+00  1.1143E+00  1.0694E+00  1.0675E+00  8.3331E-01  4.3926E-01  5.8549E-02  1.1222E+00  8.1949E-02
             3.6914E+00
 PARAMETER:  3.5164E-02  1.4844E-01  2.0818E-01  1.6710E-01  1.6532E-01 -8.2354E-02 -7.2267E-01 -2.7379E+00  2.1526E-01 -2.4017E+00
             1.4060E+00
 GRADIENT:  -6.4417E+00  1.5973E+00  9.5319E-01  4.2945E+00  2.8621E+00 -1.8984E+00  1.2361E+00  2.4632E-02  4.6028E+00  1.1126E-01
             3.4941E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1259.75437215813        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  9.3804E-01  8.5872E-01  1.0593E+00  1.1831E+00  9.5538E-01  8.3868E-01  2.3429E-01  1.8984E-02  1.0151E+00  9.2305E-02
             3.6832E+00
 PARAMETER:  3.6042E-02 -5.2313E-02  1.5765E-01  2.6813E-01  5.4349E-02 -7.5930E-02 -1.3512E+00 -3.8641E+00  1.1503E-01 -2.2827E+00
             1.4038E+00
 GRADIENT:  -1.4345E+00  3.5459E+00 -4.3916E-01  4.2898E+00 -6.6943E-01 -3.4013E-01  7.7953E-02  3.7254E-03 -1.6668E+00  1.4744E-01
             3.8711E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1260.01644963610        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      441
 NPARAMETR:  9.3581E-01  6.1735E-01  1.0581E+00  1.3304E+00  8.5910E-01  8.3876E-01  9.7075E-02  1.0000E-02  9.1769E-01  4.3242E-02
             3.6666E+00
 PARAMETER:  3.3658E-02 -3.8231E-01  1.5648E-01  3.8550E-01 -5.1869E-02 -7.5832E-02 -2.2323E+00 -5.3706E+00  1.4103E-02 -3.0409E+00
             1.3993E+00
 GRADIENT:  -1.2472E+00  3.0421E+00  1.6660E+00  7.2221E+00 -3.5812E+00 -2.1795E-01  8.7225E-03  0.0000E+00 -6.8436E-01  3.3589E-02
            -7.9200E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1260.04619835691        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      517
 NPARAMETR:  9.3545E-01  5.7813E-01  1.0504E+00  1.3507E+00  8.4354E-01  8.3880E-01  8.0964E-02  1.0000E-02  9.0402E-01  2.8959E-02
             3.6642E+00
 PARAMETER:  3.3274E-02 -4.4796E-01  1.4917E-01  4.0064E-01 -7.0154E-02 -7.5787E-02 -2.4138E+00 -5.6369E+00 -9.0110E-04 -3.4419E+00
             1.3986E+00
 GRADIENT:  -3.9051E-01  1.0414E+00  4.2552E-01  3.3044E+00 -1.1286E+00 -1.0295E-01  6.2337E-03  0.0000E+00 -3.1664E-01  1.5414E-02
            -4.8899E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1260.05676869115        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      653
 NPARAMETR:  9.3641E-01  5.5267E-01  1.0558E+00  1.3670E+00  8.3839E-01  8.3957E-01  7.1146E-02  1.0000E-02  8.9310E-01  2.2550E-02
             3.6717E+00
 PARAMETER:  3.4302E-02 -4.9299E-01  1.5428E-01  4.1263E-01 -7.6269E-02 -7.4866E-02 -2.5430E+00 -5.8326E+00 -1.3057E-02 -3.6920E+00
             1.4006E+00
 GRADIENT:   5.4175E-01  3.3023E-01  2.2883E-01  8.2473E-01 -4.3555E-01  7.1063E-02  4.1000E-03  0.0000E+00 -5.1182E-02  9.3318E-03
             1.7350E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1260.06024916252        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      828
 NPARAMETR:  9.3577E-01  5.1655E-01  1.0492E+00  1.3873E+00  8.2367E-01  8.3934E-01  5.8932E-02  1.0000E-02  8.7993E-01  1.4519E-02
             3.6691E+00
 PARAMETER:  3.3609E-02 -5.6059E-01  1.4807E-01  4.2733E-01 -9.3989E-02 -7.5143E-02 -2.7314E+00 -6.1013E+00 -2.7918E-02 -4.1323E+00
             1.3999E+00
 GRADIENT:   9.1913E-02 -8.9732E-02 -3.4399E-02 -1.7672E-01  6.5909E-02  2.1355E-02  2.5268E-03  0.0000E+00  1.5555E-02  3.9141E-03
             4.9504E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1260.06200511064        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1007
 NPARAMETR:  9.3593E-01  5.3352E-01  1.0526E+00  1.3777E+00  8.3079E-01  8.3930E-01  5.7222E-02  1.0000E-02  8.8617E-01  1.0000E-02
             3.6697E+00
 PARAMETER:  3.3788E-02 -5.2826E-01  1.5128E-01  4.2039E-01 -8.5380E-02 -7.5192E-02 -2.7608E+00 -5.9538E+00 -2.0841E-02 -4.6555E+00
             1.4001E+00
 GRADIENT:  -4.7597E-02  4.6873E-02  2.7606E-02  1.6767E-01 -5.1762E-02 -9.9970E-03  2.5173E-03  0.0000E+00 -1.2968E-02  0.0000E+00
            -1.6906E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1260.06296342473        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1184
 NPARAMETR:  9.3598E-01  5.4017E-01  1.0555E+00  1.3738E+00  8.3458E-01  8.3925E-01  2.6991E-02  1.0000E-02  8.8871E-01  1.0000E-02
             3.6699E+00
 PARAMETER:  3.3837E-02 -5.1587E-01  1.5398E-01  4.1757E-01 -8.0826E-02 -7.5252E-02 -3.5122E+00 -5.7970E+00 -1.7989E-02 -9.2491E+00
             1.4002E+00
 GRADIENT:  -7.2252E-02 -5.0806E-02 -6.0860E-02 -4.2490E-02  1.7895E-01 -1.9878E-02  5.8374E-04  0.0000E+00  1.3363E-02  0.0000E+00
            -7.1888E-03

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1260.06327776791        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1346
 NPARAMETR:  9.3599E-01  5.3825E-01  1.0533E+00  1.3748E+00  8.3272E-01  8.3932E-01  1.0000E-02  1.0000E-02  8.8823E-01  1.0000E-02
             3.6699E+00
 PARAMETER:  3.3854E-02 -5.1943E-01  1.5194E-01  4.1829E-01 -8.3062E-02 -7.5165E-02 -4.5769E+00 -5.6641E+00 -1.8528E-02 -1.5541E+01
             1.4002E+00
 GRADIENT:  -1.7377E-03  5.3734E-03  3.4110E-03  1.5116E-02 -6.9529E-03 -2.8329E-04  0.0000E+00  0.0000E+00 -1.5628E-03  0.0000E+00
            -1.8662E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1346
 NO. OF SIG. DIGITS IN FINAL EST.:  3.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6948E-03 -1.8653E-04  1.0202E-04 -1.1957E-02 -1.1044E-05
 SE:             2.8160E-02  8.8814E-05  1.0021E-04  2.3991E-02  1.6901E-04
 N:                     100         100         100         100         100

 P VAL.:         9.5201E-01  3.5705E-02  3.0864E-01  6.1820E-01  9.4790E-01

 ETASHRINKSD(%)  5.6607E+00  9.9702E+01  9.9664E+01  1.9628E+01  9.9434E+01
 ETASHRINKVR(%)  1.1001E+01  9.9999E+01  9.9999E+01  3.5403E+01  9.9997E+01
 EBVSHRINKSD(%)  5.5604E+00  9.9714E+01  9.9577E+01  1.9104E+01  9.9359E+01
 EBVSHRINKVR(%)  1.0812E+01  9.9999E+01  9.9998E+01  3.4559E+01  9.9996E+01
 RELATIVEINF(%)  8.4974E+01  2.5282E-05  1.2029E-04  3.5299E+00  1.4277E-04
 EPSSHRINKSD(%)  1.9580E+01
 EPSSHRINKVR(%)  3.5326E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1260.0632777679095     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -524.91245120417136     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.76
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.81
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1260.063       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.36E-01  5.38E-01  1.05E+00  1.37E+00  8.33E-01  8.39E-01  1.00E-02  1.00E-02  8.88E-01  1.00E-02  3.67E+00
 


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
+        1.61E+03
 
 TH 2
+       -1.01E+02  3.29E+02
 
 TH 3
+        9.32E+00  1.15E+02  1.59E+02
 
 TH 4
+       -1.17E+02  3.41E+02  3.36E+01  4.82E+02
 
 TH 5
+        2.09E+01 -3.08E+02 -3.05E+02 -1.47E+02  6.39E+02
 
 TH 6
+       -3.38E+00 -1.55E+01  7.50E+00 -2.68E+01 -9.11E+00  2.14E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.55E+01 -5.42E+01  6.53E+00 -7.56E+00  2.06E+01  4.20E+00  0.00E+00  0.00E+00  9.82E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.26E+01 -1.53E+01 -2.56E+00 -1.11E+01  7.48E+00  5.76E+00  0.00E+00  0.00E+00  1.32E+01  0.00E+00  3.29E+01
 
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
 #CPUT: Total CPU Time in Seconds,       19.630
Stop Time:
Sat Sep 25 09:04:30 CDT 2021
