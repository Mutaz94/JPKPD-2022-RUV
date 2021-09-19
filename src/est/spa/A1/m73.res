Sat Sep 18 09:24:18 CDT 2021
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
$DATA ../../../../data/spa/A1/dat73.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m73.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1204.15108461815        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.1531E+01 -1.9644E+01  3.6132E+01 -6.6092E+01 -3.1596E+01  4.8840E+01 -6.4687E+00 -1.0263E+01 -8.3647E+00  8.7850E+00
            -9.5259E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1493.19322132570        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0311E+00  1.0965E+00  1.0676E+00  1.0106E+00  1.1097E+00  8.6942E-01  1.0015E+00  9.9503E-01  9.2851E-01  7.8984E-01
             1.8451E+00
 PARAMETER:  1.3066E-01  1.9216E-01  1.6539E-01  1.1057E-01  2.0408E-01 -3.9930E-02  1.0154E-01  9.5016E-02  2.5821E-02 -1.3593E-01
             7.1254E-01
 GRADIENT:   4.2364E+01  9.8315E+00  3.6245E+00  1.7154E+01  1.2825E+01  3.7155E+00 -4.2069E+00 -6.4459E+00 -5.7870E+00  2.8943E+00
            -1.3743E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1505.14374981045        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0101E+00  7.8177E-01  1.0816E+00  1.2386E+00  9.2298E-01  8.6443E-01  1.4979E+00  1.1080E+00  7.3856E-01  4.0985E-01
             2.1238E+00
 PARAMETER:  1.1006E-01 -1.4619E-01  1.7844E-01  3.1395E-01  1.9849E-02 -4.5684E-02  5.0409E-01  2.0259E-01 -2.0306E-01 -7.9197E-01
             8.5321E-01
 GRADIENT:  -2.6735E+01  4.3606E+01  1.5002E+01  8.8092E+01 -4.4268E+01  3.6693E+00  1.3458E-01  6.8325E-01 -6.6003E+00  2.2059E+00
            -4.0063E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1510.86522039644        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0177E+00  1.1169E+00  1.2853E+00  9.9693E-01  1.1999E+00  8.5144E-01  9.8138E-01  1.5011E+00  8.9340E-01  3.3695E-01
             2.3208E+00
 PARAMETER:  1.1758E-01  2.1051E-01  3.5099E-01  9.6928E-02  2.8224E-01 -6.0821E-02  8.1208E-02  5.0619E-01 -1.2717E-02 -9.8781E-01
             9.4189E-01
 GRADIENT:  -9.9865E+00 -2.2617E-01 -2.7455E-01  1.8964E-01  1.2137E+00  7.2217E-01 -5.9883E-01 -5.5081E-01 -3.6829E-01  9.6231E-01
            -5.3920E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1511.31574135161        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  1.0212E+00  1.1483E+00  1.1991E+00  9.7085E-01  1.1758E+00  8.5004E-01  9.9528E-01  1.4774E+00  8.9791E-01  1.3207E-01
             2.3255E+00
 PARAMETER:  1.2094E-01  2.3827E-01  2.8154E-01  7.0419E-02  2.6193E-01 -6.2467E-02  9.5266E-02  4.9031E-01 -7.6810E-03 -1.9244E+00
             9.4395E-01
 GRADIENT:  -6.3688E-01 -1.8254E+00  6.9177E-01 -3.6376E+00 -9.3799E-01  1.4494E-01  1.6518E-01 -9.2156E-02  2.4765E-01  1.4224E-01
             3.1533E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1511.40826046573        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  1.0220E+00  1.1903E+00  1.1015E+00  9.4504E-01  1.1554E+00  8.5042E-01  9.8658E-01  1.4028E+00  9.0925E-01  3.2179E-02
             2.3227E+00
 PARAMETER:  1.2176E-01  2.7419E-01  1.9669E-01  4.3470E-02  2.4447E-01 -6.2022E-02  8.6494E-02  4.3848E-01  4.8694E-03 -3.3364E+00
             9.4272E-01
 GRADIENT:  -4.3993E-01  2.5750E+00  2.1098E-01  1.8845E+00 -1.0386E+00 -5.1714E-03  2.3031E-01  3.6127E-02  6.8280E-02  8.9536E-03
             3.3043E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1511.48218610972        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      528
 NPARAMETR:  1.0274E+00  1.2769E+00  1.1288E+00  8.9590E-01  1.2120E+00  8.5187E-01  9.1562E-01  1.5490E+00  9.5638E-01  1.0000E-02
             2.3344E+00
 PARAMETER:  1.2701E-01  3.4441E-01  2.2115E-01 -9.9242E-03  2.9228E-01 -6.0322E-02  1.1847E-02  5.3759E-01  5.5402E-02 -5.5165E+00
             9.4775E-01
 GRADIENT:   3.8843E+00  4.0899E+00  4.1691E-01  4.7278E+00 -8.3113E-01  3.3117E-02 -1.6908E-01 -4.4923E-02 -3.1713E-01  0.0000E+00
             1.1418E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1511.83963125496        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      706
 NPARAMETR:  1.0257E+00  1.5822E+00  7.9035E-01  6.8904E-01  1.2244E+00  8.5377E-01  7.8565E-01  1.4845E+00  1.1674E+00  1.0000E-02
             2.3225E+00
 PARAMETER:  1.2538E-01  5.5881E-01 -1.3528E-01 -2.7245E-01  3.0244E-01 -5.8097E-02 -1.4124E-01  4.9505E-01  2.5482E-01 -2.3545E+01
             9.4266E-01
 GRADIENT:  -6.0144E+00  8.3019E+00  2.7809E+00  3.3082E+00 -7.2886E+00  1.6436E-01 -1.7273E+00 -3.2356E-01 -2.6634E-01  0.0000E+00
            -4.3347E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1512.59488671356        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      882
 NPARAMETR:  1.0311E+00  1.8604E+00  3.9034E-01  4.9299E-01  1.1823E+00  8.5375E-01  7.3959E-01  8.8900E-01  1.4067E+00  1.0000E-02
             2.2992E+00
 PARAMETER:  1.3058E-01  7.2079E-01 -8.4074E-01 -6.0726E-01  2.6742E-01 -5.8122E-02 -2.0165E-01 -1.7659E-02  4.4127E-01 -5.5252E+01
             9.3255E-01
 GRADIENT:   4.7333E+00  6.7989E+00 -1.5269E+00  6.7570E+00  1.1917E+00 -7.7104E-01  5.7088E-01  2.1632E-01 -9.2163E-02  0.0000E+00
            -6.4241E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1513.46385439698        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1059
 NPARAMETR:  1.0257E+00  2.1202E+00  2.2210E-01  3.1267E-01  1.2788E+00  8.5590E-01  6.8053E-01  5.2879E-01  1.8243E+00  1.0000E-02
             2.3093E+00
 PARAMETER:  1.2537E-01  8.5152E-01 -1.4046E+00 -1.0626E+00  3.4596E-01 -5.5602E-02 -2.8488E-01 -5.3717E-01  7.0122E-01 -1.0224E+02
             9.3692E-01
 GRADIENT:  -8.2455E+00  4.0171E+00 -4.3358E-01  1.7967E+00 -1.6052E-02  4.5837E-01  1.7409E+00  3.2210E-03 -1.4637E+00  0.0000E+00
             1.1192E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1513.98410309700        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1234
 NPARAMETR:  1.0280E+00  2.3332E+00  1.2955E-01  1.7547E-01  1.4222E+00  8.5499E-01  6.3321E-01  2.9067E-01  2.5901E+00  1.0000E-02
             2.3313E+00
 PARAMETER:  1.2761E-01  9.4722E-01 -1.9437E+00 -1.6403E+00  4.5221E-01 -5.6665E-02 -3.5696E-01 -1.1356E+00  1.0517E+00 -1.6362E+02
             9.4641E-01
 GRADIENT:  -1.0221E+00  1.3315E+01  1.0289E+00  7.1812E-01 -2.8757E+00  1.6904E-01 -9.1958E-01 -4.2685E-02  3.2816E-01  0.0000E+00
             5.3411E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1514.93034920203        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1419
 NPARAMETR:  1.0245E+00  2.4721E+00  5.5008E-02  8.7114E-02  1.4797E+00  8.5449E-01  6.1907E-01  1.0208E-01  3.4729E+00  1.0000E-02
             2.3317E+00
 PARAMETER:  1.2425E-01  1.0050E+00 -2.8003E+00 -2.3405E+00  4.9184E-01 -5.7249E-02 -3.7953E-01 -2.1820E+00  1.3450E+00 -2.4033E+02
             9.4661E-01
 GRADIENT:  -8.2294E+00  4.3019E+01  5.2504E+00 -2.0060E+00 -1.8389E+01 -1.9504E-01 -1.7976E+00 -1.6025E-02  4.2620E-01  0.0000E+00
            -2.9294E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1515.05945681079        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1606
 NPARAMETR:  1.0243E+00  2.4805E+00  5.1748E-02  8.3047E-02  1.4848E+00  8.5485E-01  6.1827E-01  9.5159E-02  3.5455E+00  1.0000E-02
             2.3334E+00
 PARAMETER:  1.2405E-01  1.0085E+00 -2.8614E+00 -2.3883E+00  4.9530E-01 -5.6826E-02 -3.8083E-01 -2.2522E+00  1.3657E+00 -2.4574E+02
             9.4735E-01
 GRADIENT:   8.0615E-01  8.0124E+01  5.3434E+00 -1.1674E+00 -1.8565E+01  4.8964E-01 -1.4389E+00 -1.3356E-02  1.1199E+00  0.0000E+00
            -2.2530E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1515.05955642928        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1676
 NPARAMETR:  1.0243E+00  2.4805E+00  5.1748E-02  8.3047E-02  1.4848E+00  8.5347E-01  6.1827E-01  9.7310E-02  3.5455E+00  1.0000E-02
             2.3334E+00
 PARAMETER:  1.2405E-01  1.0085E+00 -2.8614E+00 -2.3883E+00  4.9530E-01 -5.8448E-02 -3.8083E-01 -2.2299E+00  1.3657E+00 -2.4574E+02
             9.4735E-01
 GRADIENT:   7.7703E-01  8.0168E+01  5.3444E+00 -1.1638E+00 -1.8540E+01 -9.0549E-02 -1.4374E+00 -1.3975E-02  1.1186E+00  0.0000E+00
            -2.2843E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1515.06181313131        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:     1751
 NPARAMETR:  1.0243E+00  2.4803E+00  5.1738E-02  8.3050E-02  1.4849E+00  8.4778E-01  6.1827E-01  1.4648E-01  3.5454E+00  1.0000E-02
             2.3335E+00
 PARAMETER:  1.2405E-01  1.0084E+00 -2.8616E+00 -2.3883E+00  4.9532E-01 -6.5132E-02 -3.8083E-01 -1.8208E+00  1.3657E+00 -2.4574E+02
             9.4735E-01
 GRADIENT:   6.7436E-01  8.0020E+01  5.3358E+00 -1.1487E+00 -1.8364E+01 -2.4903E+00 -1.4219E+00 -3.1615E-02  1.1224E+00  0.0000E+00
            -2.3694E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1515.06382428200        NO. OF FUNC. EVALS.:  91
 CUMULATIVE NO. OF FUNC. EVALS.:     1842
 NPARAMETR:  1.0243E+00  2.4802E+00  5.1734E-02  8.3051E-02  1.4849E+00  8.4588E-01  6.1827E-01  1.7671E-01  3.5454E+00  1.0000E-02
             2.3335E+00
 PARAMETER:  1.2405E-01  1.0083E+00 -2.8616E+00 -2.3883E+00  4.9533E-01 -6.7374E-02 -3.8083E-01 -1.6333E+00  1.3657E+00 -2.4574E+02
             9.4736E-01
 GRADIENT:   6.4381E-01  7.9916E+01  5.3281E+00 -1.1422E+00 -1.8294E+01 -3.2988E+00 -1.4141E+00 -4.5879E-02  1.1264E+00  0.0000E+00
            -2.3791E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1516.20019438741        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:     1976
 NPARAMETR:  1.0244E+00  2.4802E+00  5.1772E-02  8.3000E-02  1.4847E+00  8.5053E-01  6.1833E-01  2.0425E+00  3.5442E+00  1.0000E-02
             2.3329E+00
 PARAMETER:  1.2408E-01  1.0083E+00 -2.8609E+00 -2.3889E+00  4.9520E-01 -6.1899E-02 -3.8073E-01  8.1415E-01  1.3653E+00 -2.4574E+02
             9.4712E-01
 GRADIENT:   2.1178E+00  7.4880E+01  1.9943E+00  2.8044E-01 -2.0627E+01 -3.8111E-02 -6.5489E-01  2.7150E-02  2.8463E+00  0.0000E+00
             1.0276E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1516.20062195129        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:     2050
 NPARAMETR:  1.0244E+00  2.4801E+00  5.1771E-02  8.3000E-02  1.4847E+00  8.5078E-01  6.1833E-01  2.0038E+00  3.5442E+00  1.0000E-02
             2.3329E+00
 PARAMETER:  1.2408E-01  1.0083E+00 -2.8609E+00 -2.3889E+00  4.9520E-01 -6.1597E-02 -3.8073E-01  7.9502E-01  1.3653E+00 -2.4574E+02
             9.4712E-01
 GRADIENT:   2.0939E+00  7.4956E+01  2.0943E+00  2.1178E-01 -2.0499E+01  3.8562E-02 -6.6662E-01 -5.9958E-02  2.7569E+00  0.0000E+00
             1.0603E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1516.20424823560        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:     2125
 NPARAMETR:  1.0244E+00  2.4791E+00  5.1763E-02  8.2996E-02  1.4848E+00  8.5161E-01  6.1834E-01  1.8642E+00  3.5437E+00  1.0000E-02
             2.3329E+00
 PARAMETER:  1.2408E-01  1.0079E+00 -2.8611E+00 -2.3890E+00  4.9526E-01 -6.0625E-02 -3.8072E-01  7.2284E-01  1.3652E+00 -2.4574E+02
             9.4709E-01
 GRADIENT:   2.0987E+00  7.3520E+01  2.4408E+00 -6.2243E-02 -1.9767E+01  2.8352E-01 -6.8090E-01 -3.4485E-01  2.4724E+00  0.0000E+00
             1.1628E+00

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1516.22630239800        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:     2196
 NPARAMETR:  1.0244E+00  2.4725E+00  5.1714E-02  8.2974E-02  1.4852E+00  8.5264E-01  6.1836E-01  1.6191E+00  3.5408E+00  1.0000E-02
             2.3326E+00
 PARAMETER:  1.2409E-01  1.0052E+00 -2.8620E+00 -2.3892E+00  4.9557E-01 -5.9423E-02 -3.8068E-01  5.8190E-01  1.3643E+00 -2.4574E+02
             9.4696E-01
 GRADIENT:   2.5725E+00  6.2081E+01  2.9557E+00 -7.0441E-01 -1.6990E+01  5.9647E-01 -5.8123E-01 -7.2396E-01  2.1199E+00  0.0000E+00
             1.3389E+00

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1516.25380093239        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     2266
 NPARAMETR:  1.0244E+00  2.4563E+00  5.1592E-02  8.2917E-02  1.4864E+00  8.5249E-01  6.1842E-01  1.4594E+00  3.5335E+00  1.0000E-02
             2.3318E+00
 PARAMETER:  1.2410E-01  9.9865E-01 -2.8644E+00 -2.3899E+00  4.9636E-01 -5.9595E-02 -3.8058E-01  4.7803E-01  1.3623E+00 -2.4574E+02
             9.4663E-01
 GRADIENT:   4.0925E+00  3.1975E+01  3.0239E+00 -1.6473E+00 -1.1156E+01  6.4898E-01 -3.6417E-01 -8.5867E-01  2.0597E+00  0.0000E+00
             1.7495E+00

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1516.27855318394        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     2336
 NPARAMETR:  1.0244E+00  2.4422E+00  5.1485E-02  8.2868E-02  1.4874E+00  8.5096E-01  6.1848E-01  1.6641E+00  3.5271E+00  1.0000E-02
             2.3311E+00
 PARAMETER:  1.2411E-01  9.9289E-01 -2.8665E+00 -2.3905E+00  4.9704E-01 -6.1393E-02 -3.8049E-01  6.0929E-01  1.3605E+00 -2.4574E+02
             9.4634E-01
 GRADIENT:   5.7174E+00  4.4379E+00  2.1816E+00 -2.0718E+00 -6.7423E+00  4.2324E-01 -2.7305E-01 -5.8080E-01  2.4882E+00  0.0000E+00
             2.3923E+00

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1516.28997846030        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:     2407
 NPARAMETR:  1.0244E+00  2.4390E+00  5.1461E-02  8.2856E-02  1.4877E+00  8.4962E-01  6.1849E-01  1.8984E+00  3.5256E+00  1.0000E-02
             2.3309E+00
 PARAMETER:  1.2411E-01  9.9158E-01 -2.8669E+00 -2.3906E+00  4.9720E-01 -6.2961E-02 -3.8047E-01  7.4099E-01  1.3600E+00 -2.4574E+02
             9.4627E-01
 GRADIENT:   6.2450E+00 -2.5583E+00  1.4869E+00 -1.8871E+00 -6.2374E+00  1.2364E-01 -2.1979E-01 -1.4668E-01  2.9566E+00  0.0000E+00
             2.5070E+00

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1516.29633811571        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     2477
 NPARAMETR:  1.0244E+00  2.4393E+00  5.1463E-02  8.2857E-02  1.4876E+00  8.4929E-01  6.1849E-01  1.9515E+00  3.5257E+00  1.0000E-02
             2.3310E+00
 PARAMETER:  1.2411E-01  9.9169E-01 -2.8669E+00 -2.3906E+00  4.9719E-01 -6.3353E-02 -3.8048E-01  7.6861E-01  1.3601E+00 -2.4574E+02
             9.4628E-01
 GRADIENT:   6.2597E+00 -2.2680E+00  1.3527E+00 -1.7892E+00 -6.4775E+00  2.4300E-02 -1.9912E-01 -3.2988E-02  3.0692E+00  0.0000E+00
             2.4736E+00

0ITERATION NO.:  120    OBJECTIVE VALUE:  -1516.58330275264        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     2633            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0262E+00  2.4581E+00  4.6099E-02  8.5420E-02  1.5113E+00  8.5103E-01  6.2012E-01  1.9986E+00  3.4373E+00  1.0000E-02
             2.3278E+00
 PARAMETER:  1.2589E-01  9.9941E-01 -2.9770E+00 -2.3602E+00  5.1297E-01 -6.1312E-02 -3.7784E-01  7.9246E-01  1.3347E+00 -2.4574E+02
             9.4494E-01
 GRADIENT:   8.3028E+00  2.0496E+01 -6.5357E-01  2.3914E+00  5.7363E+00  1.9891E-01  3.5190E-01  4.7779E-01  2.5162E+00  0.0000E+00
            -1.2142E+00

0ITERATION NO.:  125    OBJECTIVE VALUE:  -1516.68282704690        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     2789
 NPARAMETR:  1.0275E+00  2.4671E+00  4.3432E-02  8.3353E-02  1.5002E+00  8.5227E-01  6.2112E-01  1.7177E+00  3.3819E+00  1.0000E-02
             2.3394E+00
 PARAMETER:  1.2713E-01  1.0030E+00 -3.0366E+00 -2.3847E+00  5.0562E-01 -5.9851E-02 -3.7623E-01  6.4096E-01  1.3184E+00 -2.4574E+02
             9.4989E-01
 GRADIENT:   1.3547E+00  4.5693E+00  1.7316E-01  4.7850E-01 -9.4254E-01  1.4820E-01  2.0513E-01  1.0864E-01 -2.0544E-01  0.0000E+00
             3.8077E-01

0ITERATION NO.:  130    OBJECTIVE VALUE:  -1516.81081033506        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2966
 NPARAMETR:  1.0283E+00  2.4804E+00  2.9260E-02  7.4893E-02  1.4846E+00  8.5220E-01  6.1982E-01  9.0099E-01  3.4599E+00  1.0000E-02
             2.3609E+00
 PARAMETER:  1.2786E-01  1.0084E+00 -3.4315E+00 -2.4917E+00  4.9514E-01 -5.9934E-02 -3.7832E-01 -4.2632E-03  1.3412E+00 -2.4574E+02
             9.5903E-01
 GRADIENT:   4.5782E+00  1.0650E+01 -1.0068E+00  1.6472E+00 -5.2554E+00  5.8494E-01  6.4744E-01  4.2889E-01 -1.0052E+00  0.0000E+00
             1.5670E+00

0ITERATION NO.:  135    OBJECTIVE VALUE:  -1516.81313053774        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     3142
 NPARAMETR:  1.0283E+00  2.4804E+00  2.9277E-02  7.4833E-02  1.4846E+00  8.5073E-01  6.1988E-01  9.0018E-01  3.4591E+00  1.0000E-02
             2.3609E+00
 PARAMETER:  1.2789E-01  1.0084E+00 -3.4309E+00 -2.4925E+00  4.9515E-01 -6.1659E-02 -3.7823E-01 -5.1562E-03  1.3410E+00 -2.4574E+02
             9.5904E-01
 GRADIENT:   4.6836E+00  1.0667E+01 -9.8117E-01  1.6002E+00 -5.3046E+00 -2.7630E-02  6.6631E-01  4.2581E-01 -1.0514E+00  0.0000E+00
             1.5300E+00

0ITERATION NO.:  140    OBJECTIVE VALUE:  -1516.86657081197        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     3269
 NPARAMETR:  1.0283E+00  2.4804E+00  2.9600E-02  7.4833E-02  1.4846E+00  8.5059E-01  6.1926E-01  7.9847E-01  3.4594E+00  1.0000E-02
             2.3609E+00
 PARAMETER:  1.2788E-01  1.0084E+00 -3.4200E+00 -2.4925E+00  4.9515E-01 -6.1820E-02 -3.7923E-01 -1.2505E-01  1.3411E+00 -2.4574E+02
             9.5904E-01
 GRADIENT:   4.5403E+00  1.1225E+01 -5.2844E-01  1.3923E+00 -5.5273E+00 -1.3441E-01  3.4812E-01  3.2039E-01 -1.1887E+00  0.0000E+00
             1.7407E+00

0ITERATION NO.:  145    OBJECTIVE VALUE:  -1516.88749681845        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     3433
 NPARAMETR:  1.0283E+00  2.4804E+00  2.9651E-02  7.4833E-02  1.4846E+00  8.5092E-01  6.1938E-01  7.4639E-01  3.4594E+00  1.0000E-02
             2.3609E+00
 PARAMETER:  1.2788E-01  1.0084E+00 -3.4183E+00 -2.4925E+00  4.9515E-01 -6.1537E-02 -3.7904E-01 -1.9254E-01  1.3411E+00 -2.4574E+02
             9.5904E-01
 GRADIENT:   7.3828E+00  1.1678E+01  1.6936E+01 -8.6464E+00 -1.5576E+00 -4.6633E-02  1.4620E-01 -7.1994E+04 -4.1808E+01  0.0000E+00
            -6.0039E+00

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3433
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -4.5908E-04 -1.9542E-02  7.5110E-04  3.4591E-02 -1.6380E-04
 SE:             2.9187E-02  2.6473E-02  1.3457E-03  1.4471E-02  1.5639E-04
 N:                     100         100         100         100         100

 P VAL.:         9.8745E-01  4.6042E-01  5.7674E-01  1.6830E-02  2.9491E-01

 ETASHRINKSD(%)  2.2196E+00  1.1311E+01  9.5492E+01  5.1520E+01  9.9476E+01
 ETASHRINKVR(%)  4.3900E+00  2.1342E+01  9.9797E+01  7.6497E+01  9.9997E+01
 EBVSHRINKSD(%)  2.5502E+00  1.2145E+01  9.4519E+01  5.5374E+01  9.9428E+01
 EBVSHRINKVR(%)  5.0354E+00  2.2815E+01  9.9700E+01  8.0085E+01  9.9997E+01
 RELATIVEINF(%)  9.2968E+01  2.9794E+01  2.2469E-01  6.3917E+00  1.0637E-03
 EPSSHRINKSD(%)  2.7976E+01
 EPSSHRINKVR(%)  4.8126E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1516.8874968184452     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -781.73667025470706     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    42.45
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.82
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1516.887       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  2.48E+00  2.97E-02  7.48E-02  1.48E+00  8.51E-01  6.19E-01  7.46E-01  3.46E+00  1.00E-02  2.36E+00
 


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
+        4.01E+08
 
 TH 2
+       -5.05E+01  1.11E+06
 
 TH 3
+       -2.69E+01  2.02E+02  6.75E+08
 
 TH 4
+       -8.46E+01  2.68E+02 -3.75E+03  1.99E+08
 
 TH 5
+       -2.11E+01 -9.91E+01 -3.74E+02  3.23E+02  1.28E+07
 
 TH 6
+       -6.76E+00 -9.43E+00 -1.45E+00 -2.26E+01 -4.56E+00  2.46E+02
 
 TH 7
+        7.59E+00 -1.41E+01  1.45E+08  4.26E+00  4.45E+00 -1.97E+00  1.26E+08
 
 TH 8
+       -8.02E-01 -2.76E-01 -1.31E+02  2.70E+01 -9.04E-01 -3.98E+03 -6.62E-01  1.68E+08
 
 TH 9
+       -5.67E+06 -2.67E+00 -1.54E+02  9.41E+01 -2.34E+00  6.49E-01  2.29E+00  4.80E+00  3.22E+05
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.76E+01 -1.09E+01  6.02E+01  2.91E+01 -4.03E+00  5.40E+00  1.95E+01  5.42E-01  2.10E+00  0.00E+00  1.35E+06
 
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
 #CPUT: Total CPU Time in Seconds,       49.331
Stop Time:
Sat Sep 18 09:25:09 CDT 2021
