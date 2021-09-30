Wed Sep 29 09:17:08 CDT 2021
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
$DATA ../../../../data/int/D/dat59.csv ignore=@
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
 (2E4.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m59.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   22440.7012285127        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.8447E+02  3.5163E+02 -2.5606E+01  2.2309E+02  4.2817E+02 -2.8612E+03 -1.1671E+03 -9.6743E+01 -1.9859E+03 -9.9443E+02
            -4.5285E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1173.86111674489        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.5242E+00  1.3156E+00  7.5264E-01  2.8045E+00  7.0509E-01  4.9711E+00  3.8740E+00  1.0823E+00  5.6650E+00  3.1603E+00
             9.7123E+00
 PARAMETER:  5.2149E-01  3.7431E-01 -1.8417E-01  1.1312E+00 -2.4944E-01  1.7037E+00  1.4543E+00  1.7908E-01  1.8343E+00  1.2507E+00
             2.3734E+00
 GRADIENT:   3.8234E+01  1.9692E+01 -3.2785E+01  6.8769E+01 -2.7272E+00  2.0898E+02  8.1672E+01  5.9926E+00  1.3096E+02  4.2807E+01
             4.3594E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1271.95868076265        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      268
 NPARAMETR:  1.8938E+00  1.8584E+00  3.7453E+00  1.1461E+00  2.3591E+00  5.7008E+00  6.9950E+00  7.9148E-01  4.7879E+00  3.2857E+00
             9.3129E+00
 PARAMETER:  7.3858E-01  7.1973E-01  1.4205E+00  2.3633E-01  9.5827E-01  1.8406E+00  2.0452E+00 -1.3386E-01  1.6661E+00  1.2896E+00
             2.3314E+00
 GRADIENT:   1.6844E+01 -5.6101E+00 -3.3212E+01  1.5620E+01  1.5722E+01  1.3371E+02  7.2451E+01  2.6962E-01  5.9447E+01  9.5783E+01
             4.1570E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1450.34657796240        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      447
 NPARAMETR:  1.7669E+00  3.6448E-01  7.8348E+01  1.8438E+00  2.9426E+00  2.5957E+00  8.8683E+00  4.9525E-02  1.6760E+00  1.5727E+00
             8.2475E+00
 PARAMETER:  6.6920E-01 -9.0928E-01  4.4612E+00  7.1181E-01  1.1793E+00  1.0538E+00  2.2825E+00 -2.9053E+00  6.1639E-01  5.5278E-01
             2.2099E+00
 GRADIENT:   8.8537E+01 -3.0851E+00 -2.4114E+00  1.9382E+01  6.9315E+01 -5.5520E+01  1.7692E+01  1.8453E-05 -2.3172E+01  2.1869E+01
             1.6262E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1491.46804925434        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      623
 NPARAMETR:  1.2621E+00  5.2306E-01  5.9234E+01  1.5661E+00  2.2395E+00  2.9019E+00  7.3738E+00  1.0000E-02  2.0603E+00  1.2738E+00
             7.0965E+00
 PARAMETER:  3.3277E-01 -5.4805E-01  4.1815E+00  5.4860E-01  9.0624E-01  1.1654E+00  2.0979E+00 -5.4751E+00  8.2287E-01  3.4199E-01
             2.0596E+00
 GRADIENT:  -1.4810E+01 -5.1470E+00  3.4727E-01 -6.0272E+00  4.5993E-01 -3.5815E+01 -5.4866E+00  0.0000E+00  4.6140E-01  1.5265E+01
             3.8559E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1497.13741533591        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      798
 NPARAMETR:  1.3374E+00  6.5973E-01  3.7430E+01  1.5355E+00  2.2246E+00  3.1949E+00  7.2511E+00  1.0000E-02  2.0636E+00  9.0512E-01
             7.1502E+00
 PARAMETER:  3.9074E-01 -3.1593E-01  3.7225E+00  5.2888E-01  8.9957E-01  1.2615E+00  2.0811E+00 -5.8598E+00  8.2444E-01  3.0916E-04
             2.0671E+00
 GRADIENT:   2.8420E-01 -8.4128E-01  2.0237E-03  8.1078E-02  2.4668E+00  5.1968E+00  6.6829E-01  0.0000E+00  4.0388E-01  2.1420E+00
             2.8781E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1497.35596980671        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      959
 NPARAMETR:  1.3386E+00  6.9870E-01  3.6385E+01  1.5246E+00  2.2171E+00  3.2109E+00  7.4876E+00  1.0000E-02  2.0552E+00  7.8038E-01
             7.1529E+00
 PARAMETER:  3.9166E-01 -2.5853E-01  3.6942E+00  5.2176E-01  8.9622E-01  1.2666E+00  2.1132E+00 -5.8831E+00  8.2036E-01 -1.4798E-01
             2.0675E+00
 GRADIENT:   4.7943E+01  6.2896E+00  3.0589E-01  3.4301E+01  3.2175E+00  1.1449E+02  3.0139E+02  0.0000E+00  1.3350E+01 -2.7061E-01
             3.2232E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1497.44643399972        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1142
 NPARAMETR:  1.3402E+00  6.3911E-01  3.4423E+01  1.5440E+00  2.2176E+00  3.2022E+00  7.4783E+00  1.0000E-02  2.0440E+00  7.9574E-01
             7.1692E+00
 PARAMETER:  3.9282E-01 -3.4768E-01  3.6387E+00  5.3438E-01  8.9643E-01  1.2638E+00  2.1120E+00 -5.8831E+00  8.1490E-01 -1.2848E-01
             2.0698E+00
 GRADIENT:   7.1437E-01 -1.3108E-01 -4.4423E-02 -9.9209E-01  1.5667E+00  6.1064E+00  5.4698E+00  0.0000E+00 -3.9343E-01  7.3223E-02
            -5.5646E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1497.52117736169        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1323
 NPARAMETR:  1.3444E+00  6.0696E-01  3.4543E+01  1.5711E+00  2.2108E+00  3.1852E+00  7.6748E+00  1.0000E-02  2.0642E+00  7.9494E-01
             7.1723E+00
 PARAMETER:  3.9596E-01 -3.9930E-01  3.6422E+00  5.5178E-01  8.9337E-01  1.2585E+00  2.1379E+00 -5.8831E+00  8.2474E-01 -1.2949E-01
             2.0702E+00
 GRADIENT:   1.3589E+00  2.9223E-01 -3.0595E-02 -9.4934E-01  8.5607E-01  3.9805E+00  8.1217E+00  0.0000E+00  3.2598E-02  2.1308E-02
             7.5297E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1497.58788813630        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1506             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3386E+00  5.6534E-01  3.4537E+01  1.5895E+00  2.2053E+00  3.2306E+00  7.8288E+00  1.0000E-02  2.0745E+00  7.9975E-01
             7.1632E+00
 PARAMETER:  3.9159E-01 -4.7033E-01  3.6420E+00  5.6344E-01  8.9085E-01  1.2727E+00  2.1578E+00 -5.8831E+00  8.2972E-01 -1.2345E-01
             2.0690E+00
 GRADIENT:   4.7712E+01  7.5380E+00  8.9129E-02  3.8525E+01  6.1144E+00  1.1739E+02  3.3461E+02  0.0000E+00  1.4318E+01  8.0399E-02
             3.7302E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1497.60774138474        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1686
 NPARAMETR:  1.3397E+00  5.5200E-01  3.5356E+01  1.5999E+00  2.2003E+00  3.2239E+00  7.9324E+00  1.0000E-02  2.0862E+00  7.9589E-01
             7.1637E+00
 PARAMETER:  3.9248E-01 -4.9422E-01  3.6655E+00  5.6995E-01  8.8857E-01  1.2706E+00  2.1710E+00 -5.8831E+00  8.3536E-01 -1.2829E-01
             2.0690E+00
 GRADIENT:   7.5137E-01  1.6080E-01  8.5411E-02 -2.3647E+00 -7.7650E-01  8.6354E+00  1.0372E+01  0.0000E+00  4.8197E-01 -8.7029E-02
             2.9481E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1497.62057683159        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1864
 NPARAMETR:  1.3406E+00  5.4030E-01  3.5001E+01  1.6131E+00  2.2020E+00  3.2117E+00  7.9475E+00  1.0000E-02  2.0866E+00  8.0197E-01
             7.1643E+00
 PARAMETER:  3.9311E-01 -5.1564E-01  3.6554E+00  5.7816E-01  8.8936E-01  1.2668E+00  2.1729E+00 -5.8831E+00  8.3554E-01 -1.2068E-01
             2.0691E+00
 GRADIENT:   8.3770E-01  3.4737E-02 -3.9061E-02 -7.0203E-01  4.2090E-01  7.1924E+00  9.4172E+00  0.0000E+00  7.9619E-02  1.8992E-02
             4.9599E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1497.62776549132        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2045
 NPARAMETR:  1.3405E+00  5.3164E-01  3.5199E+01  1.6193E+00  2.2011E+00  3.2118E+00  7.9937E+00  1.0000E-02  2.0905E+00  8.0324E-01
             7.1632E+00
 PARAMETER:  3.9304E-01 -5.3178E-01  3.6610E+00  5.8200E-01  8.8897E-01  1.2668E+00  2.1787E+00 -5.8831E+00  8.3740E-01 -1.1910E-01
             2.0690E+00
 GRADIENT:   8.3266E-01  5.3598E-02 -3.6372E-02 -6.5733E-01  3.0832E-01  7.1889E+00  9.7984E+00  0.0000E+00  1.4821E-01  2.0726E-02
            -1.0664E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1497.63418931012        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2231             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3399E+00  5.2369E-01  3.5711E+01  1.6255E+00  2.1992E+00  3.2181E+00  8.0537E+00  1.0000E-02  2.0931E+00  8.0089E-01
             7.1616E+00
 PARAMETER:  3.9257E-01 -5.4686E-01  3.6754E+00  5.8582E-01  8.8809E-01  1.2688E+00  2.1861E+00 -5.8831E+00  8.3862E-01 -1.2204E-01
             2.0687E+00
 GRADIENT:   4.8075E+01  8.7528E+00  1.5391E-01  4.2460E+01  4.6278E+00  1.1548E+02  3.5009E+02  0.0000E+00  1.5320E+01 -1.7965E-02
             3.6751E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1497.63606566222        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2412             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3399E+00  5.1758E-01  3.5960E+01  1.6295E+00  2.1989E+00  3.2181E+00  8.0793E+00  1.0000E-02  2.0950E+00  7.9610E-01
             7.1612E+00
 PARAMETER:  3.9262E-01 -5.5859E-01  3.6824E+00  5.8828E-01  8.8797E-01  1.2688E+00  2.1893E+00 -5.8831E+00  8.3955E-01 -1.2804E-01
             2.0687E+00
 GRADIENT:   4.8113E+01  8.8680E+00  1.6817E-01  4.2866E+01  4.3369E+00  1.1549E+02  3.5221E+02  0.0000E+00  1.5436E+01 -1.1014E-01
             3.6457E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1497.63863348209        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2595             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3399E+00  5.1323E-01  3.6010E+01  1.6322E+00  2.1998E+00  3.2181E+00  8.0958E+00  1.0000E-02  2.0963E+00  8.0208E-01
             7.1613E+00
 PARAMETER:  3.9258E-01 -5.6704E-01  3.6838E+00  5.8992E-01  8.8837E-01  1.2688E+00  2.1914E+00 -5.8831E+00  8.4016E-01 -1.2055E-01
             2.0687E+00
 GRADIENT:   4.8078E+01  8.9230E+00  1.3732E-01  4.3048E+01  4.8169E+00  1.1546E+02  3.5374E+02  0.0000E+00  1.5542E+01 -2.7228E-03
             3.6789E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1497.63968397648        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2774
 NPARAMETR:  1.3399E+00  5.0974E-01  3.6232E+01  1.6348E+00  2.1997E+00  3.2182E+00  8.1142E+00  1.0000E-02  2.0972E+00  8.0201E-01
             7.1610E+00
 PARAMETER:  3.9256E-01 -5.7385E-01  3.6899E+00  5.9154E-01  8.8832E-01  1.2688E+00  2.1936E+00 -5.8831E+00  8.4060E-01 -1.2064E-01
             2.0686E+00
 GRADIENT:   7.6125E-01  9.9942E-02  1.8585E-02 -5.2073E-01 -5.0520E-01  7.9432E+00  1.0723E+01  0.0000E+00  1.3165E-01 -5.8268E-02
            -5.4884E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1497.64052934459        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2955
 NPARAMETR:  1.3399E+00  5.0609E-01  3.6449E+01  1.6368E+00  2.2000E+00  3.2181E+00  8.1272E+00  1.0000E-02  2.0986E+00  8.0156E-01
             7.1614E+00
 PARAMETER:  3.9261E-01 -5.8104E-01  3.6959E+00  5.9277E-01  8.8848E-01  1.2688E+00  2.1952E+00 -5.8831E+00  8.4128E-01 -1.2119E-01
             2.0687E+00
 GRADIENT:   7.6934E-01  5.6738E-02  2.3536E-02 -5.5605E-01 -5.4915E-01  7.9395E+00  1.0717E+01  0.0000E+00  1.6986E-01 -6.5585E-02
            -4.3000E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1497.64142539066        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     3141             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3398E+00  5.0276E-01  3.6456E+01  1.6387E+00  2.2013E+00  3.2182E+00  8.1411E+00  1.0000E-02  2.0998E+00  8.0587E-01
             7.1616E+00
 PARAMETER:  3.9252E-01 -5.8765E-01  3.6961E+00  5.9390E-01  8.8906E-01  1.2688E+00  2.1969E+00 -5.8831E+00  8.4186E-01 -1.1583E-01
             2.0687E+00
 GRADIENT:   4.8041E+01  9.0925E+00  1.1631E-01  4.3503E+01  5.1722E+00  1.1544E+02  3.5762E+02  0.0000E+00  1.5814E+01  6.6393E-02
             3.7113E+01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1497.64186605493        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     3326
 NPARAMETR:  1.3398E+00  5.0058E-01  3.6615E+01  1.6405E+00  2.2016E+00  3.2182E+00  8.1527E+00  1.0000E-02  2.1006E+00  8.0545E-01
             7.1614E+00
 PARAMETER:  3.9252E-01 -5.9200E-01  3.7005E+00  5.9498E-01  8.8918E-01  1.2688E+00  2.1984E+00 -5.8831E+00  8.4224E-01 -1.1635E-01
             2.0687E+00
 GRADIENT:   7.4230E-01  1.7752E-02 -1.0527E-02 -5.7530E-01 -2.6952E-02  7.9460E+00  1.0866E+01  0.0000E+00  2.3519E-01  1.2824E-02
            -2.1214E-01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1497.64222985903        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     3514             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3398E+00  4.9787E-01  3.6837E+01  1.6425E+00  2.2020E+00  3.2182E+00  8.1662E+00  1.0000E-02  2.1014E+00  8.0398E-01
             7.1613E+00
 PARAMETER:  3.9253E-01 -5.9741E-01  3.7065E+00  5.9621E-01  8.8936E-01  1.2688E+00  2.2000E+00 -5.8831E+00  8.4262E-01 -1.1818E-01
             2.0687E+00
 GRADIENT:   4.8052E+01  9.2048E+00  1.2478E-01  4.3897E+01  5.0414E+00  1.1544E+02  3.5952E+02  0.0000E+00  1.5908E+01  2.7020E-02
             3.6937E+01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1497.64233740276        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     3695
 NPARAMETR:  1.3398E+00  4.9685E-01  3.6946E+01  1.6433E+00  2.2021E+00  3.2182E+00  8.1714E+00  1.0000E-02  2.1017E+00  8.0417E-01
             7.1612E+00
 PARAMETER:  3.9254E-01 -5.9946E-01  3.7095E+00  5.9672E-01  8.8941E-01  1.2688E+00  2.2006E+00 -5.8831E+00  8.4274E-01 -1.1795E-01
             2.0687E+00
 GRADIENT:   7.4987E-01  1.3970E-02  1.2821E-03 -4.9553E-01 -1.7056E-01  7.9429E+00  1.0950E+01  0.0000E+00  2.1597E-01 -1.5838E-02
            -3.3272E-01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1497.64240236464        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     3881             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3398E+00  4.9538E-01  3.6978E+01  1.6441E+00  2.2028E+00  3.2182E+00  8.1782E+00  1.0000E-02  2.1024E+00  8.0607E-01
             7.1616E+00
 PARAMETER:  3.9251E-01 -6.0244E-01  3.7103E+00  5.9721E-01  8.8974E-01  1.2688E+00  2.2015E+00 -5.8831E+00  8.4310E-01 -1.1559E-01
             2.0687E+00
 GRADIENT:   4.8035E+01  9.2463E+00  1.1286E-01  4.3990E+01  5.2580E+00  1.1542E+02  3.6050E+02  0.0000E+00  1.5989E+01  6.9312E-02
             3.7123E+01

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1497.64248309749        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     4062
 NPARAMETR:  1.3398E+00  4.9492E-01  3.7066E+01  1.6446E+00  2.2027E+00  3.2182E+00  8.1807E+00  1.0000E-02  2.1025E+00  8.0483E-01
             7.1614E+00
 PARAMETER:  3.9253E-01 -6.0336E-01  3.7127E+00  5.9751E-01  8.8969E-01  1.2688E+00  2.2018E+00 -5.8831E+00  8.4311E-01 -1.1713E-01
             2.0687E+00
 GRADIENT:   7.4557E-01  2.1920E-03 -5.2096E-03 -5.1047E-01 -4.9576E-02  7.9441E+00  1.1007E+01  0.0000E+00  2.4015E-01  1.2373E-04
            -2.3504E-01

0ITERATION NO.:  118    OBJECTIVE VALUE:  -1497.64252874381        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:     4163
 NPARAMETR:  1.3397E+00  4.9071E-01  3.6993E+01  1.6468E+00  2.2053E+00  3.2182E+00  8.1979E+00  1.0000E-02  2.1049E+00  8.0594E-01
             7.1616E+00
 PARAMETER:  3.9255E-01 -6.0940E-01  3.7224E+00  5.9900E-01  8.8984E-01  1.2688E+00  2.2037E+00 -5.8831E+00  8.4347E-01 -1.1691E-01
             2.0687E+00
 GRADIENT:   5.3502E-03  1.5496E-02  1.4070E-02  1.9308E-02 -2.0214E-01  5.9939E-04 -1.2039E-02  0.0000E+00 -3.4438E-02 -5.5522E-03
            -1.6623E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     4163
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.4739E-03  6.2356E-02 -2.6520E-06 -8.9425E-02 -5.0082E-03
 SE:             2.9468E-02  2.1617E-02  6.7240E-06  1.7535E-02  1.3197E-02
 N:                     100         100         100         100         100

 P VAL.:         8.7932E-01  3.9188E-03  6.9328E-01  3.4073E-07  7.0431E-01

 ETASHRINKSD(%)  1.2799E+00  2.7581E+01  9.9977E+01  4.1255E+01  5.5790E+01
 ETASHRINKVR(%)  2.5435E+00  4.7555E+01  1.0000E+02  6.5490E+01  8.0455E+01
 EBVSHRINKSD(%)  1.3576E+00  2.7914E+01  9.9973E+01  3.2651E+01  5.6295E+01
 EBVSHRINKVR(%)  2.6968E+00  4.8036E+01  1.0000E+02  5.4641E+01  8.0898E+01
 RELATIVEINF(%)  9.7227E+01  2.7696E+01  1.9081E-06  2.4313E+01  4.9703E+00
 EPSSHRINKSD(%)  9.4618E+00
 EPSSHRINKVR(%)  1.8028E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1497.6425287438133     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       156.44683102459749     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:   161.62
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    18.03
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1497.643       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.34E+00  4.92E-01  3.74E+01  1.65E+00  2.20E+00  3.22E+00  8.20E+00  1.00E-02  2.10E+00  8.05E-01  7.16E+00
 


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
+        1.28E+01
 
 TH 2
+        9.96E+00  7.74E+00
 
 TH 3
+        1.46E-02  1.14E-02  1.67E-05
 
 TH 4
+        1.17E+01  9.06E+00  1.32E-02  1.06E+01
 
 TH 5
+       -6.70E+00 -5.21E+00 -7.66E-03 -6.09E+00  3.52E+00
 
 TH 6
+        1.92E+00  1.50E+00  2.21E-03  1.73E+00 -1.01E+00  2.96E-01
 
 TH 7
+        3.94E-01  3.07E-01  4.54E-04  3.56E-01 -2.08E-01  6.09E-02  1.25E-02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.57E+00 -1.22E+00 -1.77E-03 -1.45E+00  8.14E-01 -2.27E-01 -4.64E-02  0.00E+00  2.04E-01
 
 TH10
+       -1.10E+00 -8.48E-01 -1.22E-03 -1.01E+00  5.64E-01 -1.55E-01 -3.16E-02  0.00E+00  1.46E-01  1.05E-01
 
 TH11
+       -7.82E-01 -5.86E-01 -7.29E-04 -8.11E-01  3.51E-01 -6.00E-02 -1.12E-02  0.00E+00  1.64E-01  1.37E-01  4.64E-01
 
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
+        5.75E+01
 
 TH 2
+       -1.64E-01  3.89E+01
 
 TH 3
+        3.55E-03  1.42E-02  1.32E-03
 
 TH 4
+       -6.34E-01  1.72E+01 -4.81E-02  6.35E+01
 
 TH 5
+       -8.89E-01 -6.23E+00 -2.28E-01 -5.01E+00  6.09E+01
 
 TH 6
+        3.43E-01 -5.63E-02  5.38E-04  1.82E-01 -1.81E-01  1.78E+01
 
 TH 7
+        2.55E-02  4.74E+00  6.17E-04 -3.63E+00  6.30E-01 -2.84E-02  1.54E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.22E-01 -1.03E+00 -1.69E-02 -8.54E+00  1.97E+00 -2.74E-01  3.27E-01  0.00E+00  1.49E+01
 
 TH10
+       -3.78E-01 -5.25E-01 -1.28E-02 -2.62E+00  6.63E+00  6.36E-02  1.33E-01  0.00E+00  5.65E-01  1.05E+01
 
 TH11
+       -2.68E+00 -1.85E+00  2.62E-03 -4.44E+00 -1.29E+00  1.04E+00  1.37E-01  0.00E+00  2.20E+00  3.57E+00  1.96E+01
 
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
+        5.82E+01
 
 TH 2
+        4.07E+01  4.16E+01
 
 TH 3
+        2.01E-02  2.85E-02  8.98E-04
 
 TH 4
+        2.85E+01  1.62E+01  5.09E-03  5.85E+01
 
 TH 5
+       -8.84E+00 -2.79E+00 -1.92E-01 -1.24E+01  6.81E+01
 
 TH 6
+        2.03E+01  1.86E+01 -1.70E-03 -9.31E+00  4.49E+00  2.99E+01
 
 TH 7
+        3.07E+00  4.66E+00 -1.95E-03 -3.21E+00  2.66E+00  4.38E+00  1.44E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -7.72E+00 -2.98E+00 -1.68E-02 -1.48E+01  4.99E+00  7.42E+00  7.83E-01  0.00E+00  1.96E+01
 
 TH10
+       -9.43E-01  5.04E-01 -1.54E-02 -5.70E+00  7.59E+00 -8.08E-01  8.87E-01  0.00E+00  1.41E-01  1.06E+01
 
 TH11
+       -7.96E+01 -2.82E+01 -2.03E-02 -8.33E+01  2.38E+01  3.02E+01  3.52E+00  0.00E+00  6.62E+01 -9.57E+00  9.16E+02
 
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
 #CPUT: Total CPU Time in Seconds,      179.756
Stop Time:
Wed Sep 29 09:20:09 CDT 2021
