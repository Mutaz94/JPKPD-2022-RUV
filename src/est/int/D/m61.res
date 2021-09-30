Wed Sep 29 09:21:49 CDT 2021
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
$DATA ../../../../data/int/D/dat61.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m61.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   25806.6092353077        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.2262E+02  5.4358E+02  4.6016E+01  2.0879E+02  3.5037E+02 -2.5051E+03 -9.9964E+02 -1.1157E+02 -1.6416E+03 -8.0665E+02
            -5.2757E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -952.167877813390        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.5858E+00  1.4840E+00  7.5437E-01  2.8985E+00  8.0210E-01  5.5357E+00  4.5844E+00  1.0569E+00  5.2854E+00  2.3330E+00
             1.0938E+01
 PARAMETER:  5.6111E-01  4.9473E-01 -1.8187E-01  1.1642E+00 -1.2052E-01  1.8112E+00  1.6227E+00  1.5533E-01  1.7649E+00  9.4714E-01
             2.4922E+00
 GRADIENT:   3.5535E+01  4.4811E+00 -4.1612E+01  7.1284E+01 -4.9352E+01  2.8085E+02  9.4224E+01  5.8075E+00  1.2852E+02  5.0098E+01
             4.4570E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1017.33943669017        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0807E+00  3.6342E+00  1.6125E+01  1.0650E+01  2.5242E+00  1.8519E+00  1.2688E+00  9.8684E-01  6.5473E+01  2.5890E+00
             1.1432E+01
 PARAMETER:  1.7758E-01  1.3904E+00  2.8804E+00  2.4655E+00  1.0259E+00  7.1623E-01  3.3807E-01  8.6750E-02  4.2816E+00  1.0513E+00
             2.5364E+00
 GRADIENT:  -3.3959E+01  8.0103E+01 -1.3324E+01  1.5712E+01  9.4930E+00  8.1780E-01 -5.8336E+01  1.6108E+00  1.1751E+02  5.7951E+01
             4.4667E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1166.13138133488        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.3636E+00  3.4613E+00  2.2710E+00  1.1193E+00  2.5420E+00  3.6842E+00  2.8566E+00  1.7987E-01  1.0978E+01  7.7992E-01
             1.0370E+01
 PARAMETER:  4.1013E-01  1.3416E+00  9.2024E-01  2.1270E-01  1.0330E+00  1.4041E+00  1.1496E+00 -1.6155E+00  2.4959E+00 -1.4856E-01
             2.4390E+00
 GRADIENT:   1.8501E+01  3.6114E+01 -2.2452E+01  3.3995E+01 -2.6712E+01  1.3086E+02 -2.1903E+00 -3.8479E-02  3.7994E+01  1.0135E+01
             4.7509E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1270.36680310634        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  1.1597E+00  1.9402E+00  1.4933E+00  1.9719E-01  2.1803E+00  3.1221E+00  2.7004E+00  1.0527E-02  5.6984E+00  4.6140E-01
             7.3132E+00
 PARAMETER:  2.4813E-01  7.6281E-01  5.0102E-01 -1.5236E+00  8.7946E-01  1.2385E+00  1.0934E+00 -4.4538E+00  1.8402E+00 -6.7348E-01
             2.0897E+00
 GRADIENT:   3.7533E+00 -7.3160E+01 -6.5704E+00  1.4386E+00 -3.3728E+01  1.8446E+02 -3.0512E+01  1.0745E-04  4.6070E+01  4.6685E+00
             9.7304E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1275.13172296503        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  1.1518E+00  2.1193E+00  1.8027E+00  2.0533E-01  2.2795E+00  2.6094E+00  2.7961E+00  2.0259E-02  5.4061E+00  3.7020E-01
             7.3440E+00
 PARAMETER:  2.4130E-01  8.5108E-01  6.8929E-01 -1.4831E+00  9.2396E-01  1.0591E+00  1.1282E+00 -3.7991E+00  1.7875E+00 -8.9371E-01
             2.0939E+00
 GRADIENT:  -1.1683E+01 -3.5610E+01 -3.1714E+00 -4.2871E+00 -2.9354E+00  4.4835E+01  5.8153E-01  2.9404E-04  2.3387E+01  3.1539E+00
             1.1635E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1286.80902865910        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      444
 NPARAMETR:  1.2041E+00  2.5537E+00  2.6193E+00  4.1193E-01  2.3539E+00  2.4054E+00  2.8125E+00  4.9350E-01  2.6189E+00  3.0431E-01
             7.3791E+00
 PARAMETER:  2.8573E-01  1.0375E+00  1.0629E+00 -7.8689E-01  9.5607E-01  9.7773E-01  1.1341E+00 -6.0624E-01  1.0627E+00 -1.0897E+00
             2.0987E+00
 GRADIENT:   6.0880E-01  4.0411E+01  5.7944E-01  5.3022E+00  6.3651E-01 -2.3570E+01 -9.4059E+00  8.5903E-03  5.4659E+00  1.4200E+00
            -3.0731E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1297.42548133860        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      538
 NPARAMETR:  1.1603E+00  2.1402E+00  2.4707E+00  5.3218E-01  2.2039E+00  2.4991E+00  3.1025E+00  1.6913E+00  1.1226E+00  1.8384E-01
             7.5582E+00
 PARAMETER:  2.4871E-01  8.6088E-01  1.0045E+00 -5.3078E-01  8.9022E-01  1.0159E+00  1.2322E+00  6.2547E-01  2.1563E-01 -1.5937E+00
             2.1226E+00
 GRADIENT:  -2.6843E+01 -2.6681E+01 -1.0785E+00 -3.9640E+00 -1.3963E+01 -7.8135E+01 -8.1689E+01  8.3264E-01  1.0363E+00  4.3154E-01
            -5.0454E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1328.15211732532        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      715
 NPARAMETR:  1.2674E+00  2.2249E+00  3.1713E+00  6.8928E-01  2.2251E+00  3.1198E+00  4.2962E+00  7.5388E-01  6.2043E-01  1.5734E-01
             7.7052E+00
 PARAMETER:  3.3701E-01  8.9972E-01  1.2541E+00 -2.7211E-01  8.9982E-01  1.2378E+00  1.5577E+00 -1.8252E-01 -3.7734E-01 -1.7493E+00
             2.1419E+00
 GRADIENT:  -6.7814E+00  8.6596E+00 -5.0445E+00 -4.5219E+00  6.1522E+00 -9.2725E-01  6.4600E+00  3.4016E-01 -6.0892E-01  2.8684E-01
             4.4227E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1333.30090745945        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      894
 NPARAMETR:  1.3798E+00  1.4833E+00  1.0818E+01  9.3802E-01  2.4057E+00  3.1060E+00  4.8235E+00  1.0815E-02  1.0464E+00  1.7986E-01
             7.8169E+00
 PARAMETER:  4.2196E-01  4.9429E-01  2.4812E+00  3.6017E-02  9.7786E-01  1.2333E+00  1.6735E+00 -4.4268E+00  1.4535E-01 -1.6156E+00
             2.1563E+00
 GRADIENT:   1.1504E+01 -5.8513E+00 -2.4009E-01 -1.9535E+01  3.1992E+00 -5.7447E+00 -1.2641E+01  9.1263E-06  2.8553E-01  3.6091E-01
             4.1271E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1334.09168763969        NO. OF FUNC. EVALS.: 205
 CUMULATIVE NO. OF FUNC. EVALS.:     1099             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3875E+00  1.4583E+00  1.3428E+01  9.5763E-01  2.3921E+00  3.1508E+00  4.7995E+00  1.0000E-02  1.0767E+00  1.7866E-01
             7.6164E+00
 PARAMETER:  4.2751E-01  4.7726E-01  2.6973E+00  5.6710E-02  9.7217E-01  1.2476E+00  1.6685E+00 -4.6782E+00  1.7389E-01 -1.6222E+00
             2.1303E+00
 GRADIENT:   5.5288E+01  8.4751E+00  9.7768E-01 -1.1550E+01 -4.7023E+00  1.4837E+02  1.3489E+02  0.0000E+00 -8.6305E-01  3.1218E-01
             2.5739E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1335.62366134113        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     1255
 NPARAMETR:  1.3057E+00  1.5371E+00  1.1048E+01  1.0199E+00  2.4030E+00  3.2969E+00  4.8893E+00  1.0000E-02  1.1091E+00  1.1152E-01
             7.6608E+00
 PARAMETER:  3.6675E-01  5.2989E-01  2.5022E+00  1.1970E-01  9.7670E-01  1.2930E+00  1.6870E+00 -4.6782E+00  2.0353E-01 -2.0935E+00
             2.1361E+00
 GRADIENT:   3.5432E-01  5.0491E+00 -4.9871E-01  2.7064E+00  2.3997E+00  2.1740E+01 -2.0192E+01  0.0000E+00 -2.7901E+00  1.2100E-01
            -4.1767E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1336.14869496518        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1440
 NPARAMETR:  1.3091E+00  1.3538E+00  1.2519E+01  1.0481E+00  2.3882E+00  3.2719E+00  4.9601E+00  1.0000E-02  1.2467E+00  4.0516E-02
             7.5903E+00
 PARAMETER:  3.6936E-01  4.0288E-01  2.6273E+00  1.4695E-01  9.7056E-01  1.2854E+00  1.7014E+00 -4.6782E+00  3.2052E-01 -3.1061E+00
             2.1269E+00
 GRADIENT:   1.5825E+00 -1.9751E+00  1.3026E-01 -7.8380E+00 -1.1671E+00  1.8725E+01 -2.4061E+01  0.0000E+00 -1.1683E-02  1.5610E-02
            -1.4865E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1336.39299749111        NO. OF FUNC. EVALS.: 171
 CUMULATIVE NO. OF FUNC. EVALS.:     1611
 NPARAMETR:  1.3092E+00  1.3282E+00  1.2822E+01  1.1035E+00  2.3877E+00  3.2685E+00  4.9718E+00  1.0000E-02  1.3559E+00  2.7434E-02
             7.6077E+00
 PARAMETER:  3.6939E-01  3.8389E-01  2.6511E+00  1.9855E-01  9.7050E-01  1.2841E+00  1.7041E+00 -4.6782E+00  4.0438E-01 -3.4614E+00
             2.1290E+00
 GRADIENT:  -1.3583E+03  1.0561E-01 -1.9208E-02  2.5345E+03  5.1482E+02 -1.9131E+02  2.8115E+02  0.0000E+00 -1.2445E+03  7.9517E-03
            -1.1099E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1611
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.2038E-03  2.7253E-02 -8.5961E-06 -6.6303E-02  1.7259E-05
 SE:             2.9737E-02  2.6465E-02  1.1162E-05  1.3336E-02  4.2557E-04
 N:                     100         100         100         100         100

 P VAL.:         8.8758E-01  3.0311E-01  4.4124E-01  6.6436E-07  9.6765E-01

 ETASHRINKSD(%)  3.7771E-01  1.1338E+01  9.9963E+01  5.5323E+01  9.8574E+01
 ETASHRINKVR(%)  7.5398E-01  2.1390E+01  1.0000E+02  8.0040E+01  9.9980E+01
 EBVSHRINKSD(%)  9.7497E-01  1.0415E+01  9.9955E+01  5.9609E+01  9.8471E+01
 EBVSHRINKVR(%)  1.9404E+00  1.9746E+01  1.0000E+02  8.3686E+01  9.9977E+01
 RELATIVEINF(%)  9.7978E+01  3.3313E+01  3.9098E-06  6.7799E+00  4.6405E-03
 EPSSHRINKSD(%)  7.5261E+00
 EPSSHRINKVR(%)  1.4486E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1336.3929974911139     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       317.69636227729688     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    48.57
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    16.99
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1336.393       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.31E+00  1.33E+00  1.28E+01  1.10E+00  2.39E+00  3.27E+00  4.97E+00  1.00E-02  1.36E+00  2.84E-02  7.61E+00
 


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
+        5.36E+04
 
 TH 2
+        1.23E+01  2.06E+01
 
 TH 3
+       -6.92E-02  6.93E-02  2.17E-02
 
 TH 4
+       -4.97E-01 -1.43E+00  4.67E-02  2.62E+05
 
 TH 5
+        1.90E+02 -8.04E+00 -9.06E-01  2.55E+00  2.35E+03
 
 TH 6
+       -3.41E+01  1.38E+00 -6.80E-03  3.22E-01  2.13E+01  6.81E+02
 
 TH 7
+        2.92E+01  2.23E+00 -2.66E-02 -1.18E+01  6.34E+02  3.23E+00  1.79E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.87E-03  4.50E+04 -9.36E-02 -1.05E+05  2.86E+00 -1.04E-01  2.32E+00  0.00E+00  4.19E+04
 
 TH10
+        4.21E+00 -3.31E-02 -5.39E-03 -9.40E+00 -6.13E-01  4.74E-01 -2.40E-01  0.00E+00  3.79E+00  1.01E+01
 
 TH11
+       -1.50E+03 -2.27E+00 -6.73E-03 -7.12E+00  3.11E+02  2.49E+01  8.60E+01  0.00E+00  2.60E+00  8.79E-02  1.82E+01
 
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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       65.705
Stop Time:
Wed Sep 29 09:22:57 CDT 2021
