Sat Sep 18 10:36:48 CDT 2021
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
$DATA ../../../../data/spa/A3/dat65.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m65.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   157.715532635087        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.3907E+01  6.0434E+01  5.0716E+01  5.4726E+01  3.2108E+02  5.3296E+01 -8.7387E+01 -7.2202E+01 -1.4218E+02 -2.5633E+02
            -3.1458E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1258.14531976997        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0730E+00  1.0208E+00  9.4269E-01  1.0943E+00  9.3146E-01  7.6438E-01  9.5450E-01  1.0276E+00  9.2947E-01  9.8889E-01
             3.9864E+00
 PARAMETER:  1.7042E-01  1.2063E-01  4.0987E-02  1.9009E-01  2.8993E-02 -1.6868E-01  5.3434E-02  1.2726E-01  2.6864E-02  8.8824E-02
             1.4829E+00
 GRADIENT:   5.7097E+01  5.0465E+00 -2.0913E+01  3.9311E+01  2.5102E+01 -1.8622E+01  9.2721E+00  1.4715E+00  1.4698E+01  1.6118E+01
            -5.7405E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1268.77480814140        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      162
 NPARAMETR:  1.0766E+00  1.1787E+00  1.0666E+00  1.0030E+00  1.0199E+00  7.7508E-01  5.3050E-01  4.2981E-01  9.9490E-01  7.5194E-01
             4.4753E+00
 PARAMETER:  1.7385E-01  2.6439E-01  1.6446E-01  1.0302E-01  1.1966E-01 -1.5479E-01 -5.3394E-01 -7.4441E-01  9.4887E-02 -1.8509E-01
             1.5986E+00
 GRADIENT:   4.6604E+01  1.6863E+01 -2.9941E+00  3.4539E+01 -2.8825E+00 -9.7106E+00  2.0195E+00  4.9583E-01  1.2824E+01  7.5118E+00
             2.0658E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1274.29566356095        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  1.0528E+00  1.3292E+00  1.2835E+00  8.8050E-01  1.1329E+00  8.0342E-01  6.0473E-01  6.6950E-01  7.1609E-01  3.8655E-01
             4.4415E+00
 PARAMETER:  1.5150E-01  3.8458E-01  3.4958E-01 -2.7269E-02  2.2482E-01 -1.1888E-01 -4.0298E-01 -3.0122E-01 -2.3395E-01 -8.5050E-01
             1.5910E+00
 GRADIENT:  -5.7616E+00  2.0265E+01  5.0594E-01  1.8224E+01 -8.7028E+00 -1.6616E-01 -9.9625E-02  3.0838E-01 -8.2386E-01  8.8377E-01
             3.8164E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1275.23335742793        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  1.0516E+00  1.5214E+00  1.6903E+00  7.5464E-01  1.3036E+00  8.0124E-01  5.3685E-01  2.2673E-01  7.6481E-01  1.8987E-01
             4.4288E+00
 PARAMETER:  1.5031E-01  5.1965E-01  6.2492E-01 -1.8152E-01  3.6510E-01 -1.2159E-01 -5.2204E-01 -1.3840E+00 -1.6813E-01 -1.5614E+00
             1.5881E+00
 GRADIENT:  -4.2922E+00  7.8426E+00 -6.9477E-01  6.1605E+00  1.8524E+00 -5.7760E-01 -4.3349E-01  1.0888E-02 -4.9146E-01  1.3088E-01
            -2.3123E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1275.37753731518        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  1.0517E+00  1.5772E+00  1.8002E+00  7.0976E-01  1.3135E+00  8.0052E-01  5.3022E-01  1.0729E-01  8.0518E-01  1.2492E-01
             4.4335E+00
 PARAMETER:  1.5042E-01  5.5563E-01  6.8791E-01 -2.4282E-01  3.7269E-01 -1.2250E-01 -5.3446E-01 -2.1322E+00 -1.1669E-01 -1.9801E+00
             1.5892E+00
 GRADIENT:  -3.4283E-01  5.2767E-02 -1.4991E-01 -1.0916E-01  4.8014E-01 -2.0373E-02  6.2052E-02  1.6314E-03  2.1095E-01  5.7151E-02
             3.4285E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1275.50339289667        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      513
 NPARAMETR:  1.0550E+00  1.6962E+00  2.1635E+00  6.3638E-01  1.3703E+00  8.0249E-01  5.2768E-01  3.0891E-02  7.8843E-01  5.4793E-02
             4.4508E+00
 PARAMETER:  1.5354E-01  6.2838E-01  8.7171E-01 -3.5197E-01  4.1506E-01 -1.2004E-01 -5.3927E-01 -3.3773E+00 -1.3771E-01 -2.8042E+00
             1.5931E+00
 GRADIENT:   1.9442E+00 -1.0274E+00 -2.3576E-01  3.7615E-02  8.7369E-01  1.7549E-01  1.4246E-01  3.8659E-05 -1.9276E-02  1.0685E-02
             8.2910E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1275.68568548350        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      690
 NPARAMETR:  1.0557E+00  1.9694E+00  3.7547E+00  4.5981E-01  1.4365E+00  8.0381E-01  4.9850E-01  1.0000E-02  8.4486E-01  1.0000E-02
             4.4482E+00
 PARAMETER:  1.5418E-01  7.7774E-01  1.4230E+00 -6.7694E-01  4.6222E-01 -1.1840E-01 -5.9615E-01 -7.9623E+00 -6.8590E-02 -5.5718E+00
             1.5925E+00
 GRADIENT:  -2.2025E+00  7.5221E+00 -1.1590E-01  2.6793E+00 -5.5590E-01 -1.1867E-02 -5.9479E-02  0.0000E+00 -1.7851E-01  0.0000E+00
            -1.4583E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1275.76468392990        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      865
 NPARAMETR:  1.0562E+00  2.1122E+00  6.1602E+00  3.5947E-01  1.4447E+00  8.0374E-01  4.7294E-01  1.0000E-02  1.0002E+00  1.0000E-02
             4.4495E+00
 PARAMETER:  1.5465E-01  8.4771E-01  1.9181E+00 -9.2312E-01  4.6791E-01 -1.1847E-01 -6.4878E-01 -1.2256E+01  1.0021E-01 -7.9620E+00
             1.5928E+00
 GRADIENT:  -6.3682E-01  3.7507E+00 -3.9236E-02  1.0297E+00 -8.8037E-01  4.8224E-02 -1.8697E-01  0.0000E+00  2.3013E-02  0.0000E+00
            -8.7105E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1275.79353127307        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1040
 NPARAMETR:  1.0565E+00  2.2278E+00  1.1158E+01  2.7942E-01  1.4489E+00  8.0385E-01  4.6295E-01  1.0000E-02  1.1253E+00  1.0000E-02
             4.4529E+00
 PARAMETER:  1.5498E-01  9.0101E-01  2.5122E+00 -1.1751E+00  4.7082E-01 -1.1834E-01 -6.7014E-01 -1.7004E+01  2.1809E-01 -1.0638E+01
             1.5936E+00
 GRADIENT:  -2.6145E-01  3.4489E+00 -1.3535E-02  6.7077E-01 -1.0002E+00  1.3050E-01 -1.2060E-03  0.0000E+00  4.4744E-02  0.0000E+00
             1.7343E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1275.80357797482        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1215
 NPARAMETR:  1.0565E+00  2.3052E+00  1.9902E+01  2.2483E-01  1.4522E+00  8.0364E-01  4.5801E-01  1.0000E-02  1.1676E+00  1.0000E-02
             4.4527E+00
 PARAMETER:  1.5493E-01  9.3518E-01  3.0908E+00 -1.3924E+00  4.7311E-01 -1.1861E-01 -6.8087E-01 -2.1302E+01  2.5494E-01 -1.3109E+01
             1.5935E+00
 GRADIENT:  -2.2096E-01  1.5867E+00 -5.4111E-03  2.0350E-01 -4.5185E-01  1.1396E-02 -3.7340E-02  0.0000E+00 -6.9922E-02  0.0000E+00
            -9.5204E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1275.80863144143        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1390
 NPARAMETR:  1.0564E+00  2.3626E+00  3.4034E+01  1.8480E-01  1.4555E+00  8.0324E-01  4.4997E-01  1.0000E-02  1.3589E+00  1.0000E-02
             4.4517E+00
 PARAMETER:  1.5491E-01  9.5974E-01  3.6273E+00 -1.5885E+00  4.7536E-01 -1.1911E-01 -6.9858E-01 -2.5359E+01  4.0667E-01 -1.5331E+01
             1.5933E+00
 GRADIENT:  -2.8831E-01  4.4041E-01 -2.1774E-03  6.0571E-02  1.8153E-01 -1.1917E-01 -6.3698E-03  0.0000E+00  9.3776E-03  0.0000E+00
            -1.5592E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1275.80979589354        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1569
 NPARAMETR:  1.0565E+00  2.3922E+00  4.8206E+01  1.6365E-01  1.4541E+00  8.0349E-01  4.4755E-01  1.0000E-02  1.4243E+00  1.0000E-02
             4.4523E+00
 PARAMETER:  1.5495E-01  9.7221E-01  3.9755E+00 -1.7100E+00  4.7437E-01 -1.1879E-01 -7.0397E-01 -2.7913E+01  4.5366E-01 -1.6767E+01
             1.5934E+00
 GRADIENT:  -6.7014E-03  5.5586E-02 -1.1828E-03  6.6451E-03 -1.5805E-02 -8.3754E-03 -1.1770E-02  0.0000E+00 -7.1017E-03  0.0000E+00
            -2.4054E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1275.81041631133        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1745
 NPARAMETR:  1.0565E+00  2.4319E+00  8.2551E+01  1.3621E-01  1.4542E+00  8.0348E-01  4.4373E-01  1.0000E-02  1.5771E+00  1.0000E-02
             4.4523E+00
 PARAMETER:  1.5495E-01  9.8869E-01  4.5134E+00 -1.8935E+00  4.7447E-01 -1.1880E-01 -7.1254E-01 -3.1835E+01  5.5560E-01 -1.8944E+01
             1.5934E+00
 GRADIENT:  -1.8032E-01  8.3714E-01 -4.7812E-04  7.1771E-02 -1.2922E-01 -2.0257E-02 -7.1184E-03  0.0000E+00  3.8482E-03  0.0000E+00
            -4.8134E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1275.81075515435        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1923
 NPARAMETR:  1.0565E+00  2.4577E+00  1.2815E+02  1.1787E-01  1.4543E+00  8.0352E-01  4.4144E-01  1.0000E-02  1.6793E+00  1.0000E-02
             4.4523E+00
 PARAMETER:  1.5495E-01  9.9923E-01  4.9532E+00 -2.0382E+00  4.7456E-01 -1.1875E-01 -7.1771E-01 -3.4971E+01  6.1837E-01 -2.0698E+01
             1.5934E+00
 GRADIENT:  -2.9404E-02  1.0678E-01 -2.3189E-04  7.6619E-03 -9.6757E-03  1.6391E-03 -1.4141E-02  0.0000E+00 -4.2185E-03  0.0000E+00
            -2.0944E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1275.81085521498        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2101
 NPARAMETR:  1.0565E+00  2.4807E+00  1.9995E+02  1.0197E-01  1.4543E+00  8.0349E-01  4.3939E-01  1.0000E-02  1.8175E+00  1.0000E-02
             4.4523E+00
 PARAMETER:  1.5496E-01  1.0085E+00  5.3980E+00 -2.1831E+00  4.7451E-01 -1.1879E-01 -7.2236E-01 -3.8147E+01  6.9745E-01 -2.2460E+01
             1.5934E+00
 GRADIENT:  -1.0274E-01  5.1291E-01 -1.1084E-04  3.2393E-02 -8.0265E-02 -1.4053E-02 -6.0323E-03  0.0000E+00  1.0993E-03  0.0000E+00
            -3.4583E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1275.81090680261        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2284
 NPARAMETR:  1.0565E+00  2.4904E+00  2.4856E+02  9.5096E-02  1.4543E+00  8.0350E-01  4.3849E-01  1.0000E-02  1.8877E+00  1.0000E-02
             4.4523E+00
 PARAMETER:  1.5495E-01  1.0124E+00  5.6157E+00 -2.2529E+00  4.7454E-01 -1.1878E-01 -7.2441E-01 -3.9690E+01  7.3534E-01 -2.3316E+01
             1.5934E+00
 GRADIENT:  -6.7536E-02  2.7064E-01 -7.7422E-05  1.7845E-02 -3.6191E-02 -3.7440E-03  4.5906E-04  0.0000E+00  3.3308E-03  0.0000E+00
            -1.0149E-02

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1275.81092446447        NO. OF FUNC. EVALS.: 173
 CUMULATIVE NO. OF FUNC. EVALS.:     2457
 NPARAMETR:  1.0565E+00  2.4903E+00  2.5324E+02  9.5042E-02  1.4543E+00  8.0354E-01  4.3845E-01  1.0000E-02  1.8845E+00  1.0000E-02
             4.4523E+00
 PARAMETER:  1.5495E-01  1.0124E+00  5.6344E+00 -2.2534E+00  4.7451E-01 -1.1873E-01 -7.2451E-01 -3.9690E+01  7.3366E-01 -2.3316E+01
             1.5934E+00
 GRADIENT:  -2.1336E-02  5.0562E-02 -7.5875E-05  5.2317E-03 -1.1089E-02  1.3203E-02 -5.6467E-03  0.0000E+00  1.8185E-03  0.0000E+00
            -1.5224E-02

0ITERATION NO.:   86    OBJECTIVE VALUE:  -1275.81092446447        NO. OF FUNC. EVALS.:  31
 CUMULATIVE NO. OF FUNC. EVALS.:     2488
 NPARAMETR:  1.0565E+00  2.4901E+00  2.5300E+02  9.4933E-02  1.4545E+00  8.0355E-01  4.3877E-01  1.0000E-02  1.8850E+00  1.0000E-02
             4.4525E+00
 PARAMETER:  1.5495E-01  1.0124E+00  5.6344E+00 -2.2534E+00  4.7451E-01 -1.1873E-01 -7.2451E-01 -3.9690E+01  7.3366E-01 -2.3316E+01
             1.5934E+00
 GRADIENT:  -1.0883E-02  6.4338E-02  8.5715E-05  4.1421E-03 -1.4696E-02 -3.5557E-03 -1.1094E-02  0.0000E+00 -8.3552E-04  0.0000E+00
            -1.5192E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2488
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.5709E-03 -1.2255E-02  1.5113E-09 -2.8484E-03 -2.7267E-05
 SE:             2.7996E-02  1.8387E-02  8.8036E-10  4.3077E-03  1.0805E-04
 N:                     100         100         100         100         100

 P VAL.:         9.2683E-01  5.0509E-01  8.6042E-02  5.0845E-01  8.0076E-01

 ETASHRINKSD(%)  6.2085E+00  3.8400E+01  1.0000E+02  8.5569E+01  9.9638E+01
 ETASHRINKVR(%)  1.2032E+01  6.2054E+01  1.0000E+02  9.7917E+01  9.9999E+01
 EBVSHRINKSD(%)  6.2685E+00  3.8570E+01  1.0000E+02  8.5604E+01  9.9613E+01
 EBVSHRINKVR(%)  1.2144E+01  6.2264E+01  1.0000E+02  9.7928E+01  9.9999E+01
 RELATIVEINF(%)  7.2582E+01  1.0964E-05  0.0000E+00  6.0230E-07  3.0999E-04
 EPSSHRINKSD(%)  1.5243E+01
 EPSSHRINKVR(%)  2.8163E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1275.8109244644722     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -540.66009790073406     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    31.36
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     6.38
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1275.811       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.06E+00  2.49E+00  2.53E+02  9.50E-02  1.45E+00  8.04E-01  4.38E-01  1.00E-02  1.88E+00  1.00E-02  4.45E+00
 


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
+        1.54E+03
 
 TH 2
+        3.00E+01  3.26E+02
 
 TH 3
+       -4.45E-03 -3.88E-04  1.79E-08
 
 TH 4
+        5.52E+01  4.84E+02 -5.18E-04  7.21E+02
 
 TH 5
+       -8.22E+01 -9.90E+01  2.33E-04 -1.49E+02  3.55E+01
 
 TH 6
+        7.17E+01  1.38E+01  1.35E-03  5.06E+01 -3.94E+01  5.27E+02
 
 TH 7
+       -2.45E+01 -4.18E+01  2.43E-04 -5.97E+01  1.10E+01  4.19E+01  9.44E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.68E+00  3.04E-01  8.17E-06  7.14E-01 -4.41E-01  4.52E+00  3.16E-01  0.00E+00  3.95E-02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -3.01E+01 -1.77E+01  1.38E-04 -2.58E+01  6.01E+00  9.68E+00  3.65E+00  0.00E+00  5.21E-02  0.00E+00  1.75E+00
 
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
+        1.39E+03
 
 TH 2
+       -1.64E+02  3.38E+02
 
 TH 3
+       -3.98E-03 -4.21E-05  7.35E-07
 
 TH 4
+       -2.38E+02  4.83E+02  7.20E-04  7.28E+02
 
 TH 5
+       -7.83E+00 -8.75E+01  5.11E-04 -1.22E+02  7.74E+01
 
 TH 6
+       -1.99E+01 -2.02E+01  1.03E-03 -9.08E+00 -1.17E+01  3.75E+02
 
 TH 7
+       -2.67E-01 -2.71E+01  1.66E-03 -3.32E+01  1.78E+01  1.07E+01  1.50E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.74E-01 -1.08E-01 -4.64E-05  5.50E-01 -3.28E-02  2.80E+00  2.22E+00  0.00E+00  6.06E-01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.04E+01 -1.25E+01 -9.10E-05 -1.74E+01  2.50E+00  3.60E+00  1.62E+01  0.00E+00  1.42E-01  0.00E+00  2.57E+01
 
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
+        1.35E+03
 
 TH 2
+       -3.35E+02  3.53E+02
 
 TH 3
+        5.14E-06 -2.62E-06  4.40E-13
 
 TH 4
+       -4.79E+02  5.05E+02 -3.73E-06  7.21E+02
 
 TH 5
+        6.15E+00 -7.47E+01 -2.01E-06 -1.07E+02  6.90E+01
 
 TH 6
+       -9.87E+01 -2.64E+00  4.07E-07 -3.72E+00 -1.50E+01  2.53E+02
 
 TH 7
+        5.50E+01 -9.57E+01  2.51E-06 -1.37E+02  1.71E+01  2.33E+01  1.43E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.06E-01 -1.22E+00  3.20E-08 -1.74E+00  2.18E-01  2.95E-01  1.83E+00  0.00E+00  2.33E-02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        1.40E+02 -7.91E+01  1.62E-07 -1.13E+02  1.91E+01  2.22E+01  4.42E+01  0.00E+00  5.64E-01  0.00E+00  7.47E+01
 
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
 #CPUT: Total CPU Time in Seconds,       37.808
Stop Time:
Sat Sep 18 10:37:27 CDT 2021
