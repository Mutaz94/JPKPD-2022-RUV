Wed Sep 29 16:48:53 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat74.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m74.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1615.47293024302        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1907E+02 -2.7518E+01 -1.2245E+01 -9.7219E+00  6.4533E+01  5.4442E+01 -4.3469E+00  4.6952E+00 -8.9832E+00 -1.2710E+01
            -7.1488E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1624.77469958518        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  1.0010E+00  1.0133E+00  9.9007E-01  1.0438E+00  9.6455E-01  9.5793E-01  1.0149E+00  9.8260E-01  1.0325E+00  1.0198E+00
             1.1625E+00
 PARAMETER:  1.0105E-01  1.1321E-01  9.0019E-02  1.4288E-01  6.3908E-02  5.7018E-02  1.1480E-01  8.2445E-02  1.3201E-01  1.1957E-01
             2.5055E-01
 GRADIENT:   5.4641E+00 -6.3307E+00 -6.4844E+00 -3.2856E+00  1.3251E+01  1.4807E+00 -1.3118E+00  6.2843E+00 -5.5700E-01  2.0914E-02
             1.2493E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1627.39115179267        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  1.0039E+00  8.6650E-01  8.7391E-01  1.1341E+00  8.4243E-01  9.5271E-01  1.2144E+00  5.2804E-01  9.7028E-01  9.3068E-01
             1.1605E+00
 PARAMETER:  1.0389E-01 -4.3288E-02 -3.4777E-02  2.2587E-01 -7.1469E-02  5.1553E-02  2.9425E-01 -5.3859E-01  6.9827E-02  2.8159E-02
             2.4885E-01
 GRADIENT:   1.1727E+01  2.1981E+00 -5.4618E+00  1.2379E+01  1.6244E+01 -9.6661E-01  1.7693E-01  7.9342E-01  2.3580E+00 -1.6244E+00
             5.9034E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1628.37279862626        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      532
 NPARAMETR:  9.9849E-01  7.4730E-01  7.7777E-01  1.1927E+00  7.2961E-01  9.5732E-01  1.3807E+00  3.4533E-01  8.9994E-01  8.5565E-01
             1.1551E+00
 PARAMETER:  9.8485E-02 -1.9128E-01 -1.5133E-01  2.7620E-01 -2.1525E-01  5.6384E-02  4.2260E-01 -9.6324E-01 -5.4246E-03 -5.5893E-02
             2.4422E-01
 GRADIENT:  -2.3384E+00  9.4587E+00  4.1095E+00  8.1685E+00 -8.8775E+00  8.1830E-01  2.3819E-01  3.8356E-01 -2.8328E+00  7.9861E-01
            -2.2391E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1628.85055669827        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      707
 NPARAMETR:  9.9797E-01  5.1614E-01  7.6071E-01  1.3190E+00  6.4955E-01  9.5243E-01  1.7790E+00  1.0911E-01  8.4702E-01  8.5417E-01
             1.1595E+00
 PARAMETER:  9.7968E-02 -5.6137E-01 -1.7351E-01  3.7685E-01 -3.3147E-01  5.1257E-02  6.7605E-01 -2.1154E+00 -6.6026E-02 -5.7629E-02
             2.4795E-01
 GRADIENT:   1.5653E+00  4.1503E+00 -2.8455E-01  1.0399E+01 -2.7615E+00 -3.1448E-01  1.0854E-01 -4.5093E-02 -3.1780E-01  8.8633E-02
             4.0789E-03

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1628.95464005906        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      884
 NPARAMETR:  9.9391E-01  3.6536E-01  7.7955E-01  1.4011E+00  6.2241E-01  9.5074E-01  2.2245E+00  2.9960E-02  8.1324E-01  8.7934E-01
             1.1645E+00
 PARAMETER:  9.3896E-02 -9.0688E-01 -1.4904E-01  4.3724E-01 -3.7416E-01  4.9489E-02  8.9955E-01 -3.4079E+00 -1.0673E-01 -2.8583E-02
             2.5225E-01
 GRADIENT:  -1.1314E+00  2.2746E+00  1.7986E+00  6.0540E+00 -3.1540E+00  1.4536E-01  8.3249E-02 -6.5571E-03 -1.1327E+00  1.9514E-01
             1.4409E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1628.98466113684        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1070             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9390E-01  3.3129E-01  7.7941E-01  1.4136E+00  6.1544E-01  9.4989E-01  2.3626E+00  2.2213E-02  8.1021E-01  8.8082E-01
             1.1645E+00
 PARAMETER:  9.3880E-02 -1.0047E+00 -1.4922E-01  4.4614E-01 -3.8542E-01  4.8594E-02  9.5978E-01 -3.7071E+00 -1.1046E-01 -2.6904E-02
             2.5232E-01
 GRADIENT:   2.8434E+02  3.3965E+01  7.0159E+00  4.2874E+02  3.5177E+01  2.5164E+01  1.8046E+01  2.6107E-03  8.0930E+00  7.6378E-01
             2.0888E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1628.98524243857        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1252
 NPARAMETR:  9.9398E-01  3.3182E-01  7.7907E-01  1.4131E+00  6.1508E-01  9.4982E-01  2.3632E+00  2.5726E-02  8.1029E-01  8.8011E-01
             1.1641E+00
 PARAMETER:  9.3964E-02 -1.0032E+00 -1.4965E-01  4.4578E-01 -3.8600E-01  4.8522E-02  9.6000E-01 -3.5603E+00 -1.1036E-01 -2.7711E-02
             2.5194E-01
 GRADIENT:   9.3304E-01  5.3639E-01  1.9044E+00 -4.8912E+00 -1.8377E+00  4.7244E-02  4.1267E-01 -4.9685E-03 -2.9197E-02  2.7657E-02
            -1.0591E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1628.98746192494        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1437
 NPARAMETR:  9.9399E-01  3.3089E-01  7.7796E-01  1.4130E+00  6.1520E-01  9.4983E-01  2.3616E+00  3.0651E-02  8.1033E-01  8.7977E-01
             1.1642E+00
 PARAMETER:  9.3971E-02 -1.0060E+00 -1.5109E-01  4.4570E-01 -3.8581E-01  4.8529E-02  9.5935E-01 -3.3851E+00 -1.1031E-01 -2.8097E-02
             2.5205E-01
 GRADIENT:   9.5504E-01 -3.7164E-02 -6.5379E-01 -5.3500E+00  1.4628E+00  5.0282E-02  2.2705E-01 -6.7071E-03  1.4753E-02  6.6249E-02
             7.5729E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1628.98994542869        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1623             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9402E-01  3.3209E-01  7.7781E-01  1.4125E+00  6.1445E-01  9.4986E-01  2.3599E+00  3.8664E-02  8.1045E-01  8.7882E-01
             1.1639E+00
 PARAMETER:  9.3999E-02 -1.0024E+00 -1.5127E-01  4.4535E-01 -3.8703E-01  4.8560E-02  9.5862E-01 -3.1528E+00 -1.1016E-01 -2.9171E-02
             2.5177E-01
 GRADIENT:   2.8482E+02  3.4074E+01  7.6224E+00  4.2723E+02  3.4656E+01  2.5153E+01  1.8131E+01  5.1331E-03  8.0768E+00  6.6911E-01
             1.9122E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1628.99481335029        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1809             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9397E-01  3.3277E-01  7.7633E-01  1.4131E+00  6.1447E-01  9.4979E-01  2.3610E+00  5.4249E-02  8.1039E-01  8.7762E-01
             1.1636E+00
 PARAMETER:  9.3948E-02 -1.0003E+00 -1.5318E-01  4.4577E-01 -3.8700E-01  4.8481E-02  9.5909E-01 -2.8142E+00 -1.1024E-01 -3.0539E-02
             2.5154E-01
 GRADIENT:   2.8463E+02  3.4217E+01  5.1432E+00  4.3001E+02  3.7491E+01  2.5110E+01  1.8388E+01  7.9702E-03  8.0140E+00  6.7473E-01
             1.9636E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1629.03798214743        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1988
 NPARAMETR:  9.9108E-01  3.1853E-01  7.7975E-01  1.4216E+00  6.1219E-01  9.4883E-01  2.3525E+00  2.0030E-01  8.0912E-01  8.7850E-01
             1.1641E+00
 PARAMETER:  9.1038E-02 -1.0440E+00 -1.4879E-01  4.5176E-01 -3.9071E-01  4.7472E-02  9.5546E-01 -1.5080E+00 -1.1181E-01 -2.9542E-02
             2.5193E-01
 GRADIENT:  -5.5787E+00 -7.8331E-01 -1.7290E+00 -1.5393E+00 -1.1872E+00 -2.9186E-01 -2.4661E+00 -9.5307E-02  5.0847E-01  2.0202E+00
             1.5282E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1629.15854258551        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2166
 NPARAMETR:  9.9033E-01  3.0657E-01  7.8537E-01  1.4241E+00  6.1150E-01  9.4786E-01  2.4703E+00  3.2484E-01  8.0541E-01  8.5262E-01
             1.1577E+00
 PARAMETER:  9.0286E-02 -1.0823E+00 -1.4160E-01  4.5354E-01 -3.9184E-01  4.6454E-02  1.0044E+00 -1.0244E+00 -1.1640E-01 -5.9445E-02
             2.4639E-01
 GRADIENT:  -6.2299E+00 -4.5605E-01 -1.4151E+00 -1.1730E+01 -1.4628E+00 -5.9122E-01 -8.1103E-02  4.0367E-02 -2.5052E-01  2.9162E-01
            -3.3438E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1629.17957750192        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:     2325
 NPARAMETR:  9.9275E-01  3.0631E-01  7.8729E-01  1.4289E+00  6.1247E-01  9.4939E-01  2.4706E+00  3.2564E-01  8.0612E-01  8.5271E-01
             1.1583E+00
 PARAMETER:  9.2719E-02 -1.0832E+00 -1.3916E-01  4.5691E-01 -3.9025E-01  4.8064E-02  1.0045E+00 -1.0220E+00 -1.1552E-01 -5.9341E-02
             2.4698E-01
 GRADIENT:  -2.8605E-01  3.1653E-01 -1.9708E+00 -3.0790E+00 -1.2387E+00  4.3816E-02 -2.5835E-01  1.3249E-03  1.3060E-01 -4.3329E-02
             2.5271E-03

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1629.18731966477        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2500
 NPARAMETR:  9.9190E-01  2.8619E-01  8.0412E-01  1.4467E+00  6.1762E-01  9.4878E-01  2.5543E+00  3.6072E-01  8.0214E-01  8.6179E-01
             1.1594E+00
 PARAMETER:  9.1868E-02 -1.1511E+00 -1.1800E-01  4.6930E-01 -3.8188E-01  4.7422E-02  1.0378E+00 -9.1964E-01 -1.2047E-01 -4.8741E-02
             2.4792E-01
 GRADIENT:  -7.6887E-01  5.9353E-01 -2.9704E+00  5.2708E+00 -1.3092E-01  1.0644E-02 -1.1210E+00  5.2447E-02  2.0219E-01  4.6727E-01
             4.4699E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1629.18940573061        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2675
 NPARAMETR:  9.9096E-01  2.6045E-01  8.2640E-01  1.4659E+00  6.2408E-01  9.4802E-01  2.7012E+00  4.0665E-01  7.9693E-01  8.7001E-01
             1.1598E+00
 PARAMETER:  9.0915E-02 -1.2454E+00 -9.0678E-02  4.8244E-01 -3.7147E-01  4.6623E-02  1.0937E+00 -7.9980E-01 -1.2699E-01 -3.9246E-02
             2.4824E-01
 GRADIENT:  -9.5173E-01  7.1702E-01 -2.3718E+00  8.4013E+00 -1.2128E-01 -1.5661E-02 -1.2995E+00  9.4131E-02  2.7147E-01  8.0722E-01
             6.4450E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1629.19010268221        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2851
 NPARAMETR:  9.9042E-01  2.4395E-01  8.4055E-01  1.4770E+00  6.2812E-01  9.4754E-01  2.8129E+00  4.3351E-01  7.9319E-01  8.7358E-01
             1.1597E+00
 PARAMETER:  9.0371E-02 -1.3108E+00 -7.3696E-02  4.9001E-01 -3.6503E-01  4.6118E-02  1.1342E+00 -7.3585E-01 -1.3170E-01 -3.5159E-02
             2.4814E-01
 GRADIENT:  -9.0527E-01  7.0468E-01 -1.8251E+00  8.3050E+00  1.3801E-01 -3.2993E-02 -1.1565E+00  7.3055E-02  1.3523E-01  7.2187E-01
             5.6655E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1629.23676005533        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     3035             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9093E-01  2.3997E-01  8.4369E-01  1.4727E+00  6.2872E-01  9.4752E-01  2.8683E+00  4.4827E-01  7.9360E-01  8.7148E-01
             1.1585E+00
 PARAMETER:  9.0889E-02 -1.3272E+00 -6.9967E-02  4.8711E-01 -3.6408E-01  4.6091E-02  1.1537E+00 -7.0236E-01 -1.3118E-01 -3.7566E-02
             2.4716E-01
 GRADIENT:   2.8818E+02  2.8108E+01  2.9242E+00  5.1379E+02  3.2331E+01  2.5328E+01  2.1356E+01  6.0556E-01  1.0536E+01  1.3042E+00
             2.3638E+00

0ITERATION NO.:   88    OBJECTIVE VALUE:  -1629.23791095002        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     3127
 NPARAMETR:  9.9068E-01  2.3955E-01  8.4392E-01  1.4732E+00  6.2868E-01  9.4741E-01  2.8653E+00  4.4378E-01  7.9328E-01  8.6974E-01
             1.1583E+00
 PARAMETER:  9.0639E-02 -1.3290E+00 -6.9697E-02  4.8742E-01 -3.6413E-01  4.5980E-02  1.1527E+00 -7.1243E-01 -1.3158E-01 -3.9566E-02
             2.4693E-01
 GRADIENT:  -6.8334E-01 -2.7475E-02  5.1970E-02  1.6244E+00  1.4465E-01 -6.3085E-02 -2.0455E-01  7.7931E-03  9.4846E-02  4.5191E-02
             5.3722E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     3127
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.5224E-04  2.3429E-02 -1.4114E-02 -1.6415E-02 -9.0894E-03
 SE:             2.9771E-02  1.4324E-02  9.7035E-03  2.6921E-02  2.3189E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7716E-01  1.0192E-01  1.4580E-01  5.4201E-01  6.9509E-01

 ETASHRINKSD(%)  2.6448E-01  5.2012E+01  6.7492E+01  9.8125E+00  2.2313E+01
 ETASHRINKVR(%)  5.2825E-01  7.6971E+01  8.9432E+01  1.8662E+01  3.9647E+01
 EBVSHRINKSD(%)  6.3246E-01  5.9336E+01  6.8229E+01  8.3643E+00  1.8581E+01
 EBVSHRINKVR(%)  1.2609E+00  8.3465E+01  8.9906E+01  1.6029E+01  3.3709E+01
 RELATIVEINF(%)  9.7483E+01  2.7073E+00  1.0152E+00  1.8936E+01  5.5205E+00
 EPSSHRINKSD(%)  4.3231E+01
 EPSSHRINKVR(%)  6.7773E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1629.2379109500191     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -894.08708438628094     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    44.96
 Elapsed covariance  time in seconds:     6.96
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1629.238       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.91E-01  2.40E-01  8.44E-01  1.47E+00  6.29E-01  9.47E-01  2.87E+00  4.44E-01  7.93E-01  8.70E-01  1.16E+00
 


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
 
         3.34E-02  7.61E-01  1.48E-01  4.00E-01  1.66E-01  5.84E-02  5.12E+00  3.87E-01  1.33E-01  1.22E-01  9.40E-02
 


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
+        1.11E-03
 
 TH 2
+        1.34E-02  5.80E-01
 
 TH 3
+       -9.08E-04 -2.40E-02  2.20E-02
 
 TH 4
+       -7.07E-03 -3.03E-01  1.54E-02  1.60E-01
 
 TH 5
+        2.34E-03  1.10E-01  6.98E-03 -5.58E-02  2.76E-02
 
 TH 6
+        3.43E-04  1.56E-02 -3.41E-04 -7.95E-03  3.08E-03  3.41E-03
 
 TH 7
+       -8.98E-02 -3.89E+00  1.66E-01  2.04E+00 -7.37E-01 -1.06E-01  2.62E+01
 
 TH 8
+       -3.94E-03 -1.79E-01  3.48E-02  9.67E-02 -1.92E-02 -2.97E-03  1.20E+00  1.50E-01
 
 TH 9
+        1.97E-03  8.64E-02 -6.75E-03 -4.60E-02  1.48E-02  2.07E-03 -5.75E-01 -3.10E-02  1.78E-02
 
 TH10
+        3.95E-04  1.95E-02  5.15E-03 -9.17E-03  7.39E-03 -3.86E-04 -1.20E-01 -1.16E-02  8.94E-04  1.49E-02
 
 TH11
+       -5.69E-05 -9.94E-03  2.04E-03  5.32E-03 -9.85E-04 -1.05E-03  6.94E-02 -3.84E-04 -1.79E-03  7.26E-04  8.84E-03
 
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
+        3.34E-02
 
 TH 2
+        5.26E-01  7.61E-01
 
 TH 3
+       -1.83E-01 -2.12E-01  1.48E-01
 
 TH 4
+       -5.30E-01 -9.95E-01  2.60E-01  4.00E-01
 
 TH 5
+        4.22E-01  8.70E-01  2.84E-01 -8.40E-01  1.66E-01
 
 TH 6
+        1.76E-01  3.52E-01 -3.94E-02 -3.40E-01  3.18E-01  5.84E-02
 
 TH 7
+       -5.25E-01 -9.99E-01  2.18E-01  9.94E-01 -8.67E-01 -3.56E-01  5.12E+00
 
 TH 8
+       -3.05E-01 -6.08E-01  6.07E-01  6.25E-01 -3.00E-01 -1.32E-01  6.07E-01  3.87E-01
 
 TH 9
+        4.43E-01  8.50E-01 -3.41E-01 -8.62E-01  6.66E-01  2.66E-01 -8.42E-01 -6.01E-01  1.33E-01
 
 TH10
+        9.69E-02  2.10E-01  2.84E-01 -1.88E-01  3.65E-01 -5.42E-02 -1.92E-01 -2.45E-01  5.49E-02  1.22E-01
 
 TH11
+       -1.81E-02 -1.39E-01  1.46E-01  1.41E-01 -6.31E-02 -1.91E-01  1.44E-01 -1.06E-02 -1.43E-01  6.32E-02  9.40E-02
 
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
+        1.28E+03
 
 TH 2
+       -1.29E+01  1.10E+03
 
 TH 3
+        5.56E+01 -5.83E+01  1.03E+03
 
 TH 4
+        3.08E+01  2.84E+02 -3.21E+01  8.02E+02
 
 TH 5
+        1.40E+01  3.28E+01 -1.68E+03 -1.62E+02  3.18E+03
 
 TH 6
+       -5.01E+00 -2.80E+01 -1.51E+01 -4.23E+01  1.34E+01  3.58E+02
 
 TH 7
+        1.48E+00  1.37E+02 -5.46E+01 -2.26E+01  1.11E+02  1.47E+00  2.47E+01
 
 TH 8
+       -2.47E+01 -1.23E+01 -2.47E+01  1.73E+01 -3.47E+01  2.76E-01 -4.76E+00  2.60E+01
 
 TH 9
+        1.17E+01 -2.24E+02  1.62E+02  1.19E+02 -3.10E+02  1.83E+01 -4.54E+01  1.80E+01  3.29E+02
 
 TH10
+       -3.02E+01 -1.52E+02  6.41E+01  3.42E+01 -2.52E+02  2.80E+01 -3.15E+01  3.31E+01  9.29E+01  1.57E+02
 
 TH11
+       -4.45E+01 -3.22E+01 -1.57E+01  2.11E+01 -4.05E+01  3.13E+01 -8.01E+00  1.69E+01  2.28E+01  2.36E+01  1.33E+02
 
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
 #CPUT: Total CPU Time in Seconds,       51.966
Stop Time:
Wed Sep 29 16:49:46 CDT 2021
