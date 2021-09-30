Thu Sep 30 00:44:09 CDT 2021
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
$DATA ../../../../data/spa1/A3/dat97.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m97.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -371.662446865917        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.4875E+02  3.2870E+00  2.1848E+02 -6.8779E+01  1.8999E+02  2.4652E+01 -2.9804E+01 -3.9745E+02 -3.7526E+01 -1.2285E+02
            -2.6945E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1528.33319457418        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.1151E-01  1.0415E+00  9.5992E-01  1.0864E+00  9.7187E-01  8.7055E-01  9.6005E-01  1.0273E+00  8.9451E-01  9.6629E-01
             3.6428E+00
 PARAMETER:  7.3507E-03  1.4070E-01  5.9098E-02  1.8291E-01  7.1466E-02 -3.8636E-02  5.9229E-02  1.2689E-01 -1.1479E-02  6.5712E-02
             1.3927E+00
 GRADIENT:  -4.3117E+01 -1.2432E+01 -1.5151E+01 -5.3553E+00 -9.9068E+00 -1.0174E+01  1.0205E+01  9.0753E+00  1.9303E+01  2.4792E+01
             1.6695E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1542.90961869225        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  9.1006E-01  9.5723E-01  8.5005E-01  1.1393E+00  8.8200E-01  9.1230E-01  1.0256E+00  1.6353E-01  8.6600E-01  8.2160E-01
             3.4123E+00
 PARAMETER:  5.7521E-03  5.6284E-02 -6.2458E-02  2.3041E-01 -2.5561E-02  8.2100E-03  1.2527E-01 -1.7107E+00 -4.3875E-02 -9.6497E-02
             1.3274E+00
 GRADIENT:  -3.0959E+01 -2.2802E+00 -1.7606E+01  1.4571E+01  8.2371E+00  3.1020E+00  8.5461E+00  3.0018E-01  1.7021E+01  1.9793E+01
             1.1240E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1554.47781022232        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      297
 NPARAMETR:  9.2288E-01  8.0009E-01  4.4170E-01  1.1723E+00  5.3132E-01  9.2837E-01  1.2426E+00  1.2429E-01  7.4734E-01  3.4836E-01
             3.0270E+00
 PARAMETER:  1.9743E-02 -1.2303E-01 -7.1712E-01  2.5898E-01 -5.3239E-01  2.5673E-02  3.1719E-01 -1.9851E+00 -1.9123E-01 -9.5451E-01
             1.2076E+00
 GRADIENT:  -2.3875E+01  1.9404E+01 -1.4619E+01  4.9656E+01  9.4771E+00 -1.9038E-02  4.0571E+00  4.9400E-02 -2.3994E+00  3.2864E+00
            -2.5414E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1556.56148574159        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      472
 NPARAMETR:  9.2894E-01  6.1282E-01  4.9700E-01  1.2499E+00  5.0904E-01  9.1837E-01  1.5391E+00  1.6054E-01  7.0291E-01  2.1368E-01
             3.0689E+00
 PARAMETER:  2.6284E-02 -3.8968E-01 -5.9916E-01  3.2309E-01 -5.7522E-01  1.4843E-02  5.3122E-01 -1.7292E+00 -2.5253E-01 -1.4433E+00
             1.2213E+00
 GRADIENT:   3.5489E-01  2.1495E+00  3.6735E+00 -1.7012E+00 -6.5674E+00 -3.6789E-01  3.8340E-01 -7.6979E-02  1.4817E-01  1.0724E-01
             1.2387E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1557.25314736542        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      651
 NPARAMETR:  9.2085E-01  5.1565E-01  9.0421E-01  1.3856E+00  7.4003E-01  9.0575E-01  1.8744E+00  2.6944E-01  6.2441E-01  2.3780E-01
             3.1556E+00
 PARAMETER:  1.7540E-02 -5.6233E-01 -6.9012E-04  4.2611E-01 -2.0106E-01  1.0066E-03  7.2829E-01 -1.2114E+00 -3.7095E-01 -1.3363E+00
             1.2492E+00
 GRADIENT:  -3.7082E+00  2.7837E+00  2.0396E+00 -6.0580E-01 -4.0129E+00 -6.9250E-01  1.9612E+00  4.2694E-01 -9.1757E-01  7.4797E-01
            -6.2148E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1558.23146533230        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      833
 NPARAMETR:  9.1661E-01  3.0377E-01  1.2918E+00  1.5586E+00  8.7974E-01  9.0048E-01  1.9790E+00  5.6827E-02  6.4864E-01  7.5389E-02
             3.2122E+00
 PARAMETER:  1.2924E-02 -1.0915E+00  3.5601E-01  5.4381E-01 -2.8125E-02 -4.8224E-03  7.8260E-01 -2.7677E+00 -3.3288E-01 -2.4851E+00
             1.2669E+00
 GRADIENT:  -3.3273E+00  2.5488E+00  1.8493E+00  1.3530E+01 -3.7387E+00 -8.8468E-01  5.9748E-01  1.5883E-02  1.8345E+00  4.6770E-02
             9.9580E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1558.38009439439        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1009
 NPARAMETR:  9.1584E-01  2.0658E-01  1.3515E+00  1.6154E+00  8.8191E-01  9.0114E-01  1.8916E+00  1.0000E-02  6.4262E-01  2.3658E-02
             3.2144E+00
 PARAMETER:  1.2082E-02 -1.4771E+00  4.0119E-01  5.7959E-01 -2.5665E-02 -4.0910E-03  7.3743E-01 -4.7265E+00 -3.4220E-01 -3.6440E+00
             1.2676E+00
 GRADIENT:  -4.7689E-01  3.5259E-01  9.4850E-02  1.9026E+00 -1.0885E-01 -1.1954E-01  5.0629E-02  0.0000E+00  5.7778E-01  4.0159E-03
             5.2031E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1558.44092125702        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1186
 NPARAMETR:  9.1451E-01  1.0210E-01  1.3819E+00  1.6797E+00  8.6491E-01  9.0083E-01  1.6362E+00  1.0000E-02  6.3048E-01  1.0000E-02
             3.2131E+00
 PARAMETER:  1.0636E-02 -2.1818E+00  4.2344E-01  6.1862E-01 -4.5130E-02 -4.4348E-03  5.9236E-01 -8.8253E+00 -3.6128E-01 -5.9565E+00
             1.2672E+00
 GRADIENT:   3.4564E-01  3.4619E-01  8.8516E-01  5.3938E+00 -2.0894E+00  9.1678E-02  3.0495E-02  0.0000E+00 -3.1895E-02  0.0000E+00
            -6.7608E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1558.46235697286        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1364
 NPARAMETR:  9.1336E-01  5.2565E-02  1.3898E+00  1.7086E+00  8.5661E-01  8.9985E-01  1.3217E+00  1.0000E-02  6.2132E-01  1.0000E-02
             3.2144E+00
 PARAMETER:  9.3737E-03 -2.8457E+00  4.2917E-01  6.3569E-01 -5.4778E-02 -5.5285E-03  3.7889E-01 -1.3183E+01 -3.7591E-01 -8.2837E+00
             1.2676E+00
 GRADIENT:  -4.3942E-01  9.9875E-02 -3.0269E-02  3.3467E+00 -4.2865E-01 -6.3686E-02  9.5106E-03  0.0000E+00 -1.1595E-01  0.0000E+00
            -8.6884E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1558.46528270349        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1541
 NPARAMETR:  9.1318E-01  3.5774E-02  1.3965E+00  1.7181E+00  8.5479E-01  8.9972E-01  1.1494E+00  1.0000E-02  6.1789E-01  1.0000E-02
             3.2143E+00
 PARAMETER:  9.1733E-03 -3.2305E+00  4.3397E-01  6.4119E-01 -5.6900E-02 -5.6672E-03  2.3926E-01 -1.5809E+01 -3.8144E-01 -9.6608E+00
             1.2676E+00
 GRADIENT:  -7.1370E-02  2.9442E-02  6.1912E-02  8.7650E-01 -2.0834E-01 -8.5591E-03  3.9850E-03  0.0000E+00 -3.0991E-02  0.0000E+00
            -4.9473E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1558.46609041226        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1728
 NPARAMETR:  9.1306E-01  3.2103E-02  1.3969E+00  1.7198E+00  8.5411E-01  8.9960E-01  9.0878E-01  1.0000E-02  6.1731E-01  1.0000E-02
             3.2141E+00
 PARAMETER:  9.0479E-03 -3.3388E+00  4.3424E-01  6.4219E-01 -5.7695E-02 -5.8000E-03  4.3525E-03 -1.5809E+01 -3.8238E-01 -9.6608E+00
             1.2675E+00
 GRADIENT:  -1.4569E-01  8.4859E-03  7.4869E-03 -3.2050E-01  2.4253E-02 -2.7478E-02  2.1657E-03  0.0000E+00  2.8123E-02  0.0000E+00
            -3.5351E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1558.46613717128        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1917
 NPARAMETR:  9.1314E-01  3.1162E-02  1.3969E+00  1.7198E+00  8.5394E-01  8.9969E-01  8.6811E-01  1.0000E-02  6.1712E-01  1.0000E-02
             3.2141E+00
 PARAMETER:  9.1304E-03 -3.3686E+00  4.3427E-01  6.4223E-01 -5.7889E-02 -5.7048E-03 -4.1437E-02 -1.5809E+01 -3.8268E-01 -9.6608E+00
             1.2675E+00
 GRADIENT:   1.7532E-01 -7.6997E-03 -1.2149E-02 -1.5969E+00  2.0178E-01  2.5283E-02  1.8998E-03  0.0000E+00  4.2951E-02  0.0000E+00
            -2.1663E-03

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1558.46629286419        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     2109             RESET HESSIAN, TYPE I
 NPARAMETR:  9.1313E-01  3.1464E-02  1.3961E+00  1.7197E+00  8.5369E-01  8.9967E-01  8.4357E-01  1.0000E-02  6.1719E-01  1.0000E-02
             3.2141E+00
 PARAMETER:  9.1252E-03 -3.3589E+00  4.3365E-01  6.4214E-01 -5.8186E-02 -5.7231E-03 -7.0110E-02 -1.5809E+01 -3.8259E-01 -9.6608E+00
             1.2675E+00
 GRADIENT:   3.1323E+01  8.1553E-02  8.4916E-01  9.2239E+01  7.3736E-01  2.8261E+00  1.9867E-03  0.0000E+00  1.9635E+00  0.0000E+00
             1.2671E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1558.46636582658        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     2297             RESET HESSIAN, TYPE I
 NPARAMETR:  9.1313E-01  3.1576E-02  1.3958E+00  1.7196E+00  8.5358E-01  8.9968E-01  8.2578E-01  1.0000E-02  6.1722E-01  1.0000E-02
             3.2141E+00
 PARAMETER:  9.1249E-03 -3.3554E+00  4.3344E-01  6.4209E-01 -5.8312E-02 -5.7180E-03 -9.1430E-02 -1.5809E+01 -3.8253E-01 -9.6608E+00
             1.2675E+00
 GRADIENT:   3.1318E+01  8.2585E-02  8.6143E-01  9.2267E+01  7.0855E-01  2.8266E+00  1.9140E-03  0.0000E+00  1.9593E+00  0.0000E+00
             1.2657E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1558.46644449032        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     2489             RESET HESSIAN, TYPE I
 NPARAMETR:  9.1313E-01  3.1797E-02  1.3952E+00  1.7195E+00  8.5349E-01  8.9967E-01  8.1013E-01  1.0000E-02  6.1726E-01  1.0000E-02
             3.2141E+00
 PARAMETER:  9.1241E-03 -3.3484E+00  4.3307E-01  6.4202E-01 -5.8421E-02 -5.7291E-03 -1.1056E-01 -1.5809E+01 -3.8246E-01 -9.6608E+00
             1.2675E+00
 GRADIENT:   3.1307E+01  8.4589E-02  8.2218E-01  9.2370E+01  7.4625E-01  2.8193E+00  1.8835E-03  0.0000E+00  1.9489E+00  0.0000E+00
             1.2644E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1558.46648532478        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     2677             RESET HESSIAN, TYPE I
 NPARAMETR:  9.1312E-01  3.1930E-02  1.3953E+00  1.7194E+00  8.5327E-01  8.9968E-01  8.0548E-01  1.0000E-02  6.1729E-01  1.0000E-02
             3.2140E+00
 PARAMETER:  9.1175E-03 -3.3442E+00  4.3307E-01  6.4200E-01 -5.8674E-02 -5.7174E-03 -1.1632E-01 -1.5809E+01 -3.8242E-01 -9.6608E+00
             1.2675E+00
 GRADIENT:   3.1260E+01  9.0225E-02  9.7378E-01  9.2569E+01  4.8266E-01  2.8223E+00  1.8762E-03  0.0000E+00  1.9392E+00  0.0000E+00
             1.2621E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1558.46653062763        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     2866             RESET HESSIAN, TYPE I
 NPARAMETR:  9.1313E-01  3.2012E-02  1.3949E+00  1.7193E+00  8.5323E-01  8.9968E-01  7.9801E-01  1.0000E-02  6.1732E-01  1.0000E-02
             3.2140E+00
 PARAMETER:  9.1222E-03 -3.3416E+00  4.3285E-01  6.4194E-01 -5.8727E-02 -5.7169E-03 -1.2564E-01 -1.5809E+01 -3.8236E-01 -9.6608E+00
             1.2675E+00
 GRADIENT:   3.1277E+01  8.8928E-02  9.3798E-01  9.2470E+01  5.4436E-01  2.8225E+00  1.8688E-03  0.0000E+00  1.9432E+00  0.0000E+00
             1.2634E+01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1558.46693951752        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     3056             RESET HESSIAN, TYPE I
 NPARAMETR:  9.1315E-01  3.2133E-02  1.3942E+00  1.7190E+00  8.5314E-01  8.9970E-01  5.8242E-01  1.0000E-02  6.1748E-01  1.0000E-02
             3.2141E+00
 PARAMETER:  9.1397E-03 -3.3379E+00  4.3232E-01  6.4176E-01 -5.8826E-02 -5.6930E-03 -4.4056E-01 -1.5809E+01 -3.8211E-01 -9.6608E+00
             1.2675E+00
 GRADIENT:   3.1355E+01  7.8043E-02  8.4049E-01  9.1892E+01  7.5076E-01  2.8348E+00  1.2907E-03  0.0000E+00  1.9731E+00  0.0000E+00
             1.2689E+01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1558.46704412444        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     3242             RESET HESSIAN, TYPE I
 NPARAMETR:  9.1314E-01  3.2050E-02  1.3941E+00  1.7191E+00  8.5308E-01  8.9969E-01  5.2112E-01  1.0000E-02  6.1745E-01  1.0000E-02
             3.2141E+00
 PARAMETER:  9.1328E-03 -3.3404E+00  4.3228E-01  6.4183E-01 -5.8903E-02 -5.7049E-03 -5.5178E-01 -1.5809E+01 -3.8216E-01 -9.6608E+00
             1.2675E+00
 GRADIENT:   3.1332E+01  7.9986E-02  8.5235E-01  9.2124E+01  7.0497E-01  2.8282E+00  1.1031E-03  0.0000E+00  1.9619E+00  0.0000E+00
             1.2664E+01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1558.46716031562        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     3424
 NPARAMETR:  9.1312E-01  3.2368E-02  1.3940E+00  1.7192E+00  8.5295E-01  8.9967E-01  4.5747E-01  1.0000E-02  6.1744E-01  1.0000E-02
             3.2140E+00
 PARAMETER:  9.1120E-03 -3.3306E+00  4.3219E-01  6.4187E-01 -5.9051E-02 -5.7290E-03 -6.8204E-01 -1.5809E+01 -3.8218E-01 -9.6608E+00
             1.2675E+00
 GRADIENT:   3.7363E-02  5.3091E-03  6.0152E-02 -5.8941E-01 -7.8915E-02 -1.2953E-03  5.7844E-04  0.0000E+00 -1.8904E-02  0.0000E+00
            -7.0887E-02

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1558.46744444941        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     3615
 NPARAMETR:  9.1313E-01  3.2600E-02  1.3939E+00  1.7190E+00  8.5283E-01  8.9968E-01  2.0159E-01  1.0000E-02  6.1751E-01  1.0000E-02
             3.2140E+00
 PARAMETER:  9.1260E-03 -3.3234E+00  4.3212E-01  6.4172E-01 -5.9194E-02 -5.7153E-03 -1.5015E+00 -1.5809E+01 -3.8206E-01 -9.6608E+00
             1.2675E+00
 GRADIENT:   6.6081E-02  3.7445E-03  1.5130E-01 -8.0313E-01 -2.0142E-01  5.9336E-03  1.1738E-04  0.0000E+00 -2.2025E-02  0.0000E+00
            -8.8241E-02

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1558.46753919098        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     3804             RESET HESSIAN, TYPE I
 NPARAMETR:  9.1314E-01  3.2694E-02  1.3935E+00  1.7188E+00  8.5284E-01  8.9969E-01  1.3062E-01  1.0000E-02  6.1759E-01  1.0000E-02
             3.2140E+00
 PARAMETER:  9.1374E-03 -3.3205E+00  4.3184E-01  6.4160E-01 -5.9183E-02 -5.7065E-03 -1.9354E+00 -1.5809E+01 -3.8193E-01 -9.6608E+00
             1.2675E+00
 GRADIENT:   3.1301E+01  8.4766E-02  9.4674E-01  9.2249E+01  5.1877E-01  2.8234E+00  1.3348E-04  0.0000E+00  1.9334E+00  0.0000E+00
             1.2618E+01

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1558.46761867631        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     3992
 NPARAMETR:  9.1315E-01  3.3028E-02  1.3932E+00  1.7185E+00  8.5281E-01  8.9971E-01  8.0195E-02  1.0000E-02  6.1774E-01  1.0000E-02
             3.2141E+00
 PARAMETER:  9.1410E-03 -3.3104E+00  4.3162E-01  6.4146E-01 -5.9213E-02 -5.6837E-03 -2.4233E+00 -1.5809E+01 -3.8169E-01 -9.6608E+00
             1.2675E+00
 GRADIENT:   1.1431E-01 -3.2640E-03  6.9530E-02 -1.2185E+00 -2.8440E-02  1.8800E-02  1.9889E-05  0.0000E+00  2.1241E-02  0.0000E+00
            -1.2588E-02

0ITERATION NO.:  120    OBJECTIVE VALUE:  -1558.46767260148        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     4181
 NPARAMETR:  9.1314E-01  3.3480E-02  1.3930E+00  1.7184E+00  8.5273E-01  8.9969E-01  5.6593E-02  1.0000E-02  6.1777E-01  1.0000E-02
             3.2140E+00
 PARAMETER:  9.1324E-03 -3.2968E+00  4.3146E-01  6.4140E-01 -5.9307E-02 -5.7021E-03 -2.7719E+00 -1.5809E+01 -3.8163E-01 -9.6608E+00
             1.2675E+00
 GRADIENT:   4.2775E-02  4.5845E-03  1.2696E-01 -7.1186E-01 -1.8488E-01  4.8084E-03  1.0244E-05  0.0000E+00 -1.2851E-02  0.0000E+00
            -6.3004E-02

0ITERATION NO.:  123    OBJECTIVE VALUE:  -1558.46770618965        NO. OF FUNC. EVALS.: 107
 CUMULATIVE NO. OF FUNC. EVALS.:     4288
 NPARAMETR:  9.1314E-01  3.4156E-02  1.3926E+00  1.7181E+00  8.5265E-01  8.9969E-01  5.3173E-02  1.0000E-02  6.1786E-01  1.0000E-02
             3.2139E+00
 PARAMETER:  9.1428E-03 -3.2981E+00  4.3127E-01  6.4135E-01 -5.9199E-02 -5.6982E-03 -2.8061E+00 -1.5809E+01 -3.8158E-01 -9.6608E+00
             1.2675E+00
 GRADIENT:   4.9966E-03 -4.1411E-03  6.7041E-03  1.1359E-01  4.2123E-02  9.0132E-04  1.3696E-05  0.0000E+00 -4.5133E-03  0.0000E+00
             1.0967E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     4288
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.9446E-03 -5.6168E-05  6.6226E-05 -1.0155E-02 -8.4846E-07
 SE:             2.8792E-02  3.1262E-05  1.1360E-04  2.5837E-02  1.8032E-04
 N:                     100         100         100         100         100

 P VAL.:         9.4615E-01  7.2391E-02  5.5992E-01  6.9430E-01  9.9625E-01

 ETASHRINKSD(%)  3.5435E+00  9.9895E+01  9.9619E+01  1.3444E+01  9.9396E+01
 ETASHRINKVR(%)  6.9614E+00  1.0000E+02  9.9999E+01  2.5080E+01  9.9996E+01
 EBVSHRINKSD(%)  3.4338E+00  9.9898E+01  9.9593E+01  1.3449E+01  9.9374E+01
 EBVSHRINKVR(%)  6.7496E+00  1.0000E+02  9.9998E+01  2.5090E+01  9.9996E+01
 RELATIVEINF(%)  8.5490E+01  2.3612E-06  6.9044E-05  2.3209E+00  1.1274E-04
 EPSSHRINKSD(%)  1.7265E+01
 EPSSHRINKVR(%)  3.1550E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1558.4677061896527     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -639.52917298498005     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    68.69
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     6.97
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1558.468       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.13E-01  3.34E-02  1.39E+00  1.72E+00  8.53E-01  9.00E-01  5.47E-02  1.00E-02  6.18E-01  1.00E-02  3.21E+00
 


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
+       -5.33E+01  3.83E+02
 
 TH 3
+        2.43E+01  6.32E+01  9.27E+01
 
 TH 4
+       -5.29E+01  4.81E+02 -4.48E+00  7.35E+02
 
 TH 5
+       -6.01E+01 -2.86E+02 -2.39E+02 -1.62E+02  6.61E+02
 
 TH 6
+       -1.26E+01 -5.78E+01 -1.32E+01 -5.86E+01  5.27E+01  1.93E+02
 
 TH 7
+       -5.09E-02 -4.16E-02 -2.39E-02 -2.36E-02  7.19E-02  3.79E-02  1.72E-05
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.38E+01 -1.11E+02 -1.40E+01 -2.88E+01  7.68E+01  3.52E+01  4.20E-02  0.00E+00  3.29E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -8.08E+01 -1.34E+01  6.97E-01 -1.18E+01  4.37E+00  3.27E+00  5.14E-03  0.00E+00  3.00E+01  0.00E+00  2.17E+01
 
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
+        1.54E+03
 
 TH 2
+       -7.79E+01  3.55E+02
 
 TH 3
+       -1.68E+01  4.61E+01  8.81E+01
 
 TH 4
+       -6.23E+01  4.68E+02 -1.60E+01  7.32E+02
 
 TH 5
+        4.61E+01 -2.35E+02 -2.15E+02 -1.30E+02  5.86E+02
 
 TH 6
+        7.31E+00 -1.17E+01  2.50E+00 -1.79E+01 -1.47E+00  2.12E+02
 
 TH 7
+       -3.81E-02 -5.69E-02 -1.36E-02  4.99E-03  4.16E-02  3.15E-02  6.19E-02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.03E+01 -7.91E+01  6.08E+00 -1.84E+01  1.81E+01 -3.46E+00  2.14E-02  0.00E+00  2.86E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.16E+01 -1.45E+01 -3.03E-02 -1.59E+01  4.76E+00  4.84E+00  9.15E-04  0.00E+00  1.88E+01  0.00E+00  5.93E+01
 
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
+        1.54E+03
 
 TH 2
+       -1.10E+02  3.44E+02
 
 TH 3
+       -7.40E+01  4.20E+01  8.35E+01
 
 TH 4
+       -8.00E+01  4.69E+02 -1.22E+01  7.40E+02
 
 TH 5
+        1.64E+02 -2.04E+02 -2.05E+02 -1.12E+02  5.54E+02
 
 TH 6
+        4.50E+01  3.09E+01  1.51E+01  2.32E+01 -4.37E+01  2.46E+02
 
 TH 7
+       -3.15E-04 -1.43E-03  2.26E-04 -8.14E-04  2.55E-04 -8.49E-04  7.80E-08
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -5.00E+01 -5.86E+01  2.11E+01 -1.50E+01 -3.00E+01 -4.30E+01  3.82E-03  0.00E+00  2.65E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        1.78E+02 -2.95E+01 -1.65E+01 -3.43E+01  3.96E+01  2.70E+01 -2.04E-04  0.00E+00 -1.79E+01  0.00E+00  2.17E+02
 
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
 #CPUT: Total CPU Time in Seconds,       75.718
Stop Time:
Thu Sep 30 00:45:27 CDT 2021
