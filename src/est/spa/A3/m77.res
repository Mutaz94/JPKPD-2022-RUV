Sat Sep 18 10:41:32 CDT 2021
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
$DATA ../../../../data/spa/A3/dat77.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m77.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -596.659607398177        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.8825E+01  2.9707E+01  1.4231E+02 -1.2311E+02 -1.0390E+00  8.1888E-01 -1.1675E+00 -5.5803E+01 -6.3292E+01 -5.2059E+01
            -1.9328E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1339.86821472736        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0570E+00  9.8527E-01  8.9136E-01  1.1565E+00  1.0152E+00  8.3118E-01  7.7510E-01  9.9685E-01  7.7893E-01  8.0304E-01
             3.8075E+00
 PARAMETER:  1.5541E-01  8.5158E-02 -1.5009E-02  2.4540E-01  1.1511E-01 -8.4908E-02 -1.5476E-01  9.6847E-02 -1.4984E-01 -1.1936E-01
             1.4370E+00
 GRADIENT:   4.1044E+01  2.7615E+01 -1.3352E+01  5.8703E+01 -2.0022E+00 -3.8802E+01  6.7684E+00  6.1036E+00  4.3435E+00  1.2504E+01
             8.7512E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1350.99434241994        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0535E+00  7.5044E-01  5.5043E-01  1.2661E+00  6.3507E-01  9.3648E-01  3.6797E-01  5.8469E-01  9.4923E-01  3.8245E-01
             3.5837E+00
 PARAMETER:  1.5208E-01 -1.8709E-01 -4.9705E-01  3.3591E-01 -3.5402E-01  3.4375E-02 -8.9975E-01 -4.3667E-01  4.7898E-02 -8.6116E-01
             1.3764E+00
 GRADIENT:   2.2677E+01  3.1469E+01 -5.4313E+00  9.4515E+01 -1.2895E+01 -3.8157E+00  6.9490E-01  4.0496E+00  1.7283E+01  3.9169E+00
             7.1300E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1361.16821387012        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0313E+00  6.2894E-01  4.2048E-01  1.2227E+00  4.9150E-01  9.5064E-01  9.6444E-01  2.7495E-01  7.8657E-01  1.4480E-01
             3.2323E+00
 PARAMETER:  1.3080E-01 -3.6371E-01 -7.6635E-01  3.0109E-01 -6.1029E-01  4.9377E-02  6.3789E-02 -1.1912E+00 -1.4007E-01 -1.8324E+00
             1.2732E+00
 GRADIENT:  -7.5476E+00  2.0988E+01  1.2705E+01  3.5273E+01 -2.0195E+01 -1.9132E+00  5.2226E-01 -4.2084E-01  1.1831E+00 -4.3352E-02
             5.7148E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1362.67561887802        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  1.0383E+00  5.7986E-01  2.7369E-01  1.1227E+00  3.6243E-01  9.6873E-01  1.1158E+00  1.9224E-01  7.9588E-01  6.8106E-02
             3.0912E+00
 PARAMETER:  1.3761E-01 -4.4497E-01 -1.1957E+00  2.1573E-01 -9.1493E-01  6.8230E-02  2.0957E-01 -1.5490E+00 -1.2831E-01 -2.5867E+00
             1.2286E+00
 GRADIENT:   1.4404E+01  3.3416E+01  4.9390E+01 -1.6769E+01 -7.4621E+01 -7.3385E-01 -8.7929E-01 -5.1039E-01 -9.9486E-01 -9.6375E-02
             1.0999E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1364.09675502158        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  1.0283E+00  4.9797E-01  2.0665E-01  1.0792E+00  3.0036E-01  9.8602E-01  1.1963E+00  1.6421E-01  8.5445E-01  3.3677E-02
             2.9199E+00
 PARAMETER:  1.2793E-01 -5.9722E-01 -1.4767E+00  1.7624E-01 -1.1028E+00  8.5921E-02  2.7925E-01 -1.7066E+00 -5.7301E-02 -3.2909E+00
             1.1716E+00
 GRADIENT:   5.2020E+00  4.4126E+00  8.7260E+00 -6.1285E+00 -1.2799E+01  2.4293E-01 -4.3833E-01 -4.5815E-01 -9.3490E-01 -4.2890E-02
            -1.6486E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1364.24085967811        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:      468
 NPARAMETR:  1.0179E+00  4.8258E-01  1.9382E-01  1.0754E+00  2.9062E-01  9.9088E-01  1.2068E+00  4.1433E-01  8.8192E-01  2.6066E-02
             2.8413E+00
 PARAMETER:  1.1775E-01 -6.2861E-01 -1.5408E+00  1.7273E-01 -1.1357E+00  9.0834E-02  2.8801E-01 -7.8109E-01 -2.5654E-02 -3.5471E+00
             1.1443E+00
 GRADIENT:  -1.6929E+01 -1.3701E+01 -2.5775E+01  1.1607E+01  3.5680E+01 -9.2631E-02  1.7991E+00 -9.2347E-01  1.0529E+00 -7.4982E-03
            -1.6582E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1365.05740524940        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      645
 NPARAMETR:  1.0248E+00  5.1403E-01  2.1720E-01  1.0927E+00  3.1251E-01  9.8527E-01  1.1654E+00  5.1488E-01  8.4831E-01  3.6313E-02
             2.8657E+00
 PARAMETER:  1.2445E-01 -5.6548E-01 -1.4269E+00  1.8867E-01 -1.0631E+00  8.5164E-02  2.5303E-01 -5.6382E-01 -6.4503E-02 -3.2156E+00
             1.1528E+00
 GRADIENT:  -7.3394E+00 -3.4322E+00 -7.5105E+00  2.0516E+00  1.3147E+01 -7.2128E-02  9.6567E-01  8.0745E-02  6.7239E-02  6.7490E-03
            -8.8859E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1365.15785651912        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      820
 NPARAMETR:  1.0302E+00  4.8294E-01  2.2195E-01  1.1086E+00  3.0621E-01  9.8602E-01  1.2183E+00  5.2841E-01  8.4199E-01  3.4603E-02
             2.8750E+00
 PARAMETER:  1.2975E-01 -6.2786E-01 -1.4053E+00  2.0307E-01 -1.0835E+00  8.5926E-02  2.9747E-01 -5.3788E-01 -7.1991E-02 -3.2638E+00
             1.1560E+00
 GRADIENT:   1.9878E-02 -3.2309E-03 -2.1435E-03  3.6115E-02 -3.0459E-02 -5.8259E-03 -2.0997E-02 -1.6584E-02 -3.4600E-02 -5.8054E-03
            -1.0979E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1365.15930527201        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:      940
 NPARAMETR:  1.0301E+00  4.8267E-01  2.2193E-01  1.1086E+00  3.0610E-01  9.8669E-01  1.2149E+00  5.3175E-01  8.4071E-01  5.8706E-02
             2.8837E+00
 PARAMETER:  1.2961E-01 -6.2843E-01 -1.4054E+00  2.0306E-01 -1.0838E+00  8.6598E-02  2.9463E-01 -5.3157E-01 -7.3503E-02 -2.7352E+00
             1.1591E+00
 GRADIENT:   4.3620E+00  5.8283E-01  1.2721E+00  1.3266E+00  5.9312E+00  6.6089E-01  3.0872E-01  1.9501E-01 -1.7354E-02 -5.6071E-03
             2.7019E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1365.20750856674        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:     1013
 NPARAMETR:  1.0306E+00  4.7935E-01  2.2226E-01  1.1108E+00  3.0504E-01  9.8582E-01  1.1735E+00  4.6166E-01  8.4423E-01  2.2133E-01
             2.8755E+00
 PARAMETER:  1.3014E-01 -6.3533E-01 -1.4039E+00  2.0507E-01 -1.0873E+00  8.5716E-02  2.6000E-01 -6.7294E-01 -6.9332E-02 -1.4081E+00
             1.1562E+00
 GRADIENT:   5.2570E+00  1.2199E+00 -1.2206E-01  1.1056E+00  1.0020E+01  3.8294E-01  6.6892E-01  3.1685E-01  2.3183E-01  2.3236E-01
             1.1991E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1365.38922208164        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:     1091
 NPARAMETR:  1.0370E+00  4.2844E-01  2.2941E-01  1.1399E+00  2.9564E-01  9.8711E-01  1.0742E+00  1.7073E-01  8.3907E-01  4.2855E-01
             2.9462E+00
 PARAMETER:  1.3638E-01 -7.4760E-01 -1.3722E+00  2.3097E-01 -1.1186E+00  8.7026E-02  1.7153E-01 -1.6677E+00 -7.5459E-02 -7.4736E-01
             1.1805E+00
 GRADIENT:   1.1000E+01  4.3341E+00 -1.2445E+00  3.8652E+00  1.1265E+01  1.2093E+00  1.8386E+00  1.5747E-01  2.2429E-01  2.3427E+00
             1.1157E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1366.44979052810        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:     1163
 NPARAMETR:  1.0368E+00  2.7670E-01  2.5408E-01  1.2146E+00  2.8112E-01  9.8247E-01  8.9615E-01  1.1252E-02  8.2615E-01  5.5240E-01
             2.8852E+00
 PARAMETER:  1.3615E-01 -1.1848E+00 -1.2701E+00  2.9445E-01 -1.1690E+00  8.2317E-02 -9.6434E-03 -4.3872E+00 -9.0980E-02 -4.9349E-01
             1.1596E+00
 GRADIENT:   8.2013E+00  7.2698E-01  1.2446E+00  7.1849E-01  1.5244E+01  4.4958E-01  1.6881E-01  1.4570E-03  1.1601E+00  1.0984E+00
            -3.6902E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1366.87864040988        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1233
 NPARAMETR:  1.0343E+00  1.7986E-01  2.6369E-01  1.2543E+00  2.7517E-01  9.7898E-01  6.9210E-01  1.0000E-02  8.0150E-01  5.5327E-01
             2.9425E+00
 PARAMETER:  1.3372E-01 -1.6156E+00 -1.2330E+00  3.2659E-01 -1.1904E+00  7.8755E-02 -2.6802E-01 -7.2984E+00 -1.2127E-01 -4.9191E-01
             1.1792E+00
 GRADIENT:   1.1045E+01  5.7268E-01  5.8609E+00 -7.8045E-01  5.1108E+00  8.5174E-01  1.0675E-02  0.0000E+00  8.5117E-01 -2.7822E-01
             4.2386E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1367.16075460749        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1303
 NPARAMETR:  1.0219E+00  8.1158E-02  2.7025E-01  1.2920E+00  2.7208E-01  9.7382E-01  3.9336E-01  1.0000E-02  7.7857E-01  5.7282E-01
             2.9211E+00
 PARAMETER:  1.2170E-01 -2.4114E+00 -1.2084E+00  3.5619E-01 -1.2017E+00  7.3476E-02 -8.3302E-01 -1.2953E+01 -1.5029E-01 -4.5719E-01
             1.1720E+00
 GRADIENT:   2.6787E+00  2.9921E-01  7.1670E+00  4.3173E+00  4.5552E-01  2.3032E-01  7.2404E-04  0.0000E+00 -1.8010E-01 -3.3037E-01
             1.4364E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1367.22766800681        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1373
 NPARAMETR:  1.0180E+00  2.4391E-02  2.6784E-01  1.3055E+00  2.6700E-01  9.7273E-01  1.5821E-01  1.0000E-02  7.7140E-01  5.8601E-01
             2.9087E+00
 PARAMETER:  1.1784E-01 -3.6135E+00 -1.2174E+00  3.6659E-01 -1.2205E+00  7.2349E-02 -1.7439E+00 -2.1895E+01 -1.5955E-01 -4.3441E-01
             1.1677E+00
 GRADIENT:   5.4796E+00  2.7201E-02  3.6043E+00  7.3184E+00  3.7932E+00  3.2974E-01  1.6930E-05  0.0000E+00 -3.1456E-01  6.9039E-03
             9.3450E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1367.23067880818        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:     1446
 NPARAMETR:  1.0166E+00  1.0013E-02  2.6422E-01  1.3049E+00  2.6326E-01  9.7257E-01  7.9800E-02  1.0000E-02  7.7216E-01  5.9068E-01
             2.9000E+00
 PARAMETER:  1.1645E-01 -4.5039E+00 -1.2310E+00  3.6615E-01 -1.2346E+00  7.2186E-02 -2.4282E+00 -2.8613E+01 -1.5856E-01 -4.2648E-01
             1.1647E+00
 GRADIENT:   5.0242E+00  1.8739E-02  4.9809E+00  9.8360E+00  3.5750E-01  2.6017E-01  7.6015E-07  0.0000E+00 -5.0714E-01 -5.3202E-02
             4.8708E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1367.24412767417        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:     1521
 NPARAMETR:  1.0176E+00  1.0000E-02  2.6123E-01  1.2973E+00  2.6189E-01  9.7323E-01  6.4834E-02  1.0000E-02  7.7676E-01  5.9024E-01
             2.8969E+00
 PARAMETER:  1.1742E-01 -4.7618E+00 -1.2423E+00  3.6029E-01 -1.2398E+00  7.2867E-02 -2.6359E+00 -3.0669E+01 -1.5262E-01 -4.2723E-01
             1.1636E+00
 GRADIENT:   7.9555E+00  0.0000E+00  1.2537E+00  1.7811E+00  9.2349E+00  4.8081E-01  5.9761E-07  0.0000E+00  2.4213E-01 -4.8116E-02
             6.2685E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1367.25979158086        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1703
 NPARAMETR:  1.0161E+00  1.0000E-02  2.5367E-01  1.2863E+00  2.5639E-01  9.7360E-01  2.7508E-02  1.0000E-02  7.8366E-01  5.9400E-01
             2.8893E+00
 PARAMETER:  1.1602E-01 -5.8546E+00 -1.2717E+00  3.5175E-01 -1.2610E+00  7.3249E-02 -3.4933E+00 -3.9187E+01 -1.4378E-01 -4.2087E-01
             1.1610E+00
 GRADIENT:   1.6468E-01  0.0000E+00 -4.6106E-01 -6.4962E-01  9.1152E-01  1.4279E-02  3.1601E-08  0.0000E+00  1.0551E-02  4.6732E-02
             4.3916E-01

0ITERATION NO.:   93    OBJECTIVE VALUE:  -1367.26017699501        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1795
 NPARAMETR:  1.0161E+00  1.0000E-02  2.5357E-01  1.2864E+00  2.5624E-01  9.7363E-01  2.6867E-02  1.0000E-02  7.8399E-01  5.9448E-01
             2.8866E+00
 PARAMETER:  1.1595E-01 -5.8861E+00 -1.2721E+00  3.5181E-01 -1.2616E+00  7.3278E-02 -3.5169E+00 -3.9426E+01 -1.4336E-01 -4.2007E-01
             1.1601E+00
 GRADIENT:   6.3431E-03  0.0000E+00  2.5842E-03 -2.9898E-03 -3.7159E-03  4.3073E-04  3.1809E-08  0.0000E+00 -1.7935E-04  4.6831E-04
             1.8300E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1795
 NO. OF SIG. DIGITS IN FINAL EST.:  3.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -5.4964E-04 -3.1155E-06  1.2482E-04 -1.2531E-02  6.0530E-04
 SE:             2.9010E-02  3.3503E-06  2.4353E-04  2.5410E-02  2.0453E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8488E-01  3.5240E-01  6.0825E-01  6.2190E-01  9.7639E-01

 ETASHRINKSD(%)  2.8122E+00  9.9989E+01  9.9184E+01  1.4875E+01  3.1481E+01
 ETASHRINKVR(%)  5.5453E+00  1.0000E+02  9.9993E+01  2.7537E+01  5.3052E+01
 EBVSHRINKSD(%)  2.8132E+00  9.9989E+01  9.9169E+01  1.4178E+01  3.1435E+01
 EBVSHRINKVR(%)  5.5472E+00  1.0000E+02  9.9993E+01  2.6345E+01  5.2989E+01
 RELATIVEINF(%)  6.2000E+01  1.3359E-07  2.0853E-04  7.8350E+00  1.3646E+00
 EPSSHRINKSD(%)  2.9880E+01
 EPSSHRINKVR(%)  5.0831E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1367.2601769950097     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -632.10935043127154     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.36
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.75
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1367.260       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.00E-02  2.54E-01  1.29E+00  2.56E-01  9.74E-01  2.69E-02  1.00E-02  7.84E-01  5.94E-01  2.89E+00
 


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
+        1.08E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -9.31E+01  0.00E+00  1.18E+04
 
 TH 4
+       -5.03E+01  0.00E+00 -9.01E+02  8.08E+02
 
 TH 5
+        2.98E+02  0.00E+00 -1.44E+04 -5.68E+02  2.21E+04
 
 TH 6
+       -6.71E+00  0.00E+00  2.43E+01 -1.21E+01  1.73E+01  1.85E+02
 
 TH 7
+       -3.22E+00  0.00E+00  2.40E-01  1.46E+00 -4.50E-01 -6.80E+00  6.14E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.57E+00  0.00E+00  1.65E+02 -2.98E+01  5.94E+01 -4.74E+00  5.52E+00  0.00E+00  1.76E+02
 
 TH10
+       -1.04E+01  0.00E+00 -1.74E+02  9.07E+00  2.42E+02 -2.27E+00 -1.01E-01  0.00E+00  5.79E+00  1.09E+02
 
 TH11
+       -1.30E+01  0.00E+00 -3.05E+00 -8.92E+00 -8.47E+00  2.84E+00  1.28E-01  0.00E+00  1.44E+01  2.81E+01  3.59E+01
 
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
 #CPUT: Total CPU Time in Seconds,       23.172
Stop Time:
Sat Sep 18 10:41:57 CDT 2021
