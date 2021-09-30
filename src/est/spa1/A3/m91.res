Thu Sep 30 00:40:27 CDT 2021
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
$DATA ../../../../data/spa1/A3/dat91.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m91.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -24.9335663215612        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.5786E+02  3.5717E+01  1.4692E+02 -1.0593E+01  2.7210E+02 -5.4626E+00 -7.3169E+01 -1.8585E+02 -1.0659E+02 -1.2675E+02
            -3.4980E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1497.10061062670        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.4938E-01  9.4987E-01  1.3333E+00  1.0270E+00  1.0876E+00  1.0163E+00  8.5653E-01  7.3469E-01  1.0802E+00  3.2625E-01
             2.9254E+00
 PARAMETER:  4.8057E-02  4.8572E-02  3.8767E-01  1.2668E-01  1.8394E-01  1.1614E-01 -5.4862E-02 -2.0830E-01  1.7716E-01 -1.0201E+00
             1.1734E+00
 GRADIENT:   6.6677E+01 -7.3776E+01 -7.0649E+00 -1.0241E+02  6.2808E+01  9.0039E+00  2.6819E+00  1.4987E+00  1.7976E-01  1.2566E+00
            -1.5307E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1510.21494636025        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.6422E-01  5.7701E-01  6.2264E-01  1.2918E+00  5.9085E-01  9.0515E-01  3.7143E-01  6.5492E-01  9.9117E-01  1.4251E-01
             2.9836E+00
 PARAMETER:  6.3569E-02 -4.4990E-01 -3.7379E-01  3.5607E-01 -4.2619E-01  3.4634E-04 -8.9041E-01 -3.2324E-01  9.1133E-02 -1.8484E+00
             1.1931E+00
 GRADIENT:   9.2168E+01 -1.1953E+01 -8.6988E+01  7.4387E+01  1.3561E+02 -4.0072E+01 -9.4491E-01  6.2655E+00 -1.2048E+01  5.0798E-01
            -9.3057E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1521.00779248402        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:      253
 NPARAMETR:  9.3218E-01  5.9886E-01  4.6811E-01  1.2318E+00  4.7822E-01  1.0132E+00  1.2247E-01  2.5657E-01  1.1263E+00  1.0034E-01
             3.2080E+00
 PARAMETER:  2.9770E-02 -4.1273E-01 -6.5905E-01  3.0845E-01 -6.3768E-01  1.1307E-01 -1.9999E+00 -1.2604E+00  2.1894E-01 -2.1991E+00
             1.2656E+00
 GRADIENT:  -3.8045E+01  1.5634E+01 -1.7236E+01  3.0042E+01  1.2712E+01 -7.1845E-01 -3.3030E-01 -3.9518E-02  4.1028E+00 -2.4828E-01
            -2.3489E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1525.65780041676        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      429
 NPARAMETR:  9.3209E-01  4.6293E-01  6.6211E-01  1.3375E+00  5.5368E-01  9.9748E-01  3.2705E-01  1.7612E-01  9.7352E-01  8.8934E-02
             3.3081E+00
 PARAMETER:  2.9678E-02 -6.7017E-01 -3.1232E-01  3.9080E-01 -4.9116E-01  9.7472E-02 -1.0176E+00 -1.6366E+00  7.3160E-02 -2.3199E+00
             1.2964E+00
 GRADIENT:  -3.2258E+01  9.9999E+00  4.9436E+00  1.0473E+01 -7.5653E+00 -1.8998E+00 -3.9297E-01  4.5406E-01 -6.7982E-01  1.1629E-01
             1.2010E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1531.65446007835        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      608
 NPARAMETR:  9.4742E-01  3.2184E-01  4.4761E-01  1.3265E+00  3.9489E-01  9.9172E-01  2.1393E+00  4.9176E-02  9.4087E-01  3.5247E-02
             3.1815E+00
 PARAMETER:  4.5989E-02 -1.0337E+00 -7.0383E-01  3.8255E-01 -8.2915E-01  9.1681E-02  8.6046E-01 -2.9123E+00  3.9046E-02 -3.2454E+00
             1.2573E+00
 GRADIENT:   5.7206E+00 -8.1176E-01  1.0802E+00  1.0231E+01  3.8388E+00 -3.1809E+00 -2.1244E+00 -6.5116E-03  5.8598E+00 -4.0888E-02
             9.3047E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1532.28787627648        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      785            RESET HESSIAN, TYPE II
 NPARAMETR:  9.4645E-01  3.2791E-01  3.9513E-01  1.2848E+00  3.6393E-01  1.0014E+00  2.1913E+00  5.0798E-02  9.3961E-01  4.0228E-02
             3.1223E+00
 PARAMETER:  4.4962E-02 -1.0150E+00 -8.2855E-01  3.5062E-01 -9.1080E-01  1.0139E-01  8.8450E-01 -2.8799E+00  3.7713E-02 -3.1132E+00
             1.2386E+00
 GRADIENT:   3.6610E+01  5.3633E+00  8.3233E+00  3.8444E+01  3.1134E+01  2.6315E+00  3.4056E+00 -1.5251E-02  2.3201E+00 -6.5101E-02
             1.3537E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1533.01385191572        NO. OF FUNC. EVALS.:  80
 CUMULATIVE NO. OF FUNC. EVALS.:      865
 NPARAMETR:  9.4505E-01  3.1603E-01  3.9134E-01  1.2700E+00  3.5841E-01  9.8760E-01  2.2360E+00  1.2363E-01  9.5153E-01  2.2620E-01
             3.0743E+00
 PARAMETER:  4.3480E-02 -1.0519E+00 -8.3819E-01  3.3900E-01 -9.2608E-01  8.7524E-02  9.0468E-01 -1.9905E+00  5.0320E-02 -1.3863E+00
             1.2231E+00
 GRADIENT:   3.6155E+01  3.8794E+00  1.0648E+01  1.3346E+01  4.1261E+01 -2.3505E+00  4.7083E+00  4.7548E-02  4.0671E+00 -9.0619E-01
             1.2718E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1533.19625737969        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      937
 NPARAMETR:  9.3423E-01  2.9738E-01  3.8184E-01  1.2638E+00  3.4658E-01  9.9248E-01  2.1980E+00  1.7095E-01  9.5634E-01  3.2282E-01
             2.9974E+00
 PARAMETER:  3.1972E-02 -1.1128E+00 -8.6275E-01  3.3416E-01 -9.5965E-01  9.2449E-02  8.8753E-01 -1.6664E+00  5.5359E-02 -1.0307E+00
             1.1977E+00
 GRADIENT:   1.5030E+01  3.3568E+00  2.4391E+01  6.7294E+00  3.1872E+01 -1.3674E+00  1.0196E+00  2.1100E-01  1.3642E+00 -7.5752E-01
            -1.7940E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1533.20292120213        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:     1012
 NPARAMETR:  9.3195E-01  2.8922E-01  3.7607E-01  1.2636E+00  3.4089E-01  9.9574E-01  2.1866E+00  1.6461E-01  9.5929E-01  3.5136E-01
             2.9774E+00
 PARAMETER:  2.9523E-02 -1.1406E+00 -8.7799E-01  3.3395E-01 -9.7621E-01  9.5729E-02  8.8233E-01 -1.7042E+00  5.8438E-02 -9.4593E-01
             1.1910E+00
 GRADIENT:   1.0689E+01  3.2984E+00  2.7675E+01  8.6017E+00  2.9253E+01 -5.0197E-01 -1.0874E-01  2.2367E-01  5.6524E-01 -5.9560E-01
            -5.4564E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1533.85342121336        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     1172
 NPARAMETR:  9.3519E-01  2.8687E-01  3.6902E-01  1.2807E+00  3.3758E-01  1.0118E+00  2.1912E+00  1.1194E-01  9.5765E-01  3.9843E-01
             3.0143E+00
 PARAMETER:  3.2994E-02 -1.1487E+00 -8.9690E-01  3.4742E-01 -9.8594E-01  1.1176E-01  8.8444E-01 -2.0898E+00  5.6723E-02 -8.2022E-01
             1.2034E+00
 GRADIENT:  -2.1204E+01 -5.6327E-01 -3.3835E+00 -5.8709E+00  9.1341E+00  1.4446E+00 -1.1964E+00  1.5611E-01 -2.5264E+00  1.3446E+00
            -6.3908E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1534.15060285147        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1347
 NPARAMETR:  9.4485E-01  2.6696E-01  3.6851E-01  1.2920E+00  3.3357E-01  1.0059E+00  2.3338E+00  3.3558E-02  9.6554E-01  3.8984E-01
             3.0124E+00
 PARAMETER:  4.3274E-02 -1.2207E+00 -8.9830E-01  3.5615E-01 -9.9789E-01  1.0586E-01  9.4750E-01 -3.2945E+00  6.4928E-02 -8.4201E-01
             1.2027E+00
 GRADIENT:   1.6618E+00  1.7948E-01  2.7943E+00  9.5468E-01 -1.4801E+00 -7.6575E-02 -4.1287E-01  1.0499E-02  4.1892E-01  2.8859E-02
            -3.5957E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1534.17165652153        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1523
 NPARAMETR:  9.4388E-01  2.5996E-01  3.6696E-01  1.2933E+00  3.3200E-01  1.0055E+00  2.3883E+00  1.0000E-02  9.6198E-01  3.8525E-01
             3.0213E+00
 PARAMETER:  4.2248E-02 -1.2472E+00 -9.0251E-01  3.5721E-01 -1.0026E+00  1.0544E-01  9.7059E-01 -4.6957E+00  6.1237E-02 -8.5386E-01
             1.2057E+00
 GRADIENT:  -1.8145E-01 -3.0148E-01  6.6428E-01  8.9936E-01  5.2403E-01 -8.9138E-02 -1.7303E-01  0.0000E+00 -2.7731E-01 -1.0456E-01
            -6.5643E-01

0ITERATION NO.:   64    OBJECTIVE VALUE:  -1534.18492053194        NO. OF FUNC. EVALS.: 128
 CUMULATIVE NO. OF FUNC. EVALS.:     1651
 NPARAMETR:  9.4454E-01  2.6798E-01  3.6329E-01  1.2867E+00  3.3085E-01  1.0062E+00  2.3409E+00  1.0000E-02  9.6600E-01  3.8164E-01
             3.0228E+00
 PARAMETER:  4.2945E-02 -1.2169E+00 -9.1256E-01  3.5206E-01 -1.0061E+00  1.0618E-01  9.5054E-01 -6.5437E+00  6.5406E-02 -8.6328E-01
             1.2062E+00
 GRADIENT:   3.9149E-01  1.0205E-01  6.2422E-01 -9.7643E-01  4.2595E-01  6.9111E-02  5.3232E-02  0.0000E+00  8.0725E-02 -7.9430E-02
             4.0468E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1651
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.3918E-03  3.8325E-02 -1.7595E-04 -2.3264E-02  9.1268E-03
 SE:             2.8936E-02  1.7104E-02  2.1353E-04  2.5787E-02  1.2507E-02
 N:                     100         100         100         100         100

 P VAL.:         9.3412E-01  2.5046E-02  4.0993E-01  3.6698E-01  4.6554E-01

 ETASHRINKSD(%)  3.0604E+00  4.2699E+01  9.9285E+01  1.3609E+01  5.8101E+01
 ETASHRINKVR(%)  6.0271E+00  6.7166E+01  9.9995E+01  2.5366E+01  8.2445E+01
 EBVSHRINKSD(%)  2.8175E+00  5.1836E+01  9.9186E+01  1.1750E+01  5.5588E+01
 EBVSHRINKVR(%)  5.5556E+00  7.6803E+01  9.9993E+01  2.2119E+01  8.0276E+01
 RELATIVEINF(%)  9.3612E+01  6.9325E+00  3.3281E-04  3.3508E+01  8.9683E-01
 EPSSHRINKSD(%)  2.3653E+01
 EPSSHRINKVR(%)  4.1711E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1534.1849205319440     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -615.24638732727135     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.24
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.65
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1534.185       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.45E-01  2.68E-01  3.63E-01  1.29E+00  3.31E-01  1.01E+00  2.34E+00  1.00E-02  9.66E-01  3.82E-01  3.02E+00
 


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
+        1.16E+03
 
 TH 2
+       -6.24E+01  7.16E+02
 
 TH 3
+       -6.66E+00  1.09E+03  5.72E+03
 
 TH 4
+       -2.05E+01  2.15E+02 -2.97E+02  5.63E+02
 
 TH 5
+        1.02E+02 -2.12E+03 -8.57E+03 -1.90E+02  1.42E+04
 
 TH 6
+        4.44E+00 -6.84E+00  1.80E+01 -1.14E+01 -4.90E+00  1.72E+02
 
 TH 7
+        9.33E-01  4.01E+01  4.27E+00 -8.86E-01 -1.99E+01  5.35E-01  6.16E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        9.22E+00 -4.27E+00  3.98E+01 -1.47E+01  8.02E+01  6.38E-01  4.61E+00  0.00E+00  1.34E+02
 
 TH10
+       -6.82E+00  1.75E+01 -1.54E+02 -9.79E+00  2.60E+02 -8.33E-01  4.82E+00  0.00E+00 -4.16E+00  6.26E+01
 
 TH11
+       -1.77E+01  2.50E+00 -6.30E+01 -1.21E+01  6.54E+01  3.17E+00  2.23E+00  0.00E+00  4.31E+00  1.93E+01  5.58E+01
 
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
 #CPUT: Total CPU Time in Seconds,       33.985
Stop Time:
Thu Sep 30 00:41:16 CDT 2021
