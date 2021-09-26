Sat Sep 25 08:30:49 CDT 2021
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
$DATA ../../../../data/spa/A2/dat21.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m21.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -815.153015599594        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -3.9004E+01  2.6926E+01  4.7081E+01 -6.7343E+01  7.5689E+01  9.1048E+00 -7.2915E+01 -5.2138E+00 -1.3748E+02 -7.9062E+01
            -1.4682E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1313.60752477930        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.1122E+00  8.7154E-01  9.9784E-01  1.2107E+00  9.2895E-01  8.0864E-01  1.2524E+00  8.6813E-01  1.4190E+00  9.6161E-01
             4.6643E+00
 PARAMETER:  2.0634E-01 -3.7498E-02  9.7838E-02  2.9123E-01  2.6296E-02 -1.1240E-01  3.2507E-01 -4.1417E-02  4.4994E-01  6.0853E-02
             1.6399E+00
 GRADIENT:   5.6609E+01 -1.0888E+01 -1.2621E+01  5.3802E+00 -1.2363E+01 -2.6393E+01  1.4890E+01  7.3259E+00  5.7420E+01  2.4814E+01
             2.1265E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1342.39012275505        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0710E+00  3.2629E-01  6.0753E-01  1.6662E+00  4.7386E-01  1.0046E+00  3.6083E+00  1.0905E-01  1.5768E+00  4.1076E-01
             3.3746E+00
 PARAMETER:  1.6861E-01 -1.0200E+00 -3.9836E-01  6.1057E-01 -6.4685E-01  1.0459E-01  1.3832E+00 -2.1160E+00  5.5540E-01 -7.8976E-01
             1.3163E+00
 GRADIENT:   1.4922E+01  4.0416E+01  6.4381E+01  1.3008E+02 -9.5296E+01  2.2914E+01  2.8560E+01  6.0057E-02  7.6161E+01 -3.4699E+00
             9.0342E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1407.83354823299        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0455E+00  3.2879E-01  3.3843E-01  1.2801E+00  3.3471E-01  9.6164E-01  1.2655E+00  1.0000E-02  1.1193E+00  4.2384E-01
             2.5997E+00
 PARAMETER:  1.4452E-01 -1.0124E+00 -9.8344E-01  3.4690E-01 -9.9451E-01  6.0882E-02  3.3545E-01 -5.7184E+00  2.1270E-01 -7.5840E-01
             1.0554E+00
 GRADIENT:  -3.1127E+01  1.4391E+01  2.1904E+01  4.7686E+01 -2.9419E+01  5.4396E-01 -5.2798E+00  0.0000E+00  3.1494E-01 -8.4350E+00
            -1.2307E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1410.36275912966        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  1.0606E+00  2.7057E-01  2.6670E-01  1.1933E+00  2.7928E-01  9.6049E-01  1.4645E+00  1.0000E-02  1.1301E+00  5.0926E-01
             2.5259E+00
 PARAMETER:  1.5881E-01 -1.2072E+00 -1.2216E+00  2.7669E-01 -1.1755E+00  5.9692E-02  4.8153E-01 -7.2633E+00  2.2230E-01 -5.7479E-01
             1.0266E+00
 GRADIENT:   1.7644E+00  7.9281E-01 -3.3220E-01  4.1196E+00  1.4108E+00 -8.0395E-02 -2.3482E+00  0.0000E+00 -1.2821E-01 -3.9110E+00
            -6.8831E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1410.59873460474        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      374
 NPARAMETR:  1.0566E+00  2.0074E-01  2.6699E-01  1.2095E+00  2.6780E-01  9.5503E-01  1.9407E+00  1.0000E-02  1.1117E+00  5.9884E-01
             2.4827E+00
 PARAMETER:  1.5504E-01 -1.5058E+00 -1.2206E+00  2.9018E-01 -1.2175E+00  5.3989E-02  7.6303E-01 -8.0344E+00  2.0593E-01 -4.1275E-01
             1.0094E+00
 GRADIENT:  -1.2947E+00  5.8213E+00  8.4762E+00  8.3392E-01 -1.6203E+01  1.2817E-01  2.8200E+00  0.0000E+00 -1.8079E+00  2.8782E+00
            -2.1324E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1415.94784437648        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      451
 NPARAMETR:  1.0454E+00  6.4422E-02  3.2357E-01  1.3258E+00  2.9537E-01  9.5156E-01  4.2012E+00  1.0000E-02  1.0522E+00  5.2556E-01
             2.5653E+00
 PARAMETER:  1.4441E-01 -2.6423E+00 -1.0284E+00  3.8205E-01 -1.1195E+00  5.0352E-02  1.5354E+00 -9.6719E+00  1.5093E-01 -5.4329E-01
             1.0421E+00
 GRADIENT:   1.3647E+00  1.2811E+00  4.8330E+00  1.1753E+01 -1.0577E+01  1.6733E+00 -4.9125E-01  0.0000E+00 -2.1679E+00 -4.4156E+00
            -1.8389E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1416.63227133163        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      526
 NPARAMETR:  1.0399E+00  1.6437E-02  3.1840E-01  1.3283E+00  2.8837E-01  9.4466E-01  9.5715E+00  1.0000E-02  1.0513E+00  5.8142E-01
             2.5297E+00
 PARAMETER:  1.3917E-01 -4.0082E+00 -1.0445E+00  3.8391E-01 -1.1435E+00  4.3071E-02  2.3588E+00 -1.3148E+01  1.5007E-01 -4.4229E-01
             1.0281E+00
 GRADIENT:  -2.0704E+00  4.6501E-01  2.7491E+00  3.9007E+00 -4.6266E+00 -1.2033E-01  3.4037E-01  0.0000E+00 -8.1634E-01  1.2941E+00
            -4.3938E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1417.30750132198        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      705
 NPARAMETR:  1.0429E+00  1.7700E-02  3.5804E-01  1.3730E+00  3.1419E-01  9.4514E-01  9.3333E+00  1.0000E-02  1.0393E+00  5.5963E-01
             2.5680E+00
 PARAMETER:  1.4199E-01 -3.9342E+00 -9.2712E-01  4.1700E-01 -1.0577E+00  4.3582E-02  2.3336E+00 -1.2543E+01  1.3852E-01 -4.8047E-01
             1.0431E+00
 GRADIENT:  -1.3779E+00  8.9174E-01  1.7297E-01 -7.9721E-01 -1.7499E+00  1.4546E-01  1.2444E+00  0.0000E+00  2.0084E-01  1.7106E-01
             3.1619E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1417.41678138312        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      885
 NPARAMETR:  1.0435E+00  1.0034E-02  3.7017E-01  1.3885E+00  3.2138E-01  9.4533E-01  1.2455E+01  1.0000E-02  1.0323E+00  5.6373E-01
             2.5656E+00
 PARAMETER:  1.4257E-01 -4.5017E+00 -8.9380E-01  4.2819E-01 -1.0351E+00  4.3781E-02  2.6221E+00 -1.3854E+01  1.3179E-01 -4.7317E-01
             1.0422E+00
 GRADIENT:   2.0176E+00  1.2052E+00 -1.1310E+00  6.0449E-01  2.6544E+00  3.1194E-01  2.4223E+00  0.0000E+00 -5.6502E-01  3.2125E-01
            -1.4183E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1417.43304410751        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1062
 NPARAMETR:  1.0428E+00  1.0019E-02  3.6443E-01  1.3816E+00  3.1775E-01  9.4394E-01  1.2417E+01  1.0000E-02  1.0336E+00  5.5708E-01
             2.5694E+00
 PARAMETER:  1.4193E-01 -4.5033E+00 -9.0943E-01  4.2325E-01 -1.0465E+00  4.2313E-02  2.6191E+00 -1.3887E+01  1.3307E-01 -4.8504E-01
             1.0437E+00
 GRADIENT:   8.6921E-02  3.7059E-01 -2.4095E-01 -5.5847E-01  3.4591E-01 -8.8166E-02  7.1966E-01  0.0000E+00 -2.2089E-01 -1.2665E-01
            -2.9049E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1417.43555293607        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:     1219
 NPARAMETR:  1.0427E+00  1.0000E-02  3.6458E-01  1.3817E+00  3.1769E-01  9.4418E-01  1.2269E+01  1.0000E-02  1.0341E+00  5.5798E-01
             2.5702E+00
 PARAMETER:  1.4185E-01 -4.5115E+00 -9.0901E-01  4.2335E-01 -1.0467E+00  4.2562E-02  2.6071E+00 -1.3887E+01  1.3351E-01 -4.8343E-01
             1.0440E+00
 GRADIENT:   7.2810E+00  0.0000E+00  2.5688E+00  9.1447E+00  4.9753E+00  4.8366E-01  8.1489E-02  0.0000E+00  5.8675E-01  1.6849E-01
             1.1676E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1417.43645726071        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:     1354            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0428E+00  1.0000E-02  3.6440E-01  1.3817E+00  3.1771E-01  9.4397E-01  1.2055E+01  1.0000E-02  1.0338E+00  5.5747E-01
             2.5690E+00
 PARAMETER:  1.4193E-01 -4.5186E+00 -9.0951E-01  4.2330E-01 -1.0466E+00  4.2343E-02  2.5895E+00 -1.3887E+01  1.3326E-01 -4.8434E-01
             1.0435E+00
 GRADIENT:   7.5732E+00  0.0000E+00  1.6440E+00  9.2574E+00  6.2269E+00  3.9230E-01  1.4998E-02  0.0000E+00  4.8520E-01  8.5349E-02
             8.0942E-01

0ITERATION NO.:   62    OBJECTIVE VALUE:  -1417.43645726071        NO. OF FUNC. EVALS.:  55
 CUMULATIVE NO. OF FUNC. EVALS.:     1409
 NPARAMETR:  1.0428E+00  1.0000E-02  3.6440E-01  1.3817E+00  3.1771E-01  9.4397E-01  1.2055E+01  1.0000E-02  1.0338E+00  5.5747E-01
             2.5690E+00
 PARAMETER:  1.4193E-01 -4.5186E+00 -9.0951E-01  4.2330E-01 -1.0466E+00  4.2343E-02  2.5895E+00 -1.3887E+01  1.3326E-01 -4.8434E-01
             1.0435E+00
 GRADIENT:  -1.3671E-02  0.0000E+00 -1.8580E-02 -2.4192E-02 -5.7440E-02 -2.4797E-03  1.1519E-03  0.0000E+00 -1.6543E-03  1.8890E-03
            -2.5408E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1409
 NO. OF SIG. DIGITS IN FINAL EST.:  3.8
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -4.7888E-04  6.3965E-04  7.9332E-05 -8.5150E-03 -2.2642E-03
 SE:             2.9045E-02  1.8126E-03  2.5899E-04  2.7626E-02  1.9666E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8685E-01  7.2416E-01  7.5937E-01  7.5791E-01  9.0834E-01

 ETASHRINKSD(%)  2.6945E+00  9.3928E+01  9.9132E+01  7.4510E+00  3.4118E+01
 ETASHRINKVR(%)  5.3163E+00  9.9631E+01  9.9992E+01  1.4347E+01  5.6595E+01
 EBVSHRINKSD(%)  2.6577E+00  9.4564E+01  9.9131E+01  6.6623E+00  3.3812E+01
 EBVSHRINKVR(%)  5.2447E+00  9.9705E+01  9.9992E+01  1.2881E+01  5.6192E+01
 RELATIVEINF(%)  7.7155E+01  3.3342E-02  2.4283E-04  1.7709E+01  1.1801E+00
 EPSSHRINKSD(%)  3.1974E+01
 EPSSHRINKVR(%)  5.3725E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1417.4364572607083     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -682.28563069697009     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.34
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.31
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1417.436       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.00E-02  3.64E-01  1.38E+00  3.18E-01  9.44E-01  1.21E+01  1.00E-02  1.03E+00  5.57E-01  2.57E+00
 


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
+        1.09E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -5.96E+01  0.00E+00  6.22E+03
 
 TH 4
+       -3.01E+01  0.00E+00 -3.04E+02  4.76E+02
 
 TH 5
+        1.83E+02  0.00E+00 -9.26E+03 -2.19E+02  1.55E+04
 
 TH 6
+        8.62E-01  0.00E+00  1.62E+01 -7.05E+00 -6.19E+00  2.11E+02
 
 TH 7
+        1.47E-02  0.00E+00 -5.37E-02 -3.03E-02  1.28E-01  2.12E-02  7.25E-03
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        9.35E+00  0.00E+00  4.02E+01 -8.29E+00  2.98E+01  2.02E+00  7.82E-03  0.00E+00  1.44E+02
 
 TH10
+       -9.69E+00  0.00E+00 -9.07E+01  7.15E-01  1.41E+02  1.18E+00 -6.38E-03  0.00E+00 -5.64E-01  1.17E+02
 
 TH11
+       -1.50E+01  0.00E+00 -5.93E+00 -5.65E+00 -7.26E+00  3.51E+00 -1.23E-02  0.00E+00  6.69E+00  3.16E+01  4.45E+01
 
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
 #CPUT: Total CPU Time in Seconds,       22.713
Stop Time:
Sat Sep 25 08:31:13 CDT 2021
