Sat Sep 18 15:31:23 CDT 2021
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
$DATA ../../../../data/spa/D/dat64.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m64.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   12166.3351017985        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -2.4104E+02  2.3491E+02 -5.7503E+01  1.3522E+02  6.0996E+02 -3.2325E+03 -8.9505E+02 -5.6014E+01 -1.5134E+03 -1.2662E+03
            -2.0889E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -597.426151711028        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.2526E+00  7.3970E-01  8.1624E-01  1.9062E+00  1.8275E+00  2.4985E+00  1.5406E+00  9.7901E-01  1.5526E+00  1.3716E+00
             1.3669E+01
 PARAMETER:  3.2519E-01 -2.0151E-01 -1.0305E-01  7.4513E-01  7.0296E-01  1.0157E+00  5.3218E-01  7.8791E-02  5.3990E-01  4.1595E-01
             2.7151E+00
 GRADIENT:  -4.0022E+01  2.5295E+01 -1.7941E+01  3.3143E+01 -9.7946E+00  1.9816E+01  1.7431E+00  6.8929E+00 -5.0077E+00  2.6343E+00
             8.9677E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -645.079980359391        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.3086E+00  4.3144E-01  1.0650E+00  2.7620E+00  6.2733E+00  2.4418E+00  8.6074E+00  4.2198E-01  2.8827E+00  7.5256E+00
             9.3727E+00
 PARAMETER:  3.6892E-01 -7.4062E-01  1.6300E-01  1.1160E+00  1.9363E+00  9.9273E-01  2.2526E+00 -7.6280E-01  1.1587E+00  2.1183E+00
             2.3378E+00
 GRADIENT:   3.5212E+00  1.2759E+01  1.4542E+01  8.7694E+01 -3.9315E+00 -9.5005E+01  3.5028E+00 -1.1921E+00  1.6702E+01  1.1322E+00
             8.5494E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -685.233178655662        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.2615E+00  5.2621E-01  6.2249E-01  1.6727E+00  6.9867E+00  2.8824E+00  1.6090E+00  8.2720E-02  2.2749E+00  9.0930E+00
             1.0191E+01
 PARAMETER:  3.3230E-01 -5.4205E-01 -3.7403E-01  6.1444E-01  2.0440E+00  1.1586E+00  5.7562E-01 -2.3923E+00  9.2194E-01  2.3075E+00
             2.4215E+00
 GRADIENT:   1.1762E+00  2.3509E+01  3.7477E-01  1.2185E+01 -1.5162E+01  9.7328E+00  5.3930E+00 -3.7617E-02 -1.9418E+01  1.2163E+01
             6.8638E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -757.726235981942        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  1.0801E+00  1.1329E-01  2.7802E-01  1.1939E+00  7.6710E+00  2.8827E+00  4.3272E+00  3.1035E+00  1.5357E+00  1.2111E+00
             9.2365E+00
 PARAMETER:  1.7708E-01 -2.0778E+00 -1.1801E+00  2.7719E-01  2.1374E+00  1.1587E+00  1.5649E+00  1.2325E+00  5.2898E-01  2.9149E-01
             2.3232E+00
 GRADIENT:   2.3434E+01  1.7352E+01 -4.0181E+00 -4.1122E+00 -9.1504E+00  6.9449E+01  2.4907E+00  1.6289E+01  4.5520E+01  4.4121E-01
             4.6655E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -822.636090049812        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  5.4174E-01  1.0000E-02  3.7167E-02  3.8436E-01  8.6722E+00  1.9814E+00  1.2493E+01  1.4582E+00  3.2728E-01  6.0761E-02
             8.2562E+00
 PARAMETER: -5.1297E-01 -4.8799E+00 -3.1923E+00 -8.5618E-01  2.2601E+00  7.8380E-01  2.6252E+00  4.7721E-01 -1.0169E+00 -2.7008E+00
             2.2110E+00
 GRADIENT:  -5.0766E+00  0.0000E+00  6.1920E-01  1.3446E+01 -1.9826E+00 -1.7773E+01  5.4922E-01 -2.9380E+00  1.5912E+00  3.0363E-04
            -3.9843E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -823.594543911200        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      549
 NPARAMETR:  5.5116E-01  1.0000E-02  3.8065E-02  3.8862E-01  1.0473E+01  2.0803E+00  1.0774E+01  1.4887E+00  3.1434E-01  6.1227E-02
             8.3063E+00
 PARAMETER: -4.9573E-01 -4.9042E+00 -3.1685E+00 -8.4516E-01  2.4488E+00  8.3252E-01  2.4771E+00  4.9791E-01 -1.0573E+00 -2.6932E+00
             2.2170E+00
 GRADIENT:  -3.2946E+00  0.0000E+00 -1.7179E-01 -6.0922E-01 -1.5592E-01 -1.4029E+00  1.2664E-01 -1.0427E-01  1.7014E+00  6.8712E-06
             8.2835E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -824.538379872508        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      725
 NPARAMETR:  5.9271E-01  1.0000E-02  4.6896E-02  4.4727E-01  6.5290E+01  2.0919E+00  3.4782E+00  1.6501E+00  7.7237E-02  9.7312E-02
             8.3107E+00
 PARAMETER: -4.2304E-01 -6.2540E+00 -2.9598E+00 -7.0458E-01  4.2788E+00  8.3808E-01  1.3465E+00  6.0083E-01 -2.4609E+00 -2.2298E+00
             2.2175E+00
 GRADIENT:   4.5491E-01  0.0000E+00  6.9652E-01 -1.1852E+00  2.0490E-02 -7.2007E-02  5.5690E-03  5.6421E-01  1.1384E-01  1.5856E-08
            -3.2294E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -824.579767715036        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      900
 NPARAMETR:  5.8893E-01  1.0000E-02  4.6163E-02  4.4268E-01  4.6008E+02  2.0919E+00  1.0962E+00  1.6412E+00  1.5625E-02  1.0795E-01
             8.3148E+00
 PARAMETER: -4.2944E-01 -7.8886E+00 -2.9756E+00 -7.1492E-01  6.2314E+00  8.3809E-01  1.9184E-01  5.9543E-01 -4.0589E+00 -2.1261E+00
             2.2180E+00
 GRADIENT:   1.5704E-02  0.0000E+00  2.5465E-02 -3.5671E-02  2.9341E-03  1.8351E-02  5.8903E-04  3.2061E-02  4.4952E-03  3.7431E-10
            -3.6528E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -824.580037617436        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1075
 NPARAMETR:  5.8859E-01  1.0000E-02  4.6090E-02  4.4221E-01  7.4279E+02  2.0917E+00  8.2524E-01  1.6402E+00  1.0559E-02  1.1119E-01
             8.3152E+00
 PARAMETER: -4.3003E-01 -8.2885E+00 -2.9772E+00 -7.1597E-01  6.7104E+00  8.3799E-01 -9.2083E-02  5.9484E-01 -4.4508E+00 -2.0966E+00
             2.2181E+00
 GRADIENT:  -8.3471E-03  0.0000E+00 -1.1179E-02  1.7206E-02  1.8163E-03 -6.1362E-03  3.3478E-04  5.6171E-04  2.0496E-03  1.6268E-10
            -3.3195E-03

0ITERATION NO.:   50    OBJECTIVE VALUE:  -824.691081218535        NO. OF FUNC. EVALS.: 105
 CUMULATIVE NO. OF FUNC. EVALS.:     1180
 NPARAMETR:  5.8798E-01  1.0000E-02  4.5984E-02  4.4210E-01  1.0230E+01  2.0888E+00  8.2521E-01  1.6388E+00  1.0000E-02  1.1118E-01
             8.3101E+00
 PARAMETER: -4.3107E-01 -8.2885E+00 -2.9795E+00 -7.1623E-01  2.4253E+00  8.3657E-01 -9.2123E-02  5.9398E-01 -4.7667E+00 -2.0966E+00
             2.2175E+00
 GRADIENT:   1.5415E+00  0.0000E+00  6.1294E+00  1.5505E+00 -3.1080E-02  2.0789E+00  2.9836E-04  9.9709E-01  0.0000E+00  2.9450E-05
             2.1286E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -824.700746322951        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1357
 NPARAMETR:  5.9112E-01  1.0000E-02  4.5982E-02  4.4255E-01  1.0364E+01  2.0939E+00  8.2519E-01  1.6324E+00  1.0000E-02  1.1112E-01
             8.3154E+00
 PARAMETER: -4.2574E-01 -8.2885E+00 -2.9795E+00 -7.1520E-01  2.4383E+00  8.3902E-01 -9.2137E-02  5.9003E-01 -4.7667E+00 -2.0972E+00
             2.2181E+00
 GRADIENT:  -4.3886E-02  0.0000E+00  3.9963E-02 -2.6122E-01  4.6629E-03 -1.3674E-02  2.9309E-04  6.2290E-03  0.0000E+00  2.1012E-05
             8.6254E-02

0ITERATION NO.:   58    OBJECTIVE VALUE:  -824.700777362659        NO. OF FUNC. EVALS.:  91
 CUMULATIVE NO. OF FUNC. EVALS.:     1448
 NPARAMETR:  5.9117E-01  1.0000E-02  4.5978E-02  4.4264E-01  1.0333E+01  2.0940E+00  8.2519E-01  1.6323E+00  1.0000E-02  1.1110E-01
             8.3145E+00
 PARAMETER: -4.2565E-01 -8.2885E+00 -2.9796E+00 -7.1500E-01  2.4354E+00  8.3908E-01 -9.2140E-02  5.9000E-01 -4.7667E+00 -2.0973E+00
             2.2180E+00
 GRADIENT:  -2.5536E-02  0.0000E+00 -1.7339E-01  5.5103E-02  1.3105E-03 -7.9622E-03  2.9362E-04 -1.6357E-02  0.0000E+00  2.2478E-05
            -3.7661E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1448
 NO. OF SIG. DIGITS IN FINAL EST.:  3.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.1642E-03 -8.4993E-05  7.7352E-03 -3.8622E-04  4.0559E-06
 SE:             2.9233E-02  5.7081E-05  2.4729E-02  2.2460E-04  1.4826E-05
 N:                     100         100         100         100         100

 P VAL.:         8.3299E-01  1.3649E-01  7.5443E-01  8.5511E-02  7.8441E-01

 ETASHRINKSD(%)  2.0672E+00  9.9809E+01  1.7155E+01  9.9248E+01  9.9950E+01
 ETASHRINKVR(%)  4.0916E+00  1.0000E+02  3.1367E+01  9.9994E+01  1.0000E+02
 EBVSHRINKSD(%)  1.9787E+00  9.9756E+01  1.6235E+01  9.9216E+01  9.9949E+01
 EBVSHRINKVR(%)  3.9182E+00  9.9999E+01  2.9834E+01  9.9994E+01  1.0000E+02
 RELATIVEINF(%)  2.0036E+01  5.0521E-05  1.2551E+00  7.6327E-05  1.0045E-05
 EPSSHRINKSD(%)  1.4808E+01
 EPSSHRINKVR(%)  2.7424E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -824.70077736265887     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -89.549950798920690     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    20.63
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.03
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -824.701       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         5.91E-01  1.00E-02  4.60E-02  4.43E-01  1.03E+01  2.09E+00  8.25E-01  1.63E+00  1.00E-02  1.11E-01  8.31E+00
 


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
+        6.92E+02
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -6.21E+02  0.00E+00  1.34E+05
 
 TH 4
+       -3.59E+02  0.00E+00 -1.96E+04  3.34E+03
 
 TH 5
+        1.72E-01  0.00E+00 -3.01E+00  4.17E-01  7.57E-03
 
 TH 6
+        4.07E+00  0.00E+00  9.13E+01 -2.37E+01  8.10E-03  4.07E+01
 
 TH 7
+       -9.96E-01  0.00E+00  1.94E+00  8.84E-01  2.99E-02  2.85E-01  3.60E+00
 
 TH 8
+        3.53E+00  0.00E+00 -1.12E+02 -4.80E+01 -1.38E-02  1.65E+00 -5.11E-01  3.58E+01
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+       -1.31E+00  0.00E+00 -1.74E+00 -3.51E-01  2.55E-03  1.01E-02 -1.70E+00 -1.05E-01  0.00E+00  5.34E-01
 
 TH11
+       -8.86E+00  0.00E+00  1.02E+02 -1.53E+01 -4.09E-03  9.36E-01 -3.61E-02  2.25E+00  0.00E+00 -2.67E-03  5.41E+00
 
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
 #CPUT: Total CPU Time in Seconds,       28.734
Stop Time:
Sat Sep 18 15:31:53 CDT 2021
