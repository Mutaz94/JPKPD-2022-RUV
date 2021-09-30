Thu Sep 30 09:08:51 CDT 2021
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
$DATA ../../../../data/spa2/D/dat48.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      700
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

 TOT. NO. OF OBS RECS:      600
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
 RAW OUTPUT FILE (FILE): m48.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   20203.8892546932        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.7218E+02  3.3552E+02  1.8538E+01  2.9160E+02  1.1188E+02 -1.7726E+03 -9.0245E+02 -8.4947E+01 -1.2436E+03 -4.6612E+02
            -4.0461E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -641.410312082164        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.4531E+00  1.2426E+00  9.7288E-01  1.4495E+00  1.0443E+00  1.9031E+00  1.5661E+00  9.9447E-01  1.4058E+00  1.1158E+00
             1.4288E+01
 PARAMETER:  4.7368E-01  3.1719E-01  7.2507E-02  4.7119E-01  1.4338E-01  7.4347E-01  5.4862E-01  9.4454E-02  4.4059E-01  2.0960E-01
             2.7594E+00
 GRADIENT:   4.0280E+01 -2.2511E+01 -1.8688E+01  1.4997E+01  2.0170E+01  4.6725E+01 -2.0354E+01  3.5179E+00 -4.6326E+00  1.3422E+01
             2.6833E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -723.398174488171        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.4125E+00  1.5030E+00  7.0325E+00  1.5571E+00  2.9551E+00  2.0942E+00  6.8826E+00  6.7815E-01  1.8260E+00  2.0169E+00
             1.2767E+01
 PARAMETER:  4.4535E-01  5.0748E-01  2.0505E+00  5.4279E-01  1.1835E+00  8.3916E-01  2.0290E+00 -2.8839E-01  7.0212E-01  8.0154E-01
             2.6469E+00
 GRADIENT:   5.9785E+01  1.5174E+01 -2.9418E+00 -3.1510E+01 -4.0814E+00  5.4886E+01  9.2423E+01  3.9957E-02  3.6064E+01  1.6109E+01
             2.7632E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -760.065111769870        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.2404E+00  1.0523E+00  9.3284E+00  1.6055E+00  2.6680E+00  1.7291E+00  4.9815E+00  9.4564E-01  1.3089E+00  1.6462E+00
             1.2267E+01
 PARAMETER:  3.1545E-01  1.5102E-01  2.3331E+00  5.7345E-01  1.0813E+00  6.4759E-01  1.7057E+00  4.4104E-02  3.6921E-01  5.9845E-01
             2.6069E+00
 GRADIENT:   4.7965E+00 -4.0226E+00 -1.0706E+00  1.0355E+01 -2.4641E+00  1.9338E+01  2.9478E+00  6.5808E-02  1.4104E+01  1.3095E+01
             2.6297E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -796.318767091018        NO. OF FUNC. EVALS.: 152
 CUMULATIVE NO. OF FUNC. EVALS.:      378
 NPARAMETR:  1.1107E+00  8.7972E-01  1.0943E+01  1.4640E+00  2.4463E+00  1.6268E+00  5.5288E+00  6.6466E+00  1.0572E+00  1.0589E+00
             1.0132E+01
 PARAMETER:  2.0502E-01 -2.8150E-02  2.4927E+00  4.8115E-01  9.9458E-01  5.8664E-01  1.8100E+00  1.9941E+00  1.5564E-01  1.5725E-01
             2.4157E+00
 GRADIENT:  -3.1267E+01 -2.9660E+00 -3.2222E+00  2.9985E+00 -9.9082E-02 -1.0052E+01 -1.3574E+01  5.8666E+00  6.9200E+00  7.2159E+00
             8.7909E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -808.605955204470        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      554
 NPARAMETR:  1.1484E+00  7.8272E-01  8.9135E+00  1.3865E+00  2.2048E+00  1.7442E+00  5.9602E+00  1.3003E+00  8.9445E-01  3.6392E-01
             9.2262E+00
 PARAMETER:  2.3838E-01 -1.4498E-01  2.2876E+00  4.2680E-01  8.9062E-01  6.5631E-01  1.8851E+00  3.6261E-01 -1.1547E-02 -9.1083E-01
             2.3220E+00
 GRADIENT:   4.2708E+00 -2.6708E+00 -6.4487E-01 -3.8609E+00  4.3876E-01 -1.4702E+00  1.6678E+00  1.4442E-01  8.1918E-01  8.8786E-01
            -3.3156E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -808.952083711399        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      730
 NPARAMETR:  1.1420E+00  9.5060E-01  7.2872E+00  1.3025E+00  2.1899E+00  1.7577E+00  5.5551E+00  7.6260E-01  7.7875E-01  1.3134E-01
             9.2503E+00
 PARAMETER:  2.3282E-01  4.9340E-02  2.0861E+00  3.6429E-01  8.8384E-01  6.6398E-01  1.8147E+00 -1.7102E-01 -1.5006E-01 -1.9300E+00
             2.3247E+00
 GRADIENT:  -2.8186E-01 -7.6201E-01 -7.0770E-01  4.2731E-01  1.1909E+00  4.2597E-01 -4.4261E-01  6.6036E-02 -6.4907E-01  1.1023E-01
             1.6993E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -809.180305084198        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      909
 NPARAMETR:  1.1415E+00  9.4487E-01  1.0126E+01  1.3287E+00  2.2443E+00  1.7551E+00  5.6347E+00  4.4036E-01  8.3779E-01  5.5439E-02
             9.2450E+00
 PARAMETER:  2.3238E-01  4.3289E-02  2.4151E+00  3.8424E-01  9.0838E-01  6.6253E-01  1.8289E+00 -7.2017E-01 -7.6990E-02 -2.7925E+00
             2.3241E+00
 GRADIENT:  -5.6870E-01  7.0601E-01 -1.6887E-01  7.2368E-01 -3.0987E-01  8.1151E-02  9.6349E-01  1.0500E-02  2.4625E-01  1.9525E-02
            -9.4276E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -809.258969016485        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     1069
 NPARAMETR:  1.1428E+00  9.1903E-01  1.4799E+01  1.3425E+00  2.2826E+00  1.7528E+00  5.7417E+00  2.5009E-01  8.4139E-01  1.7570E-02
             9.2461E+00
 PARAMETER:  2.3346E-01  1.5559E-02  2.7946E+00  3.9453E-01  9.2533E-01  6.6121E-01  1.8477E+00 -1.2859E+00 -7.2698E-02 -3.9416E+00
             2.3242E+00
 GRADIENT:   2.1608E-01  6.5966E-01  1.0182E-01 -3.7170E-01 -1.8203E+00 -1.0404E-01  2.6757E+00  1.4682E-03 -2.9873E-02  1.9434E-03
            -8.9347E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -809.280138038716        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1247
 NPARAMETR:  1.1436E+00  9.0474E-01  1.4837E+01  1.3481E+00  2.3041E+00  1.7554E+00  5.7538E+00  1.3517E-01  8.4072E-01  1.0000E-02
             9.2577E+00
 PARAMETER:  2.3421E-01 -1.1025E-04  2.7972E+00  3.9870E-01  9.3468E-01  6.6269E-01  1.8499E+00 -1.9012E+00 -7.3499E-02 -4.6930E+00
             2.3255E+00
 GRADIENT:   3.6209E-01  9.3074E-02  2.3362E-02 -5.1975E-01 -5.4084E-01  4.7107E-01  2.1027E+00  4.3083E-04 -1.6180E-01  0.0000E+00
             6.3295E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -809.291657789440        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1428             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1434E+00  9.0079E-01  1.5201E+01  1.3520E+00  2.3184E+00  1.7533E+00  5.8032E+00  9.9385E-02  8.5052E-01  1.0000E-02
             9.2532E+00
 PARAMETER:  2.3402E-01 -4.4784E-03  2.8214E+00  4.0160E-01  9.4086E-01  6.6149E-01  1.8584E+00 -2.2088E+00 -6.1904E-02 -4.6988E+00
             2.3250E+00
 GRADIENT:   1.5716E+01  9.2855E-01 -5.0107E-03  8.7604E+00  1.1789E+00  9.6566E+00  5.7720E+01  2.2945E-04  4.1140E-01  0.0000E+00
             2.8199E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -809.295220817992        NO. OF FUNC. EVALS.: 173
 CUMULATIVE NO. OF FUNC. EVALS.:     1601
 NPARAMETR:  1.1443E+00  8.9217E-01  1.5528E+01  1.3576E+00  2.3189E+00  1.7550E+00  5.7913E+00  9.0475E-02  8.5010E-01  1.0000E-02
             9.2625E+00
 PARAMETER:  2.3477E-01 -1.4095E-02  2.8427E+00  4.0571E-01  9.4108E-01  6.6246E-01  1.8564E+00 -2.3027E+00 -6.2407E-02 -4.6988E+00
             2.3260E+00
 GRADIENT:   5.7021E-01  6.5471E-02 -4.5727E-03 -2.9139E-01 -7.2830E-02  4.9202E-01  2.2912E+00  1.7731E-04 -1.2715E-01  0.0000E+00
             9.1754E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -809.304098245741        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1779
 NPARAMETR:  1.1440E+00  8.8498E-01  1.5809E+01  1.3616E+00  2.3210E+00  1.7545E+00  5.8245E+00  7.2259E-02  8.5366E-01  1.0000E-02
             9.2608E+00
 PARAMETER:  2.3455E-01 -2.2185E-02  2.8606E+00  4.0863E-01  9.4200E-01  6.6219E-01  1.8621E+00 -2.5275E+00 -5.8219E-02 -4.6988E+00
             2.3258E+00
 GRADIENT:   5.2863E-01  1.0557E-01 -4.5534E-03 -5.5540E-01 -4.5515E-02  4.3078E-01  2.8870E+00  1.0953E-04 -9.2086E-02  0.0000E+00
             7.0847E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -809.309712638901        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1957
 NPARAMETR:  1.1438E+00  8.7913E-01  1.6081E+01  1.3653E+00  2.3227E+00  1.7539E+00  5.8474E+00  5.4082E-02  8.5732E-01  1.0000E-02
             9.2589E+00
 PARAMETER:  2.3431E-01 -2.8823E-02  2.8776E+00  4.1136E-01  9.4273E-01  6.6182E-01  1.8660E+00 -2.8173E+00 -5.3948E-02 -4.6988E+00
             2.3256E+00
 GRADIENT:   4.5032E-01  1.4225E-01 -3.0409E-03 -5.1847E-01 -6.5098E-02  3.4056E-01  3.1778E+00  5.9650E-05 -8.5636E-02  0.0000E+00
             3.8408E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -809.317166141078        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2142
 NPARAMETR:  1.1427E+00  8.6333E-01  1.6422E+01  1.3716E+00  2.3261E+00  1.7522E+00  5.8932E+00  4.7073E-02  8.6800E-01  1.0000E-02
             9.2533E+00
 PARAMETER:  2.3339E-01 -4.6954E-02  2.8986E+00  4.1598E-01  9.4417E-01  6.6084E-01  1.8738E+00 -2.9561E+00 -4.1561E-02 -4.6988E+00
             2.3250E+00
 GRADIENT:   1.3864E-01 -8.4151E-02 -1.2275E-02 -1.6504E+00  1.8679E-01  5.9419E-02  3.7984E+00  4.3858E-05  1.8577E-01  0.0000E+00
             1.5317E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -809.320779079885        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2328             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1430E+00  8.5752E-01  1.6794E+01  1.3766E+00  2.3284E+00  1.7522E+00  5.9123E+00  4.0767E-02  8.7232E-01  1.0000E-02
             9.2557E+00
 PARAMETER:  2.3361E-01 -5.3716E-02  2.9210E+00  4.1958E-01  9.4518E-01  6.6088E-01  1.8770E+00 -3.0999E+00 -3.6603E-02 -4.6988E+00
             2.3252E+00
 GRADIENT:   1.5440E+01  5.7572E-01  2.2562E-03  1.0285E+01  1.0857E+00  9.7098E+00  5.9832E+01  3.3074E-05  3.0090E-01  0.0000E+00
             2.8091E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -809.321839162521        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2504
 NPARAMETR:  1.1433E+00  8.5706E-01  1.6992E+01  1.3792E+00  2.3286E+00  1.7528E+00  5.9129E+00  3.4735E-02  8.7150E-01  1.0000E-02
             9.2579E+00
 PARAMETER:  2.3390E-01 -5.4245E-02  2.9328E+00  4.2147E-01  9.4527E-01  6.6123E-01  1.8771E+00 -3.2600E+00 -3.7538E-02 -4.6988E+00
             2.3255E+00
 GRADIENT:   2.5362E-01  9.6827E-02 -1.9219E-03 -1.1502E-01 -9.6857E-02  2.7498E-01  3.5484E+00  2.2429E-05 -8.0653E-02  0.0000E+00
            -2.0666E-02

0ITERATION NO.:   85    OBJECTIVE VALUE:  -809.323161835951        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2688
 NPARAMETR:  1.1427E+00  8.5329E-01  1.7353E+01  1.3819E+00  2.3299E+00  1.7515E+00  5.9322E+00  2.3181E-02  8.7332E-01  1.0000E-02
             9.2522E+00
 PARAMETER:  2.3339E-01 -5.8656E-02  2.9538E+00  4.2346E-01  9.4582E-01  6.6046E-01  1.8804E+00 -3.6644E+00 -3.5457E-02 -4.6988E+00
             2.3249E+00
 GRADIENT:   7.7855E-02  1.9175E-01  3.8544E-03  1.4052E-01 -1.9148E-01  4.5941E-02  3.8135E+00  9.6036E-06 -1.3554E-01  0.0000E+00
            -9.1805E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -809.324448362686        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2871             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1427E+00  8.4835E-01  1.7600E+01  1.3842E+00  2.3326E+00  1.7512E+00  5.9433E+00  2.2980E-02  8.7672E-01  1.0000E-02
             9.2531E+00
 PARAMETER:  2.3344E-01 -6.4465E-02  2.9679E+00  4.2516E-01  9.4700E-01  6.6029E-01  1.8823E+00 -3.6731E+00 -3.1563E-02 -4.6988E+00
             2.3250E+00
 GRADIENT:   1.5335E+01  6.9483E-01  1.1259E-02  1.1614E+01  8.9020E-01  9.5904E+00  6.0297E+01  9.6968E-06  9.7911E-02  0.0000E+00
             2.7195E+01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -809.324768534398        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     3056
 NPARAMETR:  1.1428E+00  8.4682E-01  1.7667E+01  1.3852E+00  2.3331E+00  1.7520E+00  5.9458E+00  2.1744E-02  8.7764E-01  1.0000E-02
             9.2552E+00
 PARAMETER:  2.3350E-01 -6.6264E-02  2.9717E+00  4.2585E-01  9.4718E-01  6.6075E-01  1.8827E+00 -3.7284E+00 -3.0520E-02 -4.6988E+00
             2.3252E+00
 GRADIENT:   9.2530E-02  6.6557E-02  9.3999E-04 -8.7961E-02 -8.8526E-02  1.6710E-01  3.8189E+00  8.1699E-06 -6.5973E-02  0.0000E+00
            -4.2339E-01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -809.325057482707        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     3246             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1429E+00  8.4471E-01  1.7720E+01  1.3860E+00  2.3342E+00  1.7519E+00  5.9504E+00  2.1201E-02  8.7959E-01  1.0000E-02
             9.2560E+00
 PARAMETER:  2.3356E-01 -6.8763E-02  2.9747E+00  4.2641E-01  9.4768E-01  6.6069E-01  1.8835E+00 -3.7537E+00 -2.8297E-02 -4.6988E+00
             2.3253E+00
 GRADIENT:   1.5360E+01  6.0325E-01  8.2127E-03  1.1413E+01  9.7281E-01  9.7354E+00  6.0458E+01  8.2108E-06  1.7303E-01  0.0000E+00
             2.7701E+01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -809.325208397706        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     3430
 NPARAMETR:  1.1429E+00  8.4381E-01  1.7786E+01  1.3867E+00  2.3347E+00  1.7518E+00  5.9533E+00  2.0010E-02  8.8025E-01  1.0000E-02
             9.2562E+00
 PARAMETER:  2.3357E-01 -6.9824E-02  2.9784E+00  4.2689E-01  9.4789E-01  6.6067E-01  1.8840E+00 -3.8115E+00 -2.7548E-02 -4.6988E+00
             2.3253E+00
 GRADIENT:   1.2863E-01  6.9887E-03 -2.3249E-03 -3.4023E-01 -1.1370E-03  1.5182E-01  3.8917E+00  6.8386E-06  9.2804E-03  0.0000E+00
            -2.0084E-01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -809.325387526889        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     3619             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1429E+00  8.4269E-01  1.7920E+01  1.3878E+00  2.3351E+00  1.7518E+00  5.9572E+00  2.1427E-02  8.8049E-01  1.0000E-02
             9.2561E+00
 PARAMETER:  2.3358E-01 -7.1157E-02  2.9859E+00  4.2772E-01  9.4804E-01  6.6064E-01  1.8846E+00 -3.7431E+00 -2.7282E-02 -4.6988E+00
             2.3253E+00
 GRADIENT:   1.5354E+01  6.3310E-01  1.0802E-02  1.1757E+01  9.1352E-01  9.7399E+00  6.0539E+01  8.1697E-06  1.1680E-01  0.0000E+00
             2.7561E+01

0ITERATION NO.:  114    OBJECTIVE VALUE:  -809.325426195412        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:     3755
 NPARAMETR:  1.1429E+00  8.4174E-01  1.8079E+01  1.3892E+00  2.3356E+00  1.7518E+00  5.9615E+00  2.1284E-02  8.8108E-01  1.0000E-02
             9.2559E+00
 PARAMETER:  2.3358E-01 -7.1701E-02  2.9855E+00  4.2772E-01  9.4825E-01  6.6063E-01  1.8848E+00 -3.7877E+00 -2.6396E-02 -4.6988E+00
             2.3253E+00
 GRADIENT:   7.3386E-04  8.6741E-03 -1.7001E-03 -1.9135E-01 -2.2931E-04  5.2477E-04 -3.4763E-02 -1.0454E-05  2.8990E-03  0.0000E+00
             2.5293E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3755
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7342E-02  3.6511E-02  6.3248E-06 -6.7323E-02  8.2006E-06
 SE:             2.8493E-02  2.3432E-02  5.0617E-06  1.2364E-02  5.5689E-05
 N:                     100         100         100         100         100

 P VAL.:         5.4277E-01  1.1918E-01  2.1146E-01  5.1934E-08  8.8293E-01

 ETASHRINKSD(%)  4.5452E+00  2.1501E+01  9.9983E+01  5.8578E+01  9.9813E+01
 ETASHRINKVR(%)  8.8838E+00  3.8379E+01  1.0000E+02  8.2842E+01  1.0000E+02
 EBVSHRINKSD(%)  6.8395E+00  1.6083E+01  9.9974E+01  6.1413E+01  9.9747E+01
 EBVSHRINKVR(%)  1.3211E+01  2.9579E+01  1.0000E+02  8.5110E+01  9.9999E+01
 RELATIVEINF(%)  8.4798E+01  3.4744E+01  9.0671E-07  6.6353E+00  8.1151E-05
 EPSSHRINKSD(%)  8.7418E+00
 EPSSHRINKVR(%)  1.6719E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -809.32542619541243     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       293.40081365019466     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    95.16
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.35
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -809.325       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.14E+00  8.42E-01  1.79E+01  1.39E+00  2.34E+00  1.75E+00  5.96E+00  2.05E-02  8.81E-01  1.00E-02  9.26E+00
 


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
+        2.25E+02
 
 TH 2
+       -2.94E+00  3.14E+01
 
 TH 3
+       -9.88E-04  6.63E-03  8.96E-04
 
 TH 4
+       -1.22E+01  3.49E+01 -2.85E-04  1.56E+02
 
 TH 5
+       -8.76E-01 -3.99E+00 -8.24E-02 -1.11E+01  1.22E+01
 
 TH 6
+       -6.74E+00 -6.39E-01  2.46E-03  3.99E+00 -2.80E-01  4.71E+01
 
 TH 7
+        1.29E+00  3.53E+00 -9.40E-04 -8.92E+00  3.38E-01 -2.03E-01  2.85E+00
 
 TH 8
+        1.10E-01  2.47E-01 -1.64E-04  7.87E-03 -1.40E-02 -2.09E-02  8.79E-04 -2.93E-01
 
 TH 9
+        1.72E+00 -2.91E+00 -1.44E-02 -4.30E+01  3.64E+00 -2.06E+00  2.26E+00  2.57E-01  2.73E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -6.97E+00 -2.92E+00  2.38E-03 -1.03E+01  3.95E-01  2.18E+00  2.35E-01 -1.96E-03  3.09E+00  0.00E+00  7.73E+00
 
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
 #CPUT: Total CPU Time in Seconds,      107.542
Stop Time:
Thu Sep 30 09:10:40 CDT 2021
