Sat Sep 18 15:13:52 CDT 2021
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
$DATA ../../../../data/spa/D/dat26.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m26.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   8014.29815716020        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -5.3023E+01  4.5945E+00 -1.1600E+02 -7.9406E+01  2.0862E+02 -1.0444E+03 -4.1069E+02 -2.3059E+01 -7.8149E+02 -2.9428E+02
            -1.6838E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -706.199645493923        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.7382E+00  1.3082E+00  1.0622E+00  1.8239E+00  1.0460E+00  1.8715E+00  1.3193E+00  9.7705E-01  1.4879E+00  1.1549E+00
             1.3950E+01
 PARAMETER:  6.5287E-01  3.6867E-01  1.6036E-01  7.0097E-01  1.4499E-01  7.2674E-01  3.7707E-01  7.6781E-02  4.9735E-01  2.4403E-01
             2.7355E+00
 GRADIENT:   8.7549E+01  2.8974E+01 -1.5212E+00  4.6224E+01 -1.7876E+01  3.8601E+01 -3.3112E+00  3.9376E+00  4.9156E+00  5.0028E+00
             1.5916E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -731.769816836716        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.5866E+00  1.0101E+00  4.4439E+00  2.2198E+00  2.9849E+00  1.6780E+00  3.6657E+00  5.6004E-01  1.7485E+00  6.6938E+00
             1.1857E+01
 PARAMETER:  5.6161E-01  1.1010E-01  1.5915E+00  8.9742E-01  1.1936E+00  6.1760E-01  1.3990E+00 -4.7974E-01  6.5873E-01  2.0012E+00
             2.5729E+00
 GRADIENT:   8.5639E+01  2.1277E+01  4.1204E+00  5.1823E+01 -2.5109E+01 -1.2762E+01  1.1087E+01  7.1780E-03  2.5626E+01  2.8979E+01
             1.3111E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -797.732518275693        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.1726E+00  5.6023E-01  1.7254E+01  1.7661E+00  7.0396E+00  1.5665E+00  2.9036E+00  3.8265E+00  1.2537E+00  5.8195E+00
             8.4547E+00
 PARAMETER:  2.5923E-01 -4.7940E-01  2.9480E+00  6.6879E-01  2.0516E+00  5.4886E-01  1.1659E+00  1.4420E+00  3.2614E-01  1.8612E+00
             2.2347E+00
 GRADIENT:  -1.9686E+01  1.9840E+00  1.7122E+00  3.4608E+00 -4.8402E-01  4.8287E+00  2.1600E+00 -1.2266E-01  3.1058E+00 -6.0874E-01
             8.7431E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -805.978743037541        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  1.2026E+00  6.4140E-01  1.4796E+00  1.5082E+00  3.2354E+00  1.3830E+00  1.7650E+00  2.7660E-01  9.7758E-01  4.3457E+00
             8.6030E+00
 PARAMETER:  2.8451E-01 -3.4410E-01  4.9178E-01  5.1092E-01  1.2742E+00  4.2423E-01  6.6813E-01 -1.1852E+00  7.7322E-02  1.5692E+00
             2.2521E+00
 GRADIENT:   3.8948E+01  5.5956E+00  5.8307E+00 -1.2924E+01 -9.6241E+00 -9.8385E+00 -1.3119E+00 -1.0056E-02 -8.9936E+00 -2.4663E+00
            -1.7557E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -825.585694960578        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      374
 NPARAMETR:  1.0638E+00  3.0475E-01  5.7179E-01  1.5402E+00  5.3664E+00  1.3117E+00  8.6976E-01  1.3793E-02  8.3791E-01  7.0239E+00
             8.5834E+00
 PARAMETER:  1.6186E-01 -1.0883E+00 -4.5898E-01  5.3192E-01  1.7802E+00  3.7131E-01 -3.9541E-02 -4.1836E+00 -7.6850E-02  2.0493E+00
             2.2498E+00
 GRADIENT:  -2.0323E+01  2.1958E+01 -2.0383E+01  8.5726E+01 -2.2151E+01 -3.0345E+01  3.1788E-01 -6.8748E-04 -7.2228E+00  1.4071E+01
            -5.8580E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -891.777439590976        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      446
 NPARAMETR:  6.3590E-01  3.9711E-02  5.9496E-02  5.7390E-01  1.1888E+01  1.3294E+00  4.6155E-02  1.0000E-02  3.7745E-01  2.2351E+00
             8.3104E+00
 PARAMETER: -3.5271E-01 -3.1261E+00 -2.7219E+00 -4.5530E-01  2.5755E+00  3.8474E-01 -2.9758E+00 -1.2274E+01 -8.7430E-01  9.0428E-01
             2.2175E+00
 GRADIENT:   9.8987E+00 -3.5019E+00 -7.2136E+01  1.3494E+02  1.4686E+01 -1.5570E+00  3.3514E-03  0.0000E+00 -1.2076E+00 -1.3997E+01
            -3.4680E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -899.261633744401        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      517
 NPARAMETR:  4.8445E-01  1.9762E-02  2.8280E-02  3.2603E-01  1.0223E+01  1.2315E+00  1.0000E-02  1.0000E-02  2.3560E-01  1.8509E+00
             8.2950E+00
 PARAMETER: -6.2474E-01 -3.8240E+00 -3.4656E+00 -1.0208E+00  2.4246E+00  3.0827E-01 -4.5546E+00 -1.4819E+01 -1.3456E+00  7.1566E-01
             2.2157E+00
 GRADIENT:   2.4024E+01  3.9134E+00 -4.1169E+01  5.4779E+01  1.3438E+01 -8.4821E+00  0.0000E+00  0.0000E+00 -2.4397E-01 -1.1333E+01
            -5.6088E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -899.827689753598        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      590
 NPARAMETR:  4.0845E-01  1.1828E-02  1.8634E-02  2.3508E-01  9.0732E+00  1.2072E+00  1.0000E-02  1.0000E-02  1.7966E-01  1.6028E+00
             8.3517E+00
 PARAMETER: -7.9540E-01 -4.3373E+00 -3.8828E+00 -1.3478E+00  2.3053E+00  2.8833E-01 -5.5111E+00 -1.6471E+01 -1.6167E+00  5.7172E-01
             2.2225E+00
 GRADIENT:   2.2867E+01  1.4566E+00 -4.0732E+01  5.5008E+01  1.2477E+01 -8.7144E+00  0.0000E+00  0.0000E+00 -2.4007E-01 -9.8490E+00
            -2.7073E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -901.944994452828        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      769
 NPARAMETR:  4.7205E-01  1.6412E-02  2.7898E-02  3.1578E-01  9.5476E+00  1.2835E+00  1.0000E-02  1.0000E-02  2.5193E-01  1.9437E+00
             8.2961E+00
 PARAMETER: -6.5068E-01 -4.0097E+00 -3.4792E+00 -1.0527E+00  2.3563E+00  3.4959E-01 -4.7675E+00 -1.5207E+01 -1.2786E+00  7.6459E-01
             2.2158E+00
 GRADIENT:  -3.6017E+00  1.9597E+00  1.0583E+00 -7.8166E-02  1.1060E+01  5.2543E+00  0.0000E+00  0.0000E+00 -4.8867E-01 -9.6246E+00
             5.1033E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -904.737827121090        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      948
 NPARAMETR:  5.2039E-01  1.4071E-02  3.6893E-02  3.8111E-01  1.0572E+01  1.3115E+00  1.0000E-02  1.0000E-02  4.9036E-01  3.7992E+00
             7.9055E+00
 PARAMETER: -5.5318E-01 -4.1637E+00 -3.1997E+00 -8.6467E-01  2.4582E+00  3.7117E-01 -4.8965E+00 -1.4832E+01 -6.1261E-01  1.4348E+00
             2.1676E+00
 GRADIENT:  -9.1592E+00 -2.3318E-01  3.6835E+01 -4.1149E+01  4.9063E+00  1.5065E+00  0.0000E+00  0.0000E+00  7.0131E-03 -1.4115E+00
             7.0147E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -905.727850240910        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1125
 NPARAMETR:  4.6534E-01  1.0000E-02  2.6990E-02  3.0596E-01  1.0141E+01  1.2892E+00  1.0000E-02  1.0000E-02  4.9809E-01  3.8827E+00
             7.8170E+00
 PARAMETER: -6.6498E-01 -4.6140E+00 -3.5123E+00 -1.0843E+00  2.4166E+00  3.5403E-01 -5.7245E+00 -1.6115E+01 -5.9698E-01  1.4565E+00
             2.1563E+00
 GRADIENT:  -4.0169E+00  0.0000E+00  2.5191E+00 -1.8696E+00  3.5847E+00 -1.1363E-01  0.0000E+00  0.0000E+00 -1.0181E+00  1.0456E-01
            -1.5464E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -905.892589507172        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1301
 NPARAMETR:  4.7611E-01  1.0000E-02  2.8617E-02  3.1932E-01  9.7032E+00  1.2963E+00  1.0000E-02  1.0000E-02  5.3751E-01  3.5992E+00
             7.7538E+00
 PARAMETER: -6.4210E-01 -4.5908E+00 -3.4538E+00 -1.0416E+00  2.3725E+00  3.5953E-01 -5.5613E+00 -1.5788E+01 -5.2081E-01  1.3807E+00
             2.1482E+00
 GRADIENT:  -2.9298E+00  0.0000E+00  5.9400E+00 -6.1292E+00  1.3443E+00 -4.9603E-01  0.0000E+00  0.0000E+00 -7.3982E-02  3.7731E-01
            -1.7385E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -905.950975910671        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1478
 NPARAMETR:  4.7897E-01  1.0000E-02  2.8833E-02  3.2197E-01  9.3984E+00  1.3001E+00  1.0000E-02  1.0000E-02  5.4867E-01  3.2875E+00
             7.7590E+00
 PARAMETER: -6.3611E-01 -4.5667E+00 -3.4462E+00 -1.0333E+00  2.3405E+00  3.6247E-01 -5.4826E+00 -1.5659E+01 -5.0025E-01  1.2901E+00
             2.1489E+00
 GRADIENT:  -4.2758E-02  0.0000E+00  8.9970E-02  3.1507E-01  6.1727E-01 -5.6997E-02  0.0000E+00  0.0000E+00  2.5398E-01  4.9484E-05
            -3.7639E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -905.956896968582        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1653
 NPARAMETR:  4.7646E-01  1.0000E-02  2.8425E-02  3.1847E-01  9.3485E+00  1.2989E+00  1.0000E-02  1.0000E-02  5.4256E-01  3.2012E+00
             7.7731E+00
 PARAMETER: -6.4136E-01 -4.5611E+00 -3.4605E+00 -1.0442E+00  2.3352E+00  3.6152E-01 -5.4936E+00 -1.5683E+01 -5.1146E-01  1.2635E+00
             2.1507E+00
 GRADIENT:   1.2247E-01  0.0000E+00 -7.1064E-01  7.9918E-01 -1.0676E+00 -1.0581E-02  0.0000E+00  0.0000E+00 -2.8622E-01  3.4548E-01
            -6.8808E-02

0ITERATION NO.:   74    OBJECTIVE VALUE:  -905.957947483985        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1780
 NPARAMETR:  4.7500E-01  1.0000E-02  2.8197E-02  3.1648E-01  9.3383E+00  1.2984E+00  1.0000E-02  1.0000E-02  5.4607E-01  3.1886E+00
             7.7661E+00
 PARAMETER: -6.4444E-01 -4.5649E+00 -3.4685E+00 -1.0505E+00  2.3341E+00  3.6117E-01 -5.5110E+00 -1.5702E+01 -5.0501E-01  1.2596E+00
             2.1498E+00
 GRADIENT:  -8.0343E-02  0.0000E+00  4.1937E-01 -3.6862E-01  8.6012E-01  4.7406E-02  0.0000E+00  0.0000E+00  1.2402E-01 -2.7411E-01
             9.4012E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1780
 NO. OF SIG. DIGITS IN FINAL EST.:  3.8
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.9827E-03  3.6163E-06  1.3626E-04 -1.5010E-02 -9.5602E-03
 SE:             2.8852E-02  1.9127E-06  2.9229E-04  1.9643E-02  5.8588E-03
 N:                     100         100         100         100         100

 P VAL.:         9.4521E-01  5.8671E-02  6.4109E-01  4.4477E-01  1.0273E-01

 ETASHRINKSD(%)  3.3407E+00  9.9994E+01  9.9021E+01  3.4195E+01  8.0372E+01
 ETASHRINKVR(%)  6.5698E+00  1.0000E+02  9.9990E+01  5.6697E+01  9.6148E+01
 EBVSHRINKSD(%)  3.4003E+00  9.9993E+01  9.9054E+01  3.4778E+01  8.2889E+01
 EBVSHRINKVR(%)  6.6849E+00  1.0000E+02  9.9991E+01  5.7461E+01  9.7072E+01
 RELATIVEINF(%)  4.1900E+00  4.0467E-08  5.8397E-05  2.9128E-01  4.8109E-01
 EPSSHRINKSD(%)  1.3886E+01
 EPSSHRINKVR(%)  2.5844E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -905.95794748398487     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -170.80712092024669     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.11
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.61
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -905.958       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.75E-01  1.00E-02  2.82E-02  3.16E-01  9.34E+00  1.30E+00  1.00E-02  1.00E-02  5.46E-01  3.19E+00  7.77E+00
 


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
+        1.92E+06
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.13E+04  0.00E+00  1.97E+07
 
 TH 4
+       -1.04E+03  0.00E+00 -9.86E+04  1.64E+06
 
 TH 5
+        7.75E+00  0.00E+00 -4.78E+02  5.79E+01  3.77E+02
 
 TH 6
+       -1.04E+02  0.00E+00 -1.13E+02  7.04E+01 -1.43E+00  1.16E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -2.13E+06  0.00E+00 -7.50E+02  5.18E+02 -5.70E+00  2.41E+02  0.00E+00  0.00E+00  2.36E+06
 
 TH10
+       -4.89E+01  0.00E+00 -6.50E+02  1.35E+05 -3.44E+00  1.04E+01  0.00E+00  0.00E+00  4.78E+01  1.12E+04
 
 TH11
+       -1.09E+01  0.00E+00  4.27E+02 -3.13E+01  1.84E-01 -1.18E+00  0.00E+00  0.00E+00 -2.07E+00  4.70E-02  6.52E+02
 
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
 #CPUT: Total CPU Time in Seconds,       33.787
Stop Time:
Sat Sep 18 15:14:28 CDT 2021
