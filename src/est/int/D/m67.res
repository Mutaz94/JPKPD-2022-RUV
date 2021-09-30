Wed Sep 29 09:29:21 CDT 2021
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
$DATA ../../../../data/int/D/dat67.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m67.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   30118.6691446547        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.7130E+02  4.5150E+02 -3.6369E+01  7.8284E+01  4.1574E+02 -2.8163E+03 -1.0786E+03 -8.3382E+01 -2.2685E+03 -9.7688E+02
            -6.0269E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -977.831849536594        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.9910E+00  1.5524E+00  9.4307E-01  2.8372E+00  7.1556E-01  5.0068E+00  4.4164E+00  9.7907E-01  4.9562E+00  2.6197E+00
             1.0747E+01
 PARAMETER:  7.8865E-01  5.3977E-01  4.1382E-02  1.1428E+00 -2.3469E-01  1.7108E+00  1.5853E+00  7.8848E-02  1.7006E+00  1.0631E+00
             2.4747E+00
 GRADIENT:   6.7178E+01  2.9489E+01 -3.8143E+01  6.6076E+01 -1.4110E+01  1.7613E+02  7.2726E+01  4.9709E+00  1.0924E+02  1.7804E+01
             3.2589E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1080.96820937938        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  3.8061E+00  3.9822E-01  2.4084E+01  3.8379E+00  2.2358E+00  3.8770E+00  9.6877E+00  2.6807E-01  3.5932E+00  1.4747E+00
             1.1069E+01
 PARAMETER:  1.4366E+00 -8.2075E-01  3.2816E+00  1.4449E+00  9.0458E-01  1.4551E+00  2.3709E+00 -1.2165E+00  1.3790E+00  4.8847E-01
             2.5042E+00
 GRADIENT:   2.1511E+02  6.6873E+00 -1.8684E+01  1.1342E+02  9.0684E+00  2.5136E+01  2.0103E+01  1.3197E-01  5.7017E+01  3.4357E+01
             3.4777E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1193.57094152175        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.2422E+00  5.9760E-01  4.1150E+01  1.5671E+00  2.1497E+00  2.6616E+00  6.2497E+00  1.0000E-02  1.3581E+00  9.3996E-01
             1.1456E+01
 PARAMETER:  3.1691E-01 -4.1483E-01  3.8172E+00  5.4925E-01  8.6532E-01  1.0789E+00  1.9325E+00 -5.5258E+00  4.0605E-01  3.8084E-02
             2.5385E+00
 GRADIENT:  -6.5282E+00 -2.3952E+01  1.9892E+00  7.3510E+00 -4.2125E+01  1.9470E+01  9.7326E+01  0.0000E+00  4.2851E+00  1.3831E+01
             4.5529E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1238.34407523572        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  9.2036E-01  1.4976E+00  7.4463E+01  9.6280E-01  2.6228E+00  2.6212E+00  3.7933E+00  1.4175E-02  1.3648E+00  9.5227E-01
             8.4815E+00
 PARAMETER:  1.7005E-02  5.0390E-01  4.4103E+00  6.2092E-02  1.0643E+00  1.0636E+00  1.4332E+00 -4.1563E+00  4.1102E-01  5.1095E-02
             2.2379E+00
 GRADIENT:  -7.4763E+01 -2.9972E+01 -1.5464E-01  1.0032E+01  2.3245E+01  1.6660E+01 -3.1781E+01  3.8001E-07  2.0052E+00  1.5728E+01
             6.1846E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1258.45397438462        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  1.1370E+00  2.0261E+00  5.7382E+01  6.8926E-01  2.4571E+00  2.6184E+00  3.5690E+00  8.9849E-02  1.3205E+00  7.9185E-01
             7.9314E+00
 PARAMETER:  2.2836E-01  8.0613E-01  4.1497E+00 -2.7213E-01  9.9899E-01  1.0626E+00  1.3723E+00 -2.3096E+00  3.7802E-01 -1.3339E-01
             2.1708E+00
 GRADIENT:  -1.0893E+01  2.0169E-01  8.0963E-02  4.6817E+00  4.7165E+00 -7.2904E-01 -5.2148E+00  8.6507E-06  3.2914E+00  1.1212E+01
            -5.1474E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1266.52956226667        NO. OF FUNC. EVALS.: 132
 CUMULATIVE NO. OF FUNC. EVALS.:      504
 NPARAMETR:  1.1362E+00  2.0351E+00  9.6637E+00  6.8982E-01  2.4670E+00  2.6111E+00  3.6098E+00  8.6259E-02  1.3174E+00  6.2430E-02
             8.2312E+00
 PARAMETER:  2.2767E-01  8.1054E-01  2.3684E+00 -2.7132E-01  1.0030E+00  1.0598E+00  1.3837E+00 -2.3504E+00  3.7567E-01 -2.6737E+00
             2.2079E+00
 GRADIENT:  -1.5511E+01 -2.5296E+00 -8.1021E-01  1.5354E+00  2.2702E+01 -7.4434E+00 -1.5983E+00  1.9505E-04  4.1354E+00  6.9154E-02
             2.1724E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1281.26591221493        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      666
 NPARAMETR:  1.1386E+00  2.0952E+00  1.3244E+03  6.8978E-01  2.4359E+00  2.8631E+00  4.7065E+00  1.6660E-02  1.3154E+00  1.0000E-02
             8.5727E+00
 PARAMETER:  2.2981E-01  8.3964E-01  7.2888E+00 -2.7138E-01  9.9030E-01  1.1519E+00  1.6489E+00 -3.9947E+00  3.7415E-01 -8.9557E+00
             2.2486E+00
 GRADIENT:  -2.8666E+01 -9.9545E+00  7.1451E-03 -2.7825E+01 -1.3603E+01 -4.2733E+01  1.3715E+01 -2.6666E-09  7.8494E+00  0.0000E+00
             7.2363E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1287.40009358333        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:      832
 NPARAMETR:  1.1379E+00  2.0607E+00  1.1032E+03  8.1801E-01  2.4438E+00  2.8597E+00  4.7408E+00  1.7363E-02  6.9149E-01  1.0000E-02
             8.5343E+00
 PARAMETER:  2.2922E-01  8.2306E-01  7.1059E+00 -1.0089E-01  9.9357E-01  1.1507E+00  1.6562E+00 -3.9534E+00 -2.6890E-01 -8.7455E+00
             2.2441E+00
 GRADIENT:  -2.8943E+01  1.7081E+00  1.5618E-02 -1.6414E+00 -1.9448E+01 -4.2628E+01 -1.3348E+01  4.5927E-08  4.1412E-01  0.0000E+00
             3.4606E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1288.80580870130        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1009
 NPARAMETR:  1.1382E+00  2.0300E+00  2.3922E+01  8.3746E-01  2.4560E+00  2.8880E+00  4.7802E+00  1.7167E-02  7.4577E-01  1.0000E-02
             8.3736E+00
 PARAMETER:  2.2943E-01  8.0806E-01  3.2748E+00 -7.7384E-02  9.9855E-01  1.1606E+00  1.6645E+00 -3.9648E+00 -1.9333E-01 -8.7455E+00
             2.2251E+00
 GRADIENT:  -2.6952E+01  2.3982E+00  9.0203E-02  3.4350E-01 -2.6392E+00 -3.8450E+01 -1.3800E+01  3.5781E-06  2.9983E-01  0.0000E+00
            -2.6906E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1292.48336288937        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1188
 NPARAMETR:  1.1414E+00  1.7690E+00  5.3332E+01  9.3779E-01  2.5026E+00  3.2002E+00  5.1829E+00  1.1980E-02  8.6630E-01  1.0000E-02
             8.3734E+00
 PARAMETER:  2.3229E-01  6.7043E-01  4.0765E+00  3.5767E-02  1.0173E+00  1.2632E+00  1.7454E+00 -4.3245E+00 -4.3525E-02 -8.7455E+00
             2.2251E+00
 GRADIENT:  -2.0759E+01  2.7468E-01  3.8558E-02  3.5627E-02 -2.6952E-01  5.7106E+00 -8.4716E+00  4.4254E-07 -9.5092E-02  0.0000E+00
            -5.6120E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1293.02909061494        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:     1383
 NPARAMETR:  1.1641E+00  1.7129E+00  5.1045E+01  9.6264E-01  2.5127E+00  3.2347E+00  5.2921E+00  1.1412E-02  9.0099E-01  1.0000E-02
             8.4081E+00
 PARAMETER:  2.5198E-01  6.3818E-01  4.0327E+00  6.1923E-02  1.0214E+00  1.2739E+00  1.7662E+00 -4.3731E+00 -4.2644E-03 -8.7455E+00
             2.2292E+00
 GRADIENT:  -1.6512E+01 -1.8281E-01 -1.2687E-02 -1.1393E+00  2.4911E+00  1.0482E+01 -6.3321E+00  4.5787E-07  1.7513E-01  0.0000E+00
             7.1300E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1293.11891405137        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1563
 NPARAMETR:  1.1651E+00  1.6190E+00  5.2374E+01  9.9927E-01  2.5079E+00  3.2018E+00  5.3938E+00  1.1336E-02  9.8829E-01  1.0000E-02
             8.3123E+00
 PARAMETER:  2.5278E-01  5.8181E-01  4.0584E+00  9.9271E-02  1.0194E+00  1.2637E+00  1.7852E+00 -4.3798E+00  8.8223E-02 -8.7455E+00
             2.2177E+00
 GRADIENT:  -1.6017E+01 -7.8883E-01 -8.8031E-03 -5.0597E-01  2.4705E+00  6.4792E+00 -6.8745E+00  4.2980E-07  2.1084E-01  0.0000E+00
            -1.4164E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1293.18153147058        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1740
 NPARAMETR:  1.1679E+00  1.5436E+00  4.1121E+01  1.0301E+00  2.4909E+00  3.1384E+00  5.5536E+00  1.1147E-02  1.0390E+00  1.0000E-02
             8.3237E+00
 PARAMETER:  2.5519E-01  5.3412E-01  3.8165E+00  1.2966E-01  1.0126E+00  1.2437E+00  1.8144E+00 -4.3965E+00  1.3828E-01 -8.7455E+00
             2.2191E+00
 GRADIENT:  -1.6428E+01 -1.1579E+00 -4.2742E-02 -2.5721E+00  2.8341E+00 -1.5333E+00 -3.4307E+00  6.5217E-07  7.0241E-01  0.0000E+00
            -1.0311E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1293.19588521446        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1915
 NPARAMETR:  1.1694E+00  1.5292E+00  3.9422E+01  1.0454E+00  2.4816E+00  3.1122E+00  5.6241E+00  1.1043E-02  1.0331E+00  1.0000E-02
             8.3584E+00
 PARAMETER:  2.5652E-01  5.2472E-01  3.7743E+00  1.4443E-01  1.0089E+00  1.2353E+00  1.8271E+00 -4.4060E+00  1.3259E-01 -8.7455E+00
             2.2233E+00
 GRADIENT:  -1.6784E+01 -4.4789E-01 -1.7339E-02 -1.1068E+00  1.4082E+00 -4.9029E+00 -2.4507E+00  6.6831E-07  2.7115E-01  0.0000E+00
            -3.6200E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1293.20088651019        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2090
 NPARAMETR:  1.1707E+00  1.5209E+00  3.8311E+01  1.0580E+00  2.4748E+00  3.1028E+00  5.6764E+00  1.0940E-02  1.0285E+00  1.0000E-02
             8.3872E+00
 PARAMETER:  2.5763E-01  5.1927E-01  3.7457E+00  1.5635E-01  1.0062E+00  1.2323E+00  1.8363E+00 -4.4154E+00  1.2810E-01 -8.7455E+00
             2.2267E+00
 GRADIENT:  -1.6917E+01  2.1051E-01  5.9658E-03  4.5779E-01  2.4413E-01 -6.0755E+00 -1.9062E+00  7.4563E-07 -1.4180E-01  0.0000E+00
             1.7665E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1293.49697131885        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2268
 NPARAMETR:  1.1775E+00  1.3975E+00  3.3164E+01  1.1224E+00  2.4506E+00  3.1485E+00  5.9714E+00  1.0216E-02  1.1099E+00  1.0000E-02
             8.3894E+00
 PARAMETER:  2.6336E-01  4.3469E-01  3.6015E+00  2.1549E-01  9.9634E-01  1.2469E+00  1.8870E+00 -4.4838E+00  2.0428E-01 -8.7455E+00
             2.2270E+00
 GRADIENT:  -1.5104E+01  4.3300E-01  1.9377E-02  6.8426E-01 -9.2815E-01  9.0149E-02  1.5601E+00  6.3942E-07 -1.7374E-01  0.0000E+00
             3.1516E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1294.10509823334        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2447
 NPARAMETR:  1.2616E+00  1.3754E+00  3.2556E+01  1.1275E+00  2.4578E+00  3.2224E+00  5.9311E+00  1.0000E-02  1.1272E+00  1.0000E-02
             8.4135E+00
 PARAMETER:  3.3382E-01  4.2162E-01  3.5916E+00  2.2048E-01  9.9687E-01  1.2630E+00  1.8897E+00 -4.5349E+00  2.2116E-01 -8.7455E+00
             2.2243E+00
 GRADIENT:   1.4654E-01  7.5776E-02  3.8785E-03  7.5792E-02 -6.1830E-01 -1.3526E+00  1.1589E+00  5.1395E-06  2.3560E-02  0.0000E+00
            -5.6025E+00

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2447
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.4374E-03  2.5825E-02 -9.8328E-07 -6.7370E-02 -2.9261E-06
 SE:             2.9339E-02  2.5303E-02  3.9169E-06  1.2706E-02  1.2926E-04
 N:                     100         100         100         100         100

 P VAL.:         7.4770E-01  3.0744E-01  8.0179E-01  1.1466E-07  9.8194E-01

 ETASHRINKSD(%)  1.7106E+00  1.5231E+01  9.9987E+01  5.7433E+01  9.9567E+01
 ETASHRINKVR(%)  3.3919E+00  2.8142E+01  1.0000E+02  8.1880E+01  9.9998E+01
 EBVSHRINKSD(%)  2.0038E+00  9.0872E+00  9.9984E+01  6.4440E+01  9.9495E+01
 EBVSHRINKVR(%)  3.9675E+00  1.7349E+01  1.0000E+02  8.7355E+01  9.9997E+01
 RELATIVEINF(%)  9.5980E+01  4.3541E+01  6.4893E-07  6.7276E+00  6.1343E-04
 EPSSHRINKSD(%)  6.9681E+00
 EPSSHRINKVR(%)  1.3451E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1294.1050982333418     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       359.98426153506898     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    78.32
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    17.23
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1294.105       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.26E+00  1.38E+00  3.28E+01  1.13E+00  2.45E+00  3.20E+00  5.99E+00  1.00E-02  1.13E+00  1.00E-02  8.37E+00
 


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
+        6.53E+01
 
 TH 2
+       -2.68E-01  1.38E+01
 
 TH 3
+        1.59E-03  5.14E-03  4.19E-04
 
 TH 4
+       -8.49E-01  1.98E+01  1.47E-02  1.22E+02
 
 TH 5
+       -9.34E-01 -3.13E+00 -1.15E-01 -6.48E+00  4.31E+01
 
 TH 6
+        1.42E+00 -9.19E-02  6.36E-04  3.11E-01 -6.67E-01  1.85E+01
 
 TH 7
+       -3.36E-02  1.84E+00 -3.05E-03 -1.06E+01  7.19E-01 -1.37E-01  3.41E+00
 
 TH 8
+       -6.36E-02 -5.07E-02 -8.40E-04 -1.97E-01 -1.20E-02  1.20E-02 -5.94E-03  3.64E-01
 
 TH 9
+        9.62E-02 -1.75E+00 -3.90E-03 -2.85E+01  2.85E+00 -7.79E-02  2.08E+00  7.53E-02  1.32E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -3.28E+00 -1.64E+00 -3.46E-04 -8.24E+00  2.71E-01 -6.16E+01  4.61E-01 -1.60E-03  3.34E+00  0.00E+00  1.44E+01
 
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
 #CPUT: Total CPU Time in Seconds,       95.650
Stop Time:
Wed Sep 29 09:30:58 CDT 2021
