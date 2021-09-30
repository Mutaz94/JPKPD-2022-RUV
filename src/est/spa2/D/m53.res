Thu Sep 30 09:15:46 CDT 2021
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
$DATA ../../../../data/spa2/D/dat53.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m53.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   23041.6061120938        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.1596E+02  3.0807E+02  1.9890E+01  3.2344E+02  1.1536E+02 -1.9092E+03 -1.1344E+03 -7.0772E+01 -1.2069E+03 -6.0583E+02
            -4.5580E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -598.985987613961        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.3579E+00  1.2728E+00  9.4710E-01  1.4723E+00  9.9225E-01  2.0513E+00  1.4454E+00  9.9369E-01  1.2311E+00  1.1883E+00
             1.4298E+01
 PARAMETER:  4.0592E-01  3.4120E-01  4.5654E-02  4.8683E-01  9.2216E-02  8.1848E-01  4.6838E-01  9.3673E-02  3.0794E-01  2.7253E-01
             2.7601E+00
 GRADIENT:   1.6490E+01 -1.4416E+01 -8.8939E+00  1.6927E+01  1.4131E+01  5.0860E+01 -2.3349E+01  3.3159E+00 -6.8397E+00  1.3944E+01
             2.2177E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -692.610754033061        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.3386E+00  1.3674E+00  4.1898E+00  1.6064E+00  1.9839E+00  1.9886E+00  4.7354E+00  4.2195E-01  1.8407E+00  1.9922E+00
             1.3024E+01
 PARAMETER:  3.9161E-01  4.1289E-01  1.5327E+00  5.7402E-01  7.8506E-01  7.8745E-01  1.6551E+00 -7.6288E-01  7.1013E-01  7.8923E-01
             2.6668E+00
 GRADIENT:   3.1442E+01  2.1483E+00 -6.3599E+00  2.4455E+00 -1.3033E+01  1.1060E+01  2.2414E+01  8.2687E-02  3.0900E+01  2.9367E+01
             2.5900E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -739.661111995746        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.0231E+00  1.3350E+00  1.0541E+01  1.0858E+00  3.0732E+00  2.0945E+00  3.7792E+00  8.1840E+00  7.3531E-01  1.1095E+00
             9.6692E+00
 PARAMETER:  1.2280E-01  3.8892E-01  2.4553E+00  1.8229E-01  1.2227E+00  8.3929E-01  1.4295E+00  2.2022E+00 -2.0747E-01  2.0395E-01
             2.3689E+00
 GRADIENT:  -4.6310E+01 -2.3948E+01 -1.8580E+00 -1.7631E+01  1.6035E+01  2.5501E+01 -2.3106E+01  1.1265E+00  3.9800E+00  5.4379E+00
             6.9327E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -756.555444511491        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      319
 NPARAMETR:  1.1199E+00  1.5403E+00  2.5945E+00  9.8534E-01  1.8107E+00  1.9127E+00  4.0263E+00  4.1167E-01  3.8293E-01  5.2458E-01
             9.2728E+00
 PARAMETER:  2.1328E-01  5.3196E-01  1.0534E+00  8.5227E-02  6.9369E-01  7.4849E-01  1.4928E+00 -7.8754E-01 -8.5990E-01 -5.4516E-01
             2.3271E+00
 GRADIENT:  -1.0594E+00 -4.5648E+00  1.3117E+00 -4.6108E+00 -3.3354E+00 -2.2134E+01 -3.1227E+01  1.1035E-01  5.4449E-01  2.7976E+00
            -3.5848E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -762.853287060628        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      496
 NPARAMETR:  1.1246E+00  1.0596E+00  5.0560E+00  1.2515E+00  2.0588E+00  2.0253E+00  5.4643E+00  2.1889E-01  5.9047E-01  3.1040E-01
             9.5631E+00
 PARAMETER:  2.1742E-01  1.5785E-01  1.7206E+00  3.2431E-01  8.2213E-01  8.0574E-01  1.7982E+00 -1.4192E+00 -4.2684E-01 -1.0699E+00
             2.3579E+00
 GRADIENT:  -4.3930E+00 -2.2306E+00 -3.1308E-01 -5.0171E+00  2.7032E+00 -4.0641E+00  2.8432E+00  1.2730E-02 -1.3801E+00  6.7796E-01
             5.9549E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -763.316820476036        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      671
 NPARAMETR:  1.1385E+00  1.0409E+00  4.9333E+00  1.2923E+00  2.0067E+00  2.0576E+00  5.4753E+00  1.3020E-01  7.1637E-01  2.0252E-01
             9.5303E+00
 PARAMETER:  2.2975E-01  1.4008E-01  1.6960E+00  3.5645E-01  7.9650E-01  8.2153E-01  1.8003E+00 -1.9387E+00 -2.3355E-01 -1.4969E+00
             2.3545E+00
 GRADIENT:   1.5662E+00 -5.5962E-01 -3.7926E-01  5.5574E-01  2.5965E-01  1.2654E+00 -1.4666E-01  5.1755E-03  2.7763E-01  2.9993E-01
            -3.2162E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -763.401240274789        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      849
 NPARAMETR:  1.1294E+00  1.0949E+00  5.4434E+00  1.2792E+00  2.0407E+00  2.0409E+00  5.4119E+00  1.1904E-01  7.0888E-01  1.0043E-01
             9.5544E+00
 PARAMETER:  2.2168E-01  1.9064E-01  1.7944E+00  3.4623E-01  8.1328E-01  8.1340E-01  1.7886E+00 -2.0283E+00 -2.4408E-01 -2.1983E+00
             2.3570E+00
 GRADIENT:  -2.7200E+00  9.3883E-01  8.2193E-02  2.3167E-01 -4.2809E-01 -1.4502E+00  6.4105E-01  3.2494E-03  3.0242E-01  7.2378E-02
             9.0274E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -763.492298620307        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1029
 NPARAMETR:  1.1370E+00  1.0641E+00  5.8651E+00  1.2908E+00  2.0724E+00  2.0554E+00  5.4979E+00  1.0000E-02  7.0511E-01  2.6484E-02
             9.5514E+00
 PARAMETER:  2.2842E-01  1.6212E-01  1.8690E+00  3.5530E-01  8.2872E-01  8.2048E-01  1.8044E+00 -1.6930E+01 -2.4940E-01 -3.5312E+00
             2.3567E+00
 GRADIENT:   4.4284E-01  4.1920E-01  8.9459E-02 -2.0271E-01  1.2710E-01  7.7003E-01  1.7148E+00  0.0000E+00 -1.3584E-01  4.8968E-03
             1.1091E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -763.509545519257        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1209             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1362E+00  1.0492E+00  5.6222E+00  1.2956E+00  2.0523E+00  2.0548E+00  5.5461E+00  1.0000E-02  7.1619E-01  1.0000E-02
             9.5517E+00
 PARAMETER:  2.2769E-01  1.4804E-01  1.8267E+00  3.5897E-01  8.1898E-01  8.2018E-01  1.8131E+00 -1.8547E+01 -2.3381E-01 -5.5453E+00
             2.3567E+00
 GRADIENT:   1.3536E+01  1.4491E+00  5.9745E-02  6.2795E+00  1.0717E+00  2.2364E+01  6.0881E+01  0.0000E+00  4.0431E-01  0.0000E+00
             2.8866E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -763.515067969632        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1389             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1363E+00  1.0398E+00  5.6443E+00  1.2998E+00  2.0521E+00  2.0549E+00  5.5661E+00  1.0000E-02  7.1836E-01  1.0000E-02
             9.5519E+00
 PARAMETER:  2.2778E-01  1.3907E-01  1.8306E+00  3.6220E-01  8.1888E-01  8.2023E-01  1.8167E+00 -1.8547E+01 -2.3078E-01 -5.5417E+00
             2.3567E+00
 GRADIENT:   1.3597E+01  1.2774E+00  5.2696E-02  6.6527E+00  1.0286E+00  2.2410E+01  6.1330E+01  0.0000E+00  3.2605E-01  0.0000E+00
             2.8736E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -763.515535536191        NO. OF FUNC. EVALS.: 173
 CUMULATIVE NO. OF FUNC. EVALS.:     1562
 NPARAMETR:  1.1372E+00  1.0346E+00  5.6765E+00  1.3037E+00  2.0531E+00  2.0551E+00  5.5588E+00  1.0000E-02  7.1918E-01  1.0000E-02
             9.5599E+00
 PARAMETER:  2.2856E-01  1.3398E-01  1.8363E+00  3.6524E-01  8.1934E-01  8.2033E-01  1.8154E+00 -1.8547E+01 -2.2964E-01 -5.5417E+00
             2.3576E+00
 GRADIENT:   4.7143E-01  9.7840E-02 -6.8031E-03  1.1359E-01 -6.2440E-02  8.7411E-01  2.0096E+00  0.0000E+00 -7.4652E-02  0.0000E+00
             8.3811E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -763.519849584558        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1747
 NPARAMETR:  1.1360E+00  1.0265E+00  5.6966E+00  1.3049E+00  2.0544E+00  2.0551E+00  5.5885E+00  1.0000E-02  7.2280E-01  1.0000E-02
             9.5509E+00
 PARAMETER:  2.2755E-01  1.2615E-01  1.8399E+00  3.6616E-01  8.1998E-01  8.2032E-01  1.8207E+00 -1.8547E+01 -2.2462E-01 -5.5417E+00
             2.3566E+00
 GRADIENT:   2.1499E-01 -3.2227E-02 -4.2476E-02 -1.0405E+00  1.8070E-01  8.0216E-01  2.8441E+00  0.0000E+00  1.0069E-01  0.0000E+00
             3.6834E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -763.521641297854        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1934
 NPARAMETR:  1.1356E+00  1.0224E+00  5.7553E+00  1.3089E+00  2.0533E+00  2.0554E+00  5.6007E+00  1.0000E-02  7.2209E-01  1.0000E-02
             9.5467E+00
 PARAMETER:  2.2717E-01  1.2213E-01  1.8501E+00  3.6916E-01  8.1946E-01  8.2045E-01  1.8229E+00 -1.8547E+01 -2.2561E-01 -5.5417E+00
             2.3562E+00
 GRADIENT:   8.1714E-02  1.3900E-01  2.4132E-02  2.7437E-01 -2.6586E-01  8.8472E-01  2.6309E+00  0.0000E+00 -1.6978E-01  0.0000E+00
            -7.4629E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -763.523440756778        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2121             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1356E+00  1.0138E+00  5.8089E+00  1.3117E+00  2.0578E+00  2.0547E+00  5.6143E+00  1.0000E-02  7.2842E-01  1.0000E-02
             9.5481E+00
 PARAMETER:  2.2718E-01  1.1372E-01  1.8594E+00  3.7130E-01  8.2165E-01  8.2012E-01  1.8253E+00 -1.8547E+01 -2.1688E-01 -5.5417E+00
             2.3563E+00
 GRADIENT:   1.3395E+01  7.3765E-01  5.4407E-02  7.4998E+00  9.4954E-01  2.2489E+01  6.2315E+01  0.0000E+00  1.9636E-01  0.0000E+00
             2.8028E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -763.523852482729        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2308             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1358E+00  1.0132E+00  5.8278E+00  1.3128E+00  2.0592E+00  2.0550E+00  5.6196E+00  1.0000E-02  7.2958E-01  1.0000E-02
             9.5501E+00
 PARAMETER:  2.2732E-01  1.1311E-01  1.8626E+00  3.7212E-01  8.2231E-01  8.2030E-01  1.8263E+00 -1.8547E+01 -2.1529E-01 -5.5417E+00
             2.3565E+00
 GRADIENT:   1.3422E+01  7.7741E-01  5.0319E-02  7.5611E+00  9.6597E-01  2.2543E+01  6.2463E+01  0.0000E+00  2.0262E-01  0.0000E+00
             2.8226E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -763.524216498732        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2493             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1359E+00  1.0131E+00  5.8546E+00  1.3140E+00  2.0606E+00  2.0548E+00  5.6254E+00  1.0000E-02  7.3159E-01  1.0000E-02
             9.5505E+00
 PARAMETER:  2.2740E-01  1.1297E-01  1.8672E+00  3.7305E-01  8.2301E-01  8.2016E-01  1.8273E+00 -1.8547E+01 -2.1254E-01 -5.5417E+00
             2.3566E+00
 GRADIENT:   1.3437E+01  8.6342E-01  4.9768E-02  7.6546E+00  9.5938E-01  2.2483E+01  6.2634E+01  0.0000E+00  2.2621E-01  0.0000E+00
             2.8198E+01

0ITERATION NO.:   82    OBJECTIVE VALUE:  -763.524216498732        NO. OF FUNC. EVALS.:  55
 CUMULATIVE NO. OF FUNC. EVALS.:     2548
 NPARAMETR:  1.1359E+00  1.0131E+00  5.8546E+00  1.3140E+00  2.0606E+00  2.0548E+00  5.6254E+00  1.0000E-02  7.3159E-01  1.0000E-02
             9.5505E+00
 PARAMETER:  2.2740E-01  1.1297E-01  1.8672E+00  3.7305E-01  8.2301E-01  8.2016E-01  1.8273E+00 -1.8547E+01 -2.1254E-01 -5.5417E+00
             2.3566E+00
 GRADIENT:   1.3245E-01  1.7648E-02 -2.2449E-02 -4.0944E-01  3.5394E-02  8.0036E-01  3.0098E+00  0.0000E+00  3.4469E-02  0.0000E+00
            -7.7960E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2548
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4830E-02  3.0643E-02  4.9941E-06 -5.8078E-02  1.5990E-05
 SE:             2.8155E-02  2.4244E-02  8.1053E-06  1.0758E-02  5.7209E-05
 N:                     100         100         100         100         100

 P VAL.:         5.9839E-01  2.0625E-01  5.3780E-01  6.7250E-08  7.7985E-01

 ETASHRINKSD(%)  5.6756E+00  1.8778E+01  9.9973E+01  6.3960E+01  9.9808E+01
 ETASHRINKVR(%)  1.1029E+01  3.4030E+01  1.0000E+02  8.7011E+01  1.0000E+02
 EBVSHRINKSD(%)  5.8955E+00  1.3417E+01  9.9957E+01  6.8921E+01  9.9737E+01
 EBVSHRINKVR(%)  1.1443E+01  2.5033E+01  1.0000E+02  9.0341E+01  9.9999E+01
 RELATIVEINF(%)  8.7927E+01  3.7567E+01  2.4383E-06  4.6855E+00  8.8472E-05
 EPSSHRINKSD(%)  7.9008E+00
 EPSSHRINKVR(%)  1.5177E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -763.52421649873236     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       339.20202334687474     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    59.02
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    11.30
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -763.524       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.14E+00  1.01E+00  5.85E+00  1.31E+00  2.06E+00  2.05E+00  5.63E+00  1.00E-02  7.32E-01  1.00E-02  9.55E+00
 


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
+        2.34E+02
 
 TH 2
+       -1.55E+01  6.55E+00
 
 TH 3
+       -2.04E-01  5.05E-02  8.56E-04
 
 TH 4
+       -1.26E+02  3.77E+01  3.51E-01  2.29E+02
 
 TH 5
+        9.55E+00 -3.30E+00 -3.48E-02 -2.02E+01  1.85E+00
 
 TH 6
+       -4.82E+01  3.44E+00  1.71E-01  4.04E+01 -4.63E+00  4.80E+01
 
 TH 7
+        9.59E+00 -2.43E+00 -2.45E-02 -1.51E+01  1.34E+00 -3.28E+00  1.02E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.15E+01 -1.12E+01 -1.14E-01 -6.95E+01  6.20E+00 -1.50E+01  4.65E+00  0.00E+00  2.13E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -8.93E-01 -2.07E+00 -1.16E-02 -1.07E+01  9.55E-01  4.37E-01  6.41E-01  0.00E+00  2.98E+00  0.00E+00  1.15E+00
 
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
+        1.76E+02
 
 TH 2
+       -3.26E+00  2.61E+01
 
 TH 3
+        5.06E-02  1.37E-01  8.95E-02
 
 TH 4
+       -7.68E+00  3.17E+01 -3.79E-02  1.70E+02
 
 TH 5
+       -1.31E+00 -4.55E+00 -1.02E+00 -1.25E+01  1.75E+01
 
 TH 6
+       -3.37E-01 -1.17E+00  4.92E-02  2.95E+00 -1.07E+00  3.71E+01
 
 TH 7
+        2.90E-01  3.24E+00 -3.56E-02 -1.08E+01  7.87E-01 -6.65E-01  3.51E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.91E-01 -2.88E+00 -3.21E-01 -4.44E+01  4.83E+00 -7.50E-01  2.51E+00  0.00E+00  2.79E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -6.67E+00 -2.54E+00  3.03E-03 -1.17E+01  7.89E-01  8.09E-01  6.08E-01  0.00E+00  2.48E+00  0.00E+00  7.49E+00
 
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
+        1.82E+02
 
 TH 2
+        4.56E+01  2.63E+01
 
 TH 3
+        5.52E-01  2.59E-01  2.69E-02
 
 TH 4
+        7.21E+01  3.31E+01  9.95E-01  1.72E+02
 
 TH 5
+       -1.52E+01 -5.29E+00 -4.59E-01 -2.60E+01  9.08E+00
 
 TH 6
+        2.98E+01  5.98E+00 -7.79E-02 -1.19E+01  5.38E-01  3.71E+01
 
 TH 7
+       -1.58E+00  3.76E+00 -7.86E-02 -1.17E+01  2.59E+00  2.40E+00  4.39E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -6.78E+00 -3.28E+00 -2.78E-01 -4.56E+01  7.47E+00  3.99E+00  2.61E+00  0.00E+00  2.35E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -3.13E+01 -8.21E+00 -1.47E-01 -7.73E+00  2.85E+00 -3.55E+00  1.25E+00  0.00E+00 -3.22E+00  0.00E+00  1.24E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       69.482
Stop Time:
Thu Sep 30 09:17:02 CDT 2021
