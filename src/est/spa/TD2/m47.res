Wed Sep 29 19:02:54 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat47.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m47.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1645.00731704865        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3440E+02 -4.0604E+01 -9.8506E+00 -4.7602E+01  4.1130E+01  6.1974E+01 -3.8064E+00 -1.8181E+00 -2.3786E+01  1.8999E+00
             9.9966E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1652.16034737823        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.8564E-01  1.0134E+00  1.0107E+00  1.0595E+00  9.7274E-01  9.3055E-01  1.0138E+00  1.0166E+00  1.0935E+00  9.8554E-01
             9.7578E-01
 PARAMETER:  8.5539E-02  1.1333E-01  1.1061E-01  1.5777E-01  7.2356E-02  2.8019E-02  1.1371E-01  1.1643E-01  1.8941E-01  8.5431E-02
             7.5486E-02
 GRADIENT:   6.7655E+00 -2.5745E+00  2.3736E-01 -8.5819E+00 -1.7481E+00 -1.2361E+00  1.1602E+00 -2.4459E+00 -1.5952E+00  2.1943E+00
            -3.0078E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1652.40870864619        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      351
 NPARAMETR:  9.8362E-01  8.8989E-01  1.1834E+00  1.1612E+00  9.8452E-01  9.2899E-01  8.7639E-01  1.2969E+00  1.0843E+00  9.5649E-01
             9.5569E-01
 PARAMETER:  8.3481E-02 -1.6658E-02  2.6840E-01  2.4947E-01  8.4398E-02  2.6346E-02 -3.1946E-02  3.6000E-01  1.8092E-01  5.5515E-02
             5.4680E-02
 GRADIENT:   5.9351E+00  1.2214E+01  2.2646E+00  1.6635E+01 -9.1211E+00 -1.2873E+00 -6.6452E-02  2.2583E+00  1.5846E+00 -4.7007E+00
            -9.9198E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1653.01941164899        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      527
 NPARAMETR:  9.8107E-01  8.7693E-01  1.2726E+00  1.1646E+00  1.0287E+00  9.3265E-01  8.2619E-01  1.2791E+00  1.0835E+00  1.0614E+00
             9.7866E-01
 PARAMETER:  8.0891E-02 -3.1325E-02  3.4107E-01  2.5237E-01  1.2825E-01  3.0278E-02 -9.0932E-02  3.4618E-01  1.8023E-01  1.5957E-01
             7.8428E-02
 GRADIENT:   5.0843E-01  6.6952E+00  2.2592E+00  6.3438E+00 -5.6306E+00  3.6976E-01  3.4532E-01 -5.5230E-02  4.1631E-01  9.3954E-01
            -2.7836E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1653.39833569232        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      703
 NPARAMETR:  9.7760E-01  5.9103E-01  1.3995E+00  1.3474E+00  9.7786E-01  9.3077E-01  7.8646E-01  1.2849E+00  9.6553E-01  1.0439E+00
             9.8080E-01
 PARAMETER:  7.7348E-02 -4.2588E-01  4.3614E-01  3.9820E-01  7.7613E-02  2.8255E-02 -1.4021E-01  3.5067E-01  6.4920E-02  1.4299E-01
             8.0611E-02
 GRADIENT:  -8.5797E-01  3.8182E+00  2.1456E+00  5.1224E+00 -3.1632E+00  2.2146E-01 -4.9283E-01 -6.6435E-01 -1.3308E+00 -1.2990E-01
             4.7725E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1653.41025686851        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      878
 NPARAMETR:  9.7717E-01  5.2965E-01  1.4207E+00  1.3895E+00  9.6529E-01  9.2954E-01  7.9733E-01  1.2972E+00  9.4294E-01  1.0351E+00
             9.7887E-01
 PARAMETER:  7.6905E-02 -5.3554E-01  4.5112E-01  4.2897E-01  6.4678E-02  2.6938E-02 -1.2648E-01  3.6019E-01  4.1243E-02  1.3454E-01
             7.8641E-02
 GRADIENT:  -4.0014E-02  4.4721E+00  2.2094E+00  9.5925E+00 -3.4191E+00 -1.7145E-01 -4.9129E-01 -6.2562E-01 -9.3162E-01 -3.8767E-01
            -3.0506E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1653.42571023369        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1053
 NPARAMETR:  9.7653E-01  4.6092E-01  1.4426E+00  1.4348E+00  9.5186E-01  9.2853E-01  8.4497E-01  1.3155E+00  9.1721E-01  1.0253E+00
             9.7719E-01
 PARAMETER:  7.6247E-02 -6.7454E-01  4.6646E-01  4.6100E-01  5.0665E-02  2.5850E-02 -6.8454E-02  3.7420E-01  1.3576E-02  1.2495E-01
             7.6925E-02
 GRADIENT:   6.4233E-01  4.1201E+00  1.7828E+00  1.1332E+01 -2.9176E+00 -4.6133E-01 -4.0917E-01 -4.6452E-01 -4.8511E-01 -5.3047E-01
            -9.0760E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1653.43355167390        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1228
 NPARAMETR:  9.7574E-01  3.9781E-01  1.4647E+00  1.4753E+00  9.4081E-01  9.2800E-01  9.3045E-01  1.3382E+00  8.9328E-01  1.0171E+00
             9.7640E-01
 PARAMETER:  7.5443E-02 -8.2177E-01  4.8166E-01  4.8885E-01  3.8990E-02  2.5279E-02  2.7911E-02  3.9133E-01 -1.2858E-02  1.1691E-01
             7.6113E-02
 GRADIENT:   9.2834E-01  3.4472E+00  1.3322E+00  1.0935E+01 -2.3243E+00 -5.5956E-01 -3.0827E-01 -3.1096E-01 -2.0570E-01 -5.2791E-01
            -1.1250E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1653.43475049560        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1407
 NPARAMETR:  9.7542E-01  3.7664E-01  1.4730E+00  1.4886E+00  9.3753E-01  9.2795E-01  9.7021E-01  1.3472E+00  8.8524E-01  1.0147E+00
             9.7640E-01
 PARAMETER:  7.5114E-02 -8.7647E-01  4.8731E-01  4.9784E-01  3.5495E-02  2.5224E-02  6.9757E-02  3.9806E-01 -2.1900E-02  1.1461E-01
             7.6116E-02
 GRADIENT:   9.0123E-01  3.1534E+00  1.1959E+00  1.0205E+01 -2.1090E+00 -5.4085E-01 -2.7088E-01 -2.7162E-01 -1.4509E-01 -4.9160E-01
            -1.0797E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1653.51791078662        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     1597             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7546E-01  3.6586E-01  1.4731E+00  1.4850E+00  9.3670E-01  9.2938E-01  1.1478E+00  1.3481E+00  8.7827E-01  1.0175E+00
             9.7830E-01
 PARAMETER:  7.5150E-02 -9.0549E-01  4.8738E-01  4.9540E-01  3.4612E-02  2.6761E-02  2.3783E-01  3.9867E-01 -2.9804E-02  1.1735E-01
             7.8062E-02
 GRADIENT:   4.1144E+02  4.8079E+01  9.8188E+00  6.5750E+02  7.6443E+00  3.9640E+01  1.3669E+00  2.3625E+00  1.5408E+01  1.1237E+00
             8.7290E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1653.51963222082        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1784             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7546E-01  3.6636E-01  1.4717E+00  1.4851E+00  9.3598E-01  9.2938E-01  1.1341E+00  1.3463E+00  8.7825E-01  1.0159E+00
             9.7828E-01
 PARAMETER:  7.5150E-02 -9.0415E-01  4.8641E-01  4.9548E-01  3.3841E-02  2.6762E-02  2.2584E-01  3.9738E-01 -2.9822E-02  1.1580E-01
             7.8036E-02
 GRADIENT:   4.1140E+02  4.8343E+01  1.0027E+01  6.5799E+02  7.3538E+00  3.9643E+01  1.2408E+00  2.2527E+00  1.5206E+01  9.6616E-01
             8.1952E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1653.52161516893        NO. OF FUNC. EVALS.: 121
 CUMULATIVE NO. OF FUNC. EVALS.:     1905
 NPARAMETR:  9.7562E-01  3.6897E-01  1.4622E+00  1.4844E+00  9.3303E-01  9.2945E-01  1.1273E+00  1.3381E+00  8.7703E-01  1.0134E+00
             9.7807E-01
 PARAMETER:  7.5319E-02 -8.9704E-01  4.7992E-01  4.9504E-01  3.0680E-02  2.6838E-02  2.1978E-01  3.9125E-01 -3.1218E-02  1.1330E-01
             7.7822E-02
 GRADIENT:   4.1170E+02  4.9088E+01  9.9400E+00  6.5813E+02  6.8650E+00  3.9653E+01  1.1544E+00  2.1723E+00  1.4305E+01  9.7330E-01
             7.1347E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1653.52561577553        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:     2020             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7574E-01  3.6668E-01  1.4482E+00  1.4838E+00  9.2767E-01  9.2960E-01  1.1545E+00  1.3265E+00  8.7877E-01  1.0082E+00
             9.7804E-01
 PARAMETER:  7.5439E-02 -9.0327E-01  4.7035E-01  4.9460E-01  2.4918E-02  2.6999E-02  2.4371E-01  3.8252E-01 -2.9228E-02  1.0818E-01
             7.7792E-02
 GRADIENT:   4.1182E+02  4.8250E+01  9.5263E+00  6.5618E+02  7.3699E+00  3.9664E+01  1.3894E+00  2.0985E+00  1.5411E+01  9.0445E-01
             8.2997E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1653.52618113833        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     2210
 NPARAMETR:  9.7551E-01  3.6723E-01  1.4482E+00  1.4834E+00  9.2752E-01  9.2938E-01  1.1524E+00  1.3264E+00  8.7830E-01  1.0079E+00
             9.7803E-01
 PARAMETER:  7.5202E-02 -9.0176E-01  4.7032E-01  4.9432E-01  2.4759E-02  2.6762E-02  2.4186E-01  3.8250E-01 -2.9770E-02  1.0785E-01
             7.7789E-02
 GRADIENT:   1.5273E+00 -9.1656E-03 -2.0963E-02 -8.6636E+00  1.1068E-01  1.0458E-01  1.2968E-02  3.6267E-02  1.0981E-01  2.6042E-02
             1.6288E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1653.52636390518        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2395             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7552E-01  3.6855E-01  1.4478E+00  1.4830E+00  9.2735E-01  9.2938E-01  1.1404E+00  1.3256E+00  8.7827E-01  1.0076E+00
             9.7798E-01
 PARAMETER:  7.5210E-02 -8.9819E-01  4.7007E-01  4.9405E-01  2.4571E-02  2.6764E-02  2.3134E-01  3.8183E-01 -2.9805E-02  1.0758E-01
             7.7734E-02
 GRADIENT:   4.1121E+02  4.8712E+01  9.9893E+00  6.5485E+02  6.5640E+00  3.9609E+01  1.2619E+00  2.0218E+00  1.4843E+01  8.6072E-01
             7.6022E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1653.52677290109        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2580
 NPARAMETR:  9.7553E-01  3.6881E-01  1.4473E+00  1.4826E+00  9.2743E-01  9.2938E-01  1.1452E+00  1.3254E+00  8.7857E-01  1.0076E+00
             9.7800E-01
 PARAMETER:  7.5221E-02 -8.9746E-01  4.6968E-01  4.9378E-01  2.4664E-02  2.6767E-02  2.3557E-01  3.8173E-01 -2.9460E-02  1.0756E-01
             7.7758E-02
 GRADIENT:   1.5067E+00  9.9749E-02  1.0996E-01 -8.1984E+00 -1.3993E-01  1.0385E-01 -3.6562E-03  1.0088E-02 -4.9036E-02  2.3912E-03
            -1.3517E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1653.52699240046        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2765             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7555E-01  3.6936E-01  1.4462E+00  1.4818E+00  9.2761E-01  9.2939E-01  1.1518E+00  1.3251E+00  8.7909E-01  1.0075E+00
             9.7805E-01
 PARAMETER:  7.5241E-02 -8.9598E-01  4.6896E-01  4.9327E-01  2.4852E-02  2.6773E-02  2.4131E-01  3.8147E-01 -2.8866E-02  1.0752E-01
             7.7802E-02
 GRADIENT:   4.1121E+02  4.8499E+01  9.5404E+00  6.5228E+02  7.3082E+00  3.9601E+01  1.3844E+00  2.1131E+00  1.5148E+01  8.5546E-01
             8.2917E-01

0ITERATION NO.:   84    OBJECTIVE VALUE:  -1653.52720694668        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:     2900
 NPARAMETR:  9.7554E-01  3.6980E-01  1.4464E+00  1.4819E+00  9.2747E-01  9.2939E-01  1.1447E+00  1.3247E+00  8.7897E-01  1.0076E+00
             9.7801E-01
 PARAMETER:  7.5239E-02 -8.9480E-01  4.6911E-01  4.9331E-01  2.4701E-02  2.6773E-02  2.3518E-01  3.8117E-01 -2.9008E-02  1.0754E-01
             7.7766E-02
 GRADIENT:  -1.7152E-02 -3.8998E-02  1.0923E-01  4.4247E-01 -8.1899E-02 -5.0992E-04 -1.9704E-03  7.0022E-03 -7.2508E-02  1.0443E-03
            -6.2376E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2900
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.9692E-04 -9.3875E-03 -3.6071E-02 -3.1538E-03 -4.1598E-02
 SE:             2.9840E-02  6.6095E-03  1.8507E-02  2.8736E-02  2.0550E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9206E-01  1.5552E-01  5.1294E-02  9.1261E-01  4.2945E-02

 ETASHRINKSD(%)  3.1374E-02  7.7857E+01  3.7998E+01  3.7315E+00  3.1155E+01
 ETASHRINKVR(%)  6.2739E-02  9.5097E+01  6.1558E+01  7.3237E+00  5.2603E+01
 EBVSHRINKSD(%)  4.4592E-01  7.8463E+01  4.1676E+01  3.9931E+00  2.7852E+01
 EBVSHRINKVR(%)  8.8985E-01  9.5362E+01  6.5983E+01  7.8267E+00  4.7947E+01
 RELATIVEINF(%)  9.6607E+01  2.0410E-01  8.8794E+00  5.0376E+00  7.4442E+00
 EPSSHRINKSD(%)  4.5779E+01
 EPSSHRINKVR(%)  7.0601E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1653.5272069466771     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -918.37638038293892     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    37.21
 Elapsed covariance  time in seconds:     5.61
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1653.527       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.76E-01  3.70E-01  1.45E+00  1.48E+00  9.27E-01  9.29E-01  1.14E+00  1.32E+00  8.79E-01  1.01E+00  9.78E-01
 


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
 
         2.86E-02  3.36E-01  3.23E-01  2.07E-01  1.51E-01  6.54E-02  1.32E+00  2.97E-01  1.45E-01  2.29E-01  6.48E-02
 


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
+        8.20E-04
 
 TH 2
+        3.29E-03  1.13E-01
 
 TH 3
+       -5.06E-04 -2.17E-02  1.04E-01
 
 TH 4
+       -1.93E-03 -6.82E-02  1.80E-02  4.30E-02
 
 TH 5
+        5.57E-04  2.62E-02  3.36E-02 -1.38E-02  2.28E-02
 
 TH 6
+       -7.08E-05 -3.16E-03  2.94E-03  1.83E-03  3.92E-04  4.27E-03
 
 TH 7
+       -9.58E-03 -4.08E-01  6.14E-02  2.42E-01 -1.03E-01  1.00E-02  1.75E+00
 
 TH 8
+       -1.28E-03 -3.13E-02  6.52E-02  2.26E-02  1.42E-02  3.21E-03  1.18E-01  8.82E-02
 
 TH 9
+        1.49E-03  4.16E-02 -1.33E-02 -2.50E-02  7.88E-03 -1.13E-03 -1.39E-01 -1.50E-02  2.11E-02
 
 TH10
+        4.17E-04  3.22E-02  4.38E-02 -1.73E-02  2.93E-02  2.04E-03 -1.30E-01  1.09E-02  8.97E-03  5.22E-02
 
 TH11
+        1.09E-04  7.95E-04  3.43E-03 -1.37E-04  1.58E-03  1.89E-04 -5.36E-03  2.88E-04 -2.48E-04  1.72E-03  4.20E-03
 
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
+        2.86E-02
 
 TH 2
+        3.42E-01  3.36E-01
 
 TH 3
+       -5.47E-02 -2.00E-01  3.23E-01
 
 TH 4
+       -3.24E-01 -9.79E-01  2.69E-01  2.07E-01
 
 TH 5
+        1.29E-01  5.16E-01  6.89E-01 -4.41E-01  1.51E-01
 
 TH 6
+       -3.78E-02 -1.44E-01  1.39E-01  1.35E-01  3.98E-02  6.54E-02
 
 TH 7
+       -2.53E-01 -9.18E-01  1.44E-01  8.81E-01 -5.16E-01  1.16E-01  1.32E+00
 
 TH 8
+       -1.51E-01 -3.14E-01  6.80E-01  3.67E-01  3.18E-01  1.66E-01  3.02E-01  2.97E-01
 
 TH 9
+        3.58E-01  8.53E-01 -2.84E-01 -8.29E-01  3.59E-01 -1.19E-01 -7.21E-01 -3.48E-01  1.45E-01
 
 TH10
+        6.37E-02  4.19E-01  5.94E-01 -3.65E-01  8.49E-01  1.36E-01 -4.29E-01  1.61E-01  2.70E-01  2.29E-01
 
 TH11
+        5.89E-02  3.65E-02  1.64E-01 -1.02E-02  1.62E-01  4.45E-02 -6.25E-02  1.49E-02 -2.64E-02  1.16E-01  6.48E-02
 
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
+        1.58E+03
 
 TH 2
+       -2.06E+02  5.63E+02
 
 TH 3
+       -1.18E+02  6.80E+01  1.23E+02
 
 TH 4
+       -1.02E+02  5.45E+02 -3.28E+00  7.34E+02
 
 TH 5
+        1.95E+02 -1.95E+02 -2.32E+02 -5.77E+01  7.26E+02
 
 TH 6
+       -4.94E+01  6.70E+01  1.88E+01  5.08E+01  1.10E+00  2.66E+02
 
 TH 7
+       -1.32E+01  2.99E+01  8.66E-02  1.64E+01  2.88E+00  3.14E+00  4.62E+00
 
 TH 8
+        3.69E+01 -2.91E+01 -2.52E+01 -2.33E+01  1.54E+00 -1.67E+01 -2.27E+00  2.87E+01
 
 TH 9
+       -5.32E+01 -1.53E+02  2.41E+01 -8.70E+01 -5.34E+01 -1.92E+01 -1.03E+01  8.20E+00  2.18E+02
 
 TH10
+        4.17E+01 -7.12E+00 -1.41E+01  3.25E+00 -9.65E+01 -3.56E+01 -9.65E-01  1.73E+01  1.10E+01  8.33E+01
 
 TH11
+       -1.99E+01 -3.44E+01 -1.48E+01 -3.83E+01 -1.42E+01 -1.88E+01 -2.42E-01  1.30E+01  2.36E+01  1.39E+01  2.57E+02
 
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
 #CPUT: Total CPU Time in Seconds,       42.884
Stop Time:
Wed Sep 29 19:03:41 CDT 2021
