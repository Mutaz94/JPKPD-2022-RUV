Wed Sep 29 17:13:37 CDT 2021
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
$DATA ../../../../data/spa/S2/dat24.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m24.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1647.27631673078        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4565E+02 -4.1290E+01  3.2584E+01 -9.6719E+01 -3.4231E+01  3.2528E+01 -4.2851E+00  6.0757E+00 -1.0824E+01 -4.2642E+00
             2.0564E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1657.42079122701        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.7431E-01  1.0329E+00  9.8352E-01  1.0809E+00  1.0299E+00  1.0241E+00  1.0244E+00  9.6350E-01  1.0118E+00  1.0427E+00
             9.3051E-01
 PARAMETER:  7.3970E-02  1.3238E-01  8.3378E-02  1.7779E-01  1.2944E-01  1.2377E-01  1.2412E-01  6.2814E-02  1.1170E-01  1.4177E-01
             2.7976E-02
 GRADIENT:   9.0343E-01  5.4682E-01  4.8153E-01 -8.9628E-01 -2.3581E+00  2.9619E+00 -3.4591E+00  6.5576E+00  4.9868E-01 -5.2133E+00
            -7.9937E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1659.85726711577        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      351
 NPARAMETR:  9.6914E-01  8.4758E-01  8.9736E-01  1.2198E+00  9.1463E-01  1.0571E+00  1.3954E+00  4.1407E-01  8.7758E-01  1.0840E+00
             9.7190E-01
 PARAMETER:  6.8653E-02 -6.5370E-02 -8.3029E-03  2.9868E-01  1.0767E-02  1.5554E-01  4.3318E-01 -7.8173E-01 -3.0591E-02  1.8066E-01
             7.1497E-02
 GRADIENT:  -9.4046E+00  1.8626E+01 -1.7626E+01  5.5400E+01  1.7555E+01  1.5105E+01 -1.4451E+00  8.4417E-01 -5.1383E-01  8.8485E+00
             1.0981E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1662.28838276217        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      528
 NPARAMETR:  9.7469E-01  7.4794E-01  8.0844E-01  1.2578E+00  7.9132E-01  1.0159E+00  1.5652E+00  2.8646E-01  8.4197E-01  8.9880E-01
             9.4369E-01
 PARAMETER:  7.4359E-02 -1.9043E-01 -1.1265E-01  3.2934E-01 -1.3406E-01  1.1580E-01  5.4804E-01 -1.1502E+00 -7.2007E-02 -6.7000E-03
             4.2048E-02
 GRADIENT:   2.0933E+00  2.2008E+01  8.9425E+00  3.6430E+01 -9.7223E+00 -1.2741E-01 -6.0652E-01  1.0504E-01 -3.6745E-01 -4.4245E+00
            -1.6867E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1664.75210635879        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      704
 NPARAMETR:  9.6859E-01  4.0757E-01  7.6542E-01  1.4232E+00  6.5896E-01  1.0129E+00  2.3917E+00  4.4259E-02  7.6204E-01  9.0875E-01
             9.4044E-01
 PARAMETER:  6.8085E-02 -7.9753E-01 -1.6733E-01  4.5292E-01 -3.1709E-01  1.1286E-01  9.7201E-01 -3.0177E+00 -1.7176E-01  4.3103E-03
             3.8596E-02
 GRADIENT:  -2.0419E+00  6.6545E+00  5.5846E+00  1.6141E+01 -1.1569E+01  6.2455E-01  5.9385E-01 -2.3740E-03  6.0734E-01  1.7279E+00
             4.4437E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1665.01381556007        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      886
 NPARAMETR:  9.6826E-01  3.2654E-01  7.6743E-01  1.4543E+00  6.4563E-01  1.0092E+00  2.7231E+00  2.1914E-02  7.4545E-01  9.0987E-01
             9.4022E-01
 PARAMETER:  6.7744E-02 -1.0192E+00 -1.6470E-01  4.7452E-01 -3.3752E-01  1.0918E-01  1.1018E+00 -3.7206E+00 -1.9377E-01  5.5459E-03
             3.8354E-02
 GRADIENT:   1.1829E+00  4.4433E-01  5.1251E-01 -3.9416E+00 -4.9761E-02 -3.8894E-02  4.3925E-01 -2.4240E-03 -2.7864E-01 -2.8046E-01
             1.5587E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1665.02085412432        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1073             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6833E-01  3.2687E-01  7.6663E-01  1.4518E+00  6.4546E-01  1.0096E+00  2.7324E+00  2.7274E-02  7.4619E-01  9.1072E-01
             9.3963E-01
 PARAMETER:  6.7818E-02 -1.0182E+00 -1.6576E-01  4.7280E-01 -3.3780E-01  1.0959E-01  1.1052E+00 -3.5018E+00 -1.9277E-01  6.4770E-03
             3.7732E-02
 GRADIENT:   4.0735E+02  4.9425E+01  1.0360E+01  6.5580E+02  4.4349E+01  4.7749E+01  4.7229E+01  1.0665E-02  1.6923E+01  1.2974E+00
             9.1506E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1665.02203622858        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1256
 NPARAMETR:  9.6835E-01  3.2733E-01  7.6618E-01  1.4515E+00  6.4532E-01  1.0097E+00  2.7295E+00  3.1331E-02  7.4630E-01  9.0974E-01
             9.3952E-01
 PARAMETER:  6.7835E-02 -1.0168E+00 -1.6634E-01  4.7261E-01 -3.3801E-01  1.0960E-01  1.1041E+00 -3.3632E+00 -1.9263E-01  5.3999E-03
             3.7615E-02
 GRADIENT:   1.3489E+00  1.8343E-01 -1.2529E-01 -8.7744E+00  6.1994E-01  1.2130E-01  1.0585E+00 -4.1219E-03  6.6703E-02  9.7475E-02
             7.7489E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1665.02365986470        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1444             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6837E-01  3.2835E-01  7.6597E-01  1.4511E+00  6.4488E-01  1.0097E+00  2.7270E+00  3.8183E-02  7.4639E-01  9.0880E-01
             9.3928E-01
 PARAMETER:  6.7858E-02 -1.0137E+00 -1.6662E-01  4.7232E-01 -3.3869E-01  1.0963E-01  1.1032E+00 -3.1654E+00 -1.9251E-01  4.3654E-03
             3.7359E-02
 GRADIENT:   4.0765E+02  4.9770E+01  1.1427E+01  6.5508E+02  4.3229E+01  4.7800E+01  4.7224E+01  1.8288E-02  1.6860E+01  1.1140E+00
             7.4482E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1665.02465148845        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1624
 NPARAMETR:  9.6838E-01  3.2883E-01  7.6578E-01  1.4508E+00  6.4471E-01  1.0097E+00  2.7252E+00  4.3282E-02  7.4646E-01  9.0826E-01
             9.3920E-01
 PARAMETER:  6.7872E-02 -1.0122E+00 -1.6686E-01  4.7214E-01 -3.3895E-01  1.0964E-01  1.1025E+00 -3.0400E+00 -1.9241E-01  3.7756E-03
             3.7271E-02
 GRADIENT:   1.3550E+00  5.2498E-01  1.3850E+00 -8.7532E+00 -1.3729E+00  1.1951E-01  1.1883E+00 -7.9491E-03  1.2906E-02 -5.2122E-03
            -8.3262E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1665.02684751430        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1810             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6840E-01  3.2888E-01  7.6509E-01  1.4506E+00  6.4481E-01  1.0097E+00  2.7216E+00  5.0265E-02  7.4662E-01  9.0774E-01
             9.3924E-01
 PARAMETER:  6.7886E-02 -1.0121E+00 -1.6777E-01  4.7194E-01 -3.3879E-01  1.0966E-01  1.1012E+00 -2.8904E+00 -1.9220E-01  3.2010E-03
             3.7316E-02
 GRADIENT:   4.0765E+02  4.9646E+01  1.0142E+01  6.5407E+02  4.4629E+01  4.7816E+01  4.7056E+01  2.9736E-02  1.6863E+01  1.1731E+00
             8.4330E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1665.02813480941        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1996
 NPARAMETR:  9.6841E-01  3.2902E-01  7.6458E-01  1.4503E+00  6.4492E-01  1.0097E+00  2.7181E+00  5.9646E-02  7.4679E-01  9.0717E-01
             9.3926E-01
 PARAMETER:  6.7900E-02 -1.0116E+00 -1.6843E-01  4.7177E-01 -3.3863E-01  1.0967E-01  1.0999E+00 -2.7193E+00 -1.9197E-01  2.5757E-03
             3.7334E-02
 GRADIENT:   1.3539E+00  3.6095E-02 -1.6349E+00 -8.6839E+00  1.7761E+00  1.2068E-01  9.5219E-01 -1.0988E-02  1.0572E-01  1.5031E-01
             1.5124E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1665.03119060738        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2182             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6842E-01  3.2985E-01  7.6501E-01  1.4500E+00  6.4446E-01  1.0097E+00  2.7177E+00  6.8185E-02  7.4677E-01  9.0624E-01
             9.3898E-01
 PARAMETER:  6.7912E-02 -1.0091E+00 -1.6787E-01  4.7159E-01 -3.3934E-01  1.0969E-01  1.0998E+00 -2.5855E+00 -1.9199E-01  1.5514E-03
             3.7038E-02
 GRADIENT:   4.0791E+02  4.9884E+01  1.1139E+01  6.5328E+02  4.3322E+01  4.7863E+01  4.6996E+01  4.8502E-02  1.6816E+01  1.1077E+00
             7.5256E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1665.03290161085        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2365
 NPARAMETR:  9.6842E-01  3.3014E-01  7.6509E-01  1.4500E+00  6.4443E-01  1.0097E+00  2.7171E+00  7.5536E-02  7.4681E-01  9.0575E-01
             9.3890E-01
 PARAMETER:  6.7916E-02 -1.0082E+00 -1.6776E-01  4.7153E-01 -3.3939E-01  1.0969E-01  1.0995E+00 -2.4831E+00 -1.9194E-01  1.0116E-03
             3.6950E-02
 GRADIENT:   1.3672E+00  5.0993E-01  8.9821E-01 -8.9489E+00 -1.3789E+00  1.2194E-01  1.1286E+00 -1.8394E-02  2.5603E-02  3.9638E-03
            -6.5858E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1665.03522288889        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2551             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6843E-01  3.2961E-01  7.6468E-01  1.4499E+00  6.4495E-01  1.0097E+00  2.7136E+00  8.7413E-02  7.4698E-01  9.0530E-01
             9.3901E-01
 PARAMETER:  6.7919E-02 -1.0098E+00 -1.6830E-01  4.7147E-01 -3.3859E-01  1.0970E-01  1.0983E+00 -2.3371E+00 -1.9171E-01  5.0903E-04
             3.7072E-02
 GRADIENT:   4.0787E+02  4.9512E+01  8.6274E+00  6.5282E+02  4.5849E+01  4.7866E+01  4.6757E+01  7.8436E-02  1.6867E+01  1.2377E+00
             9.5135E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1665.04411895238        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2737
 NPARAMETR:  9.6814E-01  3.2949E-01  7.6543E-01  1.4504E+00  6.4478E-01  1.0096E+00  2.7088E+00  1.4395E-01  7.4688E-01  9.0357E-01
             9.3875E-01
 PARAMETER:  6.7625E-02 -1.0102E+00 -1.6732E-01  4.7187E-01 -3.3884E-01  1.0953E-01  1.0965E+00 -1.8383E+00 -1.9186E-01 -1.4063E-03
             3.6796E-02
 GRADIENT:   7.5782E-01  1.0847E-01 -2.9734E+00 -8.5653E+00  8.3896E-01  5.9602E-02  6.3206E-01  1.7514E-02  1.4738E-02  8.9588E-01
             4.8949E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1665.04814189137        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2918             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6839E-01  3.2921E-01  7.6675E-01  1.4506E+00  6.4517E-01  1.0097E+00  2.7191E+00  1.6686E-01  7.4712E-01  9.0072E-01
             9.3835E-01
 PARAMETER:  6.7880E-02 -1.0111E+00 -1.6560E-01  4.7197E-01 -3.3824E-01  1.0969E-01  1.1003E+00 -1.6906E+00 -1.9153E-01 -4.5564E-03
             3.6368E-02
 GRADIENT:   4.0860E+02  4.9825E+01  7.7043E+00  6.5481E+02  4.4308E+01  4.7947E+01  4.7013E+01  2.9690E-01  1.6920E+01  1.7636E+00
             1.1791E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1665.05118480081        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:     3081
 NPARAMETR:  9.6837E-01  3.2908E-01  7.6759E-01  1.4508E+00  6.4522E-01  1.0097E+00  2.7210E+00  1.5827E-01  7.4698E-01  8.9985E-01
             9.3816E-01
 PARAMETER:  6.7863E-02 -1.0115E+00 -1.6450E-01  4.7214E-01 -3.3817E-01  1.0967E-01  1.1010E+00 -1.7435E+00 -1.9172E-01 -5.5302E-03
             3.6171E-02
 GRADIENT:   2.9208E-02  1.6452E-01 -3.5272E-01 -3.9428E-01 -1.3461E+00  5.9584E-03  5.2265E-02 -2.1965E-03  3.8112E-02  5.1018E-02
             2.6753E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     3081
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.5276E-04  3.5269E-02 -8.9553E-03 -2.5045E-02  5.6803E-03
 SE:             2.9864E-02  1.8200E-02  4.0867E-03  2.5621E-02  2.4094E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7455E-01  5.2647E-02  2.8426E-02  3.2831E-01  8.1362E-01

 ETASHRINKSD(%)  1.0000E-10  3.9027E+01  8.6309E+01  1.4168E+01  1.9281E+01
 ETASHRINKVR(%)  1.0000E-10  6.2823E+01  9.8126E+01  2.6328E+01  3.4845E+01
 EBVSHRINKSD(%)  3.9853E-01  4.5736E+01  8.7180E+01  1.1275E+01  1.4788E+01
 EBVSHRINKVR(%)  7.9547E-01  7.0554E+01  9.8357E+01  2.1279E+01  2.7388E+01
 RELATIVEINF(%)  9.8440E+01  6.3593E+00  1.8154E-01  2.1968E+01  6.8151E+00
 EPSSHRINKSD(%)  4.4301E+01
 EPSSHRINKVR(%)  6.8976E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1665.0511848008073     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -929.90035823706910     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    46.48
 Elapsed covariance  time in seconds:     6.58
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1665.051       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.68E-01  3.29E-01  7.68E-01  1.45E+00  6.45E-01  1.01E+00  2.72E+00  1.58E-01  7.47E-01  9.00E-01  9.38E-01
 


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
 
         2.99E-02  1.50E-01  1.87E-01  8.79E-02  1.25E-01  6.05E-02  8.14E-01  7.12E-01  6.93E-02  1.37E-01  6.63E-02
 


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
+        8.97E-04
 
 TH 2
+        3.74E-04  2.25E-02
 
 TH 3
+       -9.80E-04 -1.73E-03  3.52E-02
 
 TH 4
+       -4.22E-04 -1.11E-02  6.74E-03  7.73E-03
 
 TH 5
+       -5.47E-04  3.95E-03  2.23E-02  1.74E-03  1.57E-02
 
 TH 6
+        3.37E-04  6.70E-04 -9.67E-04 -3.34E-04 -5.69E-04  3.66E-03
 
 TH 7
+       -3.48E-03 -1.12E-01  3.15E-02  6.08E-02 -5.11E-03 -5.35E-03  6.63E-01
 
 TH 8
+       -3.71E-03 -1.16E-02  1.22E-01  2.62E-02  7.57E-02 -4.18E-03  1.40E-01  5.07E-01
 
 TH 9
+        3.48E-04  5.59E-03 -2.56E-03 -2.73E-03 -4.36E-04  2.68E-04 -3.01E-02 -7.21E-03  4.80E-03
 
 TH10
+       -3.26E-04 -8.27E-04  1.36E-02  2.47E-03  8.77E-03 -7.59E-04  1.45E-02  3.03E-02 -1.97E-03  1.87E-02
 
 TH11
+        3.80E-04  1.17E-03  7.28E-04 -1.93E-04  9.48E-04 -2.92E-04 -5.94E-03  2.48E-03  7.60E-04 -8.01E-04  4.39E-03
 
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
+        2.99E-02
 
 TH 2
+        8.33E-02  1.50E-01
 
 TH 3
+       -1.74E-01 -6.15E-02  1.87E-01
 
 TH 4
+       -1.60E-01 -8.40E-01  4.09E-01  8.79E-02
 
 TH 5
+       -1.46E-01  2.11E-01  9.49E-01  1.58E-01  1.25E-01
 
 TH 6
+        1.86E-01  7.39E-02 -8.53E-02 -6.29E-02 -7.52E-02  6.05E-02
 
 TH 7
+       -1.43E-01 -9.21E-01  2.06E-01  8.49E-01 -5.01E-02 -1.09E-01  8.14E-01
 
 TH 8
+       -1.74E-01 -1.08E-01  9.15E-01  4.19E-01  8.49E-01 -9.71E-02  2.42E-01  7.12E-01
 
 TH 9
+        1.68E-01  5.39E-01 -1.97E-01 -4.49E-01 -5.03E-02  6.40E-02 -5.33E-01 -1.46E-01  6.93E-02
 
 TH10
+       -7.97E-02 -4.04E-02  5.32E-01  2.05E-01  5.12E-01 -9.18E-02  1.30E-01  3.11E-01 -2.08E-01  1.37E-01
 
 TH11
+        1.92E-01  1.18E-01  5.85E-02 -3.31E-02  1.14E-01 -7.28E-02 -1.10E-01  5.25E-02  1.66E-01 -8.84E-02  6.63E-02
 
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
+        1.30E+03
 
 TH 2
+        1.02E+02  5.82E+02
 
 TH 3
+       -4.01E+01  2.09E+02  1.51E+03
 
 TH 4
+        1.30E+02  3.08E+02 -3.27E+02  8.77E+02
 
 TH 5
+        7.58E+01 -4.72E+02 -1.73E+03  2.51E+02  2.59E+03
 
 TH 6
+       -1.20E+02 -1.96E+01 -4.17E+01 -4.93E+01  3.82E+01  2.97E+02
 
 TH 7
+        9.41E+00  5.16E+01  2.33E+00 -1.71E+01 -2.78E+00  4.17E+00  1.23E+01
 
 TH 8
+        7.99E-01  2.05E+00 -7.77E+01  6.82E+00  1.08E+01  5.69E+00 -1.16E+00  1.76E+01
 
 TH 9
+       -7.54E+01 -1.14E+02  6.81E+01 -1.02E+02  8.00E+00  4.13E+00  3.98E+00 -1.19E+01  3.39E+02
 
 TH10
+       -2.33E+01 -8.03E-01 -1.05E+02  7.97E+00 -2.26E+01  1.77E+01 -2.67E+00  2.24E+01  4.60E+00  1.07E+02
 
 TH11
+       -1.30E+02  3.31E+00  7.04E+01 -6.62E+01 -1.54E+02  3.68E+01  1.21E+00  5.22E+00 -2.58E+01  2.85E+01  2.67E+02
 
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
 #CPUT: Total CPU Time in Seconds,       53.124
Stop Time:
Wed Sep 29 17:14:32 CDT 2021
