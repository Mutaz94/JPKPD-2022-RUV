Sat Sep 25 13:26:24 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat28.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m28.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1623.36835368977        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.3205E+01 -1.1475E+02 -2.8254E+01 -1.5733E+02  5.9419E+01  7.7715E+00  2.5496E+00  8.1099E+00 -2.6098E+01  1.1553E+01
            -3.9897E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1638.31106657330        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.8358E-01  9.5312E-01  1.0463E+00  1.1322E+00  9.3342E-01  9.7511E-01  8.2741E-01  8.9872E-01  1.0830E+00  7.8970E-01
             1.0791E+00
 PARAMETER:  8.3442E-02  5.1987E-02  1.4529E-01  2.2415E-01  3.1099E-02  7.4790E-02 -8.9460E-02 -6.7885E-03  1.7975E-01 -1.3611E-01
             1.7616E-01
 GRADIENT:   2.9781E+01 -1.6730E+00  1.2435E+01  3.1688E-01  1.6578E+01 -1.6754E+00  2.2016E+00 -3.7320E+00  9.0859E+00 -1.5653E+01
            -1.1845E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1639.99274460753        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      182
 NPARAMETR:  9.8010E-01  8.5715E-01  9.9764E-01  1.1967E+00  8.7805E-01  9.7573E-01  6.7553E-01  7.3950E-01  1.0443E+00  8.5165E-01
             1.0840E+00
 PARAMETER:  7.9895E-02 -5.4142E-02  9.7637E-02  2.7959E-01 -3.0048E-02  7.5427E-02 -2.9226E-01 -2.0178E-01  1.4339E-01 -6.0582E-02
             1.8065E-01
 GRADIENT:  -1.0026E+01  9.6512E-01  3.4035E+00  2.2127E+00  1.1248E+01 -4.1238E+00 -9.5247E-01 -3.7824E+00  4.5972E+00 -6.8036E+00
            -7.4726E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1641.08219508645        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      361
 NPARAMETR:  9.8438E-01  7.1270E-01  9.1228E-01  1.2850E+00  7.7352E-01  9.8892E-01  9.6971E-01  6.7430E-01  9.3354E-01  7.9752E-01
             1.0986E+00
 PARAMETER:  8.4262E-02 -2.3869E-01  8.1884E-03  3.5074E-01 -1.5680E-01  8.8861E-02  6.9246E-02 -2.9408E-01  3.1229E-02 -1.2625E-01
             1.9405E-01
 GRADIENT:   6.3299E-03  6.1739E+00  4.2984E-01  8.5357E+00 -3.4649E+00  1.1921E+00 -6.7624E-02  3.4752E-01 -6.8909E-01  7.6216E-01
             1.0573E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1641.81768134795        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      537
 NPARAMETR:  9.7987E-01  4.6478E-01  9.7686E-01  1.4299E+00  7.3252E-01  9.8053E-01  1.0357E+00  7.1742E-01  8.6914E-01  8.1035E-01
             1.0935E+00
 PARAMETER:  7.9665E-02 -6.6620E-01  7.6583E-02  4.5763E-01 -2.1126E-01  8.0340E-02  1.3505E-01 -2.3210E-01 -4.0248E-02 -1.1028E-01
             1.8936E-01
 GRADIENT:  -1.8675E+00  1.5533E+00  3.0938E+00  6.2635E-01 -5.0898E+00 -1.0713E+00  1.5197E-01  4.3745E-02  3.4916E+00  2.1042E+00
            -4.0965E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1642.19756325036        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      712
 NPARAMETR:  9.7679E-01  2.8128E-01  9.6993E-01  1.5305E+00  6.8483E-01  9.7847E-01  1.2247E+00  7.5921E-01  8.0321E-01  7.7487E-01
             1.0944E+00
 PARAMETER:  7.6516E-02 -1.1684E+00  6.9469E-02  5.2559E-01 -2.7858E-01  7.8239E-02  3.0273E-01 -1.7548E-01 -1.1914E-01 -1.5506E-01
             1.9021E-01
 GRADIENT:  -1.4065E+00 -1.8465E-01 -1.7141E+00 -4.6012E+00  1.0227E+00 -1.0054E+00 -1.8273E-01  6.6490E-01 -7.7491E-01  6.0329E-01
            -5.6570E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1642.30653720584        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      887
 NPARAMETR:  9.7509E-01  1.7134E-01  9.8709E-01  1.5946E+00  6.6879E-01  9.7905E-01  1.7261E+00  7.9752E-01  7.7089E-01  7.5893E-01
             1.0965E+00
 PARAMETER:  7.4778E-02 -1.6641E+00  8.7009E-02  5.6664E-01 -3.0228E-01  7.8825E-02  6.4589E-01 -1.2625E-01 -1.6021E-01 -1.7585E-01
             1.9210E-01
 GRADIENT:   4.4940E-01 -2.6757E-01 -1.4792E+00 -3.3577E+00  3.5236E+00 -1.7615E-01 -3.3068E-02 -1.9772E-01 -4.9626E-01 -5.4469E-01
            -3.7847E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1642.33703555519        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1063
 NPARAMETR:  9.7343E-01  1.1115E-01  9.8535E-01  1.6289E+00  6.5400E-01  9.7848E-01  2.4987E+00  8.1659E-01  7.5529E-01  7.5208E-01
             1.0975E+00
 PARAMETER:  7.3068E-02 -2.0969E+00  8.5242E-02  5.8790E-01 -3.2465E-01  7.8244E-02  1.0158E+00 -1.0262E-01 -1.8066E-01 -1.8491E-01
             1.9305E-01
 GRADIENT:  -2.0934E-01  1.0960E-01 -3.0218E-02 -1.5102E-01 -4.8094E-01 -9.5555E-02  9.3971E-02  2.1581E-01  1.7821E-01  3.0668E-01
             9.3493E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1642.35540126883        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1240
 NPARAMETR:  9.7258E-01  6.7666E-02  1.0038E+00  1.6551E+00  6.5360E-01  9.7829E-01  3.0733E+00  8.4543E-01  7.4503E-01  7.4945E-01
             1.0978E+00
 PARAMETER:  7.2195E-02 -2.5932E+00  1.0378E-01  6.0384E-01 -3.2526E-01  7.8051E-02  1.2227E+00 -6.7908E-02 -1.9433E-01 -1.8841E-01
             1.9329E-01
 GRADIENT:   4.1978E-01  2.0453E-02 -3.9495E-02 -5.4239E-01  6.4393E-01  1.0407E-01  4.4607E-02 -3.1729E-02  3.0150E-01 -1.0524E-01
             2.8316E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1642.36642940934        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1415
 NPARAMETR:  9.7169E-01  3.7568E-02  1.0060E+00  1.6721E+00  6.4832E-01  9.7764E-01  3.1596E+00  8.5661E-01  7.3801E-01  7.4760E-01
             1.0979E+00
 PARAMETER:  7.1279E-02 -3.1816E+00  1.0597E-01  6.1406E-01 -3.3338E-01  7.7388E-02  1.2504E+00 -5.4776E-02 -2.0380E-01 -1.9089E-01
             1.9338E-01
 GRADIENT:   7.0374E-02  1.7091E-02  7.4481E-02  4.7007E-01 -2.3440E-01  2.1209E-02 -1.0270E-02  2.7943E-02 -4.1026E-02  5.9006E-03
             1.1424E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1642.36882626008        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1592
 NPARAMETR:  9.7108E-01  1.5673E-02  1.0102E+00  1.6850E+00  6.4562E-01  9.7725E-01  2.7735E+00  8.6704E-01  7.3257E-01  7.4608E-01
             1.0980E+00
 PARAMETER:  7.0655E-02 -4.0558E+00  1.1014E-01  6.2178E-01 -3.3755E-01  7.6987E-02  1.1201E+00 -4.2668E-02 -2.1120E-01 -1.9292E-01
             1.9345E-01
 GRADIENT:  -5.1880E-02  1.8879E-02  3.2835E-01  1.6216E+00 -9.0086E-01 -8.0155E-03 -4.2766E-03  5.4334E-02 -1.5497E-01  8.0696E-02
             4.0490E-03

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1642.37110216476        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1769
 NPARAMETR:  9.7094E-01  1.0000E-02  1.0124E+00  1.6878E+00  6.4571E-01  9.7719E-01  2.4677E+00  8.7011E-01  7.3133E-01  7.4569E-01
             1.0980E+00
 PARAMETER:  7.0514E-02 -4.5449E+00  1.1235E-01  6.2342E-01 -3.3740E-01  7.6921E-02  1.0033E+00 -3.9139E-02 -2.1289E-01 -1.9344E-01
             1.9345E-01
 GRADIENT:  -3.4840E-04  0.0000E+00 -6.9428E-04  1.1368E-02  1.4054E-03  4.8122E-03 -1.5155E-03  1.4213E-04 -5.6445E-03  6.8176E-04
            -8.4743E-03

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1642.37111619581        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1951
 NPARAMETR:  9.7094E-01  1.0000E-02  1.0127E+00  1.6878E+00  6.4583E-01  9.7716E-01  2.5021E+00  8.7036E-01  7.3132E-01  7.4576E-01
             1.0979E+00
 PARAMETER:  7.0513E-02 -4.5071E+00  1.1261E-01  6.2345E-01 -3.3722E-01  7.6894E-02  1.0171E+00 -3.8853E-02 -2.1290E-01 -1.9335E-01
             1.9342E-01
 GRADIENT:  -4.4796E-03  4.1466E-03  2.2106E-03  4.5224E-02  4.7405E-03 -1.3545E-02 -1.5459E-03  5.5074E-04 -4.5981E-03  3.2150E-03
            -1.7757E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1642.37125229475        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2131
 NPARAMETR:  9.7093E-01  1.0000E-02  1.0125E+00  1.6883E+00  6.4572E-01  9.7717E-01  2.9966E+00  8.7006E-01  7.3128E-01  7.4563E-01
             1.0979E+00
 PARAMETER:  7.0498E-02 -4.5116E+00  1.1244E-01  6.2370E-01 -3.3738E-01  7.6904E-02  1.1975E+00 -3.9198E-02 -2.1296E-01 -1.9352E-01
             1.9338E-01
 GRADIENT:  -6.0038E-02  0.0000E+00  2.5061E-02  1.0887E+00 -1.3663E-01 -1.0247E-02 -2.0927E-03 -1.9111E-02 -2.0602E-02 -1.1894E-02
            -5.0354E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1642.37361295738        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2313
 NPARAMETR:  9.7091E-01  1.0000E-02  1.0119E+00  1.6876E+00  6.4543E-01  9.7712E-01  6.7507E+00  8.6961E-01  7.3099E-01  7.4537E-01
             1.0980E+00
 PARAMETER:  7.0478E-02 -4.5116E+00  1.1181E-01  6.2330E-01 -3.3784E-01  7.6859E-02  2.0097E+00 -3.9716E-02 -2.1335E-01 -1.9388E-01
             1.9354E-01
 GRADIENT:  -8.7391E-02  0.0000E+00  6.7456E-02 -3.7509E-01 -5.2385E-02 -2.6479E-02 -3.4245E-04  9.2011E-03 -5.1293E-03  2.7933E-02
             3.2943E-02

0ITERATION NO.:   74    OBJECTIVE VALUE:  -1642.37363163139        NO. OF FUNC. EVALS.: 129
 CUMULATIVE NO. OF FUNC. EVALS.:     2442
 NPARAMETR:  9.7095E-01  1.0000E-02  1.0118E+00  1.6877E+00  6.4542E-01  9.7719E-01  6.8592E+00  8.6955E-01  7.3101E-01  7.4527E-01
             1.0980E+00
 PARAMETER:  7.0524E-02 -4.5116E+00  1.1176E-01  6.2336E-01 -3.3786E-01  7.6925E-02  2.0256E+00 -3.9780E-02 -2.1333E-01 -1.9400E-01
             1.9349E-01
 GRADIENT:   1.5189E-02  0.0000E+00  3.9297E-02 -1.3931E-01 -2.4735E-02 -6.4164E-04  4.1226E-05  1.9292E-03  3.0810E-03  1.3092E-02
             7.9780E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2442
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.0285E-04 -7.6061E-04 -1.7160E-02 -4.8092E-03 -2.2935E-02
 SE:             2.9812E-02  1.1317E-03  1.7752E-02  2.9088E-02  2.0591E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9457E-01  5.0153E-01  3.3370E-01  8.6868E-01  2.6535E-01

 ETASHRINKSD(%)  1.2495E-01  9.6209E+01  4.0530E+01  2.5501E+00  3.1018E+01
 ETASHRINKVR(%)  2.4975E-01  9.9856E+01  6.4633E+01  5.0351E+00  5.2414E+01
 EBVSHRINKSD(%)  5.0288E-01  9.6296E+01  4.1807E+01  2.8961E+00  2.9926E+01
 EBVSHRINKVR(%)  1.0032E+00  9.9863E+01  6.6136E+01  5.7083E+00  5.0896E+01
 RELATIVEINF(%)  9.3034E+01  5.3285E-03  3.4573E+00  5.8014E+00  2.4763E+00
 EPSSHRINKSD(%)  4.3910E+01
 EPSSHRINKVR(%)  6.8540E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1642.3736316313898     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -907.22280506765162     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    28.48
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0INVERSE COVARIANCE MATRIX SET TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S
 Elapsed covariance  time in seconds:     5.24
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1642.374       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.71E-01  1.00E-02  1.01E+00  1.69E+00  6.45E-01  9.77E-01  6.86E+00  8.70E-01  7.31E-01  7.45E-01  1.10E+00
 


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
 
         2.85E-02  0.00E+00  1.02E-01  3.86E-02  5.00E-02  6.18E-02  1.25E+01  1.06E-01  5.38E-02  1.42E-01  7.93E-02
 


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
+        8.12E-04
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.50E-04  0.00E+00  1.04E-02
 
 TH 4
+        4.86E-05  0.00E+00  8.08E-04  1.49E-03
 
 TH 5
+       -1.44E-04  0.00E+00  4.62E-03  4.19E-04  2.50E-03
 
 TH 6
+       -4.45E-04  0.00E+00 -1.07E-04 -9.11E-05 -1.35E-04  3.82E-03
 
 TH 7
+        7.58E-03  0.00E+00 -7.77E-01 -6.66E-02 -3.86E-01 -2.53E-01  1.55E+02
 
 TH 8
+        1.88E-05  0.00E+00 -2.33E-03 -2.68E-04 -1.24E-03 -1.18E-03  1.01E+00  1.12E-02
 
 TH 9
+        8.68E-05  0.00E+00 -4.45E-05 -2.49E-05  8.07E-05 -3.54E-04 -8.66E-02  1.29E-04  2.90E-03
 
 TH10
+        3.18E-04  0.00E+00  5.78E-03  6.69E-05  3.13E-03  4.76E-04 -1.30E+00 -6.12E-03 -2.23E-04  2.00E-02
 
 TH11
+       -5.36E-05  0.00E+00  1.72E-04  6.36E-04  9.28E-05 -2.37E-04 -5.09E-02 -1.30E-03  6.91E-04 -1.63E-03  6.29E-03
 
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
+        2.85E-02
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -5.18E-02  0.00E+00  1.02E-01
 
 TH 4
+        4.43E-02  0.00E+00  2.06E-01  3.86E-02
 
 TH 5
+       -1.01E-01  0.00E+00  9.08E-01  2.17E-01  5.00E-02
 
 TH 6
+       -2.53E-01  0.00E+00 -1.70E-02 -3.82E-02 -4.37E-02  6.18E-02
 
 TH 7
+        2.14E-02  0.00E+00 -6.12E-01 -1.39E-01 -6.20E-01 -3.28E-01  1.25E+01
 
 TH 8
+        6.22E-03  0.00E+00 -2.16E-01 -6.57E-02 -2.35E-01 -1.80E-01  7.68E-01  1.06E-01
 
 TH 9
+        5.66E-02  0.00E+00 -8.11E-03 -1.20E-02  3.00E-02 -1.06E-01 -1.29E-01  2.26E-02  5.38E-02
 
 TH10
+        7.88E-02  0.00E+00  4.01E-01  1.23E-02  4.42E-01  5.44E-02 -7.38E-01 -4.08E-01 -2.92E-02  1.42E-01
 
 TH11
+       -2.37E-02  0.00E+00  2.13E-02  2.08E-01  2.34E-02 -4.83E-02 -5.15E-02 -1.55E-01  1.62E-01 -1.45E-01  7.93E-02
 
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
+        1.26E+03
 
 TH 2
+       -2.75E-14 -2.02E-28
 
 TH 3
+       -1.18E+02  1.18E-13  5.15E+02
 
 TH 4
+       -6.43E+01 -4.97E-14 -2.63E+01  6.73E+02
 
 TH 5
+        3.04E+02  9.04E-14 -1.07E+03 -1.55E+02  2.48E+03
 
 TH 6
+       -7.16E+00 -1.14E-13 -3.75E+01 -1.58E+01  1.63E+01  2.61E+01
 
 TH 7
+        4.86E-03  7.51E-18  2.48E-03 -1.53E-02  1.99E-03 -4.51E-03  1.70E-05
 
 TH 8
+       -8.69E+00 -1.04E-13 -4.12E+01 -8.93E+00  2.84E+01  2.27E+01  9.09E-04  2.15E+01
 
 TH 9
+       -4.10E+01  5.60E-14  6.09E+01 -1.81E+01 -1.00E+02 -2.72E+01  6.97E-02 -4.14E+00  2.97E+02
 
 TH10
+       -1.02E+01 -1.70E-13 -7.62E+00 -1.27E+01 -8.74E+01  3.55E+01  1.15E-03  3.27E+01 -3.17E+00  5.70E+01
 
 TH11
+       -2.56E+01 -1.39E-13 -3.40E+01 -4.27E+00 -1.16E+01  2.67E+01  1.40E-02  2.92E+01  5.03E+01  4.66E+01  5.08E+01
 
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
+        1.22E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        9.79E+01  0.00E+00  4.91E+02
 
 TH 4
+        4.74E+01  0.00E+00 -7.35E+01  6.89E+02
 
 TH 5
+       -3.02E+02  0.00E+00 -9.60E+02 -1.11E+02  2.39E+03
 
 TH 6
+       -8.87E+01  0.00E+00  1.40E+01 -8.02E+00 -5.92E+01  1.49E+02
 
 TH 7
+        9.05E-03  0.00E+00 -7.40E-03 -1.90E-02  2.61E-02 -7.33E-03  1.04E-05
 
 TH 8
+        4.77E+00  0.00E+00 -5.39E+01  3.84E+00 -6.38E+00  1.16E+00  2.99E-03  4.42E+01
 
 TH 9
+        3.97E+01  0.00E+00 -2.62E+01 -1.46E+01  8.95E+01 -3.39E+01  4.61E-02  3.84E+00  3.32E+02
 
 TH10
+        6.55E+01  0.00E+00  9.80E+00 -2.34E+01 -1.72E+02  1.66E+01  8.17E-03  4.23E+01 -1.30E+01  1.42E+02
 
 TH11
+       -1.47E+01  0.00E+00 -1.90E+01  5.87E+01  1.54E+00  2.30E+00  5.07E-03  2.97E+00  5.32E+01  5.99E+00  1.88E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.26
 #CPUT: Total CPU Time in Seconds,       33.156
Stop Time:
Sat Sep 25 13:27:00 CDT 2021
