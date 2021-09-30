Thu Sep 30 02:36:34 CDT 2021
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
$DATA ../../../../data/spa1/D/dat9.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m9.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1583.36405101756        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1809E+02 -6.3999E+01 -3.8176E+01 -8.4381E+01  1.4339E+02 -2.7014E+02 -2.0744E+02 -3.2237E+01 -3.0463E+02 -8.9219E+01
            -5.1256E-01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1770.10744236833        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      171
 NPARAMETR:  9.0163E-01  1.0025E+00  1.0360E+00  1.0233E+00  9.7960E-01  1.1016E+00  1.7985E+00  1.1555E+00  1.5248E+00  1.2482E+00
             9.3510E-01
 PARAMETER: -3.5532E-03  1.0246E-01  1.3536E-01  1.2304E-01  7.9394E-02  1.9679E-01  6.8696E-01  2.4449E-01  5.2189E-01  3.2167E-01
             3.2900E-02
 GRADIENT:  -1.4114E+02 -5.9158E+01 -2.4462E+00 -1.0088E+01 -1.9618E+01 -2.7598E+02 -1.0886E+02 -7.9485E+00  1.2010E+01  1.1746E+01
            -1.7616E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1807.93847851699        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      349
 NPARAMETR:  9.1445E-01  1.0941E+00  1.1098E+00  9.9487E-01  1.0409E+00  1.1136E+00  2.9646E+00  1.3653E+00  1.1638E+00  1.4735E+00
             8.9179E-01
 PARAMETER:  1.0563E-02  1.8995E-01  2.0414E-01  9.4860E-02  1.4005E-01  2.0762E-01  1.1867E+00  4.1135E-01  2.5168E-01  4.8766E-01
            -1.4524E-02
 GRADIENT:  -1.1041E+02  6.4974E+00  9.6733E+00 -1.7128E+00 -4.3966E+01 -2.6401E+02 -3.0501E+00  4.6273E+00  1.9314E+01  2.8197E+01
            -4.9307E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1861.56200195439        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      527
 NPARAMETR:  9.6592E-01  9.1991E-01  9.0958E-01  1.0526E+00  9.6607E-01  1.5572E+00  3.7632E+00  1.0791E+00  8.7559E-01  1.1941E+00
             9.3920E-01
 PARAMETER:  6.5329E-02  1.6520E-02  5.2278E-03  1.5125E-01  6.5477E-02  5.4290E-01  1.4253E+00  1.7617E-01 -3.2854E-02  2.7736E-01
             3.7277E-02
 GRADIENT:  -5.4411E+00  4.3384E+00 -1.9852E+01  1.3236E+01  1.8129E+01 -2.3973E+01  1.8246E+01  2.7285E+00  6.2829E+00  4.2826E-01
             1.7532E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1865.58343503190        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      707
 NPARAMETR:  9.7336E-01  9.8687E-01  9.8709E-01  1.0141E+00  1.0143E+00  1.6321E+00  3.3091E+00  1.1761E+00  8.7789E-01  1.2130E+00
             9.4307E-01
 PARAMETER:  7.2998E-02  8.6782E-02  8.7002E-02  1.1402E-01  1.1418E-01  5.8984E-01  1.2967E+00  2.6217E-01 -3.0237E-02  2.9311E-01
             4.1384E-02
 GRADIENT:   1.1807E+00  1.3571E+00  7.3399E-02  1.9268E+00 -1.8625E+00 -5.8719E-01  2.1652E+00  4.5361E-01  1.5082E+00  6.8934E-01
            -6.5163E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1866.26529033414        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      868
 NPARAMETR:  9.7453E-01  1.0760E+00  9.0632E-01  9.4317E-01  1.0270E+00  1.7439E+00  2.8825E+00  1.1123E+00  9.0088E-01  1.2088E+00
             9.4521E-01
 PARAMETER:  7.4195E-02  1.7329E-01  1.6365E-03  4.1494E-02  1.2663E-01  6.5610E-01  1.1587E+00  2.0639E-01 -4.3848E-03  2.8961E-01
             4.3654E-02
 GRADIENT:   4.0775E+02  7.6719E+01  1.0776E+00  6.8982E+01  9.5077E+00  7.5172E+02  5.4448E+02  1.1072E+00  4.4451E+00  4.2381E+00
             1.3308E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1866.59045722022        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:     1026             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8399E-01  1.1261E+00  9.0329E-01  9.3797E-01  1.0319E+00  1.7074E+00  2.9035E+00  1.0916E+00  9.0026E-01  1.2024E+00
             9.4534E-01
 PARAMETER:  8.3863E-02  2.1873E-01 -1.7088E-03  3.5967E-02  1.3136E-01  6.3500E-01  1.1659E+00  1.8763E-01 -5.0748E-03  2.8431E-01
             4.3788E-02
 GRADIENT:   4.1498E+02  1.0985E+02  2.5089E+00  7.4316E+01  1.0499E+01  7.1142E+02  5.2991E+02 -1.0696E-01  3.6184E+00  2.5653E+00
             5.9713E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1866.63977132399        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:     1142
 NPARAMETR:  9.7461E-01  1.1284E+00  9.0223E-01  9.3570E-01  1.0321E+00  1.6867E+00  2.8708E+00  1.0967E+00  9.0535E-01  1.2031E+00
             9.4548E-01
 PARAMETER:  7.4278E-02  2.2076E-01 -2.8817E-03  3.3544E-02  1.3164E-01  6.2277E-01  1.1546E+00  1.9232E-01  5.6706E-04  2.8486E-01
             4.3933E-02
 GRADIENT:   1.1659E+00 -3.6746E-01  9.8590E-01  2.9257E+00 -9.4627E-01  1.4941E+01 -7.3262E+00 -2.9488E-01 -6.5690E-01 -3.9859E-01
            -3.3414E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1866.65954516432        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1323
 NPARAMETR:  9.7493E-01  1.1361E+00  8.9965E-01  9.3056E-01  1.0338E+00  1.6923E+00  2.8822E+00  1.0924E+00  9.0521E-01  1.2049E+00
             9.4563E-01
 PARAMETER:  7.4606E-02  2.2764E-01 -5.7518E-03  2.8029E-02  1.3323E-01  6.2607E-01  1.1586E+00  1.8839E-01  4.0815E-04  2.8637E-01
             4.4096E-02
 GRADIENT:   1.3959E+00  3.7152E-02  1.5480E+00  9.6297E-02 -1.5034E+00  1.6419E+01 -5.0859E+00 -4.0933E-01 -4.7790E-01 -4.0275E-01
            -2.7259E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1866.66732704856        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1512             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7495E-01  1.1387E+00  8.9675E-01  9.2847E-01  1.0346E+00  1.6923E+00  2.8753E+00  1.0938E+00  9.0815E-01  1.2068E+00
             9.4573E-01
 PARAMETER:  7.4635E-02  2.2991E-01 -8.9830E-03  2.5782E-02  1.3406E-01  6.2607E-01  1.1562E+00  1.8969E-01  3.6585E-03  2.8795E-01
             4.4202E-02
 GRADIENT:   4.0764E+02  1.1671E+02  2.4763E+00  7.0433E+01  9.6775E+00  7.0448E+02  5.2412E+02  2.4257E-01  4.0921E+00  3.2932E+00
             8.9390E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1866.66978561633        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     1702
 NPARAMETR:  9.7497E-01  1.1401E+00  8.9553E-01  9.2721E-01  1.0350E+00  1.6923E+00  2.8719E+00  1.0939E+00  9.0951E-01  1.2071E+00
             9.4576E-01
 PARAMETER:  7.4650E-02  2.3115E-01 -1.0339E-02  2.4430E-02  1.3444E-01  6.2607E-01  1.1550E+00  1.8976E-01  5.1483E-03  2.8822E-01
             4.4233E-02
 GRADIENT:   1.3934E+00 -4.7030E-01  6.9694E-01 -7.2183E-02 -8.3851E-01  1.6416E+01 -4.9093E+00 -9.7044E-02 -9.5889E-02 -2.9178E-02
            -7.6536E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1866.67190941793        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     1892             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7500E-01  1.1439E+00  8.9146E-01  9.2592E-01  1.0365E+00  1.6922E+00  2.8652E+00  1.0961E+00  9.1193E-01  1.2074E+00
             9.4585E-01
 PARAMETER:  7.4678E-02  2.3442E-01 -1.4892E-02  2.3031E-02  1.3590E-01  6.2605E-01  1.1526E+00  1.9179E-01  7.8044E-03  2.8847E-01
             4.4331E-02
 GRADIENT:   4.0751E+02  1.1959E+02  7.2125E-01  7.1529E+01  1.2044E+01  7.0420E+02  5.2152E+02  6.4802E-01  4.4310E+00  3.4322E+00
             1.0793E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1866.67298676046        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     2082
 NPARAMETR:  9.7501E-01  1.1456E+00  8.9035E-01  9.2502E-01  1.0369E+00  1.6922E+00  2.8619E+00  1.0956E+00  9.1254E-01  1.2075E+00
             9.4588E-01
 PARAMETER:  7.4692E-02  2.3596E-01 -1.6146E-02  2.2057E-02  1.3626E-01  6.2604E-01  1.1515E+00  1.9129E-01  8.4753E-03  2.8854E-01
             4.4363E-02
 GRADIENT:   1.3866E+00 -4.7329E-01 -1.0612E+00  1.9013E+00  1.4952E+00  1.6404E+01 -4.7736E+00  2.6604E-01  1.6500E-01  2.8190E-02
             9.5047E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1866.67719742796        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     2273             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7503E-01  1.1477E+00  8.9083E-01  9.2329E-01  1.0365E+00  1.6922E+00  2.8578E+00  1.0917E+00  9.1213E-01  1.2075E+00
             9.4587E-01
 PARAMETER:  7.4711E-02  2.3776E-01 -1.5605E-02  2.0190E-02  1.3583E-01  6.2604E-01  1.1501E+00  1.8773E-01  8.0232E-03  2.8855E-01
             4.4355E-02
 GRADIENT:   4.0755E+02  1.2186E+02  1.8482E+00  6.9690E+01  1.0543E+01  7.0422E+02  5.1966E+02  4.0227E-01  4.2496E+00  3.4138E+00
             9.8110E-01

0ITERATION NO.:   67    OBJECTIVE VALUE:  -1866.67719742796        NO. OF FUNC. EVALS.:  62
 CUMULATIVE NO. OF FUNC. EVALS.:     2335
 NPARAMETR:  9.7503E-01  1.1477E+00  8.9083E-01  9.2329E-01  1.0365E+00  1.6922E+00  2.8578E+00  1.0917E+00  9.1213E-01  1.2075E+00
             9.4587E-01
 PARAMETER:  7.4711E-02  2.3776E-01 -1.5605E-02  2.0190E-02  1.3583E-01  6.2604E-01  1.1501E+00  1.8773E-01  8.0232E-03  2.8855E-01
             4.4355E-02
 GRADIENT:  -9.9192E-03 -1.5454E-01  1.3284E-01  3.4197E-01 -1.7580E-01  2.3488E-03  2.6201E-01  7.5780E-03 -2.3291E-02  1.3626E-03
            -1.8624E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2335
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.3614E-04  1.1484E-02 -4.2139E-02 -1.9333E-02 -2.2264E-02
 SE:             2.9933E-02  2.6304E-02  1.3192E-02  1.8732E-02  2.1773E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8304E-01  6.6240E-01  1.4020E-03  3.0204E-01  3.0654E-01

 ETASHRINKSD(%)  1.0000E-10  1.1877E+01  5.5805E+01  3.7244E+01  2.7056E+01
 ETASHRINKVR(%)  1.0000E-10  2.2344E+01  8.0468E+01  6.0617E+01  4.6792E+01
 EBVSHRINKSD(%)  1.1532E-01  8.9186E+00  6.1708E+01  4.2367E+01  2.2704E+01
 EBVSHRINKVR(%)  2.3052E-01  1.7042E+01  8.5337E+01  6.6785E+01  4.0254E+01
 RELATIVEINF(%)  9.9688E+01  2.7137E+01  3.6231E+00  6.8979E+00  1.8627E+01
 EPSSHRINKSD(%)  3.5050E+01
 EPSSHRINKVR(%)  5.7815E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1866.6771974279554     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -947.73866422328274     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    43.71
 Elapsed covariance  time in seconds:     8.23
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1866.677       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.75E-01  1.15E+00  8.91E-01  9.23E-01  1.04E+00  1.69E+00  2.86E+00  1.09E+00  9.12E-01  1.21E+00  9.46E-01
 


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
 
         4.95E-02  2.65E-01  2.43E-01  1.30E-01  1.26E-01  1.36E-01  4.45E-01  4.14E-01  1.64E-01  1.56E-01  4.56E-02
 


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
+        2.45E-03
 
 TH 2
+        4.66E-03  7.01E-02
 
 TH 3
+        2.77E-03 -2.41E-02  5.88E-02
 
 TH 4
+       -1.50E-04 -3.00E-02  1.89E-02  1.69E-02
 
 TH 5
+        1.98E-03  1.35E-02  1.73E-02 -3.75E-03  1.60E-02
 
 TH 6
+        1.96E-03  6.67E-05  7.21E-03  2.10E-03  3.83E-03  1.84E-02
 
 TH 7
+       -1.02E-03 -1.01E-01  5.20E-02  4.81E-02 -9.33E-03  2.34E-02  1.98E-01
 
 TH 8
+        5.22E-03 -1.73E-02  8.47E-02  1.98E-02  3.17E-02  7.46E-03  4.83E-02  1.71E-01
 
 TH 9
+        1.13E-03  2.28E-02 -2.84E-04 -7.22E-03  5.50E-03 -6.60E-04 -3.73E-02  2.68E-03  2.70E-02
 
 TH10
+       -4.03E-04 -4.06E-03  1.33E-02  1.96E-03  9.18E-03  7.24E-03  2.17E-02  1.52E-02 -1.54E-03  2.42E-02
 
 TH11
+        3.10E-04  1.14E-03  8.43E-04 -1.52E-04  3.68E-04  4.35E-04 -1.86E-03  1.17E-03  1.02E-04 -1.55E-04  2.08E-03
 
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
+        4.95E-02
 
 TH 2
+        3.55E-01  2.65E-01
 
 TH 3
+        2.30E-01 -3.75E-01  2.43E-01
 
 TH 4
+       -2.33E-02 -8.72E-01  5.99E-01  1.30E-01
 
 TH 5
+        3.17E-01  4.02E-01  5.66E-01 -2.28E-01  1.26E-01
 
 TH 6
+        2.92E-01  1.86E-03  2.19E-01  1.19E-01  2.24E-01  1.36E-01
 
 TH 7
+       -4.65E-02 -8.57E-01  4.82E-01  8.33E-01 -1.66E-01  3.88E-01  4.45E-01
 
 TH 8
+        2.55E-01 -1.57E-01  8.44E-01  3.68E-01  6.06E-01  1.33E-01  2.62E-01  4.14E-01
 
 TH 9
+        1.39E-01  5.25E-01 -7.12E-03 -3.38E-01  2.65E-01 -2.97E-02 -5.10E-01  3.94E-02  1.64E-01
 
 TH10
+       -5.22E-02 -9.86E-02  3.51E-01  9.71E-02  4.67E-01  3.43E-01  3.13E-01  2.36E-01 -6.04E-02  1.56E-01
 
 TH11
+        1.37E-01  9.45E-02  7.62E-02 -2.57E-02  6.38E-02  7.04E-02 -9.17E-02  6.20E-02  1.36E-02 -2.18E-02  4.56E-02
 
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
+        1.20E+03
 
 TH 2
+       -4.99E+02  3.84E+02
 
 TH 3
+       -1.06E+02  8.85E+01  2.29E+02
 
 TH 4
+       -3.12E+02  1.97E+02 -1.92E+02  7.20E+02
 
 TH 5
+        1.23E+02 -1.65E+02 -2.52E+02  2.50E+02  5.39E+02
 
 TH 6
+        8.22E+01 -1.19E+02 -2.14E+01  1.44E+00  2.91E+01  1.39E+02
 
 TH 7
+       -1.48E+02  1.23E+02  2.38E+01 -2.29E+01 -3.89E+01 -7.32E+01  7.86E+01
 
 TH 8
+        9.37E+00 -9.34E+00 -3.79E+01  2.03E+00 -8.58E+00  5.30E+00 -2.15E+00  2.42E+01
 
 TH 9
+        6.40E+01 -4.94E+01 -3.13E+01 -4.64E+01  3.23E+01 -6.41E+00  9.40E+00  4.60E+00  6.90E+01
 
 TH10
+        7.85E+01 -1.75E+01  5.18E+00 -4.39E+00 -7.82E+01  1.89E+00 -2.45E+01  7.83E+00 -3.41E+00  8.29E+01
 
 TH11
+       -5.84E+01  1.29E+01 -4.62E+01  4.22E+00  5.42E+01 -4.04E+01  3.44E+01  4.43E+00  2.47E+01 -1.10E+01  5.26E+02
 
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
 #CPUT: Total CPU Time in Seconds,       52.000
Stop Time:
Thu Sep 30 02:37:27 CDT 2021
