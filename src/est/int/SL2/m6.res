Sat Sep 25 00:49:41 CDT 2021
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
$DATA ../../../../data/int/SL2/dat6.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      999
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      899
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
 RAW OUTPUT FILE (FILE): m6.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2041.20247247147        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4690E+01  7.2212E+01 -1.1075E+02  3.0494E+01  2.0160E+02  1.6717E+01 -1.0403E+02 -1.3731E+02 -1.2650E+02 -4.6509E+01
            -3.3278E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3043.36676656125        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0242E+00  1.0135E+00  1.3232E+00  1.0133E+00  1.0006E+00  9.4123E-01  1.0788E+00  9.9619E-01  1.0846E+00  1.0662E+00
             2.0986E+00
 PARAMETER:  1.2394E-01  1.1342E-01  3.8005E-01  1.1318E-01  1.0059E-01  3.9435E-02  1.7583E-01  9.6183E-02  1.8125E-01  1.6412E-01
             8.4126E-01
 GRADIENT:   4.6934E+01 -2.0732E+01 -9.0066E+00 -1.5370E+00 -2.2494E+00 -6.9127E+00 -2.2364E+00  1.3368E+00  4.7463E+00 -5.8507E+00
            -7.2276E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3047.70863083108        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0033E+00  1.2946E+00  1.7340E+00  8.8436E-01  1.2641E+00  1.0031E+00  9.2631E-01  9.6529E-01  1.0679E+00  1.3456E+00
             2.1033E+00
 PARAMETER:  1.0329E-01  3.5818E-01  6.5045E-01 -2.2892E-02  3.3438E-01  1.0306E-01  2.3449E-02  6.4674E-02  1.6569E-01  3.9681E-01
             8.4351E-01
 GRADIENT:  -5.1298E+00  3.6411E+01  1.0330E+00  1.2423E+01 -2.6579E+00  1.6736E+01  5.6273E+00 -1.1107E+00 -9.6245E+00 -1.3938E+00
            -1.1031E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3051.07615812620        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.0065E+00  1.1631E+00  1.8771E+00  9.4784E-01  1.2246E+00  9.5758E-01  8.8948E-01  2.0699E+00  1.0410E+00  1.2812E+00
             2.0814E+00
 PARAMETER:  1.0650E-01  2.5112E-01  7.2971E-01  4.6430E-02  3.0264E-01  5.6654E-02 -1.7121E-02  8.2752E-01  1.4016E-01  3.4781E-01
             8.3303E-01
 GRADIENT:   2.2947E+00 -1.1300E+01 -1.7746E+01 -1.2742E+01  1.0723E+01 -1.0274E-01  3.1366E-01 -1.0737E+00  2.8671E-01 -1.0290E+00
             9.5758E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3053.79443154731        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      305
 NPARAMETR:  1.0071E+00  1.1007E+00  2.1721E+00  1.0017E+00  1.2306E+00  9.6526E-01  9.0617E-01  2.5879E+00  9.6413E-01  1.2756E+00
             2.0616E+00
 PARAMETER:  1.0712E-01  1.9590E-01  8.7570E-01  1.0167E-01  3.0753E-01  6.4640E-02  1.4703E-03  1.0508E+00  6.3475E-02  3.4340E-01
             8.2349E-01
 GRADIENT:   3.7532E+00 -2.3751E+00 -9.3492E+00 -8.0898E+00  1.8070E+01  2.7425E+00 -1.3864E+00 -3.3988E+00 -3.3517E+00  5.5591E+00
             1.2263E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3054.05538808547        NO. OF FUNC. EVALS.: 128
 CUMULATIVE NO. OF FUNC. EVALS.:      433
 NPARAMETR:  1.0046E+00  1.0699E+00  2.1717E+00  1.0216E+00  1.1884E+00  9.5675E-01  9.3033E-01  2.5873E+00  9.7459E-01  1.2755E+00
             2.0619E+00
 PARAMETER:  1.0461E-01  1.6757E-01  8.7553E-01  1.2141E-01  2.7259E-01  5.5790E-02  2.7788E-02  1.0506E+00  7.4257E-02  3.4333E-01
             8.2365E-01
 GRADIENT:  -2.3105E+00  1.2298E+00 -6.8818E+00 -5.4166E+00 -1.4463E+00 -6.1679E-01  1.9336E-01 -3.7301E+00  1.8411E+00  1.0588E+01
             1.7602E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3054.11647632091        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      525
 NPARAMETR:  1.0057E+00  1.0439E+00  2.1717E+00  1.0413E+00  1.1762E+00  9.5841E-01  9.5660E-01  2.5872E+00  9.4991E-01  1.2755E+00
             2.0619E+00
 PARAMETER:  1.0568E-01  1.4297E-01  8.7551E-01  1.4045E-01  2.6226E-01  5.7525E-02  5.5635E-02  1.0506E+00  4.8613E-02  3.4331E-01
             8.2365E-01
 GRADIENT:  -1.1022E+01 -1.6936E+00 -1.0052E+01 -2.7661E+00 -3.3851E+00 -1.1179E+00 -4.5521E-02 -5.0190E+00 -9.2381E-02  1.2806E+01
             1.8323E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3054.16241407655        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      705
 NPARAMETR:  1.0104E+00  1.0545E+00  2.1718E+00  1.0366E+00  1.1857E+00  9.6207E-01  9.4781E-01  2.5873E+00  9.5564E-01  1.2754E+00
             2.0617E+00
 PARAMETER:  1.1030E-01  1.5306E-01  8.7558E-01  1.3595E-01  2.7030E-01  6.1332E-02  4.6404E-02  1.0506E+00  5.4624E-02  3.4328E-01
             8.2355E-01
 GRADIENT:  -1.6247E-01 -6.4949E-02 -1.0248E+01  4.9734E-02 -7.7972E-02  3.7095E-01 -1.1959E-01 -4.8809E+00  1.3472E-01  1.1546E+01
             1.6685E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3054.37201212841        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      886
 NPARAMETR:  1.0099E+00  1.0381E+00  2.1873E+00  1.0464E+00  1.1773E+00  9.6298E-01  9.6788E-01  2.5902E+00  9.4785E-01  1.2657E+00
             2.0338E+00
 PARAMETER:  1.0988E-01  1.3742E-01  8.8268E-01  1.4535E-01  2.6320E-01  6.2281E-02  6.7358E-02  1.0517E+00  4.6443E-02  3.3565E-01
             8.0991E-01
 GRADIENT:  -3.0273E-01  7.9892E-03 -1.0123E+01 -4.6979E-01 -5.2458E-01  6.2456E-01 -1.5193E-02 -6.6295E+00  1.5539E-01  1.1613E+01
            -1.0766E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3054.37273956323        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:     1051
 NPARAMETR:  1.0099E+00  1.0381E+00  2.1873E+00  1.0464E+00  1.1773E+00  9.6171E-01  9.6799E-01  2.5901E+00  9.4779E-01  1.2657E+00
             2.0338E+00
 PARAMETER:  1.0989E-01  1.3737E-01  8.8266E-01  1.4540E-01  2.6318E-01  6.0963E-02  6.7462E-02  1.0517E+00  4.6377E-02  3.3562E-01
             8.0992E-01
 GRADIENT:  -2.9560E-01  2.1830E-02 -1.0130E+01 -4.2956E-01 -5.1059E-01  1.2377E-01 -1.3292E-02 -6.6312E+00  1.5378E-01  1.1613E+01
            -1.0742E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3054.37465630133        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1234
 NPARAMETR:  1.0101E+00  1.0381E+00  2.1870E+00  1.0465E+00  1.1780E+00  9.6169E-01  9.6936E-01  2.5896E+00  9.4648E-01  1.2656E+00
             2.0345E+00
 PARAMETER:  1.1009E-01  1.3740E-01  8.8251E-01  1.4543E-01  2.6381E-01  6.0936E-02  6.8878E-02  1.0515E+00  4.4996E-02  3.3551E-01
             8.1025E-01
 GRADIENT:   1.4916E-01 -1.0108E-01 -1.0272E+01 -3.6288E-01  5.1791E-02  1.1687E-01  9.9753E-03 -6.6095E+00 -3.4001E-02  1.1567E+01
            -1.0114E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -3054.37833235012        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1420
 NPARAMETR:  1.0101E+00  1.0381E+00  2.1867E+00  1.0465E+00  1.1779E+00  9.6150E-01  9.6822E-01  2.5892E+00  9.4614E-01  1.2654E+00
             2.0354E+00
 PARAMETER:  1.1008E-01  1.3744E-01  8.8241E-01  1.4547E-01  2.6373E-01  6.0739E-02  6.7701E-02  1.0514E+00  4.4637E-02  3.3538E-01
             8.1068E-01
 GRADIENT:   1.1302E-01 -5.9366E-02 -1.0265E+01 -2.8761E-01 -2.0585E-02  4.5936E-02 -6.7166E-02 -6.5936E+00 -1.4600E-01  1.1544E+01
            -9.2432E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -3054.38775590443        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1599
 NPARAMETR:  1.0103E+00  1.0382E+00  2.1875E+00  1.0465E+00  1.1777E+00  9.6176E-01  9.6139E-01  2.5900E+00  9.4445E-01  1.2652E+00
             2.0366E+00
 PARAMETER:  1.1022E-01  1.3744E-01  8.8277E-01  1.4548E-01  2.6355E-01  6.1014E-02  6.0621E-02  1.0517E+00  4.2849E-02  3.3525E-01
             8.1127E-01
 GRADIENT:   3.7137E-01 -7.0102E-02 -1.0194E+01 -4.1083E-01 -2.6629E-01  1.6156E-01 -5.3944E-01 -6.6203E+00 -7.9358E-01  1.1452E+01
            -8.1401E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -3054.39101532464        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:     1765
 NPARAMETR:  1.0101E+00  1.0382E+00  2.1874E+00  1.0465E+00  1.1778E+00  9.6145E-01  9.6903E-01  2.5898E+00  9.4631E-01  1.2652E+00
             2.0368E+00
 PARAMETER:  1.1005E-01  1.3746E-01  8.8271E-01  1.4549E-01  2.6367E-01  6.0687E-02  6.8540E-02  1.0516E+00  4.4818E-02  3.3521E-01
             8.1139E-01
 GRADIENT:  -1.4834E-02 -3.5484E-02 -1.0214E+01 -3.0539E-01 -1.4672E-01  3.1507E-02  7.3607E-03 -6.5063E+00 -4.8794E-02  1.1544E+01
            -7.6567E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -3054.40741030947        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1945
 NPARAMETR:  1.0104E+00  1.0382E+00  2.1894E+00  1.0466E+00  1.1808E+00  9.6074E-01  9.7289E-01  2.5918E+00  9.5382E-01  1.2647E+00
             2.0397E+00
 PARAMETER:  1.1035E-01  1.3748E-01  8.8361E-01  1.4551E-01  2.6620E-01  5.9948E-02  7.2514E-02  1.0524E+00  5.2717E-02  3.3483E-01
             8.1279E-01
 GRADIENT:   6.5242E-01 -1.1809E+00 -1.0512E+01  2.8236E-01  2.0897E+00 -2.5377E-01  6.1125E-01 -6.2440E+00  1.7604E+00  1.1409E+01
            -4.3345E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -3054.42027844943        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     2143
 NPARAMETR:  1.0104E+00  1.0383E+00  2.1892E+00  1.0466E+00  1.1783E+00  9.6136E-01  9.7151E-01  2.5915E+00  9.4618E-01  1.2645E+00
             2.0407E+00
 PARAMETER:  1.1037E-01  1.3760E-01  8.8355E-01  1.4556E-01  2.6406E-01  6.0597E-02  7.1093E-02  1.0522E+00  4.4680E-02  3.3466E-01
             8.1329E-01
 GRADIENT:   6.4572E-01 -3.3212E-02 -1.0169E+01 -2.4276E-01 -4.0308E-02  1.2825E-02  1.9658E-01 -6.2785E+00  1.1184E-01  1.1488E+01
            -3.4966E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -3054.42470462758        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     2331             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0104E+00  1.0383E+00  2.1896E+00  1.0466E+00  1.1780E+00  9.6139E-01  9.7074E-01  2.5919E+00  9.4778E-01  1.2643E+00
             2.0411E+00
 PARAMETER:  1.1032E-01  1.3761E-01  8.8374E-01  1.4557E-01  2.6379E-01  6.0620E-02  7.0299E-02  1.0524E+00  4.6370E-02  3.3456E-01
             8.1350E-01
 GRADIENT:   1.2646E+01  2.3118E+00 -8.8548E+00  3.6174E+00  2.5523E+00  1.2535E+00  3.6385E-01 -5.5645E+00  6.1195E-01  1.2024E+01
            -1.5229E+00

0ITERATION NO.:   83    OBJECTIVE VALUE:  -3054.42474546248        NO. OF FUNC. EVALS.:  80
 CUMULATIVE NO. OF FUNC. EVALS.:     2411
 NPARAMETR:  1.0099E+00  1.0383E+00  2.1894E+00  1.0466E+00  1.1779E+00  9.6120E-01  9.7032E-01  2.5916E+00  9.4749E-01  1.2643E+00
             2.0413E+00
 PARAMETER:  1.0982E-01  1.3761E-01  8.8374E-01  1.4557E-01  2.6353E-01  6.0323E-02  6.9856E-02  1.0524E+00  4.6055E-02  3.3455E-01
             8.1350E-01
 GRADIENT:  -5.4052E+03 -8.6260E+03  1.3295E+03 -4.0774E+03 -4.6758E-01 -8.3768E-02 -1.1870E+04  1.1231E+03 -1.1870E+04  1.7851E+03
            -7.3414E+02
 NUMSIGDIG:         3.3         3.3         3.3         3.3         2.4         2.4         3.3         3.3         3.3         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2411
 NO. OF SIG. DIGITS IN FINAL EST.:  2.4
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4978E-03 -2.8641E-02 -1.6705E-02  1.2985E-02 -3.2320E-02
 SE:             2.9613E-02  1.8176E-02  2.1141E-02  2.5232E-02  2.2081E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5966E-01  1.1508E-01  4.2944E-01  6.0680E-01  1.4328E-01

 ETASHRINKSD(%)  7.9136E-01  3.9108E+01  2.9174E+01  1.5470E+01  2.6025E+01
 ETASHRINKVR(%)  1.5765E+00  6.2921E+01  4.9837E+01  2.8548E+01  4.5276E+01
 EBVSHRINKSD(%)  1.0899E+00  3.9549E+01  3.4504E+01  1.6872E+01  1.9576E+01
 EBVSHRINKVR(%)  2.1679E+00  6.3456E+01  5.7102E+01  3.0897E+01  3.5320E+01
 RELATIVEINF(%)  9.7801E+01  9.9648E+00  2.5503E+01  2.1096E+01  3.4980E+01
 EPSSHRINKSD(%)  1.8758E+01
 EPSSHRINKVR(%)  3.3997E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          899
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1652.2514827020016     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3054.4247454624792     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1402.1732627604777     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    71.05
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    15.26
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3054.425       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.04E+00  2.19E+00  1.05E+00  1.18E+00  9.61E-01  9.70E-01  2.59E+00  9.47E-01  1.26E+00  2.04E+00
 


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
+        2.41E+07
 
 TH 2
+       -6.76E+02  1.45E+07
 
 TH 3
+        4.94E+01  7.57E+01  7.88E+04
 
 TH 4
+       -6.38E+02 -3.27E+02  2.60E+01  1.28E+07
 
 TH 5
+        1.52E+04  1.17E+04 -9.09E+02  1.11E+04  3.08E+06
 
 TH 6
+        9.42E+01  6.25E+01 -5.36E+00  6.28E+01 -5.96E-01  2.10E+02
 
 TH 7
+        1.06E+02  9.54E+01 -6.45E+00  7.07E+01  1.74E+04  9.78E+01  3.15E+07
 
 TH 8
+        3.46E+01  4.33E+01  1.38E+02  2.71E+01 -6.21E+02 -3.42E+00 -1.96E+00  4.00E+04
 
 TH 9
+        7.24E+02  5.37E+02 -3.95E+01  5.57E+02 -1.01E+07  1.03E+02  3.23E+07 -2.78E+01  3.31E+07
 
 TH10
+        2.26E+02  2.49E+02 -1.99E+02  2.29E+02 -4.02E+03 -2.41E+01 -2.32E+01 -1.40E+02 -1.88E+02  1.66E+06
 
 TH11
+       -7.06E+01 -9.20E+01 -5.09E+02 -6.17E+01  1.01E+03  8.68E+00  1.28E+01 -1.55E+02  5.66E+01  2.36E+02  1.08E+05
 
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
 #CPUT: Total CPU Time in Seconds,       86.400
Stop Time:
Sat Sep 25 00:51:09 CDT 2021
