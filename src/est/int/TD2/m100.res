Wed Sep 29 07:57:41 CDT 2021
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
$DATA ../../../../data/int/TD2/dat100.csv ignore=@
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
 (2E4.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m100.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3080.49390099738        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.2951E+02 -2.6144E-02  9.5807E+01  1.7705E+02  1.2307E+02  3.3058E+01 -1.3225E+01 -1.3419E+02 -3.8713E+01 -1.7095E+01
            -1.3243E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3348.82102197840        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.8007E-01  1.2310E+00  8.0483E-01  8.3688E-01  1.0411E+00  1.0312E+00  1.0267E+00  1.2107E+00  1.0496E+00  1.0856E+00
             1.6026E+00
 PARAMETER:  7.9874E-02  3.0780E-01 -1.1713E-01 -7.8074E-02  1.4031E-01  1.3068E-01  1.2638E-01  2.9124E-01  1.4845E-01  1.8212E-01
             5.7161E-01
 GRADIENT:   2.0545E+02  1.0788E+02 -5.4832E+01 -1.0622E+01 -1.0588E+01  2.7172E+01  2.5193E+01 -6.7169E+00  2.2845E+00 -1.3207E+01
             1.0158E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3354.11437389466        NO. OF FUNC. EVALS.: 110
 CUMULATIVE NO. OF FUNC. EVALS.:      194
 NPARAMETR:  9.7886E-01  1.2427E+00  7.9916E-01  8.3745E-01  1.0592E+00  1.0194E+00  9.7501E-01  1.4560E+00  1.0692E+00  1.1912E+00
             1.5855E+00
 PARAMETER:  7.8636E-02  3.1725E-01 -1.2419E-01 -7.7388E-02  1.5753E-01  1.1918E-01  7.4693E-02  4.7567E-01  1.6692E-01  2.7497E-01
             5.6088E-01
 GRADIENT:   3.2292E+01 -3.6329E+01 -5.5637E+01 -2.3396E+01 -3.0438E+01 -5.6303E+00  2.1160E+01 -7.3434E+00  4.5177E+00 -9.4898E+00
             9.9661E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3376.12279885020        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  9.7198E-01  1.5582E+00  1.2412E+00  6.7994E-01  1.4618E+00  1.0438E+00  9.3460E-01  1.7404E+00  8.8191E-01  1.4250E+00
             1.5821E+00
 PARAMETER:  7.1579E-02  5.4352E-01  3.1608E-01 -2.8575E-01  4.7964E-01  1.4291E-01  3.2369E-02  6.5410E-01 -2.5661E-02  4.5414E-01
             5.5875E-01
 GRADIENT:   1.5655E+01 -8.6349E+00  6.3787E+00  1.2657E+01  1.2363E+01  3.8787E+00  1.3942E+01 -1.7167E+01 -8.7758E+00 -5.5041E+00
             7.6966E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3378.09821769238        NO. OF FUNC. EVALS.: 125
 CUMULATIVE NO. OF FUNC. EVALS.:      495
 NPARAMETR:  9.5955E-01  1.5572E+00  1.2412E+00  6.7910E-01  1.4610E+00  1.0293E+00  8.6908E-01  1.7634E+00  1.0206E+00  1.4313E+00
             1.5814E+00
 PARAMETER:  5.8714E-02  5.4290E-01  3.1611E-01 -2.8699E-01  4.7913E-01  1.2889E-01 -4.0316E-02  6.6723E-01  1.2037E-01  4.5860E-01
             5.5831E-01
 GRADIENT:   1.6846E+02  3.0274E+02  5.5021E+00  6.6603E+01  8.7937E+01  2.9697E+01  1.2973E+01 -9.6106E+00 -5.3327E-01  1.3517E+01
             8.7896E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3380.36645410691        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      675
 NPARAMETR:  9.6403E-01  1.5572E+00  1.2412E+00  6.7910E-01  1.4610E+00  1.0323E+00  7.7557E-01  1.7634E+00  1.1648E+00  1.4313E+00
             1.5813E+00
 PARAMETER:  6.3371E-02  5.4291E-01  3.1610E-01 -2.8699E-01  4.7913E-01  1.3181E-01 -1.5416E-01  6.6725E-01  2.5257E-01  4.5860E-01
             5.5826E-01
 GRADIENT:  -9.3251E-01 -2.5534E+01 -7.6710E+00  1.6980E+01  1.3175E+01 -2.7846E-01  2.6041E+00 -8.0579E+00  2.0563E+00  1.9752E+00
             8.2763E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3381.75272141252        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      853
 NPARAMETR:  9.6467E-01  1.5713E+00  1.3528E+00  6.6886E-01  1.4777E+00  1.0320E+00  7.6464E-01  1.7606E+00  1.1871E+00  1.4142E+00
             1.5833E+00
 PARAMETER:  6.4033E-02  5.5192E-01  4.0216E-01 -3.0218E-01  4.9049E-01  1.3145E-01 -1.6835E-01  6.6563E-01  2.7151E-01  4.4654E-01
             5.5948E-01
 GRADIENT:   5.0525E-01 -2.0349E+01 -4.7425E+00  1.0214E+01  3.1555E+00 -3.8399E-01  1.2923E+00 -9.6730E+00  1.1557E-01 -7.4479E-01
             8.3524E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3383.30228088142        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1029
 NPARAMETR:  9.6452E-01  1.7626E+00  1.1843E+00  5.4693E-01  1.5880E+00  1.0341E+00  7.3989E-01  1.7626E+00  1.2946E+00  1.4979E+00
             1.5810E+00
 PARAMETER:  6.3880E-02  6.6679E-01  2.6916E-01 -5.0344E-01  5.6245E-01  1.3353E-01 -2.0125E-01  6.6677E-01  3.5821E-01  5.0406E-01
             5.5803E-01
 GRADIENT:   1.0570E-01  1.9302E+00 -1.7826E-01  6.6478E-01  6.4989E-01  1.7025E-01 -6.2582E-02 -9.3108E+00  9.5994E-02 -2.4693E-02
             7.4699E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3385.55150235304        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:     1147
 NPARAMETR:  9.6423E-01  1.7551E+00  1.2278E+00  5.5004E-01  1.5857E+00  1.0332E+00  7.3889E-01  2.0485E+00  1.2853E+00  1.4988E+00
             1.5574E+00
 PARAMETER:  6.3574E-02  6.6251E-01  3.0524E-01 -4.9777E-01  5.6100E-01  1.3264E-01 -2.0261E-01  8.1709E-01  3.5103E-01  5.0467E-01
             5.4299E-01
 GRADIENT:   1.8507E+02  4.9128E+02 -2.0790E+00  6.3098E+01  8.0842E+01  3.2624E+01  8.0311E+00 -4.9577E+00  7.1662E+00  1.5660E+01
             5.9350E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3388.12426922441        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:     1266             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6634E-01  1.7410E+00  1.4973E+00  5.5469E-01  1.6230E+00  1.0373E+00  7.4736E-01  2.5529E+00  1.2357E+00  1.5152E+00
             1.5075E+00
 PARAMETER:  6.5765E-02  6.5449E-01  5.0367E-01 -4.8934E-01  5.8427E-01  1.3660E-01 -1.9121E-01  1.0372E+00  3.1167E-01  5.1552E-01
             5.1042E-01
 GRADIENT:   2.0325E+02  4.9765E+02  4.2514E+00  5.6399E+01  8.7343E+01  3.7007E+01  8.1218E+00 -2.7385E+00  6.7443E+00  1.6708E+01
            -2.7098E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3388.14869280731        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     1422
 NPARAMETR:  9.6382E-01  1.7410E+00  1.4610E+00  5.5469E-01  1.6230E+00  1.0346E+00  7.4972E-01  2.5528E+00  1.2209E+00  1.5152E+00
             1.5075E+00
 PARAMETER:  6.3144E-02  6.5447E-01  4.7915E-01 -4.8935E-01  5.8430E-01  1.3405E-01 -1.8805E-01  1.0372E+00  2.9956E-01  5.1554E-01
             5.1047E-01
 GRADIENT:   2.2368E-01 -1.4260E+01 -1.7217E-01 -1.3494E+01 -2.1407E+00  9.3559E-02 -1.9806E-01 -6.2354E+00 -1.1906E-02 -9.4180E-01
            -1.1462E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -3388.48598558444        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:     1562
 NPARAMETR:  9.6198E-01  1.7361E+00  1.4860E+00  5.7313E-01  1.6254E+00  1.0325E+00  7.5146E-01  2.6104E+00  1.1991E+00  1.5201E+00
             1.5157E+00
 PARAMETER:  6.1240E-02  6.5162E-01  4.9609E-01 -4.5663E-01  5.8577E-01  1.3197E-01 -1.8574E-01  1.0595E+00  2.8155E-01  5.1876E-01
             5.1585E-01
 GRADIENT:   1.9112E+02  5.0434E+02  9.1391E-01  7.1847E+01  9.4416E+01  3.3963E+01  6.3291E+00 -1.8823E+00  6.0693E+00  1.8706E+01
             1.1757E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -3388.53848034827        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1741
 NPARAMETR:  9.6391E-01  1.7361E+00  1.5174E+00  5.7204E-01  1.6254E+00  1.0348E+00  7.6662E-01  2.6106E+00  1.1686E+00  1.5187E+00
             1.5144E+00
 PARAMETER:  6.3241E-02  6.5162E-01  5.1700E-01 -4.5854E-01  5.8575E-01  1.3418E-01 -1.6577E-01  1.0596E+00  2.5580E-01  5.1787E-01
             5.1504E-01
 GRADIENT:   1.1497E-01  1.0232E+01  5.5544E-03 -2.8212E-03  8.8352E-02  9.4323E-02  3.8260E-02 -6.8977E+00 -4.4574E-02  1.3886E-02
             1.4881E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -3388.66493547776        NO. OF FUNC. EVALS.: 143
 CUMULATIVE NO. OF FUNC. EVALS.:     1884
 NPARAMETR:  9.6427E-01  1.7243E+00  1.5228E+00  5.7733E-01  1.6214E+00  1.0347E+00  7.6648E-01  2.6317E+00  1.1656E+00  1.5181E+00
             1.5134E+00
 PARAMETER:  6.3613E-02  6.4482E-01  5.2054E-01 -4.4934E-01  5.8326E-01  1.3408E-01 -1.6594E-01  1.0676E+00  2.5325E-01  5.1743E-01
             5.1434E-01
 GRADIENT:   1.9680E+02  4.8780E+02  2.3894E+00  6.7370E+01  9.2553E+01  3.5252E+01  6.6536E+00 -2.5595E+00  4.9589E+00  1.8363E+01
             9.3754E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -3388.76241903282        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:     2047
 NPARAMETR:  9.6299E-01  1.7195E+00  1.5381E+00  5.7815E-01  1.6180E+00  1.0338E+00  7.6639E-01  2.6634E+00  1.1629E+00  1.5150E+00
             1.5133E+00
 PARAMETER:  6.2289E-02  6.4205E-01  5.3056E-01 -4.4792E-01  5.8120E-01  1.3325E-01 -1.6607E-01  1.0796E+00  2.5090E-01  5.1544E-01
             5.1427E-01
 GRADIENT:  -1.7644E+00 -3.2661E-01 -3.7047E-01 -4.3684E+00 -5.7200E-01 -2.4626E-01 -2.1332E-01 -6.1794E+00  2.8434E-01  1.8474E-01
             9.2428E-01

0ITERATION NO.:   71    OBJECTIVE VALUE:  -3388.76241903282        NO. OF FUNC. EVALS.:  27
 CUMULATIVE NO. OF FUNC. EVALS.:     2074
 NPARAMETR:  9.6383E-01  1.7171E+00  1.5399E+00  5.7758E-01  1.6201E+00  1.0345E+00  7.6713E-01  2.6571E+00  1.1600E+00  1.5168E+00
             1.5150E+00
 PARAMETER:  6.2289E-02  6.4205E-01  5.3056E-01 -4.4792E-01  5.8120E-01  1.3325E-01 -1.6607E-01  1.0796E+00  2.5090E-01  5.1544E-01
             5.1427E-01
 GRADIENT:  -2.0253E+00  4.8778E+03 -2.9507E+03  6.9909E+03 -5.3955E+03 -2.7906E-01 -1.6525E-01  2.8486E+03  2.2215E-01 -6.0830E+03
            -6.1044E+03
 NUMSIGDIG:         1.7         2.3         2.3         2.3         2.3         1.9         1.9         2.3         1.6         2.3
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2074
 NO. OF SIG. DIGITS IN FINAL EST.:  1.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.9444E-03 -3.0846E-02 -2.6810E-02  3.2879E-02 -2.6341E-02
 SE:             2.9813E-02  2.3060E-02  1.6542E-02  2.1824E-02  2.6063E-02
 N:                     100         100         100         100         100

 P VAL.:         9.4800E-01  1.8101E-01  1.0508E-01  1.3192E-01  3.1217E-01

 ETASHRINKSD(%)  1.2123E-01  2.2746E+01  4.4582E+01  2.6888E+01  1.2685E+01
 ETASHRINKVR(%)  2.4230E-01  4.0318E+01  6.9289E+01  4.6547E+01  2.3760E+01
 EBVSHRINKSD(%)  5.5843E-01  2.1518E+01  5.2688E+01  3.2930E+01  8.9380E+00
 EBVSHRINKVR(%)  1.1137E+00  3.8406E+01  7.7616E+01  5.5016E+01  1.7077E+01
 RELATIVEINF(%)  9.8880E+01  1.1403E+01  1.2878E+01  7.9354E+00  4.8886E+01
 EPSSHRINKSD(%)  1.9293E+01
 EPSSHRINKVR(%)  3.4864E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3388.7624190328193     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1734.6730592644085     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    62.53
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    15.91
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3388.762       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.63E-01  1.72E+00  1.54E+00  5.78E-01  1.62E+00  1.03E+00  7.66E-01  2.66E+00  1.16E+00  1.52E+00  1.51E+00
 


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
+        8.45E+06
 
 TH 2
+        6.14E+00  6.45E+04
 
 TH 3
+        9.96E+05  2.74E+02  1.17E+05
 
 TH 4
+        3.73E+01  1.39E+02  1.07E+03  1.17E+06
 
 TH 5
+       -1.33E+01 -6.95E+01 -3.25E+02  5.26E+02  8.88E+04
 
 TH 6
+        4.87E+00  1.13E+01 -1.79E+01  5.42E+01 -1.62E+01  1.82E+02
 
 TH 7
+        1.83E+00 -3.57E+01  5.07E+01 -1.51E+02  3.95E+01  1.33E-01  1.29E+02
 
 TH 8
+        3.93E+00  2.42E+04  9.47E+01 -1.38E+02 -9.00E+00  5.01E+00 -1.44E+01  9.16E+03
 
 TH 9
+       -2.79E+06 -9.16E+01  1.23E+02 -3.65E+02  1.12E+02 -2.55E-02  4.43E+01 -3.00E+01  2.90E+01
 
 TH10
+       -1.48E+01 -1.25E+01 -3.69E+02 -3.87E+05 -6.23E+00 -1.86E+01  4.90E+01 -1.85E+00  1.35E+02  1.29E+05
 
 TH11
+       -2.44E+01 -2.72E+02  1.23E+05 -3.89E+05  2.89E+02 -1.72E+01  5.81E+01 -8.68E+01  1.36E+02  1.29E+05  1.30E+05
 
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
 #CPUT: Total CPU Time in Seconds,       78.497
Stop Time:
Wed Sep 29 07:59:02 CDT 2021
