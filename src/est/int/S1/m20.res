Fri Sep 24 23:07:46 CDT 2021
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
$DATA ../../../../data/int/S1/dat20.csv ignore=@
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
Current Date:       24 SEP 2021
Days until program expires : 205
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
 RAW OUTPUT FILE (FILE): m20.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3139.00614259868        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.2008E+02  2.0962E+01  8.7534E+01 -1.3163E+01  9.3241E+01  3.5704E-01 -3.2522E+00 -2.4999E+02 -5.2529E+01 -1.4579E+01
            -1.0835E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3556.48415818220        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  8.7534E-01  8.7611E-01  9.9085E-01  1.0809E+00  8.6646E-01  1.0484E+00  9.3043E-01  9.5251E-01  1.0624E+00  1.0447E+00
             1.6846E+00
 PARAMETER: -3.3145E-02 -3.2261E-02  9.0808E-02  1.7776E-01 -4.3338E-02  1.4725E-01  2.7897E-02  5.1342E-02  1.6054E-01  1.4371E-01
             6.2150E-01
 GRADIENT:  -1.9293E+02 -3.4885E+01  1.9959E+01  5.3344E+01 -1.9043E+01 -2.7114E+00  7.2108E+00  1.4178E+01  1.7840E+01  2.3611E+01
             5.9646E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3613.83626332657        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.3308E-01  1.0511E+00  1.0104E+00  9.4899E-01  1.0179E+00  8.5147E-01  8.3424E-01  9.9963E-02  1.1107E+00  1.0218E+00
             1.5394E+00
 PARAMETER:  3.0735E-02  1.4982E-01  1.1032E-01  4.7639E-02  1.1774E-01 -6.0787E-02 -8.1232E-02 -2.2030E+00  2.0499E-01  1.2158E-01
             5.3138E-01
 GRADIENT:  -1.0155E+02  4.1539E-01  2.6728E+01 -1.2836E+01  3.7209E+01 -8.0557E+01 -1.4234E+00 -2.0942E-01  9.2401E+00 -1.5679E+01
             4.5886E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3696.96974469747        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      234
 NPARAMETR:  9.4743E-01  9.1169E-01  8.7072E-01  1.0135E+00  8.5933E-01  9.8170E-01  9.5113E-01  4.4370E-01  1.0502E+00  9.9613E-01
             1.1597E+00
 PARAMETER:  4.5996E-02  7.5432E-03 -3.8439E-02  1.1338E-01 -5.1602E-02  8.1527E-02  4.9895E-02 -7.1260E-01  1.4901E-01  9.6121E-02
             2.4813E-01
 GRADIENT:  -1.3557E+01  2.9519E+00 -2.5966E+01 -9.8745E+00  8.0732E-02 -7.9228E+00  6.6242E-01 -6.4621E+00 -1.1103E+00  3.0150E+00
             1.4440E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3697.02340482283        NO. OF FUNC. EVALS.: 147
 CUMULATIVE NO. OF FUNC. EVALS.:      381             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8854E-01  9.1165E-01  8.7077E-01  1.0136E+00  8.5920E-01  1.0517E+00  9.5080E-01  4.4525E-01  1.0503E+00  9.9619E-01
             1.1595E+00
 PARAMETER:  8.8473E-02  7.5008E-03 -3.8382E-02  1.1353E-01 -5.1752E-02  1.5045E-01  4.9552E-02 -7.0912E-01  1.4909E-01  9.6180E-02
             2.4801E-01
 GRADIENT:   7.7646E+01  2.6957E+00 -2.5674E+01 -1.0125E+01 -3.8604E-01  2.1694E+01  6.2024E-01 -6.5046E+00 -9.4639E-01  3.1188E+00
             7.6728E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3697.72763102145        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      496
 NPARAMETR:  9.6742E-01  9.1165E-01  8.7077E-01  1.0136E+00  8.5920E-01  1.0159E+00  9.5080E-01  4.4525E-01  1.0503E+00  9.9619E-01
             1.1595E+00
 PARAMETER:  6.6882E-02  7.5009E-03 -3.8382E-02  1.1353E-01 -5.1753E-02  1.1577E-01  4.9552E-02 -7.0912E-01  1.4909E-01  9.6180E-02
             2.4801E-01
 GRADIENT:  -4.5002E-01 -3.1417E+00 -2.6363E+01 -1.8930E+01 -4.6091E+00  4.5924E-01  1.7390E-01 -6.5753E+00 -2.7217E+00  2.6932E+00
             6.0720E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3697.88174855332        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:      662
 NPARAMETR:  9.6763E-01  9.1166E-01  8.7083E-01  1.0137E+00  8.5921E-01  1.0147E+00  9.5080E-01  4.5549E-01  1.0503E+00  9.9614E-01
             1.1595E+00
 PARAMETER:  6.7098E-02  7.5091E-03 -3.8305E-02  1.1357E-01 -5.1743E-02  1.1459E-01  4.9547E-02 -6.8639E-01  1.4912E-01  9.6133E-02
             2.4801E-01
 GRADIENT:   3.9048E-03 -3.1145E+00 -2.6359E+01 -1.8967E+01 -4.7171E+00 -3.3585E-03  1.9523E-01 -6.6892E+00 -2.6409E+00  2.8846E+00
             1.7971E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3698.32772264194        NO. OF FUNC. EVALS.: 104
 CUMULATIVE NO. OF FUNC. EVALS.:      766
 NPARAMETR:  9.5895E-01  9.1686E-01  9.0310E-01  1.0229E+00  8.6894E-01  1.0052E+00  9.5074E-01  4.5533E-01  1.0549E+00  9.8363E-01
             1.1597E+00
 PARAMETER:  5.8080E-02  1.3197E-02 -1.9261E-03  1.2260E-01 -4.0487E-02  1.0520E-01  4.9486E-02 -6.8673E-01  1.5343E-01  8.3495E-02
             2.4814E-01
 GRADIENT:   1.5017E+01  3.5839E+00  7.9037E+00  1.0540E+01 -1.0380E+01  2.8355E+00 -5.6845E-01 -7.1456E+00  2.0552E+00 -2.0272E+00
             8.8581E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3698.62309762343        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      923
 NPARAMETR:  9.6803E-01  9.3156E-01  9.1068E-01  1.0143E+00  8.8744E-01  1.0146E+00  9.5921E-01  4.5533E-01  1.0565E+00  1.0085E+00
             1.1597E+00
 PARAMETER:  6.7503E-02  2.9106E-02  6.4325E-03  1.1424E-01 -1.9412E-02  1.1452E-01  5.8357E-02 -6.8673E-01  1.5494E-01  1.0847E-01
             2.4814E-01
 GRADIENT:   6.4956E-01 -2.3003E+00 -1.2257E+00 -1.1453E+00 -6.2799E-01 -2.2720E-02  1.1565E+00 -6.9751E+00  2.3620E-01  9.3500E-01
             4.1546E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3698.75220924172        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1105            RESET HESSIAN, TYPE II
 NPARAMETR:  9.6775E-01  9.3666E-01  9.1459E-01  1.0130E+00  8.9181E-01  1.0147E+00  9.4926E-01  4.6272E-01  1.0564E+00  1.0093E+00
             1.1597E+00
 PARAMETER:  6.7218E-02  3.4566E-02  1.0724E-02  1.1289E-01 -1.4507E-02  1.1459E-01  4.7927E-02 -6.7064E-01  1.5487E-01  1.0928E-01
             2.4814E-01
 GRADIENT:   3.4947E+01  6.1648E+00  3.4954E-01  8.7522E+00  4.2368E+00  7.1868E+00  4.6856E-01 -7.0233E+00  1.7971E+00  6.1565E-01
             6.6569E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3698.81612825653        NO. OF FUNC. EVALS.: 147
 CUMULATIVE NO. OF FUNC. EVALS.:     1252
 NPARAMETR:  9.6773E-01  9.3665E-01  9.1459E-01  1.0130E+00  8.9181E-01  1.0147E+00  9.4926E-01  4.6688E-01  1.0564E+00  1.0092E+00
             1.1597E+00
 PARAMETER:  6.7198E-02  3.4561E-02  1.0724E-02  1.1288E-01 -1.4510E-02  1.1457E-01  4.7927E-02 -6.6166E-01  1.5486E-01  1.0928E-01
             2.4814E-01
 GRADIENT:  -6.2482E-02  3.0840E+06  2.2949E+05 -1.3660E+06 -3.0840E+06 -1.7824E-02 -8.7525E-01  4.6596E+05  9.9564E+05  2.6236E-01
            -1.2432E+06
 NUMSIGDIG:         3.5         3.3         3.3         3.3         3.3         3.4        10.1         3.3         3.3         1.7
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1252
 NO. OF SIG. DIGITS IN FINAL EST.:  1.7

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.0130E-04 -2.3913E-02 -8.9764E-03  9.8719E-03 -1.8006E-02
 SE:             2.9863E-02  2.1755E-02  1.1856E-02  2.8025E-02  2.5792E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8928E-01  2.7169E-01  4.4898E-01  7.2465E-01  4.8510E-01

 ETASHRINKSD(%)  1.0000E-10  2.7117E+01  6.0281E+01  6.1122E+00  1.3594E+01
 ETASHRINKVR(%)  1.0000E-10  4.6881E+01  8.4224E+01  1.1851E+01  2.5339E+01
 EBVSHRINKSD(%)  3.3767E-01  2.7214E+01  6.5781E+01  6.5281E+00  1.3544E+01
 EBVSHRINKVR(%)  6.7420E-01  4.7022E+01  8.8291E+01  1.2630E+01  2.5254E+01
 RELATIVEINF(%)  9.9322E+01  2.5783E+01  6.8076E+00  6.4102E+01  2.8961E+01
 EPSSHRINKSD(%)  2.0012E+01
 EPSSHRINKVR(%)  3.6019E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3698.8161282565293     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2044.7267684881185     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    31.11
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.12
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3698.816       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.68E-01  9.37E-01  9.15E-01  1.01E+00  8.92E-01  1.01E+00  9.49E-01  4.67E-01  1.06E+00  1.01E+00  1.16E+00
 


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
+        8.23E+09
 
 TH 2
+       -8.51E+09  1.76E+10
 
 TH 3
+        9.36E+09 -4.16E+01  9.90E+09
 
 TH 4
+       -4.98E+04 -1.89E+04  2.56E+05  5.90E+09
 
 TH 5
+       -6.39E+04  1.08E+04 -9.45E+09  7.56E+09  9.69E+09
 
 TH 6
+       -1.24E+00  7.08E+09  1.09E+00 -4.51E+04 -5.79E+04  1.88E+02
 
 TH 7
+       -8.39E+09  2.54E+03  4.86E+05 -1.42E+10  9.11E+09 -6.99E+09  1.27E+09
 
 TH 8
+        1.84E+04 -7.10E+03  1.49E+05 -7.58E+04  3.40E+03  1.60E+08 -7.59E+02  8.07E+08
 
 TH 9
+        7.74E+05  1.01E+10  5.15E+09  4.12E+09 -5.37E+00  3.15E+04  4.96E+09 -2.42E+05  1.36E+02
 
 TH10
+        7.22E+09  7.46E+09  1.53E+10  5.66E+09 -7.84E+09  6.01E+09  7.36E+09 -1.42E+03  4.27E+09  6.34E+09
 
 TH11
+       -1.98E+04  2.86E+09 -1.60E+05  8.14E+04  3.00E+09 -1.79E+04  8.42E+02 -4.55E+05  1.64E+09  1.54E+03  9.32E+08
 
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
 #CPUT: Total CPU Time in Seconds,       44.309
Stop Time:
Fri Sep 24 23:08:32 CDT 2021
