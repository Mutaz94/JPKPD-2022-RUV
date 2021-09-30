Thu Sep 30 06:29:33 CDT 2021
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
$DATA ../../../../data/spa2/A3/dat49.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m49.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   346.125469060121        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.4842E+02  9.3829E+01  2.9839E+02  1.0186E+01  2.3885E+02  2.4760E+01 -2.4222E+02 -3.3118E+02 -1.2386E+02 -1.5841E+02
            -4.6174E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1699.56650243944        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.6675E-01  1.0812E+00  8.3578E-01  1.0145E+00  8.8911E-01  8.8657E-01  1.6692E+00  9.6204E-01  1.2724E+00  9.7416E-01
             2.8137E+00
 PARAMETER:  6.6189E-02  1.7810E-01 -7.9389E-02  1.1440E-01 -1.7539E-02 -2.0396E-02  6.1234E-01  6.1301E-02  3.4089E-01  7.3819E-02
             1.1345E+00
 GRADIENT:   1.2302E+02  4.2957E+01  1.5263E+01  1.3073E+01 -1.4839E+01 -2.4269E+01  2.5995E+01  1.0154E+01  3.1608E+01  1.5765E+00
            -5.0503E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1711.14051409011        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      202
 NPARAMETR:  9.4810E-01  6.6406E-01  5.1609E-01  1.2640E+00  5.0362E-01  9.1948E-01  2.5069E+00  3.7953E-01  9.7230E-01  7.0094E-01
             2.6285E+00
 PARAMETER:  4.6703E-02 -3.0938E-01 -5.6147E-01  3.3428E-01 -5.8593E-01  1.6055E-02  1.0190E+00 -8.6882E-01  7.1906E-02 -2.5533E-01
             1.0664E+00
 GRADIENT:   2.1342E+01  6.6796E+01  3.9082E+01  8.7170E+01 -3.4126E+01 -1.8254E+01  2.4274E+01  3.4658E+00 -3.9292E+01 -1.2976E+01
            -7.4239E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1732.78292214584        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      379
 NPARAMETR:  9.3074E-01  5.1461E-01  3.8429E-01  1.2692E+00  4.2153E-01  9.5428E-01  1.6459E+00  1.0000E-02  1.0844E+00  6.1582E-01
             2.8520E+00
 PARAMETER:  2.8227E-02 -5.6434E-01 -8.5636E-01  3.3838E-01 -7.6385E-01  5.3200E-02  5.9827E-01 -4.9849E+00  1.8101E-01 -3.8480E-01
             1.1480E+00
 GRADIENT:  -3.3189E+01  1.0936E-01 -1.1111E+01  8.7270E+01  5.7874E+01 -5.4111E+00 -2.6785E+01  0.0000E+00 -3.2295E+01 -1.3379E+01
             5.3634E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1760.25545026034        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      558
 NPARAMETR:  9.3475E-01  3.2525E-01  2.1313E-01  1.1522E+00  2.5697E-01  9.7909E-01  1.7516E+00  1.0000E-02  1.4231E+00  6.9652E-01
             2.5163E+00
 PARAMETER:  3.2522E-02 -1.0231E+00 -1.4459E+00  2.4167E-01 -1.2588E+00  7.8865E-02  6.6052E-01 -1.3953E+01  4.5281E-01 -2.6165E-01
             1.0228E+00
 GRADIENT:  -7.2187E+00 -9.4052E-01  1.1118E+00  8.2515E+00  1.4561E+01 -1.1321E+00  2.7732E+00  0.0000E+00  4.3666E+00 -4.0558E+00
            -6.1144E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1760.56862587417        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:      744
 NPARAMETR:  9.3814E-01  3.2504E-01  2.1107E-01  1.1485E+00  2.5538E-01  9.8160E-01  1.6696E+00  1.0000E-02  1.4265E+00  7.2915E-01
             2.5202E+00
 PARAMETER:  3.6147E-02 -1.0238E+00 -1.4556E+00  2.3844E-01 -1.2650E+00  8.1430E-02  6.1258E-01 -1.4131E+01  4.5523E-01 -2.1587E-01
             1.0244E+00
 GRADIENT:   1.9657E-01 -1.2997E+00  7.4330E-01  7.2936E+00  1.8242E+01  4.0608E-02 -4.7326E-01  0.0000E+00  3.9589E+00 -7.5096E-01
            -2.7435E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1760.76752565288        NO. OF FUNC. EVALS.: 155
 CUMULATIVE NO. OF FUNC. EVALS.:      899
 NPARAMETR:  9.3416E-01  3.2157E-01  2.0859E-01  1.1319E+00  2.4958E-01  9.7857E-01  1.6566E+00  1.0000E-02  1.3935E+00  7.3196E-01
             2.5205E+00
 PARAMETER:  3.1890E-02 -1.0345E+00 -1.4674E+00  2.2391E-01 -1.2880E+00  7.8341E-02  6.0479E-01 -1.4131E+01  4.3182E-01 -2.1203E-01
             1.0245E+00
 GRADIENT:   3.6763E+01  3.2604E+01  6.6369E+01  2.7511E+01  1.7277E+02  4.4785E+00  3.8920E+00  0.0000E+00  1.2853E+01 -6.4638E-01
             9.5029E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1760.97144930151        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:     1040
 NPARAMETR:  9.3029E-01  3.0951E-01  2.0106E-01  1.1199E+00  2.4394E-01  9.7477E-01  1.6367E+00  1.0000E-02  1.3781E+00  7.4051E-01
             2.5159E+00
 PARAMETER:  2.7736E-02 -1.0728E+00 -1.5042E+00  2.1328E-01 -1.3108E+00  7.4442E-02  5.9265E-01 -1.4131E+01  4.2074E-01 -2.0042E-01
             1.0226E+00
 GRADIENT:  -2.0276E+01 -3.5039E+00  8.4481E+00 -6.7699E+00  3.1527E+00 -2.5322E+00 -1.2013E+00  0.0000E+00 -1.1723E+01 -2.2849E+00
            -1.1517E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1762.13653993374        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1229             RESET HESSIAN, TYPE I
 NPARAMETR:  9.3565E-01  2.8975E-01  1.8180E-01  1.0902E+00  2.2670E-01  9.7636E-01  1.5889E+00  1.0000E-02  1.4441E+00  7.8462E-01
             2.4998E+00
 PARAMETER:  3.3483E-02 -1.1388E+00 -1.6048E+00  1.8636E-01 -1.3841E+00  7.6073E-02  5.6303E-01 -1.4131E+01  4.6746E-01 -1.4256E-01
             1.0162E+00
 GRADIENT:   3.9371E+01  2.4656E+01  6.3112E+01  1.3850E+01  2.2467E+02  3.0786E+00  4.8917E+00  0.0000E+00  1.1504E+01  2.3535E+00
             1.3106E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1762.21763331284        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1405
 NPARAMETR:  9.3709E-01  2.9298E-01  1.8087E-01  1.0960E+00  2.2749E-01  9.7905E-01  1.5612E+00  1.0000E-02  1.4424E+00  7.8413E-01
             2.5062E+00
 PARAMETER:  3.5027E-02 -1.1277E+00 -1.6100E+00  1.9166E-01 -1.3807E+00  7.8832E-02  5.4543E-01 -1.4131E+01  4.6630E-01 -1.4318E-01
             1.0188E+00
 GRADIENT:  -3.2457E+00 -2.1816E+00 -3.1334E-01 -4.6506E+00  1.0587E+01 -1.1621E+00 -9.0318E-02  0.0000E+00 -8.2792E+00  5.6404E-01
             3.8269E+00

0ITERATION NO.:   47    OBJECTIVE VALUE:  -1762.24193068674        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1462
 NPARAMETR:  9.3840E-01  2.9427E-01  1.8091E-01  1.1013E+00  2.2748E-01  9.8211E-01  1.5623E+00  1.0000E-02  1.4424E+00  7.8147E-01
             2.5061E+00
 PARAMETER:  3.6416E-02 -1.1232E+00 -1.6098E+00  1.9648E-01 -1.3807E+00  8.1949E-02  5.4615E-01 -1.4131E+01  4.6631E-01 -1.4658E-01
             1.0187E+00
 GRADIENT:  -6.3802E-02 -6.2952E-02 -5.8288E-01 -2.6543E-01  7.3529E+00 -1.1899E-01 -6.3758E-03  0.0000E+00 -8.3968E+00  1.6242E-01
             3.4711E+00

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1462
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5239E-03  1.7301E-02  4.3553E-06 -9.2950E-03  7.8762E-03
 SE:             2.9174E-02  2.0996E-02  2.0011E-04  2.8180E-02  2.4364E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5834E-01  4.0992E-01  9.8264E-01  7.4152E-01  7.4649E-01

 ETASHRINKSD(%)  2.2640E+00  2.9662E+01  9.9330E+01  5.5936E+00  1.8379E+01
 ETASHRINKVR(%)  4.4767E+00  5.0525E+01  9.9996E+01  1.0874E+01  3.3379E+01
 EBVSHRINKSD(%)  2.1835E+00  2.7967E+01  9.9295E+01  6.0611E+00  1.9087E+01
 EBVSHRINKVR(%)  4.3193E+00  4.8112E+01  9.9995E+01  1.1755E+01  3.4531E+01
 RELATIVEINF(%)  9.5510E+01  1.2049E+01  4.0931E-04  5.2272E+01  4.1142E+00
 EPSSHRINKSD(%)  2.7776E+01
 EPSSHRINKVR(%)  4.7838E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1762.2419306867398     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -659.51569084113271     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    29.10
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.91
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1762.242       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.38E-01  2.94E-01  1.81E-01  1.10E+00  2.27E-01  9.82E-01  1.56E+00  1.00E-02  1.44E+00  7.81E-01  2.51E+00
 


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
+        1.25E+03
 
 TH 2
+        2.73E+01  2.76E+03
 
 TH 3
+       -1.02E+02  6.54E+02  1.50E+04
 
 TH 4
+       -3.56E+00 -8.37E+00 -4.23E+02  3.91E+02
 
 TH 5
+        4.54E+06 -4.03E+03  1.44E+06 -1.97E+06  1.38E+06
 
 TH 6
+        5.23E+00 -9.26E+00  5.42E+01 -7.04E+00  4.34E+06  1.89E+02
 
 TH 7
+        1.26E+00  5.54E+01 -5.65E+01 -2.76E+00 -4.74E+01 -1.06E+00  2.27E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        6.53E+02 -2.78E+02  7.07E+03 -9.97E+02 -6.33E+05  7.29E+02  1.23E+01  0.00E+00  2.96E+05
 
 TH10
+       -7.86E+00  2.14E+01 -7.32E+00  9.66E+00  3.72E+06  3.12E+00  2.04E+01  0.00E+00  4.68E+02  1.35E+02
 
 TH11
+       -1.77E+02  6.18E+01 -1.96E+03  2.46E+02 -2.74E+02 -1.80E+02  1.59E+00  0.00E+00 -7.94E+04 -1.06E+02  2.14E+04
 
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
 #CPUT: Total CPU Time in Seconds,       39.100
Stop Time:
Thu Sep 30 06:30:13 CDT 2021
