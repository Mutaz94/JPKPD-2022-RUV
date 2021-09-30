Thu Sep 30 04:30:10 CDT 2021
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
$DATA ../../../../data/spa2/B/dat86.csv ignore=@
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m86.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2228.08190567907        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.8512E+02 -1.4977E+00  8.5919E+01 -6.6729E+00  7.8204E+01  5.2306E+01  4.7972E+00 -4.1028E+02 -6.8956E+01 -5.0538E+00
            -5.6077E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2380.04865689923        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:       88
 NPARAMETR:  1.0165E+00  9.6100E-01  9.7767E-01  1.0176E+00  9.9301E-01  1.0556E+00  9.6700E-01  1.7762E+00  8.7746E-01  9.6149E-01
             9.4972E-01
 PARAMETER:  1.1636E-01  6.0219E-02  7.7419E-02  1.1748E-01  9.2985E-02  1.5406E-01  6.6444E-02  6.7449E-01 -3.0719E-02  6.0725E-02
             4.8408E-02
 GRADIENT:   5.5841E+02 -2.2647E+01 -3.3167E+00  2.8335E+01  9.2292E+01  1.0570E+02 -8.2498E+00 -1.1882E+02 -4.0979E+00 -6.1543E+00
            -8.3748E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2382.15010996225        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:      282
 NPARAMETR:  1.0168E+00  9.6076E-01  9.7792E-01  1.0180E+00  9.4706E-01  9.8166E-01  9.6676E-01  1.7739E+00  8.7769E-01  9.6173E-01
             9.4948E-01
 PARAMETER:  1.1666E-01  5.9967E-02  7.7673E-02  1.1779E-01  4.5611E-02  8.1491E-02  6.6190E-02  6.7319E-01 -3.0463E-02  6.0981E-02
             4.8158E-02
 GRADIENT:  -2.6911E+01 -5.6373E+01  8.1900E+00 -9.5495E+01 -2.6135E-01 -3.0018E+00 -1.3327E+01 -1.4725E+02 -1.1813E+01 -3.3042E+00
            -8.5143E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2384.71623559704        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:      450
 NPARAMETR:  1.0154E+00  9.6771E-01  9.7692E-01  1.0184E+00  9.5901E-01  9.6495E-01  9.8470E-01  1.7897E+00  8.8341E-01  9.6538E-01
             9.5865E-01
 PARAMETER:  1.1527E-01  6.7178E-02  7.6647E-02  1.1828E-01  5.8149E-02  6.4320E-02  8.4578E-02  6.8203E-01 -2.3970E-02  6.4767E-02
             5.7772E-02
 GRADIENT:   5.3382E+02  6.5259E+00  8.8096E+00  2.7741E+01  3.8628E+01  4.1874E+01 -3.4573E+00 -1.1135E+02 -2.1989E+00 -6.2321E-01
            -7.2469E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2386.13973794704        NO. OF FUNC. EVALS.: 155
 CUMULATIVE NO. OF FUNC. EVALS.:      605
 NPARAMETR:  1.0138E+00  9.6760E-01  9.7566E-01  1.0184E+00  9.5901E-01  9.6129E-01  9.8249E-01  1.8076E+00  8.8270E-01  9.6436E-01
             9.6065E-01
 PARAMETER:  1.1373E-01  6.7061E-02  7.5355E-02  1.1827E-01  5.8147E-02  6.0517E-02  8.2338E-02  6.9200E-01 -2.4775E-02  6.3707E-02
             5.9853E-02
 GRADIENT:  -3.5598E+01 -5.5368E+01  1.5964E+00 -8.2664E+01  1.1093E+01 -1.1778E+01 -1.0930E+01 -1.3962E+02 -8.2661E+00 -3.2439E+00
            -7.1226E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2386.45145668586        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:      794
 NPARAMETR:  1.0137E+00  9.6776E-01  9.7579E-01  1.0183E+00  9.5888E-01  9.9563E-01  9.8263E-01  1.8094E+00  8.8258E-01  9.6440E-01
             9.6082E-01
 PARAMETER:  1.1359E-01  6.7227E-02  7.5489E-02  1.1815E-01  5.8007E-02  9.5625E-02  8.2479E-02  6.9298E-01 -2.4907E-02  6.3751E-02
             6.0027E-02
 GRADIENT:  -3.3295E+01 -5.5330E+01  1.6545E+00 -8.2947E+01  1.0525E+01  2.3280E+00 -1.0931E+01 -1.3924E+02 -8.2818E+00 -3.1754E+00
            -7.0837E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2386.60897088854        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      978
 NPARAMETR:  1.0131E+00  9.6826E-01  9.7628E-01  1.0183E+00  9.5839E-01  9.8872E-01  9.8313E-01  1.8103E+00  8.8213E-01  9.8800E-01
             9.6131E-01
 PARAMETER:  1.1302E-01  6.7744E-02  7.5996E-02  1.1815E-01  5.7499E-02  8.8652E-02  8.2987E-02  6.9348E-01 -2.5413E-02  8.7929E-02
             6.0546E-02
 GRADIENT:  -3.5143E+01 -5.5250E+01  2.3732E+00 -8.2559E+01  7.0574E+00 -4.1722E-01 -1.0229E+01 -1.3890E+02 -8.1531E+00  3.9353E-01
            -6.9238E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2394.08759747912        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1161
 NPARAMETR:  1.0142E+00  9.7095E-01  9.7625E-01  1.0198E+00  9.5801E-01  9.4311E-01  9.8367E-01  1.9673E+00  8.8228E-01  1.1843E+00
             9.6467E-01
 PARAMETER:  1.1410E-01  7.0522E-02  7.5968E-02  1.1958E-01  5.7103E-02  4.1426E-02  8.3531E-02  7.7668E-01 -2.5240E-02  2.6912E-01
             6.4033E-02
 GRADIENT:  -3.6299E+01 -5.9107E+01 -2.0494E+00 -7.3278E+01 -1.2698E+01 -1.9670E+01 -6.6829E+00 -1.0315E+02 -3.9176E+00  2.8389E+01
            -6.0767E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2399.07658505861        NO. OF FUNC. EVALS.: 122
 CUMULATIVE NO. OF FUNC. EVALS.:     1283
 NPARAMETR:  1.0115E+00  9.8646E-01  9.9453E-01  1.0199E+00  9.6437E-01  9.4769E-01  1.0107E+00  1.9937E+00  8.8497E-01  1.1115E+00
             9.8463E-01
 PARAMETER:  1.1140E-01  8.6370E-02  9.4519E-02  1.1971E-01  6.3725E-02  4.6267E-02  1.1062E-01  7.8999E-01 -2.2203E-02  2.0571E-01
             8.4507E-02
 GRADIENT:   4.7688E+02  2.0880E+01 -3.4012E-01  4.2837E+01  1.1365E+00  3.2187E+01  4.2095E+00 -6.5284E+01  5.5252E+00  2.3311E+01
            -3.7729E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2399.49817798295        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     1452
 NPARAMETR:  1.0115E+00  9.8614E-01  9.9528E-01  1.0202E+00  9.6496E-01  9.9398E-01  1.0104E+00  1.9939E+00  8.8496E-01  1.1116E+00
             9.8449E-01
 PARAMETER:  1.1146E-01  8.6045E-02  9.5271E-02  1.2001E-01  6.4330E-02  9.3965E-02  1.1038E-01  7.9009E-01 -2.2214E-02  2.0585E-01
             8.4371E-02
 GRADIENT:  -3.8538E+01 -3.8991E+01 -4.6440E+00 -6.2251E+01 -2.2560E+01  1.5003E+00 -3.5662E+00 -1.0155E+02 -2.9770E-01  1.7552E+01
            -3.7928E+01

0ITERATION NO.:   48    OBJECTIVE VALUE:  -2399.59041874737        NO. OF FUNC. EVALS.: 112
 CUMULATIVE NO. OF FUNC. EVALS.:     1564
 NPARAMETR:  1.0118E+00  9.8680E-01  9.9545E-01  1.0202E+00  9.6516E-01  9.9387E-01  1.0117E+00  1.9938E+00  8.8500E-01  1.1071E+00
             9.8546E-01
 PARAMETER:  1.1176E-01  8.6712E-02  9.5359E-02  1.2000E-01  6.4533E-02  9.3908E-02  1.1101E-01  7.9003E-01 -2.2188E-02  2.0379E-01
             8.4860E-02
 GRADIENT:  -2.2134E+04  4.9335E+04 -3.9896E+00 -8.0297E+01 -2.4715E+04  1.2582E+00 -3.4638E+00 -2.5167E+02 -2.6621E-01  1.7170E+01
            -3.7055E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1564
 NO. OF SIG. DIGITS UNREPORTABLE

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.8645E-02 -4.7779E-03 -3.8651E-02  3.3813E-02 -3.4776E-02
 SE:             2.9711E-02  2.0201E-02  3.0348E-02  2.5681E-02  2.0917E-02
 N:                     100         100         100         100         100

 P VAL.:         5.3030E-01  8.1303E-01  2.0280E-01  1.8796E-01  9.6398E-02

 ETASHRINKSD(%)  4.6435E-01  3.2324E+01  1.0000E-10  1.3964E+01  2.9926E+01
 ETASHRINKVR(%)  9.2654E-01  5.4199E+01  1.0000E-10  2.5978E+01  5.0897E+01
 EBVSHRINKSD(%)  3.1744E-01  3.3181E+01  2.8961E+01  1.5161E+01  2.3251E+01
 EBVSHRINKVR(%)  6.3386E-01  5.5352E+01  4.9535E+01  2.8024E+01  4.1096E+01
 RELATIVEINF(%)  9.9355E+01  1.4863E+01  3.3834E+01  2.9377E+01  2.4893E+01
 EPSSHRINKSD(%)  3.0331E+01
 EPSSHRINKVR(%)  5.1463E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2399.5904187473711     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1296.8641789017640     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    34.16
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    10.29
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2399.590       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  9.87E-01  9.95E-01  1.02E+00  9.65E-01  9.94E-01  1.01E+00  1.99E+00  8.85E-01  1.11E+00  9.85E-01
 


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
+        9.66E+06
 
 TH 2
+       -1.11E+07  1.27E+07
 
 TH 3
+       -1.10E+07  1.26E+07  1.25E+07
 
 TH 4
+        8.92E+06 -2.04E+07 -1.01E+07  1.65E+07
 
 TH 5
+        1.13E+07 -2.59E+07 -2.57E+07  2.09E+07  2.65E+07
 
 TH 6
+       -2.26E+02  3.02E+02  4.31E+00 -4.16E-01 -3.08E+02  1.97E+02
 
 TH 7
+       -9.73E+06  1.11E+07  1.11E+07 -8.98E+06 -1.14E+07 -3.53E-01  9.80E+06
 
 TH 8
+        2.69E+00 -8.04E+05 -1.14E+02  1.28E+06  8.22E+05 -1.65E+00 -7.07E+05  9.96E+04
 
 TH 9
+        1.23E+07 -1.41E+07 -1.40E+07  1.14E+07  1.45E+07 -5.60E-01 -1.24E+07  8.96E+05  1.58E+07
 
 TH10
+        4.83E+06 -3.13E+02 -5.49E+06 -4.46E+06  5.66E+06  1.97E-01 -4.86E+06 -1.71E+00  6.17E+06  5.64E+01
 
 TH11
+       -1.11E+07  1.27E+07  1.26E+07 -1.02E+07 -1.30E+07  5.42E+00  1.12E+07 -2.55E+01 -1.42E+07  1.90E+01  1.27E+07
 
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
 #CPUT: Total CPU Time in Seconds,       44.526
Stop Time:
Thu Sep 30 04:30:57 CDT 2021
