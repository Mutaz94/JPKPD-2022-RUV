Thu Sep 30 07:23:40 CDT 2021
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
$DATA ../../../../data/spa2/TD1/dat37.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m37.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2142.25239383564        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3235E+02  1.2383E+01  1.2714E+02 -8.3548E+00 -8.4141E-01  5.7726E+01 -2.1031E+01 -4.9990E+02 -1.1071E+02 -5.8185E+00
            -2.5434E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2354.70982621563        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.8946E-01  1.0215E+00  1.0412E+00  9.5157E-01  1.0250E+00  1.1988E+00  1.0425E+00  2.1989E+00  8.8098E-01  9.4251E-01
             9.2533E-01
 PARAMETER:  8.9408E-02  1.2131E-01  1.4035E-01  5.0357E-02  1.2465E-01  2.8135E-01  1.4166E-01  8.8795E-01 -2.6724E-02  4.0795E-02
             2.2392E-02
 GRADIENT:   5.0294E+02  1.6815E+01  1.7501E+01 -5.1690E+01 -4.5798E+01  2.5006E+02 -2.4189E+01 -3.9156E+01 -1.8370E+01  1.9144E+00
            -7.9129E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2364.74262625208        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:      271
 NPARAMETR:  9.8972E-01  1.0610E+00  1.0392E+00  9.7510E-01  1.0618E+00  1.2475E+00  1.1500E+00  2.2051E+00  8.9157E-01  9.4546E-01
             9.9735E-01
 PARAMETER:  8.9663E-02  1.5924E-01  1.3842E-01  7.4788E-02  1.5995E-01  3.2115E-01  2.3978E-01  8.9079E-01 -1.4776E-02  4.3912E-02
             9.7349E-02
 GRADIENT:  -2.1098E+01 -3.1637E+01 -1.6757E+00 -4.3756E+01 -4.8086E+01  7.6455E+01 -1.9574E+01 -1.0929E+02 -4.9071E+00  2.0148E+00
             1.8106E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2369.36663339390        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      452
 NPARAMETR:  9.8860E-01  9.2562E-01  1.1195E+00  1.0678E+00  1.0109E+00  1.2435E+00  1.3261E+00  2.3247E+00  8.1958E-01  8.0994E-01
             9.9028E-01
 PARAMETER:  8.8538E-02  2.2707E-02  2.1287E-01  1.6563E-01  1.1087E-01  3.1792E-01  3.8221E-01  9.4360E-01 -9.8966E-02 -1.1079E-01
             9.0229E-02
 GRADIENT:  -2.2143E+01 -2.4958E+01 -2.5044E+00 -2.0816E+01 -4.6212E+01  7.5735E+01 -1.3533E+01 -9.6796E+01 -3.3738E+00 -3.3622E+00
            -3.7030E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2372.05827887805        NO. OF FUNC. EVALS.: 172
 CUMULATIVE NO. OF FUNC. EVALS.:      624
 NPARAMETR:  9.8807E-01  9.3003E-01  1.1231E+00  1.0665E+00  1.0140E+00  1.2097E+00  1.3266E+00  2.3373E+00  8.1956E-01  8.1336E-01
             9.9044E-01
 PARAMETER:  8.8003E-02  2.7461E-02  2.1611E-01  1.6442E-01  1.1394E-01  2.9037E-01  3.8261E-01  9.4900E-01 -9.8991E-02 -1.0659E-01
             9.0394E-02
 GRADIENT:   4.3223E+02  3.4271E+01  8.9041E+00  1.3440E+02 -1.9938E+01  2.3318E+02  3.0733E+01 -2.2879E+01  2.6245E+00 -1.5553E+00
            -2.3150E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2379.95561715214        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:      814
 NPARAMETR:  1.0067E+00  9.3065E-01  1.1225E+00  1.0670E+00  1.0143E+00  9.9819E-01  1.3254E+00  2.3434E+00  8.2024E-01  8.3477E-01
             9.9076E-01
 PARAMETER:  1.0663E-01  2.8129E-02  2.1558E-01  1.6482E-01  1.1422E-01  9.8192E-02  3.8168E-01  9.5159E-01 -9.8161E-02 -8.0596E-02
             9.0722E-02
 GRADIENT:   4.8267E+00 -2.3455E+01 -2.6619E+00 -1.6882E+01 -4.7234E+01  5.9203E+00 -1.2745E+01 -9.3757E+01 -2.7620E+00 -5.0159E-02
            -3.2481E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2380.97180092869        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      991
 NPARAMETR:  9.8814E-01  9.3069E-01  1.1226E+00  1.0670E+00  1.0144E+00  9.5779E-01  1.3257E+00  2.3821E+00  8.2023E-01  8.3290E-01
             9.9075E-01
 PARAMETER:  8.8065E-02  2.8171E-02  2.1564E-01  1.6488E-01  1.1433E-01  5.6869E-02  3.8194E-01  9.6798E-01 -9.8174E-02 -8.2841E-02
             9.0710E-02
 GRADIENT:  -3.9890E+01 -2.2873E+01 -3.7409E+00 -1.6874E+01 -4.8544E+01 -1.1259E+01 -1.2695E+01 -8.7796E+01 -2.4630E+00 -3.2579E-02
            -2.9302E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2381.37470993085        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1177
 NPARAMETR:  9.8786E-01  9.3801E-01  1.1251E+00  1.0675E+00  1.0144E+00  9.6511E-01  1.3259E+00  2.3865E+00  8.2438E-01  8.3335E-01
             9.9153E-01
 PARAMETER:  8.7785E-02  3.6004E-02  2.1785E-01  1.6530E-01  1.1431E-01  6.4489E-02  3.8207E-01  9.6982E-01 -9.3119E-02 -8.2296E-02
             9.1494E-02
 GRADIENT:  -3.9994E+01 -1.6712E+01 -2.6908E+00 -1.0629E+01 -5.2306E+01 -8.1541E+00 -1.1603E+01 -8.7063E+01 -1.8266E+00 -1.4881E-01
            -2.1551E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2382.38500205695        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:     1317
 NPARAMETR:  1.0049E+00  9.3749E-01  1.1237E+00  1.0685E+00  1.0151E+00  9.7447E-01  1.3231E+00  2.4031E+00  8.2485E-01  8.3601E-01
             9.9209E-01
 PARAMETER:  1.0486E-01  3.5446E-02  2.1663E-01  1.6623E-01  1.1496E-01  7.4141E-02  3.7994E-01  9.7674E-01 -9.2559E-02 -7.9110E-02
             9.2055E-02
 GRADIENT:   6.7084E-01 -1.7274E+01 -3.6861E+00 -9.0180E+00 -5.1388E+01 -3.5069E+00 -1.1942E+01 -8.4494E+01 -1.6586E+00  3.3043E-01
            -1.5307E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2386.59862310286        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1502             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0057E+00  9.5933E-01  1.1243E+00  1.0683E+00  1.0153E+00  9.8526E-01  1.3250E+00  2.5372E+00  8.2831E-01  8.2777E-01
             9.9235E-01
 PARAMETER:  1.0568E-01  5.8477E-02  2.1716E-01  1.6605E-01  1.1520E-01  8.5151E-02  3.8138E-01  1.0311E+00 -8.8369E-02 -8.9014E-02
             9.2324E-02
 GRADIENT:   4.8377E+02  5.9203E+01  6.3068E+00  1.5869E+02 -3.9778E+01  5.4701E+01  3.1974E+01  1.0898E+01  5.0838E+00  6.6095E-01
            -8.5931E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2386.65628714587        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1686
 NPARAMETR:  1.0044E+00  9.5947E-01  1.1244E+00  1.0682E+00  1.0153E+00  9.8232E-01  1.3252E+00  2.5396E+00  8.2795E-01  8.2604E-01
             9.9328E-01
 PARAMETER:  1.0441E-01  5.8621E-02  2.1723E-01  1.6600E-01  1.1518E-01  8.2166E-02  3.8154E-01  1.0320E+00 -8.8801E-02 -9.1112E-02
             9.3258E-02
 GRADIENT:  -5.2581E-01  7.9010E-01 -4.7071E+00  5.7972E+00 -6.6269E+01 -3.4531E-01 -9.6457E+00 -6.5041E+01 -8.1801E-01 -1.0296E+00
            -3.9733E-01

0ITERATION NO.:   51    OBJECTIVE VALUE:  -2386.65628714587        NO. OF FUNC. EVALS.:  28
 CUMULATIVE NO. OF FUNC. EVALS.:     1714
 NPARAMETR:  1.0045E+00  9.5944E-01  1.1244E+00  1.0682E+00  1.0153E+00  9.8244E-01  1.3250E+00  2.5389E+00  8.2793E-01  8.2687E-01
             9.9331E-01
 PARAMETER:  1.0441E-01  5.8621E-02  2.1723E-01  1.6600E-01  1.1518E-01  8.2166E-02  3.8154E-01  1.0320E+00 -8.8801E-02 -9.1112E-02
             9.3258E-02
 GRADIENT:  -1.2678E+04  1.3234E+04 -1.4093E+03  7.9743E+03 -1.1558E+04 -4.2879E-01  6.9281E+03  2.4562E+03  2.6472E+04 -1.0044E+00
            -2.6477E+04
 NUMSIGDIG:         2.3         2.3         3.5         2.3         2.3         1.7         2.3         2.3         2.3         0.7
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1714
 NO. OF SIG. DIGITS IN FINAL EST.:  0.7
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.8275E-03 -6.9131E-03 -2.9943E-02 -5.6764E-04 -2.5225E-02
 SE:             2.9922E-02  2.3695E-02  2.9287E-02  2.4230E-02  1.9113E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5130E-01  7.7048E-01  3.0659E-01  9.8131E-01  1.8690E-01

 ETASHRINKSD(%)  1.0000E-10  2.0619E+01  1.8849E+00  1.8827E+01  3.5969E+01
 ETASHRINKVR(%)  1.0000E-10  3.6987E+01  3.7343E+00  3.4110E+01  5.9001E+01
 EBVSHRINKSD(%)  3.3291E-01  2.3292E+01  1.8803E+01  2.1072E+01  3.6833E+01
 EBVSHRINKVR(%)  6.6470E-01  4.1159E+01  3.4070E+01  3.7704E+01  6.0099E+01
 RELATIVEINF(%)  9.9325E+01  2.2108E+01  4.5328E+01  2.6363E+01  1.8284E+01
 EPSSHRINKSD(%)  3.2325E+01
 EPSSHRINKVR(%)  5.4200E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2386.6562871458709     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1283.9300473002638     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    36.73
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    10.36
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2386.656       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  9.59E-01  1.12E+00  1.07E+00  1.02E+00  9.82E-01  1.33E+00  2.54E+00  8.28E-01  8.26E-01  9.93E-01
 


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
+        6.02E+06
 
 TH 2
+        2.94E+03  7.19E+06
 
 TH 3
+        2.28E+06  2.50E+06  3.07E+06
 
 TH 4
+       -3.56E+06 -1.39E+03 -2.70E+06  3.96E+06
 
 TH 5
+        4.74E+02 -7.29E+02 -4.09E+06  2.82E+06  4.84E+06
 
 TH 6
+       -1.55E+02  1.71E+02  6.11E+01  9.24E+01 -1.43E+02  2.05E+02
 
 TH 7
+        1.68E+02  1.36E+06  9.47E+05 -6.52E+05 -9.88E+05  3.24E+01  4.87E+05
 
 TH 8
+        9.66E+00 -1.16E+02 -1.82E+05  2.51E+05  1.94E+05  6.33E+00 -4.48E+04  1.76E+04
 
 TH 9
+       -1.88E+02  1.88E+02 -2.89E+06  1.13E+02 -6.84E+06  2.00E+02  6.90E+01  1.19E+01  9.65E+06
 
 TH10
+       -5.15E+01  3.07E+01 -3.58E+01  4.78E+01 -6.85E+06 -3.95E-01  1.32E+01  5.36E+00  7.11E+01  9.70E+06
 
 TH11
+        2.14E+03  3.09E+03 -2.41E+06 -1.26E+03  4.97E+02 -1.64E+02  1.79E+02 -8.18E+01 -1.86E+02 -3.59E+01  6.71E+06
 
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
 #CPUT: Total CPU Time in Seconds,       47.170
Stop Time:
Thu Sep 30 07:24:28 CDT 2021
