Thu Sep 30 02:38:57 CDT 2021
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
$DATA ../../../../data/spa1/D/dat13.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m13.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1638.45047642222        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.3183E+02 -1.3692E+02 -6.9082E+01 -1.9637E+02  1.3785E+02 -2.0053E+02 -2.2879E+02 -4.2228E+01 -3.8094E+02 -7.3745E+01
            -2.5085E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1885.87570590676        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:       90
 NPARAMETR:  8.3627E-01  1.3178E+00  1.2687E+00  1.0846E+00  1.1405E+00  1.1776E+00  2.8954E+00  1.4118E+00  1.3544E+00  1.1673E+00
             9.6985E-01
 PARAMETER: -7.8801E-02  3.7594E-01  3.3800E-01  1.8121E-01  2.3147E-01  2.6346E-01  1.1631E+00  4.4487E-01  4.0337E-01  2.5465E-01
             6.9386E-02
 GRADIENT:   3.2763E+01  2.3086E+02 -5.7521E+00  2.1649E+02 -2.0849E+01  2.4556E+02  3.7166E+02 -4.8344E+00  6.2140E+01  3.1565E+01
            -2.5477E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1898.83020504948        NO. OF FUNC. EVALS.: 122
 CUMULATIVE NO. OF FUNC. EVALS.:      212
 NPARAMETR:  8.9698E-01  1.4864E+00  1.6183E+00  1.0457E+00  1.2344E+00  1.1097E+00  2.7520E+00  1.8511E+00  1.3387E+00  1.1465E+00
             9.8590E-01
 PARAMETER: -8.7220E-03  4.9639E-01  5.8140E-01  1.4471E-01  3.1062E-01  2.0413E-01  1.1123E+00  7.1577E-01  3.9168E-01  2.3668E-01
             8.5799E-02
 GRADIENT:  -3.4430E+02  4.8947E+01  8.2895E+00  5.0679E+01 -5.2409E+01 -2.8374E+02 -2.8296E+00 -4.9036E+00  6.8041E+00  1.9317E+01
            -1.6108E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1951.87290280950        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:      401
 NPARAMETR:  9.0038E-01  1.4611E+00  1.6343E+00  1.0576E+00  1.2651E+00  1.4929E+00  2.7578E+00  1.9125E+00  1.4336E+00  9.9208E-01
             1.0131E+00
 PARAMETER: -4.9394E-03  4.7917E-01  5.9120E-01  1.5599E-01  3.3511E-01  5.0073E-01  1.1144E+00  7.4840E-01  4.6019E-01  9.2044E-02
             1.1297E-01
 GRADIENT:  -1.8654E+02  4.2325E+01  7.5228E-01  5.6054E+01 -9.3005E+00 -4.6046E+01  1.6366E+00 -4.1000E+00  1.1693E+01  2.3451E+00
             4.9747E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1976.13206169703        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      563
 NPARAMETR:  1.0482E+00  1.2285E+00  1.6896E+00  9.7080E-01  1.2804E+00  1.6684E+00  2.7140E+00  2.0123E+00  1.2867E+00  9.7694E-01
             1.0128E+00
 PARAMETER:  1.4707E-01  3.0581E-01  6.2450E-01  7.0366E-02  3.4714E-01  6.1186E-01  1.0984E+00  7.9928E-01  3.5212E-01  7.6667E-02
             1.1273E-01
 GRADIENT:   5.8876E+02  1.5817E+02  1.4057E+01  5.7911E+01  4.9300E+00  8.1508E+02  3.8789E+02  2.0882E+00  3.7696E+01 -3.1941E-01
             4.1580E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1977.74332365650        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      740
 NPARAMETR:  1.0767E+00  1.3453E+00  1.6891E+00  9.5698E-01  1.3311E+00  1.5867E+00  2.7826E+00  2.1108E+00  1.1011E+00  1.0458E+00
             1.0128E+00
 PARAMETER:  1.7392E-01  3.9661E-01  6.2420E-01  5.6031E-02  3.8601E-01  5.6164E-01  1.1234E+00  8.4707E-01  1.9632E-01  1.4482E-01
             1.1268E-01
 GRADIENT:  -4.3236E+00  7.2332E+00  3.4244E+00  6.4354E+00  5.7645E+00  1.1233E+01  4.2945E+00 -7.1588E-01  4.0608E-02 -1.0568E+00
             2.0539E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1977.89134546007        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:      880
 NPARAMETR:  1.0860E+00  1.3039E+00  1.6036E+00  9.4340E-01  1.3334E+00  1.6240E+00  2.7101E+00  2.1242E+00  1.0949E+00  1.0534E+00
             1.0134E+00
 PARAMETER:  1.8251E-01  3.6534E-01  5.7222E-01  4.1730E-02  3.8773E-01  5.8492E-01  1.0970E+00  8.5339E-01  1.9063E-01  1.5198E-01
             1.1327E-01
 GRADIENT:   3.7940E+00 -6.0117E+00 -4.3352E+00  3.3558E+00  1.2308E+01  2.2367E+01 -2.7001E+00  5.7950E+00  8.4336E-01  8.4605E-01
            -1.5429E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1978.16200920832        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:     1037
 NPARAMETR:  1.0882E+00  1.3303E+00  1.5959E+00  9.3066E-01  1.3167E+00  1.6176E+00  2.7187E+00  2.0375E+00  1.0826E+00  1.0443E+00
             1.0141E+00
 PARAMETER:  1.8451E-01  3.8538E-01  5.6743E-01  2.8139E-02  3.7509E-01  5.8092E-01  1.1002E+00  8.1174E-01  1.7935E-01  1.4335E-01
             1.1402E-01
 GRADIENT:   7.7646E+02  2.2615E+02  8.4844E+00  6.3353E+01  1.5678E+01  7.3472E+02  4.0575E+02  3.8735E+00  1.0028E+01  1.3720E+00
             2.2877E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1978.23611564747        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1217
 NPARAMETR:  1.0915E+00  1.3468E+00  1.5427E+00  9.3411E-01  1.3130E+00  1.6078E+00  2.7413E+00  2.0291E+00  1.0910E+00  1.0320E+00
             1.0127E+00
 PARAMETER:  1.8757E-01  3.9776E-01  5.3352E-01  3.1836E-02  3.7231E-01  5.7487E-01  1.1084E+00  8.0758E-01  1.8713E-01  1.3153E-01
             1.1266E-01
 GRADIENT:   8.2742E+00  1.2963E+00 -3.7124E-01  3.3488E+00  6.8286E+00  1.6827E+01  3.2989E+00  2.0193E+00  4.9465E-01  2.9923E-01
            -2.8546E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1978.32531511912        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1403
 NPARAMETR:  1.0895E+00  1.3408E+00  1.5165E+00  9.3276E-01  1.3045E+00  1.6102E+00  2.7355E+00  1.9971E+00  1.0854E+00  1.0247E+00
             1.0127E+00
 PARAMETER:  1.8574E-01  3.9328E-01  5.1640E-01  3.0391E-02  3.6584E-01  5.7634E-01  1.1063E+00  7.9171E-01  1.8198E-01  1.2438E-01
             1.1261E-01
 GRADIENT:   6.5723E+00  1.1304E-01  1.5402E-01  2.1901E+00  5.7868E+00  1.7472E+01  2.4524E+00  1.3450E+00  2.1064E-01  4.2192E-01
             2.1752E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1978.39214879694        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1587
 NPARAMETR:  1.0878E+00  1.3375E+00  1.4750E+00  9.3336E-01  1.2947E+00  1.6136E+00  2.7266E+00  1.9685E+00  1.0890E+00  1.0124E+00
             1.0121E+00
 PARAMETER:  1.8415E-01  3.9077E-01  4.8865E-01  3.1036E-02  3.5824E-01  5.7846E-01  1.1031E+00  7.7729E-01  1.8523E-01  1.1231E-01
             1.1205E-01
 GRADIENT:   5.0711E+00 -5.7155E-01 -1.2158E+00  4.1506E+00  7.9677E+00  1.8323E+01  1.7384E+00  1.5421E+00  4.4745E-01  4.8976E-01
             5.1847E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1978.43301350756        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1768             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0877E+00  1.3385E+00  1.4338E+00  9.3112E-01  1.2867E+00  1.6185E+00  2.7233E+00  1.9465E+00  1.0869E+00  9.9872E-01
             1.0114E+00
 PARAMETER:  1.8403E-01  3.9154E-01  4.6031E-01  2.8637E-02  3.5205E-01  5.8152E-01  1.1018E+00  7.6601E-01  1.8329E-01  9.8718E-02
             1.1130E-01
 GRADIENT:   7.7537E+02  2.3199E+02  1.1173E+00  7.3183E+01  2.7313E+01  7.3301E+02  3.9928E+02  5.4900E+00  1.1237E+01  5.6522E-01
             1.7096E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1978.49225879554        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1949
 NPARAMETR:  1.0874E+00  1.3426E+00  1.4417E+00  9.2441E-01  1.2774E+00  1.6070E+00  2.7217E+00  1.9193E+00  1.0810E+00  9.9673E-01
             1.0110E+00
 PARAMETER:  1.8376E-01  3.9464E-01  4.6584E-01  2.1402E-02  3.4482E-01  5.7435E-01  1.1013E+00  7.5194E-01  1.7786E-01  9.6729E-02
             1.1094E-01
 GRADIENT:   4.6432E+00 -8.7402E-01  2.6557E+00 -1.8069E+00 -5.8750E-01  1.6251E+01  2.1615E+00 -1.5622E-01 -2.3174E-01  6.2930E-01
             2.8287E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1978.54016975173        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2130
 NPARAMETR:  1.0889E+00  1.3454E+00  1.4175E+00  9.2374E-01  1.2720E+00  1.6120E+00  2.7223E+00  1.9065E+00  1.0825E+00  9.8864E-01
             1.0104E+00
 PARAMETER:  1.8515E-01  3.9667E-01  4.4890E-01  2.0675E-02  3.4063E-01  5.7745E-01  1.1015E+00  7.4527E-01  1.7927E-01  8.8574E-02
             1.1037E-01
 GRADIENT:   5.8811E+00 -5.7923E-01  1.7813E+00 -4.9593E-01  8.6413E-01  1.7487E+01  2.6431E+00  5.6923E-02 -4.1952E-02  5.0146E-01
             2.6019E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1978.57651350023        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2313
 NPARAMETR:  1.0889E+00  1.3469E+00  1.3945E+00  9.2005E-01  1.2647E+00  1.6147E+00  2.7155E+00  1.8898E+00  1.0801E+00  9.8069E-01
             1.0098E+00
 PARAMETER:  1.8520E-01  3.9780E-01  4.3251E-01  1.6673E-02  3.3482E-01  5.7912E-01  1.0990E+00  7.3645E-01  1.7709E-01  8.0501E-02
             1.0978E-01
 GRADIENT:   5.8803E+00 -1.0996E+00  2.2627E+00 -1.7054E+00 -8.2139E-01  1.8150E+01  2.4987E+00 -1.1746E-01 -1.8173E-01  5.5017E-01
             2.1791E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1978.60924510585        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     2501             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0886E+00  1.3493E+00  1.3631E+00  9.2071E-01  1.2609E+00  1.6165E+00  2.7091E+00  1.8811E+00  1.0839E+00  9.7145E-01
             1.0093E+00
 PARAMETER:  1.8486E-01  3.9957E-01  4.0974E-01  1.7388E-02  3.3186E-01  5.8027E-01  1.0966E+00  7.3184E-01  1.8057E-01  7.1030E-02
             1.0928E-01
 GRADIENT:   7.8134E+02  2.4015E+02  3.0692E+00  6.9028E+01  2.0275E+01  7.2848E+02  3.9405E+02  4.1526E+00  1.0601E+01  6.6302E-01
             1.4544E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1978.62656744841        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2686
 NPARAMETR:  1.0886E+00  1.3512E+00  1.3507E+00  9.1903E-01  1.2573E+00  1.6168E+00  2.7058E+00  1.8715E+00  1.0839E+00  9.6721E-01
             1.0090E+00
 PARAMETER:  1.8485E-01  4.0098E-01  4.0062E-01  1.5564E-02  3.2897E-01  5.8042E-01  1.0954E+00  7.2676E-01  1.8053E-01  6.6655E-02
             1.0896E-01
 GRADIENT:   5.4973E+00 -1.0813E+00 -4.3193E-01  2.2642E+00  4.0025E+00  1.8566E+01  2.4076E+00  5.6494E-01  2.1097E-01  1.5211E-01
             3.1865E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1978.64770593803        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     2874             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0886E+00  1.3541E+00  1.3412E+00  9.1438E-01  1.2513E+00  1.6169E+00  2.7009E+00  1.8538E+00  1.0812E+00  9.6269E-01
             1.0086E+00
 PARAMETER:  1.8488E-01  4.0311E-01  3.9354E-01  1.0486E-02  3.2418E-01  5.8052E-01  1.0936E+00  7.1723E-01  1.7805E-01  6.1974E-02
             1.0852E-01
 GRADIENT:   7.8230E+02  2.4397E+02  4.9779E+00  6.5170E+01  1.5132E+01  7.2894E+02  3.9269E+02  3.3310E+00  1.0003E+01  7.6384E-01
             1.1773E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1978.65843051788        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     3063
 NPARAMETR:  1.0886E+00  1.3563E+00  1.3314E+00  9.1306E-01  1.2494E+00  1.6170E+00  2.6975E+00  1.8485E+00  1.0817E+00  9.5979E-01
             1.0084E+00
 PARAMETER:  1.8491E-01  4.0478E-01  3.8621E-01  9.0437E-03  3.2263E-01  5.8055E-01  1.0923E+00  7.1435E-01  1.7851E-01  5.8961E-02
             1.0832E-01
 GRADIENT:   5.5003E+00 -1.3343E+00  1.2447E+00 -5.7500E-01 -1.4954E-01  1.8565E+01  2.4398E+00 -2.6250E-02 -8.9801E-02  2.1984E-01
             7.6791E-02

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1978.66976235847        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     3251             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0887E+00  1.3612E+00  1.3100E+00  9.1220E-01  1.2474E+00  1.6170E+00  2.6908E+00  1.8442E+00  1.0851E+00  9.5521E-01
             1.0081E+00
 PARAMETER:  1.8496E-01  4.0835E-01  3.7002E-01  8.0991E-03  3.2104E-01  5.8057E-01  1.0898E+00  7.1202E-01  1.8163E-01  5.4180E-02
             1.0808E-01
 GRADIENT:   7.8275E+02  2.4896E+02  2.6009E+00  6.7959E+01  1.8803E+01  7.2801E+02  3.8908E+02  3.9074E+00  1.0556E+01  5.8250E-01
             1.3083E+00

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1978.67802246690        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     3435
 NPARAMETR:  1.0887E+00  1.3632E+00  1.3029E+00  9.1045E-01  1.2455E+00  1.6170E+00  2.6876E+00  1.8373E+00  1.0851E+00  9.5334E-01
             1.0080E+00
 PARAMETER:  1.8498E-01  4.0986E-01  3.6462E-01  6.1812E-03  3.1957E-01  5.8058E-01  1.0887E+00  7.0828E-01  1.8166E-01  5.2213E-02
             1.0795E-01
 GRADIENT:   5.5000E+00 -1.1795E+00 -5.6270E-01  2.0200E+00  2.6957E+00  1.8475E+01  2.4845E+00  4.6906E-01  2.0961E-01  1.0127E-01
             2.0581E-01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1978.68854956478        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     3622             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0887E+00  1.3663E+00  1.2983E+00  9.0649E-01  1.2427E+00  1.6171E+00  2.6827E+00  1.8259E+00  1.0830E+00  9.5122E-01
             1.0078E+00
 PARAMETER:  1.8501E-01  4.1212E-01  3.6102E-01  1.8194E-03  3.1731E-01  5.8061E-01  1.0868E+00  7.0208E-01  1.7976E-01  4.9991E-02
             1.0772E-01
 GRADIENT:   7.8349E+02  2.5297E+02  4.1520E+00  6.4676E+01  1.5236E+01  7.2832E+02  3.8815E+02  3.2253E+00  1.0035E+01  5.9571E-01
             1.0944E+00

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1978.69227453912        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     3803
 NPARAMETR:  1.0888E+00  1.3697E+00  1.2919E+00  9.0383E-01  1.2409E+00  1.6171E+00  2.6779E+00  1.8187E+00  1.0828E+00  9.4957E-01
             1.0076E+00
 PARAMETER:  1.8504E-01  4.1457E-01  3.5614E-01 -1.1181E-03  3.1581E-01  5.8062E-01  1.0850E+00  6.9811E-01  1.7954E-01  4.8250E-02
             1.0759E-01
 GRADIENT:   5.5138E+00 -1.4132E+00  1.4309E+00 -1.2383E+00 -1.5159E+00  1.8493E+01  2.4465E+00 -2.2854E-01 -1.8137E-01  1.3383E-01
            -4.1275E-02

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1978.70576518556        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     3989             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0889E+00  1.3799E+00  1.2753E+00  8.9980E-01  1.2384E+00  1.6171E+00  2.6633E+00  1.8108E+00  1.0849E+00  9.4468E-01
             1.0073E+00
 PARAMETER:  1.8513E-01  4.2201E-01  3.4315E-01 -5.5878E-03  3.1380E-01  5.8065E-01  1.0796E+00  6.9378E-01  1.8152E-01  4.3092E-02
             1.0725E-01
 GRADIENT:   7.8453E+02  2.6285E+02  4.3025E+00  6.3962E+01  1.4136E+01  7.2820E+02  3.8289E+02  2.9929E+00  9.8930E+00  4.3812E-01
             9.1959E-01

0ITERATION NO.:  120    OBJECTIVE VALUE:  -1978.71032346943        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     4178
 NPARAMETR:  1.0889E+00  1.3803E+00  1.2690E+00  8.9820E-01  1.2384E+00  1.6171E+00  2.6623E+00  1.8089E+00  1.0864E+00  9.4471E-01
             1.0073E+00
 PARAMETER:  1.8514E-01  4.2229E-01  3.3820E-01 -7.3583E-03  3.1382E-01  5.8064E-01  1.0792E+00  6.9273E-01  1.8285E-01  4.3124E-02
             1.0728E-01
 GRADIENT:   5.5257E+00 -1.3291E+00  6.6172E-01 -2.1317E-01 -4.3262E-01  1.8443E+01  2.3692E+00 -1.7017E-02 -4.9806E-02  4.9017E-03
            -3.2639E-02

0ITERATION NO.:  125    OBJECTIVE VALUE:  -1978.71432810740        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     4363             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0889E+00  1.3839E+00  1.2598E+00  8.9676E-01  1.2383E+00  1.6171E+00  2.6579E+00  1.8064E+00  1.0878E+00  9.4441E-01
             1.0073E+00
 PARAMETER:  1.8516E-01  4.2491E-01  3.3097E-01 -8.9687E-03  3.1371E-01  5.8064E-01  1.0776E+00  6.9135E-01  1.8414E-01  4.2805E-02
             1.0727E-01
 GRADIENT:   7.8443E+02  2.6586E+02  2.8817E+00  6.4575E+01  1.6381E+01  7.2770E+02  3.8234E+02  3.3241E+00  1.0316E+01  5.1677E-01
             1.1473E+00

0ITERATION NO.:  130    OBJECTIVE VALUE:  -1978.71623369974        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     4547
 NPARAMETR:  1.0889E+00  1.3865E+00  1.2547E+00  8.9550E-01  1.2379E+00  1.6171E+00  2.6543E+00  1.8044E+00  1.0887E+00  9.4378E-01
             1.0072E+00
 PARAMETER:  1.8518E-01  4.2681E-01  3.2687E-01 -1.0377E-02  3.1343E-01  5.8064E-01  1.0762E+00  6.9022E-01  1.8501E-01  4.2135E-02
             1.0722E-01
 GRADIENT:   5.5230E+00 -1.2289E+00 -3.9667E-01  1.2524E+00  1.2712E+00  1.8407E+01  2.4079E+00  2.6963E-01  1.3055E-01  3.7824E-02
             8.8760E-02

0ITERATION NO.:  135    OBJECTIVE VALUE:  -1978.71948306981        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     4732             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0889E+00  1.3896E+00  1.2530E+00  8.9213E-01  1.2369E+00  1.6171E+00  2.6498E+00  1.7982E+00  1.0875E+00  9.4316E-01
             1.0071E+00
 PARAMETER:  1.8520E-01  4.2903E-01  3.2555E-01 -1.4144E-02  3.1259E-01  5.8064E-01  1.0745E+00  6.8676E-01  1.8388E-01  4.1481E-02
             1.0712E-01
 GRADIENT:   7.8489E+02  2.7025E+02  3.7199E+00  6.2555E+01  1.4584E+01  7.2792E+02  3.8116E+02  2.9704E+00  1.0050E+01  5.0801E-01
             1.0294E+00

0ITERATION NO.:  140    OBJECTIVE VALUE:  -1978.72161545256        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     4919
 NPARAMETR:  1.0890E+00  1.3913E+00  1.2495E+00  8.9165E-01  1.2369E+00  1.6171E+00  2.6476E+00  1.7980E+00  1.0886E+00  9.4290E-01
             1.0071E+00
 PARAMETER:  1.8522E-01  4.3021E-01  3.2273E-01 -1.4683E-02  3.1257E-01  5.8064E-01  1.0737E+00  6.8665E-01  1.8485E-01  4.1202E-02
             1.0711E-01
 GRADIENT:   5.5292E+00 -1.3296E+00  3.7526E-01 -6.0449E-02 -2.7143E-01  1.8416E+01  2.3821E+00  1.4885E-02 -1.7106E-02  4.2464E-02
             7.0935E-04

0ITERATION NO.:  141    OBJECTIVE VALUE:  -1978.72161545256        NO. OF FUNC. EVALS.:  24
 CUMULATIVE NO. OF FUNC. EVALS.:     4943
 NPARAMETR:  1.0890E+00  1.3913E+00  1.2495E+00  8.9165E-01  1.2369E+00  1.6171E+00  2.6476E+00  1.7980E+00  1.0886E+00  9.4290E-01
             1.0071E+00
 PARAMETER:  1.8522E-01  4.3021E-01  3.2273E-01 -1.4683E-02  3.1257E-01  5.8064E-01  1.0737E+00  6.8665E-01  1.8485E-01  4.1202E-02
             1.0711E-01
 GRADIENT:  -9.0897E-03 -3.1720E-01  2.3288E-01  1.8240E-02 -2.0137E-01  1.3982E-03  2.1590E-01 -1.6649E-02 -3.8280E-02  4.1707E-02
             1.2027E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     4943
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.2425E-04  8.2822E-03 -5.4209E-02 -1.5232E-02 -4.0856E-02
 SE:             2.9954E-02  2.6659E-02  1.6414E-02  1.7832E-02  1.8985E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8337E-01  7.5605E-01  9.5792E-04  3.9298E-01  3.1400E-02

 ETASHRINKSD(%)  1.0000E-10  1.0690E+01  4.5011E+01  4.0261E+01  3.6396E+01
 ETASHRINKVR(%)  1.0000E-10  2.0237E+01  6.9762E+01  6.4313E+01  5.9546E+01
 EBVSHRINKSD(%)  1.3657E-01  8.4915E+00  5.1629E+01  4.4608E+01  3.2232E+01
 EBVSHRINKVR(%)  2.7295E-01  1.6262E+01  7.6602E+01  6.9317E+01  5.4074E+01
 RELATIVEINF(%)  9.9642E+01  2.3199E+01  5.4459E+00  5.4749E+00  1.6971E+01
 EPSSHRINKSD(%)  3.4859E+01
 EPSSHRINKVR(%)  5.7566E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1978.7216154525604     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1059.7830822478877     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    97.39
 Elapsed covariance  time in seconds:     8.39
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1978.722       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.09E+00  1.39E+00  1.25E+00  8.92E-01  1.24E+00  1.62E+00  2.65E+00  1.80E+00  1.09E+00  9.43E-01  1.01E+00
 


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
 
         5.32E-02  3.00E-01  4.59E-01  1.77E-01  1.16E-01  1.14E-01  3.92E-01  4.22E-01  2.23E-01  1.95E-01  4.98E-02
 


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
+        2.83E-03
 
 TH 2
+        5.18E-03  9.01E-02
 
 TH 3
+       -1.13E-04 -1.00E-01  2.11E-01
 
 TH 4
+       -1.75E-04 -4.79E-02  6.72E-02  3.13E-02
 
 TH 5
+        9.30E-05  1.41E-03  2.63E-02  7.38E-04  1.34E-02
 
 TH 6
+        1.33E-03  3.73E-03 -4.52E-03 -1.19E-03  1.20E-03  1.30E-02
 
 TH 7
+       -1.90E-03 -9.62E-02  1.29E-01  5.74E-02  5.21E-03  6.50E-03  1.54E-01
 
 TH 8
+        1.80E-03  1.18E-02  3.34E-02 -8.04E-03  1.40E-02 -9.17E-03 -7.30E-03  1.78E-01
 
 TH 9
+        1.30E-03 -5.79E-03  1.89E-02  9.01E-03 -2.47E-03  3.95E-03  1.55E-03 -3.10E-02  4.97E-02
 
 TH10
+       -1.43E-03 -6.20E-03  3.70E-02  4.48E-03  1.65E-02  1.99E-03  1.85E-02 -2.16E-03 -3.36E-03  3.79E-02
 
 TH11
+        3.74E-04 -3.42E-04  3.77E-03  1.07E-03  3.70E-04  2.02E-04  2.83E-03 -8.24E-04  2.18E-03 -1.63E-04  2.48E-03
 
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
+        5.32E-02
 
 TH 2
+        3.24E-01  3.00E-01
 
 TH 3
+       -4.62E-03 -7.27E-01  4.59E-01
 
 TH 4
+       -1.86E-02 -9.02E-01  8.26E-01  1.77E-01
 
 TH 5
+        1.51E-02  4.05E-02  4.94E-01  3.60E-02  1.16E-01
 
 TH 6
+        2.19E-01  1.09E-01 -8.63E-02 -5.90E-02  9.05E-02  1.14E-01
 
 TH 7
+       -9.12E-02 -8.17E-01  7.16E-01  8.26E-01  1.15E-01  1.45E-01  3.92E-01
 
 TH 8
+        8.01E-02  9.34E-02  1.72E-01 -1.08E-01  2.87E-01 -1.91E-01 -4.41E-02  4.22E-01
 
 TH 9
+        1.10E-01 -8.66E-02  1.85E-01  2.28E-01 -9.54E-02  1.55E-01  1.77E-02 -3.30E-01  2.23E-01
 
 TH10
+       -1.38E-01 -1.06E-01  4.13E-01  1.30E-01  7.30E-01  8.94E-02  2.41E-01 -2.63E-02 -7.75E-02  1.95E-01
 
 TH11
+        1.41E-01 -2.29E-02  1.65E-01  1.22E-01  6.42E-02  3.56E-02  1.45E-01 -3.93E-02  1.96E-01 -1.68E-02  4.98E-02
 
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
+        8.65E+02
 
 TH 2
+       -2.72E+02  1.68E+02
 
 TH 3
+       -4.78E+01  2.92E+01  7.35E+01
 
 TH 4
+       -3.24E+02  1.70E+02 -9.98E+01  5.68E+02
 
 TH 5
+        1.10E+02 -8.04E+01 -1.28E+02  1.22E+02  4.43E+02
 
 TH 6
+       -8.31E+01  1.06E+01  2.47E+01  2.71E+01 -6.45E+01  1.22E+02
 
 TH 7
+       -3.65E+00  1.77E+01 -2.96E+00 -3.11E+01  1.13E+01 -2.88E+01  3.39E+01
 
 TH 8
+       -2.84E+00 -2.51E+00 -1.22E+01  2.41E+01  1.05E+00  3.35E+00 -1.41E+00  1.06E+01
 
 TH 9
+        3.48E+01 -2.07E+01 -2.03E+01 -1.72E+01  3.57E+01 -2.12E+01  1.04E+01  6.82E+00  3.59E+01
 
 TH10
+        3.44E+01 -7.82E+00 -3.71E+00  6.58E+00 -8.99E+01  5.33E+00 -9.86E+00  9.81E+00  3.83E+00  7.38E+01
 
 TH11
+        1.02E+01 -4.60E+01 -2.33E+01  1.66E+01  3.73E+00  1.73E+01 -2.73E+01  7.56E+00 -1.43E+01  2.55E+01  4.70E+02
 
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
 #CPUT: Total CPU Time in Seconds,      105.820
Stop Time:
Thu Sep 30 02:40:44 CDT 2021
