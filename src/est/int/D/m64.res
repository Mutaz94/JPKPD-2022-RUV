Wed Sep 29 09:26:12 CDT 2021
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
$DATA ../../../../data/int/D/dat64.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m64.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   22898.4568951153        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.0479E+02  3.8892E+02 -3.7443E+01  3.2803E+02  5.9228E+02 -3.1984E+03 -1.3153E+03 -9.8037E+01 -1.9436E+03 -1.2286E+03
            -4.5498E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1095.80683587286        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:      113
 NPARAMETR:  4.3578E+00  1.8942E+00  1.0553E+00  2.3900E+00  4.2183E-01  8.0171E+00  4.0266E+00  1.0044E+00  3.7458E+00  3.8485E+00
             8.7351E+00
 PARAMETER:  1.5720E+00  7.3879E-01  1.5383E-01  9.7130E-01 -7.6316E-01  2.1816E+00  1.4929E+00  1.0436E-01  1.4206E+00  1.4477E+00
             2.2673E+00
 GRADIENT:   4.4349E+01  7.8613E+01 -3.6247E+00  8.8250E+01 -7.6965E+01  9.0219E+01  1.4920E+01  4.3742E+00  8.4665E+00  3.0589E+01
             2.3355E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1130.36889796112        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  8.0748E+00  1.4437E+00  3.4796E+00  2.2071E+00  7.8547E-01  7.6113E+00  2.2271E+00  8.0585E-01  5.8702E+00  3.8608E+00
             8.7528E+00
 PARAMETER:  2.1887E+00  4.6723E-01  1.3469E+00  8.9166E-01 -1.4148E-01  2.1296E+00  9.0071E-01 -1.1586E-01  1.8699E+00  1.4509E+00
             2.2694E+00
 GRADIENT:   7.9891E+01  1.6053E+01 -5.7853E+00  4.4479E+01 -5.8606E+01 -3.9787E+01  2.2449E+01  1.4134E+00  4.2974E+01  4.5707E+01
             2.7926E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1131.62440639056        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:      490
 NPARAMETR:  8.0508E+00  1.4474E+00  3.6586E+00  2.2260E+00  8.0705E-01  7.6282E+00  2.2086E+00  8.0043E-01  5.9544E+00  3.8681E+00
             8.7477E+00
 PARAMETER:  2.1858E+00  4.6973E-01  1.3971E+00  9.0020E-01 -1.1437E-01  2.1318E+00  8.9234E-01 -1.2260E-01  1.8841E+00  1.4528E+00
             2.2688E+00
 GRADIENT:   7.9490E+01  1.4528E+01 -7.0722E+00  4.4230E+01 -5.6969E+01 -3.8576E+01  2.2427E+01  1.3656E+00  4.2726E+01  4.8035E+01
             2.8019E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1132.61659604696        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      670
 NPARAMETR:  8.0472E+00  1.4444E+00  3.6873E+00  2.2289E+00  8.1565E-01  7.5574E+00  2.2057E+00  7.9444E-01  5.9666E+00  3.8692E+00
             8.7456E+00
 PARAMETER:  2.1853E+00  4.6772E-01  1.4049E+00  9.0152E-01 -1.0377E-01  2.1225E+00  8.9105E-01 -1.3012E-01  1.8862E+00  1.4530E+00
             2.2685E+00
 GRADIENT:   8.0706E+01  1.3815E+01 -7.4943E+00  4.4119E+01 -5.5833E+01 -4.0874E+01  2.2465E+01  1.3405E+00  4.2939E+01  4.8977E+01
             2.7995E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1156.56308474477        NO. OF FUNC. EVALS.: 149
 CUMULATIVE NO. OF FUNC. EVALS.:      819
 NPARAMETR:  8.5155E+00  1.4684E+00  3.8285E+00  2.1499E+00  8.1267E-01  6.8033E+00  2.2745E+00  1.0000E-02  6.2230E+00  3.7051E+00
             7.4280E+00
 PARAMETER:  2.2419E+00  4.8417E-01  1.4425E+00  8.6542E-01 -1.0743E-01  2.0174E+00  9.2174E-01 -7.3836E+00  1.9283E+00  1.4097E+00
             2.1053E+00
 GRADIENT:   1.0663E+02  1.8584E+01 -8.0995E+00  4.2278E+01 -6.6549E+01 -1.2765E+02  2.0533E+01  0.0000E+00  3.6243E+01  4.2138E+01
             4.8798E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1156.73479571211        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      995
 NPARAMETR:  8.4180E+00  1.4637E+00  3.8017E+00  2.1632E+00  8.1325E-01  7.1827E+00  2.2608E+00  1.0000E-02  6.1687E+00  3.7329E+00
             7.5701E+00
 PARAMETER:  2.2304E+00  4.8099E-01  1.4355E+00  8.7160E-01 -1.0672E-01  2.0717E+00  9.1570E-01 -5.9989E+00  1.9195E+00  1.4172E+00
             2.1242E+00
 GRADIENT:   9.6120E+01  1.8489E+01 -7.7512E+00  4.2649E+01 -6.4870E+01 -1.0152E+02  2.0503E+01  0.0000E+00  3.5850E+01  4.4002E+01
             8.1840E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1272.37906802901        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1178             RESET HESSIAN, TYPE I
 NPARAMETR:  2.5187E+00  1.3042E+00  3.8148E+00  1.2229E+00  1.1187E+00  7.2645E+00  1.5496E+00  1.0000E-02  5.2714E+00  3.2570E+00
             7.3456E+00
 PARAMETER:  1.0237E+00  3.6561E-01  1.4389E+00  3.0124E-01  2.1216E-01  2.0830E+00  5.3800E-01 -4.7613E+00  1.7623E+00  1.2808E+00
             2.0941E+00
 GRADIENT:   1.3743E+02  1.2256E+01 -1.6297E+01  2.4168E+01 -1.4482E+01  3.7105E+02  1.1452E+01  0.0000E+00  7.4161E+01  6.0629E+01
             9.1574E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1348.50881618546        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:     1293
 NPARAMETR:  1.2669E+00  1.1171E+00  3.9892E+00  1.0261E+00  1.1744E+00  3.7580E+00  1.3489E+00  1.0000E-02  4.2969E+00  2.9372E+00
             7.1203E+00
 PARAMETER:  3.3660E-01  2.1072E-01  1.4836E+00  1.2578E-01  2.6074E-01  1.4239E+00  3.9929E-01 -4.7613E+00  1.5579E+00  1.1775E+00
             2.0629E+00
 GRADIENT:  -9.4305E+00 -2.0305E+00 -2.5276E+01 -3.7027E+00  5.7925E-01  2.2496E+00  3.4707E+00  0.0000E+00 -4.5889E+00 -9.2332E+00
             8.9119E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1349.69116089222        NO. OF FUNC. EVALS.: 145
 CUMULATIVE NO. OF FUNC. EVALS.:     1438
 NPARAMETR:  1.3234E+00  1.1028E+00  4.0503E+00  1.0512E+00  1.1730E+00  3.7277E+00  1.0995E+00  1.0000E-02  4.3198E+00  2.9857E+00
             7.1596E+00
 PARAMETER:  3.8019E-01  1.9781E-01  1.4988E+00  1.4994E-01  2.5954E-01  1.4158E+00  1.9486E-01 -4.7613E+00  1.5632E+00  1.1938E+00
             2.0685E+00
 GRADIENT:   3.8711E+01  7.6828E+00 -2.2715E+01  3.0392E+00  3.4453E+00  1.4101E+02  8.6546E-01  0.0000E+00  4.5552E+01  3.2963E+01
             3.8078E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1350.03969296334        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1618
 NPARAMETR:  1.3548E+00  1.0208E+00  4.0790E+00  1.1183E+00  1.1545E+00  3.7304E+00  9.1137E-01  1.0000E-02  4.3243E+00  2.9894E+00
             7.1639E+00
 PARAMETER:  4.0368E-01  1.2059E-01  1.5058E+00  2.1179E-01  2.4369E-01  1.4165E+00  7.1976E-03 -4.7613E+00  1.5643E+00  1.1951E+00
             2.0691E+00
 GRADIENT:   9.4661E-01  2.0249E-01 -2.6106E+01  8.4468E-01 -1.4587E-01 -9.1543E-01 -1.6943E-01  0.0000E+00  2.2147E+00 -2.1581E+00
            -2.2448E+00

0ITERATION NO.:   52    OBJECTIVE VALUE:  -1350.05537166973        NO. OF FUNC. EVALS.:  68
 CUMULATIVE NO. OF FUNC. EVALS.:     1686
 NPARAMETR:  1.3538E+00  1.0194E+00  4.0808E+00  1.1193E+00  1.1548E+00  3.7296E+00  9.0835E-01  1.0000E-02  4.3231E+00  2.9904E+00
             7.1670E+00
 PARAMETER:  4.0296E-01  1.1929E-01  1.5065E+00  2.1276E-01  2.4398E-01  1.4165E+00  2.8777E-03 -4.7613E+00  1.5642E+00  1.1952E+00
             2.0692E+00
 GRADIENT:   6.5447E+02  1.0494E-01  1.5397E+02  6.2132E+02  1.0827E+03  1.8065E+02 -1.7780E-01  0.0000E+00  8.6178E+01 -2.2419E+02
            -1.3365E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1686
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.8798E-03 -4.5389E-02  8.0767E-05  5.7694E-03 -1.7657E-02
 SE:             3.0074E-02  1.0297E-02  2.8881E-05  2.6715E-02  2.6700E-02
 N:                     100         100         100         100         100

 P VAL.:         8.9735E-01  1.0447E-05  5.1655E-03  8.2902E-01  5.0841E-01

 ETASHRINKSD(%)  1.0000E-10  6.5503E+01  9.9903E+01  1.0502E+01  1.0550E+01
 ETASHRINKVR(%)  1.0000E-10  8.8099E+01  1.0000E+02  1.9902E+01  1.9987E+01
 EBVSHRINKSD(%)  7.4151E-01  7.5682E+01  9.9860E+01  6.6522E+00  8.1647E+00
 EBVSHRINKVR(%)  1.4775E+00  9.4086E+01  1.0000E+02  1.2862E+01  1.5663E+01
 RELATIVEINF(%)  9.8502E+01  2.7298E+00  1.0179E-04  4.3406E+01  4.5872E+01
 EPSSHRINKSD(%)  1.1706E+01
 EPSSHRINKVR(%)  2.2041E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1350.0553716697270     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       304.03398809868372     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    69.52
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    17.16
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1350.055       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.35E+00  1.02E+00  4.08E+00  1.12E+00  1.15E+00  3.73E+00  9.07E-01  1.00E-02  4.32E+00  2.99E+00  7.17E+00
 


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
+        2.22E+04
 
 TH 2
+        5.12E+01  4.45E+05
 
 TH 3
+        7.79E-01  3.77E+00  1.82E+02
 
 TH 4
+        7.42E+01  9.29E+01 -9.93E+00  1.17E+05
 
 TH 5
+        7.12E+01  3.13E+01 -6.65E+00  4.15E+02  8.34E+04
 
 TH 6
+        1.76E+00  1.53E+00  2.06E+00 -4.63E+00  4.09E+00  2.55E+02
 
 TH 7
+       -4.49E-01 -1.32E+01 -9.32E-02 -3.50E+00 -4.14E-01  6.81E-02  6.21E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        6.14E-01 -1.47E+01  3.93E+00 -7.36E+00 -6.44E-01  1.32E+00  2.19E+00  0.00E+00  1.49E+02
 
 TH10
+       -1.44E+00 -6.58E+00  1.13E+00  1.45E+01 -5.40E+00 -3.81E+00  9.36E-01  0.00E+00  2.40E+00  5.38E+02
 
 TH11
+       -2.33E+00 -7.70E+00 -2.13E+00  1.24E+00  2.15E+00 -1.03E-01  3.05E+00  0.00E+00 -1.75E+00 -5.57E-01  5.07E+01
 
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
 
 Elapsed finaloutput time in seconds:     0.05
 #CPUT: Total CPU Time in Seconds,       86.796
Stop Time:
Wed Sep 29 09:27:41 CDT 2021
