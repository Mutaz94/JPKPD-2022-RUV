Thu Sep 30 06:21:27 CDT 2021
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
$DATA ../../../../data/spa2/A3/dat32.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m32.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -21.3052911132931        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.8759E+02  1.8087E+02  2.6057E+02 -1.7066E+01  2.2149E+02  8.5668E+01 -1.6954E+02 -2.7321E+02 -1.4411E+02 -1.6447E+02
            -4.1539E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1653.14079626472        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1661E+00  9.5151E-01  8.4302E-01  1.1946E+00  8.2200E-01  7.5502E-01  1.0204E+00  1.0730E+00  9.8019E-01  9.6865E-01
             5.2633E+00
 PARAMETER:  2.5366E-01  5.0293E-02 -7.0763E-02  2.7783E-01 -9.6015E-02 -1.8101E-01  1.2023E-01  1.7047E-01  7.9994E-02  6.8147E-02
             1.7608E+00
 GRADIENT:   2.6653E+02  2.2084E+01 -6.2092E+00  3.6939E+01 -2.0396E+01 -1.4985E+01  1.3722E+01  1.0802E+01  2.2983E+01  2.9578E+01
             4.0150E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1692.36522201105        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.1193E+00  5.7179E-01  3.3937E-01  1.3391E+00  3.9139E-01  7.7496E-01  1.6124E+00  6.7323E-01  1.2464E+00  1.7334E-01
             4.7355E+00
 PARAMETER:  2.1274E-01 -4.5899E-01 -9.8066E-01  3.9200E-01 -8.3805E-01 -1.5494E-01  5.7772E-01 -2.9567E-01  3.2027E-01 -1.6525E+00
             1.6551E+00
 GRADIENT:   1.5530E+02  7.1868E+01 -3.6717E+01  1.3602E+02  2.3686E+01 -9.0853E+00  2.6607E+01  1.0333E+01  3.7255E+01  1.8948E+00
             3.7803E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1794.92143606148        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      293
 NPARAMETR:  9.8571E-01  4.4276E-01  3.4213E-01  1.2632E+00  3.6602E-01  7.5987E-01  1.3695E+00  5.1699E-01  1.2685E+00  1.8392E-01
             3.4324E+00
 PARAMETER:  8.5610E-02 -7.1472E-01 -9.7257E-01  3.3361E-01 -9.0507E-01 -1.7461E-01  4.1448E-01 -5.5974E-01  3.3787E-01 -1.5933E+00
             1.3333E+00
 GRADIENT:  -2.0558E+02  1.3426E+01  5.9653E+00  7.8239E+01  1.8056E+01 -3.2120E+01 -2.3007E+01 -1.7918E-01  3.8057E+01 -1.6514E+00
             1.5294E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1838.61045691028        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      469
 NPARAMETR:  1.0226E+00  2.9604E-01  2.0601E-01  1.1436E+00  2.4805E-01  8.0834E-01  1.8280E+00  8.9197E-01  1.2190E+00  8.0336E-02
             2.7042E+00
 PARAMETER:  1.2235E-01 -1.1172E+00 -1.4798E+00  2.3422E-01 -1.2941E+00 -1.1277E-01  7.0320E-01 -1.4324E-02  2.9805E-01 -2.4215E+00
             1.0948E+00
 GRADIENT:  -4.0865E+01  1.6334E+01 -1.9960E+01  6.4977E+01 -2.5833E+00 -1.0314E+01  6.7989E+00  2.7790E+00 -1.9005E+01 -4.0312E-02
            -1.1497E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1841.71829739069        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      644
 NPARAMETR:  1.0340E+00  2.7726E-01  1.9783E-01  1.0760E+00  2.3849E-01  8.2929E-01  1.7836E+00  8.8726E-01  1.3099E+00  4.5308E-02
             2.7178E+00
 PARAMETER:  1.3348E-01 -1.1828E+00 -1.5203E+00  1.7324E-01 -1.3334E+00 -8.7183E-02  6.7861E-01 -1.9612E-02  3.6995E-01 -2.9943E+00
             1.0998E+00
 GRADIENT:  -1.3098E+00  2.0640E+00  1.0258E+00 -2.7964E-01 -5.1355E+00  4.0044E-01 -1.8821E-01 -1.8128E-01 -1.1687E+00 -7.6936E-02
            -1.3152E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1841.83054208586        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      806
 NPARAMETR:  1.0351E+00  2.7871E-01  1.9938E-01  1.0698E+00  2.3811E-01  8.3083E-01  1.8047E+00  8.9006E-01  1.3166E+00  1.5939E-01
             2.7205E+00
 PARAMETER:  1.3448E-01 -1.1776E+00 -1.5125E+00  1.6748E-01 -1.3350E+00 -8.5333E-02  6.9037E-01 -1.6469E-02  3.7503E-01 -1.7364E+00
             1.1008E+00
 GRADIENT:   5.9795E+01  3.0061E+01  4.9771E+01  9.8661E+00  1.6837E+02  3.6838E+00  9.0231E+00  2.9420E+00  9.3375E+00  3.6452E-01
             1.5577E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1842.42367782020        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      924
 NPARAMETR:  1.0335E+00  2.6776E-01  1.9249E-01  1.0651E+00  2.3199E-01  8.2974E-01  1.7031E+00  7.6851E-01  1.3126E+00  3.0545E-01
             2.7103E+00
 PARAMETER:  1.3300E-01 -1.2177E+00 -1.5477E+00  1.6311E-01 -1.3610E+00 -8.6647E-02  6.3247E-01 -1.6331E-01  3.7198E-01 -1.0860E+00
             1.0970E+00
 GRADIENT:  -3.5437E+00 -9.1672E-01  5.8882E+00 -3.0784E+00 -2.2335E+00  3.0982E-01 -1.1235E+00 -7.8131E-01 -1.7735E+00 -1.1663E-01
            -4.0094E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1842.47427895624        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:     1041
 NPARAMETR:  1.0325E+00  2.6678E-01  1.9017E-01  1.0649E+00  2.3008E-01  8.2834E-01  1.7086E+00  7.8191E-01  1.3199E+00  3.0066E-01
             2.7137E+00
 PARAMETER:  1.3196E-01 -1.2213E+00 -1.5598E+00  1.6290E-01 -1.3693E+00 -8.8331E-02  6.3567E-01 -1.4602E-01  3.7754E-01 -1.1018E+00
             1.0983E+00
 GRADIENT:   4.9714E+01  2.2791E+01  4.6376E+01  1.7946E+01  1.8543E+02  2.1104E+00  2.9680E+00 -1.3250E-01  5.5442E+00  1.0633E+00
             1.2227E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1842.50891541776        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1223
 NPARAMETR:  1.0350E+00  2.6724E-01  1.8936E-01  1.0649E+00  2.3030E-01  8.2907E-01  1.7121E+00  7.8692E-01  1.3332E+00  2.9998E-01
             2.7135E+00
 PARAMETER:  1.3438E-01 -1.2196E+00 -1.5641E+00  1.6289E-01 -1.3684E+00 -8.7451E-02  6.3770E-01 -1.3962E-01  3.8757E-01 -1.1040E+00
             1.0982E+00
 GRADIENT:   8.3439E-01  8.2324E-02  3.3558E-01  4.2572E-01  2.4999E+00  6.4112E-02  2.1377E-01  8.9678E-02  7.3778E-02  7.9207E-02
             4.1738E-01

0ITERATION NO.:   49    OBJECTIVE VALUE:  -1842.51206947279        NO. OF FUNC. EVALS.: 128
 CUMULATIVE NO. OF FUNC. EVALS.:     1351
 NPARAMETR:  1.0347E+00  2.6714E-01  1.8918E-01  1.0643E+00  2.3014E-01  8.2893E-01  1.7102E+00  7.8654E-01  1.3334E+00  2.9929E-01
             2.7127E+00
 PARAMETER:  1.3413E-01 -1.2200E+00 -1.5651E+00  1.6233E-01 -1.3691E+00 -8.7622E-02  6.3663E-01 -1.4011E-01  3.8770E-01 -1.1063E+00
             1.0980E+00
 GRADIENT:   1.3402E-01  1.8850E-01  6.8236E-01  1.8120E-01  1.6149E+00 -1.1477E-03  8.8061E-03 -9.8914E-03 -6.2801E-02 -3.2667E-03
             5.3740E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1351
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.4118E-05  2.3464E-02 -1.2556E-02 -1.0695E-02  1.7494E-02
 SE:             2.8988E-02  2.2840E-02  1.7027E-02  2.7530E-02  1.2044E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9851E-01  3.0429E-01  4.6088E-01  6.9765E-01  1.4636E-01

 ETASHRINKSD(%)  2.8881E+00  2.3482E+01  4.2956E+01  7.7725E+00  5.9651E+01
 ETASHRINKVR(%)  5.6928E+00  4.1449E+01  6.7460E+01  1.4941E+01  8.3719E+01
 EBVSHRINKSD(%)  3.0324E+00  2.2382E+01  4.1790E+01  6.9055E+00  6.0898E+01
 EBVSHRINKVR(%)  5.9728E+00  3.9755E+01  6.6115E+01  1.3334E+01  8.4710E+01
 RELATIVEINF(%)  9.3891E+01  1.8104E+01  3.1742E+00  4.7279E+01  9.4673E-01
 EPSSHRINKSD(%)  2.6371E+01
 EPSSHRINKVR(%)  4.5787E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1842.5120694727875     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -739.78582962718042     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.86
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.77
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1842.512       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  2.67E-01  1.89E-01  1.06E+00  2.30E-01  8.29E-01  1.71E+00  7.87E-01  1.33E+00  2.99E-01  2.71E+00
 


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
+        1.42E+03
 
 TH 2
+        3.27E+00  3.29E+03
 
 TH 3
+       -5.44E+01  2.73E+03  1.61E+04
 
 TH 4
+       -2.17E+01  7.11E+01 -3.46E+02  4.77E+02
 
 TH 5
+        1.25E+02 -7.10E+03 -2.09E+04 -6.12E+02  3.77E+04
 
 TH 6
+        2.53E+00 -6.00E+00  3.83E+01 -6.16E+00 -1.89E+01  2.55E+02
 
 TH 7
+        2.85E-01  8.10E+01 -2.92E+01 -9.29E-01 -4.48E+01  3.26E-02  2.67E+01
 
 TH 8
+       -2.92E+00  7.79E+00 -9.94E+01 -4.05E+00  1.67E+02  2.47E+00  4.49E+00  3.82E+01
 
 TH 9
+        1.14E+01 -2.69E+01  1.56E+02 -1.26E+01  1.80E+02  2.39E+00  1.39E+00 -2.52E+00  7.10E+01
 
 TH10
+       -4.41E+00  1.60E+01 -1.79E+02 -6.43E+00 -7.01E+05  6.61E-01  1.76E+01  3.20E+01  4.70E+00  6.97E+01
 
 TH11
+       -2.36E+01 -1.74E+01 -1.27E+02 -1.73E+00  1.35E+02  2.84E+00  3.97E+00  1.26E+01  9.31E+00  9.74E+00  7.43E+01
 
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
 #CPUT: Total CPU Time in Seconds,       35.705
Stop Time:
Thu Sep 30 06:22:04 CDT 2021
