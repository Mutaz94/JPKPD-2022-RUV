Sat Sep 25 00:02:12 CDT 2021
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
$DATA ../../../../data/int/SL1/dat17.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m17.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2811.76508220913        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -3.6210E+01 -7.8730E+01  1.9765E+02 -2.7522E+01  5.4915E+00 -5.3820E+00 -8.0849E+01 -1.5286E+02 -5.0139E+01 -5.8500E+01
            -1.8810E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3277.65106526705        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0961E+00  1.4416E+00  9.3575E-01  8.0499E-01  1.3489E+00  1.0123E+00  1.2106E+00  1.2762E+00  9.7071E-01  1.5357E+00
             1.7888E+00
 PARAMETER:  1.9172E-01  4.6574E-01  3.3589E-02 -1.1693E-01  3.9928E-01  1.1226E-01  2.9111E-01  3.4388E-01  7.0275E-02  5.2899E-01
             6.8157E-01
 GRADIENT:   1.3203E+02 -3.6890E+00  1.9654E+01  1.8075E+01  1.4229E+01 -3.7339E+00  3.3981E+01 -1.1421E+01  3.1314E+00  1.8398E+01
             9.7026E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3284.45239735698        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.1013E+00  1.8019E+00  1.0340E+00  6.6695E-01  1.4978E+00  9.8995E-01  1.0241E+00  2.6775E+00  8.6902E-01  1.6193E+00
             1.7160E+00
 PARAMETER:  1.9646E-01  6.8886E-01  1.3346E-01 -3.0504E-01  5.0402E-01  8.9897E-02  1.2379E-01  1.0849E+00 -4.0393E-02  5.8200E-01
             6.4001E-01
 GRADIENT:   1.5137E+02  1.5649E+02  1.7267E+01  5.8625E+01 -3.9851E+01 -1.4535E+01  1.0463E+01  6.2109E+00 -9.0881E+00 -3.6785E+00
             4.9148E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3292.80484028587        NO. OF FUNC. EVALS.:  83
 CUMULATIVE NO. OF FUNC. EVALS.:      239
 NPARAMETR:  1.0728E+00  1.6761E+00  9.8697E-01  7.0268E-01  1.4786E+00  9.9242E-01  9.6495E-01  2.2718E+00  9.9897E-01  1.5666E+00
             1.7058E+00
 PARAMETER:  1.7028E-01  6.1648E-01  8.6885E-02 -2.5285E-01  4.9109E-01  9.2392E-02  6.4324E-02  9.2055E-01  9.8970E-02  5.4889E-01
             6.3405E-01
 GRADIENT:   8.9162E+01  7.6730E+01  1.4053E+01  4.2074E+01 -1.4979E+01 -8.7680E+00  4.2099E+00  2.6109E+00 -3.6285E+00 -3.4642E+00
             3.7013E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3296.05325734499        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      377
 NPARAMETR:  1.0727E+00  1.6764E+00  8.7059E-01  6.4852E-01  1.5308E+00  1.0111E+00  9.3067E-01  2.2721E+00  9.9897E-01  1.6076E+00
             1.7062E+00
 PARAMETER:  1.7021E-01  6.1666E-01 -3.8581E-02 -3.3307E-01  5.2580E-01  1.1106E-01  2.8148E-02  9.2071E-01  9.8970E-02  5.7474E-01
             6.3425E-01
 GRADIENT:   8.7532E+01 -6.3957E+00  4.6938E+00 -1.6478E+00  1.1726E+00 -1.0138E+00 -8.5063E-01  5.3092E+00 -3.7215E+00 -2.7103E+00
             4.1450E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3296.14820005330        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      492
 NPARAMETR:  1.0727E+00  1.6764E+00  8.6545E-01  6.4905E-01  1.5379E+00  1.0160E+00  9.3446E-01  2.2721E+00  9.9897E-01  1.6305E+00
             1.7062E+00
 PARAMETER:  1.7021E-01  6.1666E-01 -4.4503E-02 -3.3224E-01  5.3042E-01  1.1588E-01  3.2213E-02  9.2070E-01  9.8971E-02  5.8886E-01
             6.3425E-01
 GRADIENT:   6.0201E+01 -4.1107E+01  2.5457E+00 -3.2263E+00 -2.3595E+00 -9.4013E-01 -2.2507E-01  5.0258E+00 -3.4731E+00 -7.8759E-01
             4.1040E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3298.02365850822        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      668            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0427E+00  1.7125E+00  8.6281E-01  6.5025E-01  1.5430E+00  1.0185E+00  9.3502E-01  2.1880E+00  1.0356E+00  1.6332E+00
             1.6752E+00
 PARAMETER:  1.4183E-01  6.3796E-01 -4.7562E-02 -3.3040E-01  5.3373E-01  1.1834E-01  3.2814E-02  8.8300E-01  1.3500E-01  5.9055E-01
             6.1592E-01
 GRADIENT:   2.3325E+01  3.7962E+01  6.0826E+00  2.1086E+01  4.7793E+00  3.2797E+00  3.4217E+00  8.9433E-01 -1.8445E+00  1.3650E+00
             1.9453E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3298.06919622480        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      739
 NPARAMETR:  1.0331E+00  1.7023E+00  8.5718E-01  6.4337E-01  1.5378E+00  1.0140E+00  9.2527E-01  2.1614E+00  1.0555E+00  1.6287E+00
             1.6705E+00
 PARAMETER:  1.3261E-01  6.3199E-01 -5.4106E-02 -3.4103E-01  5.3037E-01  1.1392E-01  2.2334E-02  8.7073E-01  1.5402E-01  5.8779E-01
             6.1311E-01
 GRADIENT:   2.5754E+00  1.8113E+01  6.4869E+00  8.0517E+00  2.1587E+00  1.3867E+00  1.9501E+00  3.9718E-01 -9.5301E-01  8.2564E-01
            -3.9451E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3298.07104768721        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      810
 NPARAMETR:  1.0327E+00  1.6987E+00  8.4915E-01  6.4243E-01  1.5342E+00  1.0129E+00  9.2212E-01  2.1411E+00  1.0610E+00  1.6263E+00
             1.6713E+00
 PARAMETER:  1.3215E-01  6.2985E-01 -6.3521E-02 -3.4250E-01  5.2798E-01  1.1285E-01  1.8922E-02  8.6131E-01  1.5920E-01  5.8631E-01
             6.1362E-01
 GRADIENT:   1.4668E+00  1.2602E+01  5.5455E+00  5.5048E+00  1.6326E+00  9.5158E-01  1.3349E+00  4.0427E-01 -6.3883E-01  5.8481E-01
            -2.9489E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3298.07270376173        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      880
 NPARAMETR:  1.0324E+00  1.6959E+00  8.3838E-01  6.4157E-01  1.5299E+00  1.0121E+00  9.2002E-01  2.1142E+00  1.0646E+00  1.6242E+00
             1.6722E+00
 PARAMETER:  1.3187E-01  6.2820E-01 -7.6289E-02 -3.4384E-01  5.2522E-01  1.1204E-01  1.6637E-02  8.4866E-01  1.6257E-01  5.8501E-01
             6.1416E-01
 GRADIENT:   7.7290E-01  8.2085E+00  4.2942E+00  3.6480E+00  1.1135E+00  6.1927E-01  8.6990E-01  3.3472E-01 -4.0684E-01  3.8169E-01
            -2.0705E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3298.07319572658        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      952
 NPARAMETR:  1.0322E+00  1.6946E+00  8.2987E-01  6.4078E-01  1.5271E+00  1.0116E+00  9.1887E-01  2.0938E+00  1.0664E+00  1.6231E+00
             1.6729E+00
 PARAMETER:  1.3174E-01  6.2744E-01 -8.6484E-02 -3.4506E-01  5.2337E-01  1.1158E-01  1.5394E-02  8.3900E-01  1.6427E-01  5.8431E-01
             6.1455E-01
 GRADIENT:   4.3725E-01  5.5848E+00  3.3264E+00  2.5882E+00  7.6574E-01  4.2607E-01  6.0079E-01  2.5605E-01 -2.7503E-01  2.5575E-01
            -1.4980E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -3298.07440802881        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1022
 NPARAMETR:  1.0321E+00  1.6939E+00  8.1937E-01  6.3951E-01  1.5243E+00  1.0112E+00  9.1772E-01  2.0707E+00  1.0680E+00  1.6222E+00
             1.6736E+00
 PARAMETER:  1.3163E-01  6.2705E-01 -9.9221E-02 -3.4706E-01  5.2152E-01  1.1111E-01  1.4137E-02  8.2787E-01  1.6576E-01  5.8381E-01
             6.1501E-01
 GRADIENT:   1.5251E-01  2.8943E+00  2.1418E+00  1.5239E+00  3.8382E-01  2.2949E-01  3.2976E-01  1.4870E-01 -1.4596E-01  1.2451E-01
            -8.7053E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -3298.08575869075        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:     1093
 NPARAMETR:  1.0320E+00  1.6952E+00  8.0124E-01  6.3616E-01  1.5214E+00  1.0105E+00  9.1583E-01  2.0366E+00  1.0702E+00  1.6222E+00
             1.6748E+00
 PARAMETER:  1.3153E-01  6.2782E-01 -1.2159E-01 -3.5231E-01  5.1965E-01  1.1046E-01  1.2072E-02  8.1128E-01  1.6785E-01  5.8381E-01
             6.1572E-01
 GRADIENT:  -1.4488E-01 -8.8568E-01  1.1795E-01  2.9484E-02 -1.8884E-01 -4.3945E-02 -4.2690E-02 -3.2811E-02  2.0999E-02 -5.4344E-02
             7.7953E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -3299.20893644206        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     1261
 NPARAMETR:  1.0404E+00  1.7581E+00  7.3557E-01  5.9619E-01  1.5611E+00  1.0102E+00  8.9240E-01  2.0085E+00  1.1080E+00  1.6491E+00
             1.6811E+00
 PARAMETER:  1.3958E-01  6.6421E-01 -2.0712E-01 -4.1719E-01  5.4538E-01  1.1012E-01 -1.3841E-02  7.9738E-01  2.0260E-01  6.0024E-01
             6.1944E-01
 GRADIENT:  -4.3539E+00 -3.6391E+01 -4.6761E+00 -4.0161E+00 -6.8416E+00 -2.0263E+00 -1.4870E+00 -7.9646E-01 -3.0614E-02 -1.2039E+00
             4.4274E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -3299.21294976746        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1445
 NPARAMETR:  1.0404E+00  1.7579E+00  7.3572E-01  5.9617E-01  1.5612E+00  1.0102E+00  8.9830E-01  2.0086E+00  1.1080E+00  1.6491E+00
             1.6812E+00
 PARAMETER:  1.3957E-01  6.6415E-01 -2.0691E-01 -4.1723E-01  5.4542E-01  1.1015E-01 -7.2510E-03  7.9745E-01  2.0258E-01  6.0022E-01
             6.1950E-01
 GRADIENT:  -4.3665E+00 -3.5957E+01 -4.6279E+00 -4.0583E+00 -6.8054E+00 -2.0263E+00 -8.1337E-02 -7.8593E-01  2.9112E-01 -1.1422E+00
             4.5911E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -3299.26165537581        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:     1562
 NPARAMETR:  1.0342E+00  1.7579E+00  7.5505E-01  5.9725E-01  1.5640E+00  1.0106E+00  8.9758E-01  2.0486E+00  1.1062E+00  1.6444E+00
             1.6812E+00
 PARAMETER:  1.3362E-01  6.6413E-01 -1.8097E-01 -4.1542E-01  5.4727E-01  1.1051E-01 -8.0486E-03  8.1715E-01  2.0089E-01  5.9737E-01
             6.1952E-01
 GRADIENT:   4.4626E+00  7.3316E+00 -1.6053E+00 -7.2504E-02 -1.8204E-01 -1.7608E-02  1.7242E-01 -4.0529E-01  9.9747E-02  2.5029E-01
             6.4565E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -3299.27254849811        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:     1656
 NPARAMETR:  1.0320E+00  1.7579E+00  7.7273E-01  5.9850E-01  1.5683E+00  1.0108E+00  8.9598E-01  2.0931E+00  1.1072E+00  1.6428E+00
             1.6812E+00
 PARAMETER:  1.3149E-01  6.6413E-01 -1.5783E-01 -4.1333E-01  5.4999E-01  1.1075E-01 -9.8351E-03  8.3866E-01  2.0181E-01  5.9642E-01
             6.1952E-01
 GRADIENT:  -2.1720E+01 -3.1299E+01  3.7894E-01 -5.2376E+00 -7.1543E+00 -2.0052E+00 -4.8582E-01 -1.7749E-01 -3.2743E-01 -2.0049E+00
             6.1186E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -3299.61342684564        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1834            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0425E+00  1.7708E+00  7.7987E-01  6.0171E-01  1.5848E+00  1.0154E+00  8.9736E-01  2.1201E+00  1.1067E+00  1.6565E+00
             1.6786E+00
 PARAMETER:  1.4164E-01  6.7143E-01 -1.4862E-01 -4.0798E-01  5.6046E-01  1.1527E-01 -8.2980E-03  8.5145E-01  2.0143E-01  6.0471E-01
             6.1796E-01
 GRADIENT:   2.2936E+01  2.6704E+01  7.0996E-01  1.0881E+01  5.6715E+00  1.9908E+00  9.0838E-01 -1.5095E-01 -2.6250E-02  1.8796E+00
             3.8834E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -3299.95461536677        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     1994
 NPARAMETR:  1.0427E+00  1.7928E+00  7.6427E-01  5.8730E-01  1.6013E+00  1.0162E+00  8.9049E-01  2.1474E+00  1.1271E+00  1.6631E+00
             1.6777E+00
 PARAMETER:  1.4180E-01  6.8377E-01 -1.6884E-01 -4.3222E-01  5.7084E-01  1.1609E-01 -1.5988E-02  8.6424E-01  2.1969E-01  6.0871E-01
             6.1745E-01
 GRADIENT:   3.5539E-01 -1.4266E+01 -2.1168E-01  5.0990E+00 -1.7082E+00  2.2122E-01  6.6946E-01  1.3995E-01  1.0227E-01 -4.6395E-01
             1.4776E+00

0ITERATION NO.:   95    OBJECTIVE VALUE:  -3300.00643188464        NO. OF FUNC. EVALS.: 143
 CUMULATIVE NO. OF FUNC. EVALS.:     2137
 NPARAMETR:  1.0426E+00  1.8048E+00  7.6420E-01  5.8717E-01  1.6018E+00  1.0162E+00  8.9053E-01  2.1483E+00  1.1270E+00  1.6626E+00
             1.6783E+00
 PARAMETER:  1.4173E-01  6.9047E-01 -1.6892E-01 -4.3244E-01  5.7113E-01  1.1603E-01 -1.5938E-02  8.6467E-01  2.1958E-01  6.0840E-01
             6.1776E-01
 GRADIENT:   9.0803E-02 -7.0868E-01  2.9117E-01  1.0320E+01 -3.0136E+00  1.5978E-01  9.8116E-01 -9.0304E-02 -1.8053E-01 -6.8435E-01
             1.6028E+00

0ITERATION NO.:  100    OBJECTIVE VALUE:  -3300.20511130600        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:     2276
 NPARAMETR:  1.0370E+00  1.8142E+00  7.4640E-01  5.7300E-01  1.6160E+00  1.0130E+00  8.8302E-01  2.1760E+00  1.1481E+00  1.6666E+00
             1.6770E+00
 PARAMETER:  1.3636E-01  6.9567E-01 -1.9249E-01 -4.5687E-01  5.7996E-01  1.1288E-01 -2.4408E-02  8.7747E-01  2.3813E-01  6.1076E-01
             6.1700E-01
 GRADIENT:  -1.1240E+01 -1.4324E+01 -1.4203E+00  4.4488E+00 -2.2583E+00 -1.1412E+00  6.0707E-01  8.8905E-01  4.4843E-01 -1.2461E+00
             3.6284E-01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -3300.47531143105        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2458
 NPARAMETR:  1.0335E+00  1.8468E+00  7.2634E-01  5.5422E-01  1.6276E+00  1.0141E+00  8.7151E-01  2.1353E+00  1.1751E+00  1.6722E+00
             1.6789E+00
 PARAMETER:  1.3297E-01  7.1345E-01 -2.1974E-01 -4.9019E-01  5.8709E-01  1.1402E-01 -3.7524E-02  8.5860E-01  2.6136E-01  6.1416E-01
             6.1811E-01
 GRADIENT:  -1.8506E+01 -6.5286E+00  7.3975E-04  2.8814E+00 -6.0269E+00 -8.3535E-01 -5.6305E-02 -8.1514E-01  2.8881E-02 -1.2329E+00
             6.3182E-01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -3300.48879260772        NO. OF FUNC. EVALS.: 172
 CUMULATIVE NO. OF FUNC. EVALS.:     2630
 NPARAMETR:  1.0335E+00  1.8481E+00  7.2581E-01  5.5390E-01  1.6287E+00  1.0153E+00  8.7048E-01  2.1366E+00  1.1756E+00  1.6788E+00
             1.6794E+00
 PARAMETER:  1.3296E-01  7.1417E-01 -2.2047E-01 -4.9078E-01  5.8776E-01  1.1521E-01 -3.8706E-02  8.5919E-01  2.6181E-01  6.1807E-01
             6.1844E-01
 GRADIENT:  -1.8429E+01 -5.8668E+00 -8.7690E-03  3.4438E+00 -5.8361E+00 -3.9328E-01 -2.2335E-01 -8.6146E-01  2.6845E-02 -1.5858E-01
             1.4771E+00

0ITERATION NO.:  115    OBJECTIVE VALUE:  -3300.55220657567        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     2790
 NPARAMETR:  1.0365E+00  1.8463E+00  7.2565E-01  5.5197E-01  1.6332E+00  1.0157E+00  8.7093E-01  2.1495E+00  1.1746E+00  1.6787E+00
             1.6788E+00
 PARAMETER:  1.3584E-01  7.1317E-01 -2.2068E-01 -4.9427E-01  5.9057E-01  1.1554E-01 -3.8191E-02  8.6524E-01  2.6089E-01  6.1799E-01
             6.1807E-01
 GRADIENT:   9.7962E+00  3.5356E+01  1.3318E-02  6.3890E+00  2.8688E+00  1.9528E+00  3.2371E-01 -4.4129E-01  2.8414E-01  1.6930E+00
             1.9429E+00

0ITERATION NO.:  120    OBJECTIVE VALUE:  -3300.70911628224        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:     2948
 NPARAMETR:  1.0340E+00  1.8704E+00  7.3796E-01  5.3395E-01  1.6692E+00  1.0142E+00  8.5647E-01  2.2792E+00  1.1950E+00  1.6939E+00
             1.6755E+00
 PARAMETER:  1.3339E-01  7.2615E-01 -2.0386E-01 -5.2744E-01  6.1234E-01  1.1411E-01 -5.4939E-02  9.2382E-01  2.7811E-01  6.2702E-01
             6.1611E-01
 GRADIENT:  -1.7389E+01 -1.5079E+01  2.6807E+00 -2.9008E+00 -1.9951E+00 -7.7186E-01 -1.8046E+00  3.1383E-01 -8.8161E-01  9.7045E-02
            -2.7382E+00

0ITERATION NO.:  125    OBJECTIVE VALUE:  -3301.25372173053        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     3124
 NPARAMETR:  1.0428E+00  1.9724E+00  6.4987E-01  4.7404E-01  1.7411E+00  1.0155E+00  8.3624E-01  2.2823E+00  1.3128E+00  1.7205E+00
             1.6808E+00
 PARAMETER:  1.4191E-01  7.7924E-01 -3.3098E-01 -6.4647E-01  6.5450E-01  1.1537E-01 -7.8836E-02  9.2518E-01  3.7214E-01  6.4262E-01
             6.1926E-01
 GRADIENT:   6.9093E-01  9.5408E-02 -2.0029E-01  4.6062E-01  1.6101E-01 -2.5932E-01  3.2786E-01  1.3213E-01 -3.2172E-02  1.2897E-01
            -4.9034E-01

0ITERATION NO.:  128    OBJECTIVE VALUE:  -3301.25567040806        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     3216
 NPARAMETR:  1.0425E+00  1.9775E+00  6.4656E-01  4.7045E-01  1.7446E+00  1.0162E+00  8.3337E-01  2.2786E+00  1.3229E+00  1.7207E+00
             1.6814E+00
 PARAMETER:  1.4157E-01  7.8181E-01 -3.3609E-01 -6.5408E-01  6.5655E-01  1.1605E-01 -8.2282E-02  9.2358E-01  3.7979E-01  6.4271E-01
             6.1963E-01
 GRADIENT:  -2.6160E-02  2.8375E-02 -7.5591E-03  7.5746E-03  2.0218E-03  1.7500E-03 -2.8753E-03  4.3492E-03 -1.2199E-02  8.6392E-03
            -2.0091E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     3216
 NO. OF SIG. DIGITS IN FINAL EST.:  3.4

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.4231E-04 -2.5381E-02 -2.5563E-02  3.5150E-02 -2.2855E-02
 SE:             2.9747E-02  2.5615E-02  1.4246E-02  1.9137E-02  2.6880E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7473E-01  3.2175E-01  7.2755E-02  6.6242E-02  3.9519E-01

 ETASHRINKSD(%)  3.4472E-01  1.4186E+01  5.2274E+01  3.5889E+01  9.9493E+00
 ETASHRINKVR(%)  6.8826E-01  2.6360E+01  7.7222E+01  5.8898E+01  1.8909E+01
 EBVSHRINKSD(%)  6.7744E-01  1.4329E+01  5.3813E+01  4.0441E+01  7.5872E+00
 EBVSHRINKVR(%)  1.3503E+00  2.6604E+01  7.8668E+01  6.4527E+01  1.4599E+01
 RELATIVEINF(%)  9.8641E+01  1.3372E+01  1.1237E+01  5.5445E+00  4.0308E+01
 EPSSHRINKSD(%)  1.9271E+01
 EPSSHRINKVR(%)  3.4827E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3301.2556704080580     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1647.1663106396472     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    86.37
 Elapsed covariance  time in seconds:    14.51
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3301.256       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.98E+00  6.47E-01  4.70E-01  1.74E+00  1.02E+00  8.33E-01  2.28E+00  1.32E+00  1.72E+00  1.68E+00
 


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
 
         3.19E-02  2.36E-01  2.35E-01  1.60E-01  1.79E-01  6.77E-02  1.01E-01  3.65E-01  4.18E-01  1.65E-01  1.20E-01
 


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
+        1.02E-03
 
 TH 2
+       -4.80E-04  5.58E-02
 
 TH 3
+       -2.75E-05 -4.87E-02  5.51E-02
 
 TH 4
+        4.19E-04 -3.66E-02  3.37E-02  2.56E-02
 
 TH 5
+       -5.86E-05  3.42E-02 -2.75E-02 -2.36E-02  3.20E-02
 
 TH 6
+       -2.39E-04  5.33E-04  6.36E-04 -1.49E-04  5.68E-04  4.59E-03
 
 TH 7
+       -8.40E-05 -1.58E-02  1.40E-02  1.09E-02 -9.88E-03 -6.48E-04  1.02E-02
 
 TH 8
+       -1.56E-03 -1.74E-02  4.37E-02  1.26E-02 -2.27E-03 -5.10E-04  3.27E-03  1.33E-01
 
 TH 9
+       -9.88E-04  8.60E-02 -7.60E-02 -5.88E-02  5.79E-02  7.36E-04 -2.64E-02 -5.16E-02  1.75E-01
 
 TH10
+       -1.87E-05  1.09E-02 -1.09E-02 -7.94E-03  1.03E-02 -1.80E-03 -4.94E-03  2.93E-03  2.02E-02  2.71E-02
 
 TH11
+       -3.46E-04  1.11E-02 -1.42E-02 -7.35E-03  3.09E-03 -6.19E-04 -3.51E-03 -7.91E-03  1.35E-02 -8.96E-05  1.44E-02
 
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
+        3.19E-02
 
 TH 2
+       -6.38E-02  2.36E-01
 
 TH 3
+       -3.67E-03 -8.79E-01  2.35E-01
 
 TH 4
+        8.21E-02 -9.70E-01  8.96E-01  1.60E-01
 
 TH 5
+       -1.03E-02  8.11E-01 -6.56E-01 -8.26E-01  1.79E-01
 
 TH 6
+       -1.11E-01  3.33E-02  4.00E-02 -1.38E-02  4.69E-02  6.77E-02
 
 TH 7
+       -2.61E-02 -6.64E-01  5.93E-01  6.77E-01 -5.48E-01 -9.50E-02  1.01E-01
 
 TH 8
+       -1.34E-01 -2.02E-01  5.10E-01  2.16E-01 -3.47E-02 -2.06E-02  8.87E-02  3.65E-01
 
 TH 9
+       -7.41E-02  8.71E-01 -7.74E-01 -8.79E-01  7.74E-01  2.60E-02 -6.26E-01 -3.38E-01  4.18E-01
 
 TH10
+       -3.56E-03  2.81E-01 -2.81E-01 -3.02E-01  3.50E-01 -1.61E-01 -2.98E-01  4.87E-02  2.94E-01  1.65E-01
 
 TH11
+       -9.03E-02  3.91E-01 -5.03E-01 -3.83E-01  1.44E-01 -7.61E-02 -2.90E-01 -1.80E-01  2.69E-01 -4.54E-03  1.20E-01
 
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
+        1.12E+03
 
 TH 2
+       -4.90E+01  3.60E+02
 
 TH 3
+        4.98E+01  1.05E+02  3.84E+02
 
 TH 4
+       -1.95E+02  3.04E+02 -4.23E+02  1.30E+03
 
 TH 5
+       -6.76E+01 -1.15E+01 -9.69E+00  1.02E+02  1.32E+02
 
 TH 6
+        8.00E+01 -3.90E+01 -4.69E+01  1.37E+00 -1.70E+01  2.49E+02
 
 TH 7
+        8.36E+01 -5.60E+00 -9.04E+00 -6.52E+01 -1.98E+01  4.13E+01  2.02E+02
 
 TH 8
+        1.65E+01 -3.20E+01 -9.21E+01  9.05E+01 -1.01E+01  1.57E+01  1.03E+01  3.31E+01
 
 TH 9
+        2.16E+01 -3.77E+01 -6.40E+01  9.57E+01 -1.47E+01  1.26E+01  1.48E+01  2.37E+01  4.43E+01
 
 TH10
+        1.19E+01  1.86E+01  4.47E+01 -4.39E+01 -1.22E+01  1.67E+01  1.47E+01 -1.20E+01 -8.95E+00  5.10E+01
 
 TH11
+        4.17E+01 -2.66E+00  9.00E+01 -7.10E+01  2.39E+01  7.76E+00  1.14E+01 -1.81E+01 -5.84E+00  1.66E+01  1.18E+02
 
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
 #CPUT: Total CPU Time in Seconds,      101.024
Stop Time:
Sat Sep 25 00:03:55 CDT 2021
