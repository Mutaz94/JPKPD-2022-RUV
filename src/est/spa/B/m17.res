Wed Sep 29 11:00:33 CDT 2021
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
$DATA ../../../../data/spa/B/dat17.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      500
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

 TOT. NO. OF OBS RECS:      400
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1652.80965336836        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.5157E+02 -3.1870E+01  2.4094E+01 -4.4598E+01 -2.7213E+01  4.5050E+01 -1.3656E+01  7.1683E-01  4.8879E+00 -2.9380E+01
            -2.9561E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1664.12862871455        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      194
 NPARAMETR:  1.0128E+00  1.1405E+00  1.0077E+00  1.0109E+00  1.1343E+00  1.0061E+00  1.1602E+00  9.6847E-01  9.2669E-01  1.3365E+00
             1.1594E+00
 PARAMETER:  1.1274E-01  2.3151E-01  1.0772E-01  1.1079E-01  2.2604E-01  1.0608E-01  2.4857E-01  6.7959E-02  2.3861E-02  3.9005E-01
             2.4790E-01
 GRADIENT:  -5.4823E+01  6.0636E+00 -1.7834E+01  2.4794E+01  3.9505E-01  9.8134E-02  4.6516E-02  4.7194E+00  7.7254E-01  8.6071E+00
             3.5203E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1665.45370401768        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  1.0304E+00  1.3260E+00  1.0680E+00  9.0193E-01  1.3359E+00  1.0094E+00  1.0877E+00  1.1106E+00  1.0648E+00  1.4787E+00
             1.0837E+00
 PARAMETER:  1.2996E-01  3.8216E-01  1.6576E-01 -3.2133E-03  3.8961E-01  1.0934E-01  1.8403E-01  2.0490E-01  1.6282E-01  4.9116E-01
             1.8035E-01
 GRADIENT:  -1.4898E+01  1.0264E+01 -1.1653E+01  2.6670E+01  2.6623E+01  2.2021E+00  9.3754E+00  9.9524E-02  6.0653E+00  3.7553E+00
             6.1839E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1667.24351696538        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      546
 NPARAMETR:  1.0395E+00  1.3556E+00  9.5443E-01  8.5881E-01  1.2254E+00  1.0051E+00  9.8482E-01  9.4910E-01  1.0822E+00  1.3435E+00
             1.0639E+00
 PARAMETER:  1.3875E-01  4.0422E-01  5.3360E-02 -5.2207E-02  3.0324E-01  1.0512E-01  8.4701E-02  4.7761E-02  1.7900E-01  3.9526E-01
             1.6197E-01
 GRADIENT:   2.9331E+00 -4.8771E-01  2.8720E+00 -7.9339E-01 -2.2867E+00  2.3874E-01 -3.6136E-01 -5.3441E-01  5.9885E-01 -4.0304E-01
            -1.8393E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1667.30525712359        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      723
 NPARAMETR:  1.0402E+00  1.5047E+00  8.5186E-01  7.6882E-01  1.2559E+00  1.0069E+00  9.0497E-01  8.7456E-01  1.1785E+00  1.3582E+00
             1.0666E+00
 PARAMETER:  1.3938E-01  5.0860E-01 -6.0328E-02 -1.6290E-01  3.2783E-01  1.0685E-01  1.4621E-04 -3.4031E-02  2.6426E-01  4.0616E-01
             1.6451E-01
 GRADIENT:   2.3262E+00  1.2137E+01  5.3634E+00  6.4972E+00 -6.0189E+00  4.9707E-01 -1.7135E+00 -8.6577E-01 -1.8991E-01 -9.2151E-02
            -1.7488E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1667.49328429127        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      906
 NPARAMETR:  1.0384E+00  1.7078E+00  6.4987E-01  6.3814E-01  1.2829E+00  1.0070E+00  8.5118E-01  6.3772E-01  1.2970E+00  1.3499E+00
             1.0763E+00
 PARAMETER:  1.3773E-01  6.3519E-01 -3.3098E-01 -3.4920E-01  3.4916E-01  1.0696E-01 -6.1137E-02 -3.4986E-01  3.6008E-01  4.0003E-01
             1.7355E-01
 GRADIENT:  -4.5087E+00  1.6755E+01 -1.0367E+00  1.2629E+01 -1.0546E+00 -1.1329E-01 -6.9362E-01  1.7134E-01 -1.8732E+00  7.4872E-01
             2.5270E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1667.78153744301        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1088
 NPARAMETR:  1.0420E+00  1.8393E+00  5.4800E-01  5.3205E-01  1.3242E+00  1.0079E+00  8.1153E-01  4.1149E-01  1.4468E+00  1.3592E+00
             1.0696E+00
 PARAMETER:  1.4114E-01  7.0938E-01 -5.0148E-01 -5.3102E-01  3.8081E-01  1.0790E-01 -1.0884E-01 -7.8797E-01  4.6936E-01  4.0690E-01
             1.6725E-01
 GRADIENT:   3.0023E+00 -6.1848E+00  3.4621E+00 -5.3421E+00 -5.3515E+00  1.6079E-01  9.4832E-02  2.1352E-03 -2.7800E+00 -1.9790E-01
            -8.6567E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1667.84312958105        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1265
 NPARAMETR:  1.0398E+00  1.8830E+00  5.1982E-01  5.0900E-01  1.3468E+00  1.0077E+00  7.9706E-01  2.1877E-01  1.5219E+00  1.3735E+00
             1.0730E+00
 PARAMETER:  1.3901E-01  7.3288E-01 -5.5428E-01 -5.7531E-01  3.9773E-01  1.0768E-01 -1.2683E-01 -1.4198E+00  5.1999E-01  4.1738E-01
             1.7048E-01
 GRADIENT:  -1.9929E+00 -8.2111E-01  1.3830E+00  2.1137E+00 -1.1229E+00 -1.5443E-02 -3.2027E-01  9.9192E-04  5.1988E-02 -6.6359E-03
             5.0732E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1667.84350715540        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     1458
 NPARAMETR:  1.0397E+00  1.8984E+00  5.1040E-01  4.9922E-01  1.3535E+00  1.0077E+00  7.9291E-01  1.7356E-01  1.5419E+00  1.3777E+00
             1.0733E+00
 PARAMETER:  1.3894E-01  7.4099E-01 -5.7255E-01 -5.9471E-01  4.0270E-01  1.0770E-01 -1.3205E-01 -1.6512E+00  5.3299E-01  4.2039E-01
             1.7077E-01
 GRADIENT:  -2.1938E+00 -2.9851E-01  1.4424E+00  2.3084E+00 -1.1514E+00 -1.9321E-02 -3.7255E-01  1.4360E-03  1.0826E-01  5.3915E-02
             6.5035E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1667.87719205896        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     1652
 NPARAMETR:  1.0426E+00  1.8981E+00  5.0311E-01  4.9473E-01  1.3553E+00  1.0081E+00  7.9455E-01  1.0016E-01  1.5430E+00  1.3773E+00
             1.0728E+00
 PARAMETER:  1.4169E-01  7.4086E-01 -5.8694E-01 -6.0375E-01  4.0400E-01  1.0804E-01 -1.2997E-01 -2.2010E+00  5.3373E-01  4.2009E-01
             1.7025E-01
 GRADIENT:   3.9057E+00 -9.5647E+00  2.9343E-01  1.2095E-01  1.0113E+00  1.1949E-01  1.7709E-01  2.7610E-03  1.0224E-01  2.1704E-01
             2.2231E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1667.88338776097        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1837             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0422E+00  1.8978E+00  4.9992E-01  4.9399E-01  1.3492E+00  1.0082E+00  7.9355E-01  1.0000E-02  1.5388E+00  1.3714E+00
             1.0716E+00
 PARAMETER:  1.4130E-01  7.4068E-01 -5.9330E-01 -6.0524E-01  3.9953E-01  1.0816E-01 -1.3124E-01 -5.4080E+00  5.3101E-01  4.1580E-01
             1.6918E-01
 GRADIENT:   5.4540E+02  7.3247E+02  3.5621E+00  9.2121E+01  1.6862E+01  4.3906E+01  7.0588E+00  0.0000E+00  1.6615E+01  5.2413E+00
             1.1755E+00

0ITERATION NO.:   52    OBJECTIVE VALUE:  -1667.88338776097        NO. OF FUNC. EVALS.:  64
 CUMULATIVE NO. OF FUNC. EVALS.:     1901
 NPARAMETR:  1.0421E+00  1.8974E+00  4.9697E-01  4.9472E-01  1.3501E+00  1.0082E+00  7.9433E-01  1.0000E-02  1.5423E+00  1.3708E+00
             1.0720E+00
 PARAMETER:  1.4130E-01  7.4068E-01 -5.9330E-01 -6.0524E-01  3.9953E-01  1.0816E-01 -1.3124E-01 -5.4080E+00  5.3101E-01  4.1580E-01
             1.6918E-01
 GRADIENT:   4.3720E-02  5.2250E-01  5.5607E-01 -7.7330E-01 -5.3676E-01  1.5721E-02 -9.9708E-02  0.0000E+00 -3.7677E-01  7.3665E-02
            -1.7716E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1901
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3110E-05 -3.7148E-02 -3.0318E-04  3.4369E-02 -4.7613E-02
 SE:             2.9829E-02  2.3934E-02  9.6394E-05  2.1974E-02  2.2456E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9965E-01  1.2064E-01  1.6597E-03  1.1780E-01  3.3982E-02

 ETASHRINKSD(%)  6.7528E-02  1.9818E+01  9.9677E+01  2.6386E+01  2.4770E+01
 ETASHRINKVR(%)  1.3501E-01  3.5708E+01  9.9999E+01  4.5810E+01  4.3405E+01
 EBVSHRINKSD(%)  4.8659E-01  1.8649E+01  9.9743E+01  2.9915E+01  2.1322E+01
 EBVSHRINKVR(%)  9.7081E-01  3.3820E+01  9.9999E+01  5.0880E+01  3.8098E+01
 RELATIVEINF(%)  9.8946E+01  5.0902E+00  8.9461E-05  3.5542E+00  1.9684E+01
 EPSSHRINKSD(%)  4.4169E+01
 EPSSHRINKVR(%)  6.8829E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1667.8833877609663     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -932.73256119722816     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    26.97
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.70
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1667.883       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.90E+00  5.00E-01  4.94E-01  1.35E+00  1.01E+00  7.94E-01  1.00E-02  1.54E+00  1.37E+00  1.07E+00
 


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
+        9.97E+02
 
 TH 2
+       -5.44E+00  3.48E+02
 
 TH 3
+        1.09E+01  1.27E+02  3.40E+02
 
 TH 4
+       -1.63E+01  3.10E+02 -2.92E+02  9.74E+02
 
 TH 5
+       -1.20E+00 -8.16E+01 -1.55E+02  1.51E+02  2.05E+02
 
 TH 6
+       -8.32E-02 -9.52E-01  2.66E+00 -4.46E+00 -3.59E-01  1.92E+02
 
 TH 7
+        4.17E-01  4.61E+00 -1.26E+01 -1.53E+01 -7.11E+00 -7.41E-01  1.47E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        6.52E-01 -1.67E+01 -3.52E+01  6.35E+01  5.68E-01 -5.72E-01  1.52E+01  0.00E+00  3.16E+01
 
 TH10
+        5.43E-01 -1.03E+01 -2.42E+01  3.78E+00 -3.53E+01  2.58E-01  6.93E+00  0.00E+00  4.10E+00  4.77E+01
 
 TH11
+       -7.29E+00 -2.11E+01 -3.85E+01  8.21E+00 -2.52E+00  2.95E+00  1.43E+01  0.00E+00  5.30E+00  7.45E+00  1.86E+02
 
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
 #CPUT: Total CPU Time in Seconds,       33.714
Stop Time:
Wed Sep 29 11:01:08 CDT 2021
