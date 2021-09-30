Wed Sep 29 13:05:12 CDT 2021
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
$DATA ../../../../data/spa/A2/dat83.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m83.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -974.151633453325        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5527E+02  1.1957E+02  1.1219E+02  4.5970E+01  2.6887E+01  2.7512E+01 -1.5889E+01 -4.1302E+01 -3.1727E+01 -9.0021E+01
            -1.1550E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1359.81424352224        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.1092E+00  8.8550E-01  8.3662E-01  1.1324E+00  8.7016E-01  1.1610E+00  9.2820E-01  1.0180E+00  9.8693E-01  1.1115E+00
             3.1253E+00
 PARAMETER:  2.0364E-01 -2.1608E-02 -7.8381E-02  2.2431E-01 -3.9082E-02  2.4927E-01  2.5497E-02  1.1786E-01  8.6842E-02  2.0573E-01
             1.2395E+00
 GRADIENT:   2.2341E+02  2.2811E+01 -5.2693E+00  5.3047E+01 -1.2217E+01  4.2544E+01  8.0793E+00  9.2742E+00  1.2838E+01  2.3668E+01
             1.0126E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1371.47705690616        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0969E+00  5.1483E-01  3.7204E-01  1.3451E+00  4.0935E-01  1.2909E+00  3.3916E-01  5.8512E-01  9.8607E-01  5.0707E-01
             2.7356E+00
 PARAMETER:  1.9248E-01 -5.6392E-01 -8.8875E-01  3.9644E-01 -7.9319E-01  3.5530E-01 -9.8129E-01 -4.3593E-01  8.5968E-02 -5.7911E-01
             1.1064E+00
 GRADIENT:   1.8315E+02  8.3213E+01  1.1826E+01  2.9552E+02 -2.1831E+01  7.5271E+01 -4.4833E-01  5.4344E+00 -1.4229E+01  1.5593E+00
             5.0305E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1397.40691993006        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      314
 NPARAMETR:  1.0115E+00  4.9115E-01  4.3234E-01  1.2897E+00  4.5967E-01  1.0635E+00  2.4826E-01  1.3180E-01  1.0864E+00  7.5159E-01
             2.2340E+00
 PARAMETER:  1.1144E-01 -6.1100E-01 -7.3854E-01  3.5443E-01 -6.7725E-01  1.6161E-01 -1.2933E+00 -1.9265E+00  1.8287E-01 -1.8556E-01
             9.0378E-01
 GRADIENT:   2.7034E+00  2.3086E+01 -2.2960E+01  1.1552E+02  3.3352E+01  1.0842E+01  3.2377E-02  1.5574E-01  2.2117E+01  5.8507E+00
            -3.4068E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1412.44461399543        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      489
 NPARAMETR:  1.0055E+00  3.3596E-01  2.7486E-01  1.1591E+00  3.0328E-01  1.0374E+00  7.0397E-02  1.0000E-02  9.5778E-01  6.0768E-01
             2.3902E+00
 PARAMETER:  1.0551E-01 -9.9075E-01 -1.1915E+00  2.4766E-01 -1.0931E+00  1.3671E-01 -2.5536E+00 -4.7719E+00  5.6861E-02 -3.9811E-01
             9.7136E-01
 GRADIENT:  -1.6492E+01  9.6153E+00  2.1253E+01  1.2993E+00 -3.7304E+01 -6.6415E-01  2.0658E-04  0.0000E+00 -2.8676E+00  1.3124E+00
             1.9944E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1413.15399097375        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      645
 NPARAMETR:  1.0119E+00  3.3076E-01  2.6892E-01  1.1534E+00  3.0140E-01  1.0398E+00  3.9695E-02  1.0000E-02  9.7655E-01  6.1442E-01
             2.2863E+00
 PARAMETER:  1.1188E-01 -1.0064E+00 -1.2133E+00  2.4268E-01 -1.0993E+00  1.3902E-01 -3.1265E+00 -4.7796E+00  7.6269E-02 -3.8708E-01
             9.2694E-01
 GRADIENT:   7.2006E+01  8.7943E+00  2.8270E+01  4.3666E+01  9.3695E+01  8.4967E+00  4.3003E-03  0.0000E+00  2.2758E+00  1.1584E+00
             7.2908E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1413.16657177157        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      828
 NPARAMETR:  1.0127E+00  3.1871E-01  2.7133E-01  1.1616E+00  3.0023E-01  1.0398E+00  1.0000E-02  1.0000E-02  9.7228E-01  6.1378E-01
             2.2909E+00
 PARAMETER:  1.1260E-01 -1.0435E+00 -1.2044E+00  2.4976E-01 -1.1032E+00  1.3898E-01 -7.5263E+00 -4.7796E+00  7.1886E-02 -3.8812E-01
             9.2894E-01
 GRADIENT:   1.1384E+00  1.2074E+00  1.4631E+00  1.6624E+00 -4.6009E+00  1.4375E-01  0.0000E+00  0.0000E+00  4.6236E-02 -1.5870E-01
             4.0722E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1413.19583336822        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1005
 NPARAMETR:  1.0111E+00  2.5647E-01  2.8866E-01  1.1999E+00  3.0142E-01  1.0359E+00  1.0000E-02  1.0000E-02  9.4779E-01  6.2107E-01
             2.2930E+00
 PARAMETER:  1.1105E-01 -1.2607E+00 -1.1425E+00  2.8227E-01 -1.0992E+00  1.3528E-01 -2.6134E+01 -4.7796E+00  4.6374E-02 -3.7631E-01
             9.2985E-01
 GRADIENT:   1.7755E+00  7.9559E-01  3.5377E+00  6.9683E-01 -1.0210E+01  1.6302E-01  0.0000E+00  0.0000E+00  4.2575E-02  5.6953E-02
            -7.6546E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1413.37225915393        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1184
 NPARAMETR:  1.0058E+00  2.0085E-01  3.2711E-01  1.2557E+00  3.2239E-01  1.0300E+00  1.0000E-02  1.0000E-02  9.1317E-01  6.3669E-01
             2.3113E+00
 PARAMETER:  1.0583E-01 -1.5052E+00 -1.0174E+00  3.2770E-01 -1.0320E+00  1.2958E-01 -4.5946E+01 -4.7796E+00  9.1722E-03 -3.5147E-01
             9.3782E-01
 GRADIENT:   3.8752E-02  3.1259E-01 -2.8655E-01 -4.9572E-01  1.0295E+00  3.1719E-02  0.0000E+00  0.0000E+00  1.9165E-01  3.0748E-01
            -4.0870E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1413.37412221166        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1359
 NPARAMETR:  1.0049E+00  1.9058E-01  3.3246E-01  1.2647E+00  3.2463E-01  1.0289E+00  1.0000E-02  1.0000E-02  9.0687E-01  6.3642E-01
             2.3179E+00
 PARAMETER:  1.0485E-01 -1.5577E+00 -1.0012E+00  3.3483E-01 -1.0251E+00  1.2851E-01 -5.0359E+01 -4.7796E+00  2.2386E-03 -3.5190E-01
             9.4068E-01
 GRADIENT:  -4.9995E-01  6.6211E-01  2.7162E+00  7.9649E-01 -2.8672E+00 -8.1029E-02  0.0000E+00  0.0000E+00 -1.4145E-01 -1.5259E-01
             1.0030E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1413.40262535707        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1540
 NPARAMETR:  1.0020E+00  1.4661E-01  3.3720E-01  1.2843E+00  3.2351E-01  1.0268E+00  1.0000E-02  1.0000E-02  8.9392E-01  6.3772E-01
             2.3247E+00
 PARAMETER:  1.0199E-01 -1.8200E+00 -9.8709E-01  3.5023E-01 -1.0285E+00  1.2643E-01 -7.4732E+01 -4.7796E+00 -1.2140E-02 -3.4985E-01
             9.4359E-01
 GRADIENT:  -2.2901E-01  2.0751E-01  1.6755E+00  8.7369E-01 -2.5069E+00 -7.2078E-02  0.0000E+00  0.0000E+00 -3.1201E-01 -5.2081E-01
             7.0471E-01

0ITERATION NO.:   52    OBJECTIVE VALUE:  -1413.40473223256        NO. OF FUNC. EVALS.:  66
 CUMULATIVE NO. OF FUNC. EVALS.:     1606
 NPARAMETR:  1.0019E+00  1.4057E-01  3.3744E-01  1.2848E+00  3.2358E-01  1.0270E+00  1.0000E-02  1.0000E-02  8.9418E-01  6.4260E-01
             2.3184E+00
 PARAMETER:  1.0162E-01 -1.8560E+00 -9.8581E-01  3.5168E-01 -1.0286E+00  1.2632E-01 -7.8100E+01 -4.7796E+00 -1.2662E-02 -3.4569E-01
             9.4218E-01
 GRADIENT:  -2.7738E-01  2.9205E-02  4.1657E-01  1.1044E+00 -3.6511E-01 -4.3781E-02  0.0000E+00  0.0000E+00 -9.8656E-02 -1.9624E-01
             3.0496E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1606
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5469E-04 -2.4300E-05  1.4406E-05 -8.4843E-03 -2.5366E-03
 SE:             2.9314E-02  2.0954E-05  2.7169E-04  2.7209E-02  2.2530E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9579E-01  2.4618E-01  9.5771E-01  7.5518E-01  9.1036E-01

 ETASHRINKSD(%)  1.7933E+00  9.9930E+01  9.9090E+01  8.8463E+00  2.4522E+01
 ETASHRINKVR(%)  3.5544E+00  1.0000E+02  9.9992E+01  1.6910E+01  4.3031E+01
 EBVSHRINKSD(%)  1.8515E+00  9.9929E+01  9.9103E+01  8.1015E+00  2.4079E+01
 EBVSHRINKVR(%)  3.6688E+00  1.0000E+02  9.9992E+01  1.5547E+01  4.2360E+01
 RELATIVEINF(%)  8.3100E+01  5.1611E-06  2.7991E-04  1.6605E+01  1.5465E+00
 EPSSHRINKSD(%)  3.4398E+01
 EPSSHRINKVR(%)  5.6964E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1413.4047322325578     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -678.25390566881958     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    20.05
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.15
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1413.405       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.41E-01  3.38E-01  1.29E+00  3.23E-01  1.03E+00  1.00E-02  1.00E-02  8.93E-01  6.40E-01  2.32E+00
 


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
+        1.01E+03
 
 TH 2
+       -9.58E+01  2.82E+02
 
 TH 3
+       -1.79E+01  8.11E+02  7.83E+03
 
 TH 4
+       -3.17E+01  2.72E+02 -4.92E+02  7.13E+02
 
 TH 5
+        1.31E+02 -1.50E+03 -1.02E+04 -2.16E+02  1.50E+04
 
 TH 6
+        1.89E-01 -1.14E+01  1.47E+01 -8.73E+00  6.18E+00  1.75E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.28E+00 -4.62E+01  8.24E+01 -1.33E+01  2.78E+01  2.56E+00  0.00E+00  0.00E+00  1.77E+02
 
 TH10
+       -5.61E+00  7.38E+00 -1.19E+02  5.69E+00  4.47E+01  1.25E-01  0.00E+00  0.00E+00 -3.06E+00  1.63E+02
 
 TH11
+       -1.19E+01 -2.26E+00 -2.45E+01 -8.18E+00  6.39E+00  1.87E+00  0.00E+00  0.00E+00  1.01E+01  3.09E+01  5.00E+01
 
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
 #CPUT: Total CPU Time in Seconds,       26.257
Stop Time:
Wed Sep 29 13:05:40 CDT 2021
