Thu Sep 30 09:22:52 CDT 2021
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
$DATA ../../../../data/spa2/D/dat59.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m59.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   16229.7437378849        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.4618E+02  2.5567E+02  2.5355E+01  2.1632E+02  3.5841E+02 -2.7727E+03 -1.2243E+03 -5.4277E+01 -1.5370E+03 -9.7882E+02
            -3.0644E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -717.725020244149        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0563E+00  2.3333E+00  8.8998E-01  3.9304E+00  9.2346E-01  3.8962E+00  5.6244E+00  9.7645E-01  4.8000E+00  1.7053E+00
             9.7454E+00
 PARAMETER:  1.5473E-01  9.4730E-01 -1.6558E-02  1.4687E+00  2.0369E-02  1.4600E+00  1.8271E+00  7.6164E-02  1.6686E+00  6.3373E-01
             2.3768E+00
 GRADIENT:  -2.7897E+01  4.2013E+01 -4.6997E+01  1.1898E+02 -4.3233E+01  1.5766E+02  7.1299E+01  3.4204E+00  6.2749E+01  2.2878E+01
             1.4995E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -729.872780006220        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:      272             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0576E+00  1.4982E+00  8.9242E-01  3.9049E+00  1.2061E+00  3.8941E+00  5.6058E+00  6.0094E-01  4.7878E+00  1.0996E+00
             9.7273E+00
 PARAMETER:  1.5600E-01  5.0426E-01 -1.3818E-02  1.4622E+00  2.8739E-01  1.4595E+00  1.8238E+00 -4.0926E-01  1.6661E+00  1.9495E-01
             2.3749E+00
 GRADIENT:  -3.1011E+01 -2.6379E+00 -7.1091E+01  1.2059E+02  2.1633E+01  1.4320E+02  9.3277E+01  1.1590E+00  3.8473E+01  1.0464E+01
             1.7635E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -768.977619800624        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      345
 NPARAMETR:  1.0594E+00  1.5525E+00  8.9195E-01  3.2074E+00  1.2111E+00  3.6263E+00  4.6430E+00  1.0000E-02  4.3622E+00  1.0000E-02
             7.4682E+00
 PARAMETER:  1.5769E-01  5.3987E-01 -1.4347E-02  1.2655E+00  2.9155E-01  1.3882E+00  1.6354E+00 -5.3354E+00  1.5730E+00 -5.3884E+00
             2.1107E+00
 GRADIENT:  -1.8697E+01  3.4557E+00 -5.2719E+01  1.4128E+02  1.6546E+01  1.6092E+02  6.2256E+01  0.0000E+00  1.7842E+01  0.0000E+00
            -9.3749E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -791.675189619236        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      484
 NPARAMETR:  1.0602E+00  1.5679E+00  8.9177E-01  2.8486E+00  1.1317E+00  3.4740E+00  4.2424E+00  1.0000E-02  4.1672E+00  1.0000E-02
             7.0016E+00
 PARAMETER:  1.5848E-01  5.4973E-01 -1.4550E-02  1.1468E+00  2.2372E-01  1.3453E+00  1.5451E+00 -7.2191E+00  1.5272E+00 -7.4914E+00
             2.0461E+00
 GRADIENT:  -4.9175E+01 -3.6653E+00 -5.6177E+01  1.2166E+02 -4.6366E+00  1.0305E+01 -5.6886E+01  0.0000E+00 -1.6113E+01  0.0000E+00
            -1.5600E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -894.301659161321        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      662
 NPARAMETR:  1.0672E+00  2.5781E+00  8.9256E-01  6.2336E-01  1.6636E+00  2.3882E+00  3.3155E+00  1.0000E-02  3.2159E+00  1.0000E-02
             9.0825E+00
 PARAMETER:  1.6504E-01  1.0471E+00 -1.3659E-02 -3.7262E-01  6.0901E-01  9.7053E-01  1.2986E+00 -1.6772E+01  1.2681E+00 -1.8821E+01
             2.3063E+00
 GRADIENT:  -9.5501E+01 -1.2732E+01 -1.1705E+01  2.3894E+01 -8.8405E+00 -8.8258E+01 -8.4066E+01  0.0000E+00  2.3864E+01  0.0000E+00
             1.5100E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -923.423755599538        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      839
 NPARAMETR:  1.0698E+00  2.8045E+00  8.9273E-01  5.1445E-01  1.8383E+00  2.4726E+00  4.1951E+00  1.0000E-02  3.0276E+00  1.0000E-02
             7.8045E+00
 PARAMETER:  1.6749E-01  1.1312E+00 -1.3472E-02 -5.6466E-01  7.0882E-01  1.0053E+00  1.5339E+00 -1.7787E+01  1.2078E+00 -2.0073E+01
             2.1547E+00
 GRADIENT:  -7.5381E+01  4.2997E+00 -1.3846E+01  1.3577E+01  4.9301E+00 -7.9767E+01  8.9569E+00  0.0000E+00  1.9505E+01  0.0000E+00
            -2.7689E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -936.381104372027        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     1003
 NPARAMETR:  1.0726E+00  2.7966E+00  1.7489E+00  5.1332E-01  1.8381E+00  2.4724E+00  4.1748E+00  1.0000E-02  5.3197E-01  1.0000E-02
             7.7813E+00
 PARAMETER:  1.7006E-01  1.1284E+00  6.5898E-01 -5.6685E-01  7.0874E-01  1.0052E+00  1.5291E+00 -1.7787E+01 -5.3118E-01 -2.0073E+01
             2.1517E+00
 GRADIENT:  -5.4988E+01  5.0866E+01  2.1397E+00 -6.3572E+00 -1.5931E+01 -4.4808E+00  5.7845E+01  0.0000E+00  9.2632E-01  0.0000E+00
             1.2666E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -937.710519298544        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:     1141
 NPARAMETR:  1.0728E+00  2.7861E+00  1.6284E+00  5.1365E-01  1.8400E+00  2.4937E+00  4.1707E+00  1.0000E-02  2.0593E-01  1.0000E-02
             7.8228E+00
 PARAMETER:  1.7026E-01  1.1247E+00  5.8758E-01 -5.6621E-01  7.0974E-01  1.0138E+00  1.5281E+00 -1.7787E+01 -1.4802E+00 -2.0073E+01
             2.1570E+00
 GRADIENT:  -7.1680E+01  1.3481E+01  1.9082E+00 -1.0254E+01 -1.4944E+01 -7.2846E+01 -1.4511E+01  0.0000E+00  1.3626E-01  0.0000E+00
            -8.2298E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -941.718646975751        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1324
 NPARAMETR:  1.0750E+00  2.7001E+00  2.1040E+00  5.1668E-01  1.8566E+00  2.7274E+00  4.2108E+00  1.0000E-02  1.0000E-02  1.0000E-02
             8.2443E+00
 PARAMETER:  1.7231E-01  1.0933E+00  8.4382E-01 -5.6033E-01  7.1875E-01  1.1034E+00  1.5376E+00 -1.7787E+01 -8.1527E+00 -2.0073E+01
             2.2095E+00
 GRADIENT:  -6.3207E+01  7.4359E+00  2.2809E+00 -1.6060E+01 -1.8761E+01 -3.0944E+01 -9.8823E+00  0.0000E+00  0.0000E+00  0.0000E+00
             6.2689E+01

0ITERATION NO.:   49    OBJECTIVE VALUE:  -941.810047858730        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     1492
 NPARAMETR:  1.0733E+00  2.7263E+00  2.0562E+00  5.1954E-01  1.8701E+00  2.7042E+00  4.2753E+00  1.0000E-02  1.0000E-02  1.0000E-02
             8.3289E+00
 PARAMETER:  1.7231E-01  1.0931E+00  8.1337E-01 -5.6028E-01  7.1883E-01  1.1036E+00  1.5375E+00 -1.7787E+01 -8.1945E+00 -2.0073E+01
             2.2093E+00
 GRADIENT:   1.8539E+03 -2.9492E+02 -4.0504E+02 -3.1152E+02 -2.4678E+02  2.4513E+02 -1.1690E+02  0.0000E+00  0.0000E+00  0.0000E+00
            -7.8557E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1492
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.4545E-02 -1.5991E-02 -7.5259E-06 -3.2041E-04  4.2671E-05
 SE:             3.0664E-02  2.9441E-02  8.9704E-06  8.1323E-05  7.3312E-05
 N:                     100         100         100         100         100

 P VAL.:         5.8307E-03  5.8701E-01  4.0149E-01  8.1504E-05  5.6053E-01

 ETASHRINKSD(%)  1.0000E-10  1.3690E+00  9.9970E+01  9.9728E+01  9.9754E+01
 ETASHRINKVR(%)  1.0000E-10  2.7192E+00  1.0000E+02  9.9999E+01  9.9999E+01
 EBVSHRINKSD(%)  2.7665E+00  2.8655E+00  9.9961E+01  9.9823E+01  9.9668E+01
 EBVSHRINKVR(%)  5.4565E+00  5.6490E+00  1.0000E+02  1.0000E+02  9.9999E+01
 RELATIVEINF(%)  9.4191E+01  5.2001E+01  4.0941E-06  1.6577E-04  3.0132E-04
 EPSSHRINKSD(%)  1.4047E+01
 EPSSHRINKVR(%)  2.6122E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -941.81004785873040     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       160.91619198687670     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    35.74
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    10.16
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -941.810       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.07E+00  2.70E+00  2.04E+00  5.17E-01  1.86E+00  2.73E+00  4.21E+00  1.00E-02  1.00E-02  1.00E-02  8.24E+00
 


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
+        2.40E+05
 
 TH 2
+        8.05E+01  9.58E+02
 
 TH 3
+       -3.42E-01  2.12E-01  3.01E+03
 
 TH 4
+        8.23E+02  1.82E+01  2.96E-01  9.90E+04
 
 TH 5
+        3.34E+02 -2.41E+01 -3.32E+00 -2.03E+02  4.64E+03
 
 TH 6
+       -5.71E+01  1.29E+00 -2.16E-03  1.30E+00  1.98E+01  8.33E+02
 
 TH 7
+        3.74E+01 -7.15E+00 -1.62E-01 -2.46E+01 -8.88E+00  6.38E+00  2.04E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        9.04E+00 -1.46E+00 -6.85E-02 -8.60E+00 -2.71E+00  1.69E+02 -1.82E-01  0.00E+00  0.00E+00  0.00E+00  3.09E+01
 
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
 #CPUT: Total CPU Time in Seconds,       46.000
Stop Time:
Thu Sep 30 09:23:40 CDT 2021
