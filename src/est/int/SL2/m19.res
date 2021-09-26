Sat Sep 25 00:59:43 CDT 2021
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
$DATA ../../../../data/int/SL2/dat19.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      999
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

 TOT. NO. OF OBS RECS:      899
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
 RAW OUTPUT FILE (FILE): m19.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1839.48436949548        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.5190E+01 -4.1262E+01  1.2067E+02  1.0366E+02  7.2391E+01  1.5691E+01 -9.2376E+01 -1.4433E+02 -9.0909E+01 -4.4245E+01
            -3.7926E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3019.14539838284        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0623E+00  1.6531E+00  8.2544E-01  6.1825E-01  1.5058E+00  8.3530E-01  9.6939E-01  8.6427E-01  7.1486E-01  1.4953E+00
             1.9784E+00
 PARAMETER:  1.6042E-01  6.0267E-01 -9.1843E-02 -3.8086E-01  5.0933E-01 -7.9962E-02  6.8915E-02 -4.5868E-02 -2.3567E-01  5.0230E-01
             7.8227E-01
 GRADIENT:   9.0636E+01  1.3222E+00  4.6694E-01  1.0976E+01  4.0533E+01 -5.9323E+01  4.2556E+01 -3.5827E+00 -1.4335E+01 -2.0584E+01
            -2.6922E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3029.23589462258        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      163
 NPARAMETR:  1.0574E+00  1.7618E+00  7.5207E-01  5.6284E-01  1.5195E+00  9.5355E-01  8.0737E-01  6.7199E-01  5.0716E-01  1.7809E+00
             2.0150E+00
 PARAMETER:  1.5578E-01  6.6633E-01 -1.8493E-01 -4.7475E-01  5.1839E-01  5.2437E-02 -1.1397E-01 -2.9752E-01 -5.7894E-01  6.7711E-01
             8.0062E-01
 GRADIENT:   6.1611E+01  3.8425E+01  7.9582E+00  3.5210E+01  4.3707E+00 -1.0731E+00 -2.1490E+00 -3.0475E+00 -1.3382E+01  1.6505E+01
            -2.2961E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3048.65728904597        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      235
 NPARAMETR:  1.0323E+00  1.8219E+00  7.0450E-01  5.1935E-01  1.5369E+00  9.4928E-01  6.9326E-01  1.2777E-01  1.2655E+00  1.5844E+00
             2.1875E+00
 PARAMETER:  1.3183E-01  6.9988E-01 -2.5027E-01 -5.5518E-01  5.2980E-01  4.7944E-02 -2.6635E-01 -1.9575E+00  3.3546E-01  5.6021E-01
             8.8277E-01
 GRADIENT:  -5.4764E+00  5.9819E+00 -2.5227E+00  1.1402E+01  1.0171E+00 -2.0097E+00  3.4583E+00 -3.2106E-02  7.7002E+00 -2.3502E+00
             7.4455E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3052.49787524335        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:      336
 NPARAMETR:  1.0383E+00  2.0434E+00  5.1055E-01  3.7982E-01  1.6862E+00  9.3769E-01  7.2092E-01  1.5042E+00  1.1443E+00  1.7111E+00
             2.1519E+00
 PARAMETER:  1.3762E-01  8.1460E-01 -5.7227E-01 -8.6806E-01  6.2247E-01  3.5661E-02 -2.2723E-01  5.0824E-01  2.3480E-01  6.3711E-01
             8.6636E-01
 GRADIENT:  -4.8703E+00 -2.2932E+00 -2.6382E+00  3.6614E+00 -5.4437E+00 -7.9933E+00  1.6533E+01 -9.6812E-01 -1.4093E+00 -3.6116E+00
            -1.4877E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3055.06720198796        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      512
 NPARAMETR:  1.0396E+00  2.2015E+00  4.3815E-01  2.8382E-01  1.8397E+00  9.6082E-01  6.2573E-01  2.2723E+00  1.4551E+00  1.8057E+00
             2.1647E+00
 PARAMETER:  1.3884E-01  8.8912E-01 -7.2520E-01 -1.1594E+00  7.0958E-01  6.0028E-02 -3.6883E-01  9.2081E-01  4.7509E-01  6.9097E-01
             8.7229E-01
 GRADIENT:  -2.6142E+00  1.7024E+01 -4.2487E-01  6.7416E+00  1.4561E+00  1.1733E+00 -3.9234E+00  2.8888E-01 -1.2183E-01  7.6802E-01
             4.2370E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3055.44795100863        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      687
 NPARAMETR:  1.0406E+00  2.2926E+00  3.0752E-01  2.1496E-01  1.8994E+00  9.5856E-01  6.2031E-01  2.0216E+00  1.7445E+00  1.8376E+00
             2.1584E+00
 PARAMETER:  1.3984E-01  9.2968E-01 -1.0792E+00 -1.4373E+00  7.4151E-01  5.7682E-02 -3.7754E-01  8.0387E-01  6.5645E-01  7.0848E-01
             8.6937E-01
 GRADIENT:   3.9409E-01  8.1038E-02 -5.2693E-01  5.4077E-01  9.5755E-02  2.0968E-01 -4.5899E-01  1.8512E-01 -3.5171E-01 -3.0723E-01
            -7.8287E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3055.57721674781        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      864
 NPARAMETR:  1.0420E+00  2.3838E+00  2.0781E-01  1.5601E-01  1.9528E+00  9.5706E-01  6.1332E-01  1.0955E+00  2.3682E+00  1.8560E+00
             2.1568E+00
 PARAMETER:  1.4115E-01  9.6871E-01 -1.4712E+00 -1.7578E+00  7.6929E-01  5.6112E-02 -3.8887E-01  1.9120E-01  9.6212E-01  7.1842E-01
             8.6860E-01
 GRADIENT:   3.5219E+00  1.1933E+01 -1.2071E+00  3.5137E+00 -2.4120E+00 -7.1885E-01  1.6443E+00  5.5130E-01  2.7231E+00 -3.8208E-01
            -3.9700E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3055.83455400624        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1040
 NPARAMETR:  1.0395E+00  2.4100E+00  1.7943E-01  1.3087E-01  1.9879E+00  9.5822E-01  6.0549E-01  4.0331E-01  2.5898E+00  1.8708E+00
             2.1541E+00
 PARAMETER:  1.3874E-01  9.7964E-01 -1.6179E+00 -1.9336E+00  7.8710E-01  5.7324E-02 -4.0171E-01 -8.0804E-01  1.0516E+00  7.2638E-01
             8.6737E-01
 GRADIENT:  -1.9186E+00 -4.2765E+00  4.6300E-02 -8.5662E-01 -5.5248E-01 -2.8117E-01 -8.5769E-02  8.5607E-02  1.1584E-01 -2.0918E-01
            -1.4766E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3055.89408474998        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1215
 NPARAMETR:  1.0403E+00  2.3982E+00  1.8961E-01  1.4067E-01  1.9774E+00  9.5873E-01  6.0769E-01  9.1729E-02  2.4910E+00  1.8647E+00
             2.1558E+00
 PARAMETER:  1.3951E-01  9.7471E-01 -1.5628E+00 -1.8614E+00  7.8179E-01  5.7850E-02 -3.9809E-01 -2.2889E+00  1.0127E+00  7.2311E-01
             8.6816E-01
 GRADIENT:  -2.4257E-01  2.2367E-01 -3.7865E-02  1.0523E-01  1.9443E-01 -6.9891E-02 -1.0037E-02  3.8512E-03  2.6708E-02  5.9226E-02
            -3.9274E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3055.89609117048        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1390
 NPARAMETR:  1.0404E+00  2.3986E+00  1.8933E-01  1.4025E-01  1.9773E+00  9.5889E-01  6.0765E-01  1.0000E-02  2.4955E+00  1.8645E+00
             2.1558E+00
 PARAMETER:  1.3960E-01  9.7489E-01 -1.5643E+00 -1.8643E+00  7.8173E-01  5.8023E-02 -3.9816E-01 -4.6635E+00  1.0145E+00  7.2298E-01
             8.6817E-01
 GRADIENT:  -1.1211E-02  6.1322E-03 -3.7146E-03  7.0047E-03  1.2515E-02 -2.5636E-03  3.0231E-03  0.0000E+00  2.4039E-03  4.4003E-03
             1.0958E-02

0ITERATION NO.:   52    OBJECTIVE VALUE:  -3055.89609139700        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1447
 NPARAMETR:  1.0404E+00  2.3986E+00  1.8934E-01  1.4024E-01  1.9773E+00  9.5890E-01  6.0764E-01  1.0000E-02  2.4957E+00  1.8645E+00
             2.1558E+00
 PARAMETER:  1.3961E-01  9.7489E-01 -1.5642E+00 -1.8644E+00  7.8172E-01  5.8028E-02 -3.9817E-01 -4.6575E+00  1.0146E+00  7.2297E-01
             8.6817E-01
 GRADIENT:  -1.6664E-03 -1.5546E-04 -5.5574E-04  8.0800E-04  2.1053E-03 -6.6233E-04  4.6991E-04  0.0000E+00  1.3575E-04  6.7296E-04
             5.0621E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1447
 NO. OF SIG. DIGITS IN FINAL EST.:  3.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0248E-03 -1.9230E-02 -2.5080E-05  2.9148E-02 -1.5433E-02
 SE:             2.9556E-02  2.7004E-02  2.5135E-05  1.5248E-02  2.7168E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7234E-01  4.7641E-01  3.1837E-01  5.5939E-02  5.6999E-01

 ETASHRINKSD(%)  9.8279E-01  9.5315E+00  9.9916E+01  4.8916E+01  8.9846E+00
 ETASHRINKVR(%)  1.9559E+00  1.8155E+01  1.0000E+02  7.3904E+01  1.7162E+01
 EBVSHRINKSD(%)  1.1632E+00  9.8552E+00  9.9907E+01  5.5639E+01  5.8410E+00
 EBVSHRINKVR(%)  2.3129E+00  1.8739E+01  1.0000E+02  8.0321E+01  1.1341E+01
 RELATIVEINF(%)  9.7648E+01  2.1597E+01  3.7953E-05  4.1957E+00  6.3327E+01
 EPSSHRINKSD(%)  1.7053E+01
 EPSSHRINKVR(%)  3.1198E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          899
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1652.2514827020016     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3055.8960913970013     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1403.6446086949998     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    34.80
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.75
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3055.896       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  2.40E+00  1.89E-01  1.40E-01  1.98E+00  9.59E-01  6.08E-01  1.00E-02  2.50E+00  1.86E+00  2.16E+00
 


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
+        1.09E+03
 
 TH 2
+       -1.18E+01  4.36E+02
 
 TH 3
+        2.04E+00  6.73E+01  3.70E+02
 
 TH 4
+       -2.25E+01  4.72E+02 -6.00E+02  2.55E+03
 
 TH 5
+       -1.55E+00 -1.69E+01 -1.72E+01  8.50E+01  7.30E+01
 
 TH 6
+        5.56E+00 -3.60E+00  4.29E-01 -8.84E+00 -7.91E-01  2.04E+02
 
 TH 7
+        4.88E+00  2.83E+00 -5.21E+00 -5.91E+01  2.98E-01  1.48E+00  3.81E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -2.23E-01 -5.45E+00 -2.32E+01  8.47E+01 -4.33E-01 -1.23E-01  7.33E-01  0.00E+00  6.55E+00
 
 TH10
+        3.30E-01 -3.32E+00  1.15E+01  1.83E+01 -5.01E+00 -4.85E-02  2.30E+00  0.00E+00  1.14E+00  4.20E+01
 
 TH11
+       -1.39E+01 -1.51E+01  2.69E+00 -2.16E+01  1.16E+00  2.87E+00  1.32E+01  0.00E+00  1.25E+00  2.41E+00  2.53E+02
 
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
 #CPUT: Total CPU Time in Seconds,       47.651
Stop Time:
Sat Sep 25 01:00:32 CDT 2021
