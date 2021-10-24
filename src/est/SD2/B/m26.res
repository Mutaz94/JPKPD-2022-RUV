Sat Oct 23 16:58:55 CDT 2021
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
$DATA ../../../../data/SD2/B/dat26.csv ignore=@
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

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       23 OCT 2021
Days until program expires : 176
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
 NO. OF DATA RECS IN DATA SET:      800
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

 TOT. NO. OF OBS RECS:      700
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

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

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
 RAW OUTPUT FILE (FILE): m26.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2495.25312019207        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2278E+02  5.7533E+01  3.1027E+01  3.3657E+01  2.6451E+01  5.7186E+01 -2.7102E+01 -2.0493E+02 -1.2238E+01 -1.5758E+01
            -6.4878E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2893.30910899782        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0364E+00  1.0715E+00  9.0030E-01  1.0477E+00  1.1966E+00  6.7462E-01  1.0192E+00  3.7237E+00  9.3031E-01  8.9057E-01
             8.6122E-01
 PARAMETER:  1.3573E-01  1.6902E-01 -5.0237E-03  1.4664E-01  2.7947E-01 -2.9360E-01  1.1901E-01  1.4147E+00  2.7767E-02 -1.5892E-02
            -4.9404E-02
 GRADIENT:   1.0020E+03  1.8282E+02 -5.7237E+01  2.6596E+02  1.3351E+02 -9.9338E+00  1.2280E+01  1.1598E+02 -3.8949E+00 -3.4366E+01
            -1.0461E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2900.10626209551        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      199
 NPARAMETR:  1.0275E+00  1.0751E+00  9.2586E-01  1.0246E+00  1.1978E+00  8.7515E-01  1.0453E+00  4.5262E+00  9.2734E-01  9.2089E-01
             9.6168E-01
 PARAMETER:  1.2715E-01  1.7240E-01  2.2967E-02  1.2430E-01  2.8045E-01 -3.3365E-02  1.4432E-01  1.6099E+00  2.4567E-02  1.7583E-02
             6.0922E-02
 GRADIENT:   1.0122E+02 -1.4692E+00 -4.3146E+01  7.2622E+01 -1.2438E+01 -3.3139E+01 -2.6243E+01  3.7813E+01 -1.8880E+00 -1.8551E+01
             5.1048E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2953.47252089991        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  1.0063E+00  1.0794E+00  1.9153E+00  1.0081E+00  1.2866E+00  7.9637E-01  1.1774E+00  3.8888E+00  9.3357E-01  1.0040E+00
             9.4682E-01
 PARAMETER:  1.0625E-01  1.7642E-01  7.4987E-01  1.0810E-01  3.5204E-01 -1.2770E-01  2.6334E-01  1.4581E+00  3.1260E-02  1.0401E-01
             4.5354E-02
 GRADIENT:   4.3915E+01  1.7639E+01  6.5246E+00 -1.4469E+00 -3.7320E+01 -7.1290E+01  1.3365E+01  2.0980E+01  8.8266E+00 -4.3578E+00
             2.6418E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2963.06540002771        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      550
 NPARAMETR:  9.9270E-01  1.1555E+00  1.6893E+00  9.5294E-01  1.3496E+00  9.3046E-01  1.0275E+00  3.5033E+00  9.3100E-01  1.1525E+00
             9.2633E-01
 PARAMETER:  9.2669E-02  2.4455E-01  6.2433E-01  5.1799E-02  3.9977E-01  2.7921E-02  1.2714E-01  1.3537E+00  2.8502E-02  2.4191E-01
             2.3476E-02
 GRADIENT:  -1.9987E+00  3.2217E-01  4.6537E-01  4.3849E+00 -2.5476E+00  4.3118E-01  4.4973E-01  2.1294E-02 -4.0880E-01  7.5485E-01
             7.6392E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2963.25970678773        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      727
 NPARAMETR:  1.0007E+00  1.3016E+00  1.5971E+00  8.5820E-01  1.4362E+00  9.2930E-01  9.4155E-01  3.5369E+00  9.9525E-01  1.2442E+00
             9.2994E-01
 PARAMETER:  1.0068E-01  3.6363E-01  5.6816E-01 -5.2917E-02  4.6197E-01  2.6672E-02  3.9777E-02  1.3633E+00  9.5234E-02  3.1848E-01
             2.7361E-02
 GRADIENT:   1.7835E+01  1.2888E+00  1.8577E+00  6.4260E+00  5.8776E+00 -2.9581E-01 -6.7591E-01 -2.8759E+00  1.0794E+00 -7.1665E-01
            -4.3898E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2963.31256275146        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      902
 NPARAMETR:  1.0019E+00  1.4423E+00  1.4568E+00  7.7040E-01  1.5021E+00  9.3138E-01  9.1782E-01  3.6436E+00  1.0178E+00  1.3364E+00
             9.3470E-01
 PARAMETER:  1.0192E-01  4.6623E-01  4.7627E-01 -1.6085E-01  5.0686E-01  2.8912E-02  1.4248E-02  1.3930E+00  1.1761E-01  3.8998E-01
             3.2472E-02
 GRADIENT:   2.0383E+01  1.1219E+01  1.6254E+00  1.2302E+01  2.5363E+00  4.1048E-01  3.3071E-01 -1.0384E+00  7.0360E-01 -7.1262E-03
            -6.4636E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2963.35062880291        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1078
 NPARAMETR:  9.9890E-01  1.5447E+00  1.3217E+00  7.0253E-01  1.5441E+00  9.3224E-01  9.0530E-01  3.7499E+00  1.0253E+00  1.4020E+00
             9.3998E-01
 PARAMETER:  9.8901E-02  5.3481E-01  3.7892E-01 -2.5307E-01  5.3447E-01  2.9836E-02  5.0620E-04  1.4217E+00  1.2502E-01  4.3791E-01
             3.8106E-02
 GRADIENT:   1.2168E+01  1.3911E+01  6.5741E-01  1.2149E+01 -1.6454E+00  7.8597E-01  8.9006E-01  8.9649E-01  3.3609E-02  5.6732E-01
            -4.9865E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2963.55913237150        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1262             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9491E-01  1.5668E+00  1.2644E+00  6.7459E-01  1.5583E+00  9.3078E-01  8.9374E-01  3.6417E+00  1.0400E+00  1.4165E+00
             9.4447E-01
 PARAMETER:  9.4897E-02  5.4901E-01  3.3459E-01 -2.9365E-01  5.4362E-01  2.8271E-02 -1.2338E-02  1.3925E+00  1.3926E-01  4.4820E-01
             4.2866E-02
 GRADIENT:   4.9786E+02  7.2902E+02  6.0105E+00  1.1992E+02  9.1693E+01  5.5158E+01  9.6296E+00  7.4333E+01  5.5785E+00  1.6458E+01
             1.5933E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2963.57355823652        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1442             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9492E-01  1.5673E+00  1.2552E+00  6.7576E-01  1.5554E+00  9.3077E-01  8.9018E-01  3.6640E+00  1.0438E+00  1.4128E+00
             9.4434E-01
 PARAMETER:  9.4902E-02  5.4936E-01  3.2729E-01 -2.9191E-01  5.4174E-01  2.8257E-02 -1.6336E-02  1.3985E+00  1.4289E-01  4.4555E-01
             4.2732E-02
 GRADIENT:   4.9797E+02  7.3180E+02  5.1183E+00  1.2154E+02  9.1020E+01  5.5190E+01  9.1783E+00  7.6020E+01  5.7997E+00  1.5997E+01
             1.4398E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2963.57463440585        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1626
 NPARAMETR:  9.9492E-01  1.5667E+00  1.2543E+00  6.7619E-01  1.5543E+00  9.3077E-01  8.8920E-01  3.6607E+00  1.0452E+00  1.4110E+00
             9.4431E-01
 PARAMETER:  9.4905E-02  5.4898E-01  3.2660E-01 -2.9128E-01  5.4105E-01  2.8252E-02 -1.7433E-02  1.3977E+00  1.4419E-01  4.4432E-01
             4.2704E-02
 GRADIENT:   1.8580E+00 -4.1091E+00 -1.2335E-02 -4.3876E-01 -3.5076E-01  2.1431E-01  2.5945E-01 -4.7994E+00  1.6206E-01  2.0296E-01
             3.8133E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2963.57554805199        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1812             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9492E-01  1.5656E+00  1.2547E+00  6.7691E-01  1.5533E+00  9.3077E-01  8.8742E-01  3.6582E+00  1.0432E+00  1.4097E+00
             9.4425E-01
 PARAMETER:  9.4905E-02  5.4829E-01  3.2693E-01 -2.9021E-01  5.4038E-01  2.8252E-02 -1.9438E-02  1.3970E+00  1.4232E-01  4.4337E-01
             4.2641E-02
 GRADIENT:   4.9803E+02  7.2942E+02  5.0720E+00  1.2115E+02  9.0645E+01  5.5190E+01  8.5387E+00  7.5920E+01  5.4591E+00  1.5691E+01
             1.4228E+00

0ITERATION NO.:   59    OBJECTIVE VALUE:  -2963.57593961784        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:     1946
 NPARAMETR:  9.9492E-01  1.5654E+00  1.2548E+00  6.7706E-01  1.5532E+00  9.3076E-01  8.8858E-01  3.6580E+00  1.0450E+00  1.4096E+00
             9.4425E-01
 PARAMETER:  9.4906E-02  5.4814E-01  3.2701E-01 -2.8999E-01  5.4029E-01  2.8249E-02 -1.8130E-02  1.3969E+00  1.4406E-01  4.4331E-01
             4.2633E-02
 GRADIENT:   1.9262E-03  5.1219E-01 -2.9431E-02 -3.5792E-01  1.5070E-01 -5.6784E-04  5.2104E-02  3.8293E-02  4.4399E-02  5.0457E-02
             4.2591E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1946
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.0590E-04 -2.9185E-02 -4.1039E-02  3.0023E-02 -5.2098E-02
 SE:             2.9902E-02  2.3715E-02  2.0704E-02  2.1014E-02  2.3017E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8383E-01  2.1845E-01  4.7457E-02  1.5309E-01  2.3608E-02

 ETASHRINKSD(%)  1.0000E-10  2.0553E+01  3.0640E+01  2.9600E+01  2.2890E+01
 ETASHRINKVR(%)  1.0000E-10  3.6881E+01  5.1892E+01  5.0439E+01  4.0540E+01
 EBVSHRINKSD(%)  2.7183E-01  2.0290E+01  3.2145E+01  3.3463E+01  1.6878E+01
 EBVSHRINKVR(%)  5.4292E-01  3.6463E+01  5.3956E+01  5.5728E+01  3.0907E+01
 RELATIVEINF(%)  9.9446E+01  1.0157E+01  2.4258E+01  6.7044E+00  3.9567E+01
 EPSSHRINKSD(%)  2.6600E+01
 EPSSHRINKVR(%)  4.6125E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          700
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1286.5139464865417     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2963.5759396178405     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1677.0619931312988     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.11
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2963.576       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.95E-01  1.57E+00  1.25E+00  6.77E-01  1.55E+00  9.31E-01  8.89E-01  3.66E+00  1.05E+00  1.41E+00  9.44E-01
 


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
 
 Elapsed finaloutput time in seconds:     0.00
 #CPUT: Total CPU Time in Seconds,      136.691
Stop Time:
Sat Oct 23 16:59:16 CDT 2021
