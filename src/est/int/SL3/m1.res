Wed Sep 29 03:48:18 CDT 2021
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
$DATA ../../../../data/int/SL3/dat1.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      986
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

 TOT. NO. OF OBS RECS:      886
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
 RAW OUTPUT FILE (FILE): m1.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   1046.39263229769        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3419E+02  1.4100E+02  9.1176E+01  9.4078E+01  6.0896E+01  5.3455E+01 -1.2616E+02 -1.9107E+02 -5.5527E+01 -1.8108E+01
            -9.2919E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2292.00373718885        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.1028E+00  1.1593E+00  1.0020E+00  1.0656E+00  1.0898E+00  9.0849E-01  1.0983E+00  1.0265E+00  7.8921E-01  9.4909E-01
             5.3017E+00
 PARAMETER:  1.9786E-01  2.4779E-01  1.0198E-01  1.6349E-01  1.8600E-01  4.0237E-03  1.9378E-01  1.2611E-01 -1.3672E-01  4.7751E-02
             1.7680E+00
 GRADIENT:   9.0014E+01  1.9444E+01 -2.1487E+01  5.2492E+01  1.2902E+01 -1.0474E+01  9.3425E+00  7.3378E+00  1.5695E+01  9.8272E+00
             7.5638E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2433.71821288432        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0003E+00  1.1048E+00  1.9813E+00  1.0232E+00  1.4314E+00  8.9900E-01  1.3030E+00  5.0410E+00  8.1278E-01  3.0134E-01
             3.8927E+00
 PARAMETER:  1.0030E-01  1.9962E-01  7.8377E-01  1.2289E-01  4.5868E-01 -6.4748E-03  3.6463E-01  1.7176E+00 -1.0730E-01 -1.0995E+00
             1.4591E+00
 GRADIENT:  -7.6017E+01 -3.6193E+01 -2.4391E+01 -1.4206E+01  6.7787E+01 -2.4239E+01  2.4672E+01  5.7151E+01  1.3517E+01 -3.5307E+00
             3.9573E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2510.16531154415        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.9442E-01  8.9678E-01  1.2875E+00  1.1112E+00  1.0286E+00  9.7398E-01  1.3638E+00  1.6070E+00  6.9478E-01  9.0979E-01
             3.1470E+00
 PARAMETER:  9.4404E-02 -8.9459E-03  3.5267E-01  2.0540E-01  1.2824E-01  7.3635E-02  4.1029E-01  5.7439E-01 -2.6415E-01  5.4630E-03
             1.2464E+00
 GRADIENT:  -3.1539E+01 -1.8133E+01  2.3297E+00  4.5344E+00 -1.6835E+01  6.8914E+00  1.1403E+01  1.2758E+01 -5.3446E+00  3.9633E+00
             1.4130E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2518.06130237634        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  1.0052E+00  9.7854E-01  9.7054E-01  1.0510E+00  9.7577E-01  9.4883E-01  1.1800E+00  4.8107E-01  7.8181E-01  9.6408E-01
             3.1356E+00
 PARAMETER:  1.0522E-01  7.8303E-02  7.0102E-02  1.4974E-01  7.5476E-02  4.7474E-02  2.6550E-01 -6.3174E-01 -1.4614E-01  6.3421E-02
             1.2428E+00
 GRADIENT:  -8.6867E+00 -7.1479E+00  5.9848E+00 -8.5372E+00 -1.2375E+01 -2.2449E+00 -1.0979E+00  1.1015E+00 -1.1823E+00 -1.7693E+00
            -2.1440E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2519.79191255065        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      395
 NPARAMETR:  1.0080E+00  1.0741E+00  9.8334E-01  1.0051E+00  1.0568E+00  9.5295E-01  1.0925E+00  2.6136E-01  8.0682E-01  1.0555E+00
             3.1399E+00
 PARAMETER:  1.0797E-01  1.7144E-01  8.3199E-02  1.0505E-01  1.5520E-01  5.1809E-02  1.8846E-01 -1.2419E+00 -1.1466E-01  1.5399E-01
             1.2442E+00
 GRADIENT:  -5.1796E+01 -1.3164E+01 -3.3922E+00 -8.8062E+00 -3.4452E+00 -5.7312E+00 -2.1066E+00  1.3640E-01 -4.2245E-01 -3.3013E-01
            -2.1536E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2523.61362331018        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      576
 NPARAMETR:  1.0310E+00  1.5511E+00  1.1030E+00  7.4126E-01  1.4937E+00  9.6613E-01  8.2561E-01  3.4691E-02  9.7862E-01  1.3356E+00
             3.1649E+00
 PARAMETER:  1.3057E-01  5.3894E-01  1.9800E-01 -1.9940E-01  5.0126E-01  6.5545E-02 -9.1630E-02 -3.2613E+00  7.8391E-02  3.8939E-01
             1.2521E+00
 GRADIENT:  -1.5674E+00  1.2666E+01 -2.8025E-01  1.7476E+01  7.0695E+00  2.4834E-01 -1.2123E+00 -5.0011E-03 -1.3945E+00 -2.6930E+00
            -9.0541E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2525.10253578282        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      753
 NPARAMETR:  1.0368E+00  1.9301E+00  8.6250E-01  5.0384E-01  1.7494E+00  9.6724E-01  6.9708E-01  1.0000E-02  1.3502E+00  1.4801E+00
             3.1551E+00
 PARAMETER:  1.3618E-01  7.5757E-01 -4.7919E-02 -5.8550E-01  6.5925E-01  6.6692E-02 -2.6086E-01 -5.4663E+00  4.0026E-01  4.9210E-01
             1.2490E+00
 GRADIENT:   1.0385E+01  4.1022E+01  1.8301E+00  1.8679E+01  1.4415E-01  1.9916E-01 -2.1587E+00  0.0000E+00 -1.4088E+00 -3.6030E+00
             1.9579E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2528.98000413389        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      930
 NPARAMETR:  1.0284E+00  2.3244E+00  3.6319E-01  2.2045E-01  2.0263E+00  9.6011E-01  6.2654E-01  1.0000E-02  2.4359E+00  1.6403E+00
             3.1411E+00
 PARAMETER:  1.2799E-01  9.4345E-01 -9.1282E-01 -1.4121E+00  8.0622E-01  5.9288E-02 -3.6754E-01 -1.0417E+01  9.9034E-01  5.9488E-01
             1.2446E+00
 GRADIENT:  -7.6533E+00  1.1101E+00 -3.4033E+00  7.5415E+00  2.3573E+00 -2.6775E+00  1.9689E+00  0.0000E+00  1.1309E+00 -2.2463E+00
             9.5372E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2530.18981120928        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1105
 NPARAMETR:  1.0301E+00  2.4692E+00  2.2455E-01  1.1366E-01  2.1740E+00  9.6565E-01  6.0485E-01  1.0000E-02  3.3406E+00  1.7335E+00
             3.1351E+00
 PARAMETER:  1.2966E-01  1.0039E+00 -1.3937E+00 -2.0746E+00  8.7658E-01  6.5051E-02 -4.0277E-01 -1.2824E+01  1.3061E+00  6.5015E-01
             1.2427E+00
 GRADIENT:  -3.0594E+00 -1.8235E+01 -3.3663E+00  1.3213E+00  4.4918E+00 -4.5892E-01 -3.6591E-01  0.0000E+00 -8.9060E-01  2.8215E-01
             6.2918E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2530.75699763167        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1280
 NPARAMETR:  1.0311E+00  2.5153E+00  2.4358E-01  8.2416E-02  2.2019E+00  9.6640E-01  5.9448E-01  1.0000E-02  4.2248E+00  1.7559E+00
             3.1214E+00
 PARAMETER:  1.3067E-01  1.0224E+00 -1.3123E+00 -2.3960E+00  8.8931E-01  6.5821E-02 -4.2006E-01 -1.2061E+01  1.5410E+00  6.6301E-01
             1.2383E+00
 GRADIENT:  -1.4357E-01 -6.2180E-01  3.2490E-01 -2.1254E-01 -6.3144E-01 -1.7958E-01  6.2870E-02  0.0000E+00 -1.5456E-01  2.0318E-01
             4.2999E-01

0ITERATION NO.:   51    OBJECTIVE VALUE:  -2530.75699763167        NO. OF FUNC. EVALS.:  28
 CUMULATIVE NO. OF FUNC. EVALS.:     1308
 NPARAMETR:  1.0315E+00  2.5249E+00  2.4290E-01  8.0916E-02  2.2047E+00  9.6736E-01  5.9441E-01  1.0000E-02  4.1734E+00  1.7540E+00
             3.1205E+00
 PARAMETER:  1.3067E-01  1.0224E+00 -1.3123E+00 -2.3960E+00  8.8931E-01  6.5821E-02 -4.2006E-01 -1.2061E+01  1.5410E+00  6.6301E-01
             1.2383E+00
 GRADIENT:  -4.5962E-01 -2.3577E+01  1.3552E-01  1.0742E+01 -4.1716E-01 -1.0382E-01  2.1031E-02  0.0000E+00  1.9114E+01  1.5362E-01
             4.0750E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1308
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.2738E-03 -3.4616E-02  1.7890E-06  4.3238E-02 -2.6699E-02
 SE:             2.9230E-02  2.4611E-02  1.8653E-05  1.5947E-02  2.5715E-02
 N:                     100         100         100         100         100

 P VAL.:         9.3799E-01  1.5958E-01  9.2359E-01  6.7015E-03  2.9914E-01

 ETASHRINKSD(%)  2.0758E+00  1.7549E+01  9.9938E+01  4.6575E+01  1.3852E+01
 ETASHRINKVR(%)  4.1085E+00  3.2018E+01  1.0000E+02  7.1458E+01  2.5785E+01
 EBVSHRINKSD(%)  2.2281E+00  1.4259E+01  9.9913E+01  5.7372E+01  1.1288E+01
 EBVSHRINKVR(%)  4.4066E+00  2.6485E+01  1.0000E+02  8.1829E+01  2.1302E+01
 RELATIVEINF(%)  9.5455E+01  2.9570E+01  6.5142E-05  7.2738E+00  6.3202E+01
 EPSSHRINKSD(%)  1.5504E+01
 EPSSHRINKVR(%)  2.8604E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          886
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1628.3590808386800     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2530.7569976316727     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -902.39791679299265     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    31.52
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.69
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2530.757       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  2.52E+00  2.44E-01  8.24E-02  2.20E+00  9.66E-01  5.94E-01  1.00E-02  4.22E+00  1.76E+00  3.12E+00
 


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
+        1.07E+03
 
 TH 2
+       -2.40E+01  5.09E+02
 
 TH 3
+       -3.27E+01  4.11E+01  7.14E+03
 
 TH 4
+        2.68E+01 -1.12E+03 -5.78E+02  2.13E+04
 
 TH 5
+       -4.01E+00 -1.69E+01  6.50E+01 -3.82E+01  5.98E+01
 
 TH 6
+        7.24E+00 -6.80E+00 -6.54E-01 -8.39E+00 -1.14E+00  1.95E+02
 
 TH 7
+       -2.57E+00 -2.52E+01  3.61E+02 -1.10E+02  6.75E+00 -4.96E-01  4.05E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.68E+00 -1.07E+01 -3.67E+02  6.75E+02 -2.04E+00  1.19E-01  3.06E+00  0.00E+00  2.27E+01
 
 TH10
+       -5.74E-01 -2.71E+00 -4.63E+00  4.81E+01 -4.91E+00 -3.02E-01  2.37E+00  0.00E+00  5.19E-01  3.84E+01
 
 TH11
+       -1.87E+01 -9.68E+00  6.24E+02 -1.88E+02  3.07E+00  2.61E+00  2.35E+01  0.00E+00 -3.94E+00  4.19E+00  1.75E+02
 
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
 #CPUT: Total CPU Time in Seconds,       46.321
Stop Time:
Wed Sep 29 03:49:06 CDT 2021
