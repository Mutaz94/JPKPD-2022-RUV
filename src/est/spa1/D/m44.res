Thu Sep 30 03:05:49 CDT 2021
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
$DATA ../../../../data/spa1/D/dat44.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m44.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   15470.9637340281        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.2695E+02  3.6631E+02 -1.8542E+01  4.4010E+02  1.1775E+02 -1.7828E+03 -9.8121E+02 -6.8385E+01 -1.2242E+03 -2.7691E+02
            -3.0225E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -525.967640107133        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.6765E-01  1.0624E+00  9.3286E-01  1.7779E+00  1.3609E+00  2.7039E+00  1.7226E+00  9.5628E-01  1.9376E+00  1.0482E+00
             1.3329E+01
 PARAMETER:  6.7117E-02  1.6053E-01  3.0505E-02  6.7543E-01  4.0812E-01  1.0947E+00  6.4383E-01  5.5295E-02  7.6143E-01  1.4704E-01
             2.6900E+00
 GRADIENT:  -5.8287E+01  2.8740E+01 -2.5863E+01  7.7166E+01  1.2234E+00  9.3761E+01  2.2308E+00  5.1504E+00  6.1099E+00  2.1936E+00
             2.3423E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -554.241096759939        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0845E+00  1.0675E+00  1.6869E+00  1.9766E+00  5.5927E+00  1.8720E+00  6.2725E+00  6.8873E-01  1.8936E+00  1.8292E+00
             1.2964E+01
 PARAMETER:  1.8116E-01  1.6535E-01  6.2291E-01  7.8140E-01  1.8215E+00  7.2700E-01  1.9362E+00 -2.7290E-01  7.3850E-01  7.0386E-01
             2.6622E+00
 GRADIENT:  -8.9175E+00  3.2067E+01  1.7419E+00  6.6855E+01 -6.0277E-01  1.8923E+01  2.6263E+01  1.2976E-01  3.0970E+01  6.5611E-02
             2.1116E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -599.915257744387        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0577E+00  4.5304E-01  1.7751E+00  1.4974E+00  2.9959E+00  1.8044E+00  2.5556E+00  9.1652E-01  1.3380E+00  4.1623E+00
             1.1050E+01
 PARAMETER:  1.5610E-01 -6.9177E-01  6.7386E-01  5.0376E-01  1.1972E+00  6.9023E-01  1.0383E+00  1.2830E-02  3.9118E-01  1.5261E+00
             2.5024E+00
 GRADIENT:   8.6008E+00 -3.3688E+00  1.5414E+01 -5.1577E+01 -9.9274E+00  1.8202E+01  3.1469E+00 -4.4959E-01  2.2631E+00 -3.0282E-02
             1.4596E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -646.421197715813        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  7.7227E-01  1.9904E-01  2.8680E-01  1.1406E+00  2.4276E+01  1.4994E+00  1.4062E+00  2.1111E-01  9.9597E-01  6.4951E+00
             7.8360E+00
 PARAMETER: -1.5842E-01 -1.5142E+00 -1.1490E+00  2.3158E-01  3.2895E+00  5.0507E-01  4.4093E-01 -1.4554E+00  9.5957E-02  1.9710E+00
             2.1587E+00
 GRADIENT:  -3.3311E+01  4.7783E+01  1.6154E+01  5.1946E+01 -1.5779E-01 -7.7131E+01  3.4524E+00 -1.0923E+00 -4.1149E+01 -3.9939E+00
            -1.6356E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -669.218103038435        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  6.0092E-01  5.7933E-02  9.0841E-02  7.5209E-01  1.2619E+02  1.7179E+00  3.6821E-01  1.0516E+00  9.8539E-01  1.3952E+01
             6.8190E+00
 PARAMETER: -4.0929E-01 -2.7485E+00 -2.2986E+00 -1.8490E-01  4.9378E+00  6.4112E-01 -8.9911E-01  1.5032E-01  8.5287E-02  2.7356E+00
             2.0197E+00
 GRADIENT:   4.7514E+01  2.2069E+01 -6.0221E+01  1.8750E+02  4.9891E-01  1.8698E+01  3.6355E-01 -7.3406E-01  1.0975E+01  1.6440E+00
            -2.7734E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -686.994400561652        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      488
 NPARAMETR:  5.9566E-01  4.7423E-02  7.9678E-02  6.8834E-01  1.6834E+02  1.7922E+00  2.6755E-01  1.1353E+00  9.1486E-01  1.8164E+01
             7.1594E+00
 PARAMETER: -4.1809E-01 -2.9486E+00 -2.4298E+00 -2.7348E-01  5.2260E+00  6.8347E-01 -1.2184E+00  2.2691E-01  1.1014E-02  2.9995E+00
             2.0684E+00
 GRADIENT:   1.5653E+01  7.1884E+00 -9.4475E+01  1.5652E+02  2.6023E-01  2.1808E+01  1.0308E-01  1.9352E+00  1.5057E+01  6.8840E-01
            -2.2227E+02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -723.276114053563        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      664
 NPARAMETR:  5.7778E-01  3.0693E-02  7.3225E-02  5.8968E-01  1.7608E+02  1.5371E+00  1.6876E-01  1.1759E+00  5.6381E-01  1.9416E+01
             8.7362E+00
 PARAMETER: -4.4856E-01 -3.3837E+00 -2.5142E+00 -4.2818E-01  5.2709E+00  5.2990E-01 -1.6793E+00  2.6207E-01 -4.7304E-01  3.0661E+00
             2.2675E+00
 GRADIENT:  -4.7783E+00 -4.5425E-01 -3.7141E+00  1.1336E+01  2.4081E-01 -3.0345E+00  1.3459E-03 -1.7378E+00  6.2757E+00  8.4343E-01
             1.7316E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -727.026014562650        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      840
 NPARAMETR:  6.0532E-01  3.5802E-02  8.5802E-02  6.3607E-01  3.5240E+01  1.5708E+00  2.6931E-01  1.5273E+00  2.1987E-01  5.8342E+00
             8.5995E+00
 PARAMETER: -4.0200E-01 -3.2298E+00 -2.3557E+00 -3.5245E-01  3.6622E+00  5.5156E-01 -1.2119E+00  5.2350E-01 -1.4147E+00  1.8637E+00
             2.2517E+00
 GRADIENT:   2.7440E+00  2.0330E-01  9.4230E-02  2.1735E+00  1.4359E+00  2.4254E+00  3.4240E-02  1.8884E-02  8.6364E-01 -1.5066E+00
             6.5380E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -729.925679452428        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1017
 NPARAMETR:  5.7520E-01  2.7677E-02  7.2913E-02  5.8026E-01  8.9726E+00  1.5485E+00  2.7653E-01  1.4783E+00  7.1952E-02  2.1473E+00
             8.5950E+00
 PARAMETER: -4.5304E-01 -3.4871E+00 -2.5185E+00 -4.4427E-01  2.2942E+00  5.3732E-01 -1.1854E+00  4.9086E-01 -2.5318E+00  8.6420E-01
             2.2512E+00
 GRADIENT:  -3.8306E+00  2.2118E+00 -3.3754E+00  1.4394E+00 -3.4506E+00 -1.7345E+00  3.0486E-02  2.0671E+00  8.9284E-02  2.7047E-01
             3.8248E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -730.292033677882        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:     1157
 NPARAMETR:  5.7631E-01  2.8013E-02  7.3035E-02  5.7814E-01  1.1620E+01  1.5554E+00  1.4895E-01  1.4564E+00  5.1028E-02  2.5334E+00
             8.6187E+00
 PARAMETER: -4.5111E-01 -3.4751E+00 -2.5168E+00 -4.4793E-01  2.5528E+00  5.4170E-01 -1.8041E+00  4.7600E-01 -2.8754E+00  1.0296E+00
             2.2539E+00
 GRADIENT:   3.8434E+01  9.6072E-01  3.8248E+01  9.6701E+00 -3.8868E-01  1.0616E+01  6.9715E-03  1.1248E+00  6.7297E-02  6.0406E-02
             2.9212E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -730.325442137241        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:     1276
 NPARAMETR:  5.7683E-01  2.7901E-02  7.2973E-02  5.7802E-01  1.1648E+01  1.5494E+00  8.6579E-02  1.4600E+00  2.9748E-02  2.5318E+00
             8.5946E+00
 PARAMETER: -4.5022E-01 -3.4791E+00 -2.5177E+00 -4.4814E-01  2.5551E+00  5.3786E-01 -2.3467E+00  4.7846E-01 -3.4150E+00  1.0289E+00
             2.2511E+00
 GRADIENT:   6.7314E-01  4.6591E-01  2.0717E-01 -3.0845E+00 -4.3211E-01 -1.2280E+00  2.2616E-03 -2.8810E-01  1.4205E-02  5.7068E-02
             3.5421E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -730.358338944628        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1453
 NPARAMETR:  5.7558E-01  2.6995E-02  7.2945E-02  5.7839E-01  1.1852E+01  1.5575E+00  2.0041E-02  1.4700E+00  1.0000E-02  2.5306E+00
             8.5617E+00
 PARAMETER: -4.5237E-01 -3.5121E+00 -2.5181E+00 -4.4751E-01  2.5725E+00  5.4306E-01 -3.8100E+00  4.8526E-01 -4.9328E+00  1.0285E+00
             2.2473E+00
 GRADIENT:  -8.2968E-02  5.7359E-02 -1.1399E+00 -4.0942E-01 -2.7643E-01  1.5674E-01  8.8366E-05  5.0622E-01  0.0000E+00  4.1764E-02
             6.1481E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -730.407740549889        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1615
 NPARAMETR:  5.7536E-01  2.9382E-02  7.2960E-02  5.7855E-01  2.4245E+01  1.5566E+00  1.0000E-02  1.4649E+00  1.0000E-02  1.5956E-02
             8.5407E+00
 PARAMETER: -4.5276E-01 -3.4274E+00 -2.5178E+00 -4.4722E-01  3.2882E+00  5.4252E-01 -5.8482E+00  4.8178E-01 -5.5073E+00 -4.0379E+00
             2.2448E+00
 GRADIENT:   6.3490E-01  1.3865E-02 -3.2276E+00  4.0447E+00 -4.0887E-03  6.0404E-03  0.0000E+00 -8.5857E-01  0.0000E+00  1.7635E-07
            -3.2060E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -730.422719402282        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1791
 NPARAMETR:  5.7332E-01  3.0224E-02  7.2881E-02  5.7614E-01  5.7537E+01  1.5532E+00  1.0000E-02  1.4726E+00  1.0000E-02  1.0000E-02
             8.5542E+00
 PARAMETER: -4.5632E-01 -3.3991E+00 -2.5189E+00 -4.5141E-01  4.1524E+00  5.4031E-01 -7.1748E+00  4.8703E-01 -5.5073E+00 -8.5273E+00
             2.2464E+00
 GRADIENT:  -7.8959E-01  6.5021E-02 -3.8406E-01 -8.3545E-01 -6.5530E-03 -2.1379E-01  0.0000E+00 -6.0231E-02  0.0000E+00  0.0000E+00
            -3.6117E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -730.422890070993        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1967
 NPARAMETR:  5.7355E-01  2.9966E-02  7.2912E-02  5.7648E-01  5.6830E+01  1.5535E+00  1.0000E-02  1.4729E+00  1.0000E-02  1.0000E-02
             8.5560E+00
 PARAMETER: -4.5590E-01 -3.4077E+00 -2.5185E+00 -4.5081E-01  4.1401E+00  5.4052E-01 -7.0029E+00  4.8725E-01 -5.5073E+00 -7.9464E+00
             2.2466E+00
 GRADIENT:  -5.7520E-01 -3.8807E-02 -6.9013E-01 -3.0728E-01 -4.5302E-03 -1.6248E-01  0.0000E+00 -9.0349E-03  0.0000E+00  0.0000E+00
            -1.8961E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -730.422967309957        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2142
 NPARAMETR:  5.7390E-01  2.9856E-02  7.3036E-02  5.7716E-01  6.2972E+01  1.5537E+00  1.0000E-02  1.4741E+00  1.0000E-02  1.0000E-02
             8.5573E+00
 PARAMETER: -4.5530E-01 -3.4114E+00 -2.5168E+00 -4.4964E-01  4.2427E+00  5.4067E-01 -6.7057E+00  4.8804E-01 -5.5073E+00 -6.9425E+00
             2.2468E+00
 GRADIENT:  -4.2591E-01 -1.0822E-01 -8.7443E-01  7.3261E-02 -2.7281E-03 -1.2643E-01  0.0000E+00  2.2936E-02  0.0000E+00  0.0000E+00
            -9.6607E-02

0ITERATION NO.:   85    OBJECTIVE VALUE:  -730.423100331481        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2317
 NPARAMETR:  5.7504E-01  2.9908E-02  7.3531E-02  5.7948E-01  1.0432E+02  1.5541E+00  1.0000E-02  1.4786E+00  1.0000E-02  2.4852E-02
             8.5604E+00
 PARAMETER: -4.5332E-01 -3.4096E+00 -2.5101E+00 -4.4563E-01  4.7475E+00  5.4089E-01 -5.7144E+00  4.9113E-01 -5.5073E+00 -3.5948E+00
             2.2471E+00
 GRADIENT:  -1.8048E-01 -1.9480E-01 -1.0265E+00  5.7557E-01 -5.4337E-04 -6.2127E-02  0.0000E+00  7.3883E-02  0.0000E+00  1.1616E-08
             8.2542E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -730.423374432695        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2492
 NPARAMETR:  5.7663E-01  3.0353E-02  7.4275E-02  5.8272E-01  3.2936E+02  1.5543E+00  1.0000E-02  1.4854E+00  1.0000E-02  1.2323E+00
             8.5632E+00
 PARAMETER: -4.5056E-01 -3.3949E+00 -2.5000E+00 -4.4005E-01  5.8971E+00  5.4101E-01 -4.5577E+00  4.9567E-01 -5.5073E+00  3.0892E-01
             2.2475E+00
 GRADIENT:  -8.4373E-02 -1.8607E-01 -8.2712E-01  5.8096E-01 -1.3846E-04 -2.6680E-02  0.0000E+00  9.7110E-02  0.0000E+00  2.6357E-06
             2.1417E-01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -730.426588634915        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     2682
 NPARAMETR:  5.7706E-01  3.1289E-02  7.4316E-02  5.8276E-01  7.7504E+02  1.5550E+00  1.0000E-02  1.4850E+00  1.0000E-02  1.6107E+00
             8.5603E+00
 PARAMETER: -4.4981E-01 -3.3645E+00 -2.4994E+00 -4.3997E-01  6.7529E+00  5.4146E-01 -4.6193E+00  4.9544E-01 -5.5073E+00  5.7669E-01
             2.2471E+00
 GRADIENT:   2.5591E-01  1.2283E-01 -7.8745E-01  1.4756E-01 -3.4316E-04  1.6101E-01  0.0000E+00 -1.2056E-01  0.0000E+00  1.0586E-06
            -4.2466E-01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -730.427124946898        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     2874
 NPARAMETR:  5.7703E-01  3.1075E-02  7.4296E-02  5.8258E-01  2.6685E+04  1.5547E+00  1.0000E-02  1.4859E+00  1.0000E-02  1.6221E+00
             8.5620E+00
 PARAMETER: -4.4986E-01 -3.3714E+00 -2.4997E+00 -4.4029E-01  1.0292E+01  5.4130E-01 -4.6193E+00  4.9603E-01 -5.5073E+00  5.8372E-01
             2.2473E+00
 GRADIENT:   4.2506E-01  2.1143E-02 -6.4844E-01 -1.7119E-01 -7.4166E-06  1.4660E-01  0.0000E+00  2.2779E-02  0.0000E+00  8.1800E-10
            -9.5403E-02

0ITERATION NO.:  105    OBJECTIVE VALUE:  -730.427337230170        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     3072             RESET HESSIAN, TYPE I
 NPARAMETR:  5.7692E-01  3.1066E-02  7.4228E-02  5.8244E-01  7.3709E+08  1.5546E+00  1.0000E-02  1.4855E+00  1.0000E-02  1.6180E+00
             8.5614E+00
 PARAMETER: -4.5006E-01 -3.3716E+00 -2.5006E+00 -4.4053E-01  2.0518E+01  5.4120E-01 -4.6193E+00  4.9577E-01 -5.5073E+00  5.8120E-01
             2.2473E+00
 GRADIENT:   3.9758E+01  3.7471E-01  3.5416E+01  1.4289E+01 -2.2385E-10  1.0611E+01  0.0000E+00  1.5425E+00  0.0000E+00  0.0000E+00
             2.2234E+01

0ITERATION NO.:  107    OBJECTIVE VALUE:  -730.427337230170        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     3129
 NPARAMETR:  5.7692E-01  3.1066E-02  7.4228E-02  5.8244E-01  7.3709E+08  1.5546E+00  1.0000E-02  1.4855E+00  1.0000E-02  1.6180E+00
             8.5614E+00
 PARAMETER: -4.5006E-01 -3.3716E+00 -2.5006E+00 -4.4053E-01  2.0518E+01  5.4120E-01 -4.6193E+00  4.9577E-01 -5.5073E+00  5.8120E-01
             2.2473E+00
 GRADIENT:   3.9711E-01  2.4640E-02 -9.8428E-01  3.9569E-01 -2.4933E-10  9.5505E-02  0.0000E+00  1.9162E-02  0.0000E+00  0.0000E+00
            -2.4043E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     3129
 NO. OF SIG. DIGITS IN FINAL EST.:  3.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.0236E-03 -3.6291E-05  3.7479E-03 -4.5274E-04 -1.2573E-13
 SE:             2.8663E-02  6.5004E-06  2.3311E-02  2.8034E-04  3.4326E-12
 N:                     100         100         100         100         100

 P VAL.:         7.7953E-01  2.3709E-08  8.7227E-01  1.0631E-01  9.7078E-01

 ETASHRINKSD(%)  3.9746E+00  9.9978E+01  2.1906E+01  9.9061E+01  1.0000E+02
 ETASHRINKVR(%)  7.7912E+00  1.0000E+02  3.9013E+01  9.9991E+01  1.0000E+02
 EBVSHRINKSD(%)  3.5211E+00  9.9960E+01  2.1105E+01  9.9057E+01  1.0000E+02
 EBVSHRINKVR(%)  6.9182E+00  1.0000E+02  3.7755E+01  9.9991E+01  1.0000E+02
 RELATIVEINF(%)  2.9378E+01  1.4975E-05  3.1561E+00  4.2909E-04  0.0000E+00
 EPSSHRINKSD(%)  1.0183E+01
 EPSSHRINKVR(%)  1.9329E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -730.42733723016966     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       188.51119597450304     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    54.28
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.92
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -730.427       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         5.77E-01  3.11E-02  7.42E-02  5.82E-01  7.37E+08  1.55E+00  1.00E-02  1.49E+00  1.00E-02  1.62E+00  8.56E+00
 


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
+        1.27E+03
 
 TH 2
+       -3.45E+02  6.26E+03
 
 TH 3
+       -7.77E+02 -7.19E+02  5.62E+04
 
 TH 4
+       -4.30E+02 -8.51E+00 -1.14E+04  2.81E+03
 
 TH 5
+       -2.47E-14  3.08E-14  1.03E-13  3.88E-14  1.52E-24
 
 TH 6
+        1.58E+01  1.04E+01  1.08E+02 -4.48E+01 -8.25E-15  6.78E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        2.47E+00 -4.55E+01 -1.57E+02 -2.38E+01 -1.91E-14  2.48E+00  0.00E+00  3.35E+01
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+       -1.29E-03 -1.22E-03 -5.71E-04 -1.23E-03  8.54E-15  8.98E-05  0.00E+00  6.16E-04  0.00E+00  3.92E-04
 
 TH11
+       -1.63E+01 -1.30E+01  1.08E+02 -2.37E+01  2.31E-17  2.42E+00  0.00E+00  2.94E+00  0.00E+00 -2.02E-06  7.30E+00
 
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
 #CPUT: Total CPU Time in Seconds,       63.259
Stop Time:
Thu Sep 30 03:06:54 CDT 2021
