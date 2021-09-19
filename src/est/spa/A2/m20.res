Sat Sep 18 09:41:57 CDT 2021
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
$DATA ../../../../data/spa/A2/dat20.csv ignore=@
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
Current Date:       18 SEP 2021
Days until program expires : 211
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m20.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -887.310015969539        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.1107E+02 -2.2691E+01  3.2905E+01 -4.8323E+01  7.7907E+01  7.3379E+00 -9.2210E+00 -4.1167E+01 -2.8934E+01 -1.7807E+01
            -1.3907E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1342.04862760752        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.8019E-01  9.4371E-01  5.8796E-01  1.1823E+00  7.7082E-01  8.6312E-01  8.1314E-01  1.9309E+00  1.1909E+00  3.5683E-01
             1.8330E+00
 PARAMETER:  7.9993E-02  4.2065E-02 -4.3109E-01  2.6745E-01 -1.6030E-01 -4.7204E-02 -1.0686E-01  7.5798E-01  2.7475E-01 -9.3048E-01
             7.0596E-01
 GRADIENT:   2.3641E+01 -2.6324E+01 -9.6165E+01  1.6698E+02  1.6041E+02 -3.6514E+01  5.5813E+00 -4.0711E-01 -1.3132E+00 -1.4746E-02
            -1.7635E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1379.58519412505        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  9.8231E-01  1.1270E+00  5.4054E-01  9.4113E-01  6.8148E-01  9.0018E-01  7.4817E-01  2.4848E+00  1.1164E+00  2.0104E-02
             2.0127E+00
 PARAMETER:  8.2155E-02  2.1953E-01 -5.1519E-01  3.9329E-02 -2.8349E-01 -5.1610E-03 -1.9013E-01  1.0102E+00  2.1014E-01 -3.8069E+00
             7.9948E-01
 GRADIENT:   2.8887E+01  4.0018E+01  3.4874E+01  4.9522E+01 -7.6010E+01 -1.3839E+01  2.0095E+00  1.8261E+01 -3.6155E+00  1.9725E-02
            -6.5809E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1397.86900350296        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  9.8060E-01  1.4916E+00  5.7354E-01  7.0884E-01  9.4666E-01  9.2067E-01  5.4623E-01  2.3833E+00  1.3873E+00  1.1649E-02
             2.5171E+00
 PARAMETER:  8.0406E-02  4.9982E-01 -4.5592E-01 -2.4413E-01  4.5181E-02  1.7350E-02 -5.0472E-01  9.6848E-01  4.2734E-01 -4.3526E+00
             1.0231E+00
 GRADIENT:   3.4995E+00  1.3287E+00  2.2901E+00  5.4303E+00 -5.6822E+00 -4.4724E+00 -4.7789E+00 -3.3821E-02 -2.1900E+00  2.2344E-03
             1.2374E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1412.17203903481        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  9.7783E-01  1.3126E+00  1.8620E-01  6.8131E-01  5.6187E-01  9.4083E-01  7.5563E-01  1.5570E+00  1.1543E+00  5.2899E-02
             2.2970E+00
 PARAMETER:  7.7581E-02  3.7198E-01 -1.5810E+00 -2.8374E-01 -4.7648E-01  3.9009E-02 -1.8020E-01  5.4275E-01  2.4353E-01 -2.8394E+00
             9.3161E-01
 GRADIENT:   4.0216E+01  4.6350E+01  4.5011E+01 -3.9594E+01 -9.2201E+01  1.8863E+00 -1.9283E+00 -1.9943E+01  5.9209E+00  1.4971E-01
             5.1541E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1426.72963993794        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  9.5575E-01  1.5243E+00  1.5256E-01  5.7835E-01  6.6952E-01  9.3009E-01  7.1754E-01  2.1104E+00  1.3166E+00  7.3263E-02
             2.0074E+00
 PARAMETER:  5.4737E-02  5.2156E-01 -1.7802E+00 -4.4758E-01 -3.0119E-01  2.7531E-02 -2.3193E-01  8.4686E-01  3.7503E-01 -2.5137E+00
             7.9684E-01
 GRADIENT:  -6.5324E-01  1.1785E+01  6.1453E+00  4.1873E+00 -1.3325E+01 -8.1766E-01 -4.7719E+00  1.4891E+00 -4.4360E+00  1.8689E-01
             4.4151E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1433.29384169107        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      447
 NPARAMETR:  9.5629E-01  1.9939E+00  5.7367E-02  2.9405E-01  8.9714E-01  9.3774E-01  6.6041E-01  1.5262E+00  2.4245E+00  6.3331E-02
             2.1635E+00
 PARAMETER:  5.5309E-02  7.9008E-01 -2.7583E+00 -1.1240E+00 -8.5433E-03  3.5720E-02 -3.1490E-01  5.2277E-01  9.8563E-01 -2.6594E+00
             8.7173E-01
 GRADIENT:   1.3637E+00  1.2608E+01 -2.8379E+00  5.5696E+00 -2.0258E+01  2.3004E+00  8.6297E-01  3.9487E+00  7.8217E-02  6.8114E-03
             7.6423E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1433.85006502213        NO. OF FUNC. EVALS.: 107
 CUMULATIVE NO. OF FUNC. EVALS.:      554
 NPARAMETR:  9.5594E-01  2.0232E+00  5.4616E-02  2.8304E-01  9.1535E-01  9.3578E-01  6.6031E-01  1.4896E+00  2.5018E+00  6.5047E-02
             2.1646E+00
 PARAMETER:  5.4935E-02  8.0467E-01 -2.8074E+00 -1.1622E+00  1.1551E-02  3.3626E-02 -3.1504E-01  4.9854E-01  1.0170E+00 -2.6327E+00
             8.7223E-01
 GRADIENT:  -9.1802E+00  1.4510E+00 -3.5225E+00  6.0781E+00 -1.6402E+01  8.7626E-01  1.2504E+00  3.5516E+00 -1.8670E+00 -5.9886E-03
             4.6712E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1434.21596829500        NO. OF FUNC. EVALS.: 103
 CUMULATIVE NO. OF FUNC. EVALS.:      657
 NPARAMETR:  9.5686E-01  2.0223E+00  5.4694E-02  2.8287E-01  9.2133E-01  9.3307E-01  6.5580E-01  1.2508E+00  2.5004E+00  1.3156E-01
             2.1655E+00
 PARAMETER:  5.5900E-02  8.0425E-01 -2.8060E+00 -1.1628E+00  1.8063E-02  3.0728E-02 -3.2189E-01  3.2380E-01  1.0165E+00 -1.9283E+00
             8.7267E-01
 GRADIENT:   1.6970E-02  1.2228E+01  7.9544E-01  8.4942E+00 -2.4340E-02  1.7974E-02  1.7970E-01 -3.3254E-02 -1.8766E+00  1.7558E-03
             5.5808E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1434.22983337306        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      834
 NPARAMETR:  9.5999E-01  2.0223E+00  5.4694E-02  2.8287E-01  9.2161E-01  9.3439E-01  6.5607E-01  1.2543E+00  2.5005E+00  1.3904E-01
             2.1655E+00
 PARAMETER:  5.9171E-02  8.0425E-01 -2.8060E+00 -1.1628E+00  1.8365E-02  3.2135E-02 -3.2149E-01  3.2661E-01  1.0165E+00 -1.8730E+00
             8.7267E-01
 GRADIENT:   5.6045E-03 -1.0742E+01 -8.2562E-01  7.2662E+00 -1.1989E-01  2.5998E-03 -2.5367E-02 -8.6611E-04 -2.6514E+00 -1.4178E-03
             4.9290E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1434.48623255967        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1021
 NPARAMETR:  9.6163E-01  2.0518E+00  5.6156E-02  2.7181E-01  9.4054E-01  9.3449E-01  6.5457E-01  1.3389E+00  2.5148E+00  2.1078E-01
             2.1487E+00
 PARAMETER:  6.0871E-02  8.1870E-01 -2.7796E+00 -1.2026E+00  3.8697E-02  3.2241E-02 -3.2378E-01  3.9185E-01  1.0222E+00 -1.4569E+00
             8.6485E-01
 GRADIENT:  -2.6821E-01  1.0356E+01  3.3193E+00  4.1411E+00 -1.5703E-01 -3.2900E-02  1.8631E-02  5.3614E-03 -5.2572E+00 -6.9805E-02
             1.9513E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1434.49223585822        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1200
 NPARAMETR:  9.6185E-01  2.0522E+00  5.6169E-02  2.7165E-01  9.4103E-01  9.3462E-01  6.5379E-01  1.3388E+00  2.5150E+00  2.4790E-01
             2.1484E+00
 PARAMETER:  6.1099E-02  8.1890E-01 -2.7794E+00 -1.2033E+00  3.9215E-02  3.2381E-02 -3.2497E-01  3.9179E-01  1.0223E+00 -1.2947E+00
             8.6473E-01
 GRADIENT:   7.1262E-02  9.1991E+00  3.3528E+00  4.1274E+00 -2.5676E-01  1.5830E-02  9.7797E-02  3.7027E-03 -5.2963E+00  2.2221E-02
             3.7422E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1434.77837082037        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1382
 NPARAMETR:  9.6189E-01  2.0501E+00  5.2613E-02  2.6132E-01  9.4096E-01  9.3514E-01  6.5736E-01  1.2568E+00  2.5441E+00  2.9994E-01
             2.1307E+00
 PARAMETER:  6.1145E-02  8.1789E-01 -2.8448E+00 -1.2420E+00  3.9141E-02  3.2939E-02 -3.1952E-01  3.2860E-01  1.0338E+00 -1.1042E+00
             8.5646E-01
 GRADIENT:   3.2442E+00 -1.0080E+01  8.5506E-01 -4.6094E-01 -4.0612E+00  4.3803E-01  1.3906E+00  3.3043E-02 -7.4114E+00  1.9496E-01
             1.7484E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1434.82127760455        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1544
 NPARAMETR:  9.6056E-01  2.0529E+00  5.2339E-02  2.5886E-01  9.4372E-01  9.3412E-01  6.5479E-01  1.2556E+00  2.5492E+00  3.0469E-01
             2.1268E+00
 PARAMETER:  5.9761E-02  8.1927E-01 -2.8500E+00 -1.2515E+00  4.2075E-02  3.1850E-02 -3.2345E-01  3.2758E-01  1.0358E+00 -1.0885E+00
             8.5462E-01
 GRADIENT:  -2.9988E-01 -1.0933E+01  1.0216E+00 -1.2361E+00 -3.1066E+00  4.1703E-02  4.7144E-01 -7.4745E-03 -8.1468E+00  1.5245E-01
             6.1612E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1434.82186810591        NO. OF FUNC. EVALS.: 211
 CUMULATIVE NO. OF FUNC. EVALS.:     1755
 NPARAMETR:  9.6067E-01  2.0529E+00  5.2340E-02  2.5885E-01  9.4372E-01  9.3401E-01  6.5348E-01  1.2561E+00  2.5492E+00  3.0418E-01
             2.1268E+00
 PARAMETER:  5.9879E-02  8.1927E-01 -2.8500E+00 -1.2515E+00  4.2076E-02  3.1733E-02 -3.2545E-01  3.2802E-01  1.0358E+00 -1.0901E+00
             8.5463E-01
 GRADIENT:  -1.7451E-02 -1.0977E+01  1.0128E+00 -1.2487E+00 -3.1353E+00 -7.7248E-03  5.3065E-03 -5.9344E-03 -8.1733E+00  1.4282E-01
             4.8404E-01

0ITERATION NO.:   71    OBJECTIVE VALUE:  -1434.82186810591        NO. OF FUNC. EVALS.:  29
 CUMULATIVE NO. OF FUNC. EVALS.:     1784
 NPARAMETR:  9.6070E-01  2.0525E+00  5.2378E-02  2.5877E-01  9.4375E-01  9.3402E-01  6.5347E-01  1.2565E+00  2.5485E+00  3.0426E-01
             2.1273E+00
 PARAMETER:  5.9879E-02  8.1927E-01 -2.8500E+00 -1.2515E+00  4.2076E-02  3.1733E-02 -3.2545E-01  3.2802E-01  1.0358E+00 -1.0901E+00
             8.5463E-01
 GRADIENT:  -1.2228E+05  1.4914E+04 -4.2800E+03  9.7636E+03 -6.1141E+04 -8.6553E-03  4.9420E-03 -4.8324E-03  1.1788E+04 -1.1217E+04
            -1.4312E+04
 NUMSIGDIG:         3.3         3.3         3.3         3.3         3.3         3.6         4.2         2.7         3.3         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1784
 NO. OF SIG. DIGITS IN FINAL EST.:  2.7

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.8514E-04 -1.9979E-02  1.3772E-02  2.4399E-02 -1.3686E-02
 SE:             2.9340E-02  2.6085E-02  7.6791E-03  2.2473E-02  7.7284E-03
 N:                     100         100         100         100         100

 P VAL.:         9.7321E-01  4.4372E-01  7.2899E-02  2.7762E-01  7.6579E-02

 ETASHRINKSD(%)  1.7057E+00  1.2612E+01  7.4274E+01  2.4712E+01  7.4109E+01
 ETASHRINKVR(%)  3.3823E+00  2.3633E+01  9.3382E+01  4.3317E+01  9.3296E+01
 EBVSHRINKSD(%)  1.9637E+00  1.3594E+01  7.3650E+01  2.4141E+01  7.3788E+01
 EBVSHRINKVR(%)  3.8888E+00  2.5339E+01  9.3057E+01  4.2454E+01  9.3129E+01
 RELATIVEINF(%)  9.3989E+01  2.1129E+01  4.6883E+00  2.0445E+01  1.8159E+00
 EPSSHRINKSD(%)  3.5336E+01
 EPSSHRINKVR(%)  5.8185E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1434.8218681059061     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -699.67104154216793     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.49
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.11
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1434.822       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.61E-01  2.05E+00  5.23E-02  2.59E-01  9.44E-01  9.34E-01  6.53E-01  1.26E+00  2.55E+00  3.04E-01  2.13E+00
 


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
+        3.31E+08
 
 TH 2
+        2.66E+02  1.08E+06
 
 TH 3
+       -3.83E+03  5.43E+04  1.37E+08
 
 TH 4
+        1.60E+03 -6.45E+03  7.24E+04  2.91E+07
 
 TH 5
+       -5.54E+03  2.76E+02 -7.76E+03  3.85E+03  3.43E+08
 
 TH 6
+       -4.20E+03  2.33E+02 -2.70E+03  1.23E+03 -4.27E+03  2.10E+02
 
 TH 7
+       -1.74E+03  1.03E+02 -1.12E+03  5.31E+02 -1.76E+03  3.80E+00  2.71E+02
 
 TH 8
+        4.46E+04 -2.55E+03  2.86E+04 -1.32E+04  4.54E+04  1.29E+00  1.51E+00  5.47E+00
 
 TH 9
+        2.00E+02 -1.14E+03  3.43E+04 -4.16E+03  4.03E+02  1.54E+02  6.73E+01 -1.62E+03  4.38E+05
 
 TH10
+        9.60E+07 -6.58E+01  6.42E+02 -2.97E+02  1.01E+03 -1.21E+03 -4.97E+02  1.29E+04 -3.61E+01  2.78E+07
 
 TH11
+       -3.02E+02 -6.36E+02  7.16E+03  6.09E+03 -6.26E+02 -2.20E+02 -7.49E+01  2.36E+03 -3.92E+02  6.91E+01  9.26E+05
 
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
 #CPUT: Total CPU Time in Seconds,       30.660
Stop Time:
Sat Sep 18 09:42:29 CDT 2021
