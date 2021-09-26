Sat Sep 25 09:28:36 CDT 2021
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
$DATA ../../../../data/spa/A3/dat75.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m75.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -264.281769160912        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.0218E+01 -2.0592E-01  3.7118E+01 -9.9231E+01  1.7264E+02 -2.3365E-01 -6.4962E+01 -1.3410E+01 -1.8838E+02 -8.8378E+01
            -2.4330E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1111.70230559526        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.1073E+00  9.3332E-01  8.8867E-01  1.3310E+00  8.1372E-01  8.4547E-01  1.1225E+00  9.5237E-01  1.7080E+00  9.5554E-01
             7.4302E+00
 PARAMETER:  2.0196E-01  3.0990E-02 -1.8034E-02  3.8591E-01 -1.0614E-01 -6.7868E-02  2.1553E-01  5.1196E-02  6.3534E-01  5.4517E-02
             2.1056E+00
 GRADIENT:  -1.2599E+01 -8.5685E+00 -1.1554E+01  2.5901E+01 -9.2833E+00 -7.7511E+00  1.0720E+01  7.2121E+00  5.6930E+01  1.9521E+01
             3.3732E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1221.59775697011        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0231E+00  5.1249E-01  2.0072E-01  1.3437E+00  2.6425E-01  1.0162E+00  7.6638E-01  9.2057E-02  1.5850E+00  1.7783E-01
             5.0096E+00
 PARAMETER:  1.2280E-01 -5.6846E-01 -1.5058E+00  3.9542E-01 -1.2309E+00  1.1607E-01 -1.6607E-01 -2.2853E+00  5.6056E-01 -1.6269E+00
             1.7113E+00
 GRADIENT:  -7.6360E+01  7.4323E+01  8.5328E+00  1.3522E+02 -1.1326E+02 -9.5909E-02  4.7013E+00  1.5819E-01  1.2220E+01  1.3054E+00
             2.1703E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1301.72205213497        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.9996E-01  4.6463E-01  2.6728E-01  1.2056E+00  3.0629E-01  1.0008E+00  1.6727E+00  2.4121E-02  1.1954E+00  2.1052E-01
             3.1315E+00
 PARAMETER:  9.9957E-02 -6.6651E-01 -1.2194E+00  2.8694E-01 -1.0832E+00  1.0077E-01  6.1441E-01 -3.6247E+00  2.7846E-01 -1.4582E+00
             1.2415E+00
 GRADIENT:  -4.8544E+01  3.0131E+01 -6.9480E+00  4.3181E+01 -1.6163E+01 -4.7578E-01  1.2638E+01  5.7251E-03 -2.3748E+00  3.6680E-01
             1.7612E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1306.20461446133        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  1.0186E+00  2.9636E-01  2.4824E-01  1.2006E+00  2.5944E-01  1.0211E+00  1.6888E+00  1.0000E-02  1.2862E+00  4.1886E-01
             3.1424E+00
 PARAMETER:  1.1844E-01 -1.1162E+00 -1.2934E+00  2.8283E-01 -1.2492E+00  1.2087E-01  6.2402E-01 -6.5714E+00  3.5169E-01 -7.7022E-01
             1.2450E+00
 GRADIENT:  -1.1036E+01  1.2386E+01  7.4831E+00  4.5726E+00 -1.7497E+01  7.9561E+00  1.6096E+00  0.0000E+00  1.0698E+01  1.1234E+00
             2.4918E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1307.96367107153        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  1.0223E+00  1.9014E-01  2.3015E-01  1.1812E+00  2.3315E-01  9.7734E-01  2.1436E+00  1.0000E-02  1.2298E+00  5.1514E-01
             2.8915E+00
 PARAMETER:  1.2207E-01 -1.5600E+00 -1.3690E+00  2.6651E-01 -1.3561E+00  7.7081E-02  8.6249E-01 -9.2100E+00  3.0689E-01 -5.6331E-01
             1.1618E+00
 GRADIENT:   7.1011E+00  2.9103E+00  1.0908E+01 -7.8259E+00 -1.6736E+01 -6.8267E+00 -4.0736E+00  0.0000E+00 -7.0578E+00 -9.9182E-01
            -1.6396E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1310.64933549352        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      446
 NPARAMETR:  1.0158E+00  1.5244E-01  2.5343E-01  1.2408E+00  2.4656E-01  9.9088E-01  2.9312E+00  1.0000E-02  1.2257E+00  4.8808E-01
             2.9564E+00
 PARAMETER:  1.1572E-01 -1.7810E+00 -1.2727E+00  3.1579E-01 -1.3001E+00  9.0841E-02  1.1754E+00 -1.0550E+01  3.0348E-01 -6.1728E-01
             1.1840E+00
 GRADIENT:   3.3243E-01  2.6893E+00  6.0190E+00  5.2604E+00 -9.9974E+00 -1.0118E-01 -1.4199E+00  0.0000E+00  2.1986E+00  9.9794E-01
            -3.2775E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1311.69305333426        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      516
 NPARAMETR:  1.0112E+00  1.0153E-01  2.6145E-01  1.2608E+00  2.5046E-01  9.8894E-01  4.0844E+00  1.0000E-02  1.1762E+00  4.5493E-01
             2.9876E+00
 PARAMETER:  1.1116E-01 -2.1874E+00 -1.2415E+00  3.3176E-01 -1.2844E+00  8.8874E-02  1.5072E+00 -1.3140E+01  2.6227E-01 -6.8760E-01
             1.1945E+00
 GRADIENT:   1.8555E-01 -4.6402E-01  5.2052E-01  5.6980E-02 -6.5655E-01 -7.0568E-02 -8.9566E-01  0.0000E+00 -5.5593E-02  6.1175E-01
             4.9467E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1313.32919302554        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      593
 NPARAMETR:  1.0098E+00  9.1943E-02  2.7248E-01  1.2804E+00  2.5776E-01  9.8519E-01  4.3579E+00  1.0000E-02  1.1545E+00  4.3013E-01
             3.0048E+00
 PARAMETER:  1.0979E-01 -2.2866E+00 -1.2002E+00  3.4721E-01 -1.2557E+00  8.5076E-02  1.5720E+00 -1.3789E+01  2.4366E-01 -7.4368E-01
             1.2002E+00
 GRADIENT:  -1.5158E+00 -5.7125E+00  1.7606E+01  8.0058E+00 -1.9793E+01 -1.1952E-01 -1.1616E+01  0.0000E+00 -2.1027E+00  2.9202E+00
             5.7519E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1316.75250979352        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      667
 NPARAMETR:  1.0074E+00  6.1946E-02  2.6386E-01  1.2813E+00  2.5324E-01  9.7632E-01  5.4423E+00  1.0000E-02  1.1779E+00  3.8506E-01
             3.0173E+00
 PARAMETER:  1.0734E-01 -2.6815E+00 -1.2323E+00  3.4788E-01 -1.2734E+00  7.6034E-02  1.7942E+00 -1.6525E+01  2.6377E-01 -8.5435E-01
             1.2044E+00
 GRADIENT:   1.0348E+00 -3.3742E+00 -1.4427E+01  8.7097E+00  2.0398E+01 -3.2332E+00 -7.7278E+00  0.0000E+00  1.0211E-01 -5.0279E-01
            -1.2489E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1317.10106813012        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      740
 NPARAMETR:  1.0061E+00  5.5512E-02  2.6432E-01  1.2803E+00  2.5190E-01  9.8783E-01  5.9444E+00  1.0000E-02  1.1671E+00  3.8638E-01
             3.0283E+00
 PARAMETER:  1.0603E-01 -2.7912E+00 -1.2306E+00  3.4709E-01 -1.2787E+00  8.7752E-02  1.8824E+00 -1.7215E+01  2.5451E-01 -8.5093E-01
             1.2080E+00
 GRADIENT:   3.1971E+00  9.7027E+00 -2.3133E+01 -8.0998E+00  2.8981E+01 -2.5907E+00  1.4729E+01  0.0000E+00 -3.4247E-02 -4.1545E+00
            -1.1950E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1317.65425642870        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      897
 NPARAMETR:  1.0077E+00  5.5077E-02  2.8996E-01  1.3153E+00  2.6859E-01  9.8828E-01  6.2020E+00  1.0000E-02  1.1452E+00  3.7978E-01
             3.0684E+00
 PARAMETER:  1.0771E-01 -2.7990E+00 -1.1380E+00  3.7404E-01 -1.2146E+00  8.8216E-02  1.9249E+00 -1.7344E+01  2.3560E-01 -8.6817E-01
             1.2211E+00
 GRADIENT:   2.2635E+01  5.0603E+01 -9.9157E+01 -7.2437E+01  1.2644E+02 -1.5284E+01  8.9264E+01  0.0000E+00  1.0508E+01 -1.4587E+01
            -5.7145E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1317.75608235176        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1076
 NPARAMETR:  1.0071E+00  5.2022E-02  3.0340E-01  1.3356E+00  2.7729E-01  9.8260E-01  6.4517E+00  1.0000E-02  1.1333E+00  3.5501E-01
             3.0792E+00
 PARAMETER:  1.0711E-01 -2.8561E+00 -1.0927E+00  3.8938E-01 -1.1827E+00  8.2442E-02  1.9643E+00 -1.7742E+01  2.2517E-01 -9.3562E-01
             1.2247E+00
 GRADIENT:   4.7536E+00  1.0043E+01 -1.9678E+01 -1.2645E+01  2.5698E+01 -3.4841E+00  1.7442E+01  0.0000E+00  2.4836E+00 -2.7915E+00
            -1.2004E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1317.77817334523        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1254
 NPARAMETR:  1.0064E+00  4.9122E-02  3.0041E-01  1.3313E+00  2.7502E-01  9.8552E-01  6.6419E+00  1.0000E-02  1.1311E+00  3.4925E-01
             3.0881E+00
 PARAMETER:  1.0634E-01 -2.9135E+00 -1.1026E+00  3.8618E-01 -1.1909E+00  8.5413E-02  1.9934E+00 -1.8103E+01  2.2317E-01 -9.5197E-01
             1.2276E+00
 GRADIENT:   1.2505E+00  5.6918E+00 -1.1611E+01 -8.1848E+00  1.3847E+01 -1.1465E+00  9.9357E+00  0.0000E+00  1.5212E+00 -2.1426E+00
            -5.6696E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1317.85045945808        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1432
 NPARAMETR:  1.0067E+00  4.4017E-02  2.9974E-01  1.3303E+00  2.7438E-01  9.8200E-01  7.0482E+00  1.0000E-02  1.1316E+00  3.7121E-01
             3.0688E+00
 PARAMETER:  1.0664E-01 -3.0232E+00 -1.1048E+00  3.8539E-01 -1.1932E+00  8.1839E-02  2.0528E+00 -1.8926E+01  2.2367E-01 -8.9099E-01
             1.2213E+00
 GRADIENT:   1.4645E+00  2.3620E+00 -8.0591E+00 -5.2639E+00  1.1252E+01 -1.2941E+00  4.3561E+00  0.0000E+00  9.9733E-01 -6.0123E-01
            -2.0539E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1317.87272383483        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1608
 NPARAMETR:  1.0066E+00  4.0507E-02  3.0214E-01  1.3355E+00  2.7518E-01  9.8355E-01  7.3943E+00  1.0000E-02  1.1278E+00  3.7281E-01
             3.0636E+00
 PARAMETER:  1.0663E-01 -3.1063E+00 -1.0969E+00  3.8928E-01 -1.1903E+00  8.3415E-02  2.1007E+00 -1.9499E+01  2.2031E-01 -8.8670E-01
             1.2196E+00
 GRADIENT:   1.9249E+00  2.5814E+00 -3.9010E+00 -3.5799E+00  5.2536E+00 -7.3919E-01  5.1080E+00  0.0000E+00  6.0215E-01 -7.0705E-01
            -3.3891E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1317.87998695310        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1770
 NPARAMETR:  1.0061E+00  3.7908E-02  3.0266E-01  1.3361E+00  2.7550E-01  9.8263E-01  7.6199E+00  1.0000E-02  1.1263E+00  3.7041E-01
             3.0645E+00
 PARAMETER:  1.0607E-01 -3.1726E+00 -1.0951E+00  3.8976E-01 -1.1892E+00  8.2474E-02  2.1308E+00 -1.9969E+01  2.1893E-01 -8.9314E-01
             1.2199E+00
 GRADIENT:  -3.9993E-02  2.7724E-02 -5.0888E-02 -3.1911E-02  5.6558E-02 -1.4750E-02  5.8514E-02  0.0000E+00 -1.5273E-02 -7.8391E-03
            -3.9143E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1770
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.9378E-04  1.7587E-02  9.3224E-05 -1.5460E-02  4.6319E-03
 SE:             2.8695E-02  9.0555E-03  2.4217E-04  2.6847E-02  1.3697E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8349E-01  5.2116E-02  7.0027E-01  5.6473E-01  7.3523E-01

 ETASHRINKSD(%)  3.8676E+00  6.9663E+01  9.9189E+01  1.0058E+01  5.4114E+01
 ETASHRINKVR(%)  7.5857E+00  9.0797E+01  9.9993E+01  1.9104E+01  7.8945E+01
 EBVSHRINKSD(%)  3.3800E+00  7.6086E+01  9.9121E+01  9.0003E+00  5.3894E+01
 EBVSHRINKVR(%)  6.6458E+00  9.4281E+01  9.9992E+01  1.7191E+01  7.8742E+01
 RELATIVEINF(%)  9.1825E+01  4.1263E+00  3.0461E-04  3.1290E+01  8.2938E-01
 EPSSHRINKSD(%)  2.8071E+01
 EPSSHRINKVR(%)  4.8262E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1317.8799869531031     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -582.72916038936489     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    20.98
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.14
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1317.880       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  3.79E-02  3.03E-01  1.34E+00  2.75E-01  9.83E-01  7.62E+00  1.00E-02  1.13E+00  3.70E-01  3.06E+00
 


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
+        1.11E+03
 
 TH 2
+        2.99E+03  1.60E+07
 
 TH 3
+       -6.97E+02 -5.82E+06  1.62E+04
 
 TH 4
+       -1.43E+02 -7.31E+03  1.06E+03  6.77E+02
 
 TH 5
+        1.38E+03  5.89E+06 -3.23E+04 -2.92E+03  2.19E+06
 
 TH 6
+       -3.62E+01 -2.31E+03  4.57E+02  7.21E+01 -8.10E+02  1.94E+02
 
 TH 7
+        5.41E+01  1.19E+05 -4.31E+04 -1.60E+02  4.37E+04 -4.02E+01  8.80E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.06E+01  1.08E+03 -1.65E+02 -5.18E+01  5.23E+02 -1.87E+01  2.01E+01  0.00E+00  1.14E+02
 
 TH10
+       -8.45E+01 -4.61E+03  6.76E+02  1.61E+02 -1.27E+03  5.32E+01 -8.20E+01  0.00E+00 -2.25E+01  1.80E+02
 
 TH11
+       -5.96E+01 -3.30E+03  5.17E+02  9.89E+01 -1.15E+03  3.52E+01 -3.82E+03  0.00E+00 -1.03E+01  9.03E+01  7.76E+01
 
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
 #CPUT: Total CPU Time in Seconds,       28.199
Stop Time:
Sat Sep 25 09:29:06 CDT 2021
