Sat Sep 25 01:58:10 CDT 2021
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
$DATA ../../../../data/int/SL3/dat16.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      989
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

 TOT. NO. OF OBS RECS:      889
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
 RAW OUTPUT FILE (FILE): m16.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   543.507561059846        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -3.1624E+01 -1.3452E+02  3.7707E+01  1.7488E+02  1.9073E+02  4.5035E+01 -1.6949E+02 -1.6435E+02 -1.4865E+02 -4.9394E+01
            -8.3817E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2389.35141218736        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1531E+00  1.4636E+00  9.9485E-01  8.3920E-01  1.0778E+00  7.0947E-01  1.2952E+00  9.9317E-01  8.3602E-01  1.1824E+00
             5.1400E+00
 PARAMETER:  2.4247E-01  4.8090E-01  9.4840E-02 -7.5309E-02  1.7492E-01 -2.4323E-01  3.5864E-01  9.3146E-02 -7.9106E-02  2.6751E-01
             1.7371E+00
 GRADIENT:   1.7892E+02  1.2710E+01  6.2303E-01 -3.4389E+01 -3.2037E+01 -6.1676E+01  4.4863E+01  3.2679E+00  1.3319E+01  3.1176E+00
             7.5050E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2467.30193322989        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.1215E+00  1.6127E+00  1.0143E+01  7.9100E-01  2.0072E+00  8.6022E-01  8.2038E-01  8.2446E-01  1.3771E-01  2.4837E+00
             4.6616E+00
 PARAMETER:  2.1468E-01  5.7789E-01  2.4168E+00 -1.3446E-01  7.9672E-01 -5.0562E-02 -9.7987E-02 -9.3032E-02 -1.8826E+00  1.0097E+00
             1.6394E+00
 GRADIENT:   7.2513E+01  2.0340E+01 -4.1081E+00  5.0213E+01  2.6818E+01  4.2891E+00 -3.5027E+01  6.4367E-03 -1.3810E+00  4.8897E+01
             6.8694E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2672.10633318095        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.0529E+00  1.6094E+00  6.3364E+00  7.0535E-01  1.7704E+00  8.7519E-01  9.4561E-01  2.9020E+00  6.0978E-01  2.0404E+00
             2.9539E+00
 PARAMETER:  1.5157E-01  5.7584E-01  1.9463E+00 -2.4906E-01  6.7120E-01 -3.3313E-02  4.4071E-02  1.1654E+00 -3.9465E-01  8.1316E-01
             1.1831E+00
 GRADIENT:  -1.7490E+00 -1.6683E+01 -3.3570E+00 -5.7177E+00  3.6306E+00  4.7784E+00  6.0374E-01 -5.5880E-01 -1.3707E+00  9.1673E+00
             2.7304E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2675.09023945727        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      296
 NPARAMETR:  1.0522E+00  1.5759E+00  2.3354E+01  7.7078E-01  1.9424E+00  8.6632E-01  8.3560E-01  1.5472E+00  1.0105E+00  2.1267E+00
             2.9249E+00
 PARAMETER:  1.5089E-01  5.5485E-01  3.2508E+00 -1.6035E-01  7.6390E-01 -4.3505E-02 -7.9602E-02  5.3643E-01  1.1045E-01  8.5458E-01
             1.1733E+00
 GRADIENT:  -6.1158E+00  2.5194E+01 -2.5412E+00  1.8713E+01  1.6658E+01  1.1111E-01  4.2129E+00 -1.2796E-01  3.6977E+00  4.0082E+00
             9.0155E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2677.15293655041        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      366
 NPARAMETR:  1.0544E+00  1.3266E+00  1.9403E+02  9.2925E-01  1.9058E+00  8.5938E-01  6.5103E-01  1.0264E+00  1.0608E+00  2.1027E+00
             2.9032E+00
 PARAMETER:  1.5301E-01  3.8259E-01  5.3680E+00  2.6621E-02  7.4489E-01 -5.1547E-02 -3.2921E-01  1.2604E-01  1.5901E-01  8.4323E-01
             1.1658E+00
 GRADIENT:   1.4347E+00  7.3189E+00 -4.4615E-01  5.7494E+00  1.9277E+00 -2.6828E+00 -1.9611E+00 -2.1132E-03  8.6728E-01 -4.8651E-01
            -5.2639E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2677.33557411165        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      440
 NPARAMETR:  1.0534E+00  1.0861E+00  1.2735E+03  1.0937E+00  1.8979E+00  8.6836E-01  6.1324E-01  7.3691E-01  9.5803E-01  2.1056E+00
             2.9130E+00
 PARAMETER:  1.5207E-01  1.8256E-01  7.2495E+00  1.8960E-01  7.4072E-01 -4.1145E-02 -3.8900E-01 -2.0529E-01  5.7119E-02  8.4460E-01
             1.1692E+00
 GRADIENT:  -1.2894E+00  1.6760E+00 -9.0515E-02  2.2919E+00 -1.0660E+00  1.2224E+00 -2.3335E+00 -4.3143E-05 -1.4480E+00 -1.1437E-01
             3.7458E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2677.72033086543        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      512
 NPARAMETR:  1.0540E+00  1.0429E+00  8.4000E+02  1.1284E+00  1.9042E+00  8.6479E-01  8.6021E-01  8.9898E-01  8.6864E-01  2.1044E+00
             2.9075E+00
 PARAMETER:  1.5259E-01  1.4196E-01  6.8334E+00  2.2079E-01  7.4408E-01 -4.5269E-02 -5.0579E-02 -6.4945E-03 -4.0824E-02  8.4404E-01
             1.1673E+00
 GRADIENT:   3.5424E-01 -6.6210E-01 -1.4304E-01  1.5343E+00  6.5209E-01 -2.9070E-01 -2.7959E-01 -1.6860E-04  4.2678E-01 -1.7718E-01
            -2.1991E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2677.85070944817        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      691
 NPARAMETR:  1.0588E+00  8.8486E-01  2.5847E+03  1.2458E+00  1.9174E+00  8.6786E-01  1.0639E+00  8.1948E-01  7.7242E-01  2.1144E+00
             2.9104E+00
 PARAMETER:  1.5711E-01 -2.2329E-02  7.9574E+00  3.1981E-01  7.5096E-01 -4.1726E-02  1.6193E-01 -9.9086E-02 -1.5822E-01  8.4878E-01
             1.1683E+00
 GRADIENT:   4.0511E+00  3.1080E+00 -5.7210E-02  6.3876E+00  2.1440E-01  5.5234E-01 -3.1129E-01 -2.5418E-05 -5.4520E-01 -9.8933E-02
            -6.3495E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2677.90355011207        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      867
 NPARAMETR:  1.0574E+00  7.6554E-01  9.7660E+03  1.3250E+00  1.9171E+00  8.6672E-01  1.1624E+00  7.2081E-01  7.4597E-01  2.1193E+00
             2.9110E+00
 PARAMETER:  1.5582E-01 -1.6717E-01  9.2867E+00  3.8144E-01  7.5081E-01 -4.3035E-02  2.5048E-01 -2.2738E-01 -1.9307E-01  8.5108E-01
             1.1685E+00
 GRADIENT:   4.6035E-01 -5.4269E-01 -1.7051E-02 -1.2910E+00  1.8686E-01  1.2945E-01  1.1010E-01 -6.7857E-06  3.0655E-01  5.3815E-01
             8.5222E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2677.91480669687        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1042
 NPARAMETR:  1.0577E+00  6.7628E-01  3.7219E+04  1.3881E+00  1.9179E+00  8.6581E-01  1.1833E+00  6.2167E-01  7.3732E-01  2.1172E+00
             2.9107E+00
 PARAMETER:  1.5608E-01 -2.9115E-01  1.0625E+01  4.2791E-01  7.5121E-01 -4.4090E-02  2.6831E-01 -3.7535E-01 -2.0473E-01  8.5011E-01
             1.1684E+00
 GRADIENT:  -2.7222E-01  7.4461E-02 -4.9326E-03  3.1206E-01  2.8407E-01 -3.4883E-01 -8.5614E-02 -7.4800E-06 -6.3079E-02  2.3972E-01
             4.8451E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2677.91735993033        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     1233             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0572E+00  6.5422E-01  5.3475E+04  1.4040E+00  1.9164E+00  8.6663E-01  1.2134E+00  5.9835E-01  7.3072E-01  2.1153E+00
             2.9102E+00
 PARAMETER:  1.5561E-01 -3.2432E-01  1.0987E+01  4.3934E-01  7.5044E-01 -4.3140E-02  2.9344E-01 -4.1358E-01 -2.1372E-01  8.4921E-01
             1.1682E+00
 GRADIENT:   4.8638E+00  1.5765E+00 -3.4737E-03  1.2545E+01  2.4910E+00  3.3723E-01 -5.6232E-03  1.5429E-04  7.2705E-02  1.0954E+00
             2.1734E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2677.91772919461        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1414
 NPARAMETR:  1.0573E+00  6.5365E-01  5.5678E+04  1.4039E+00  1.9164E+00  8.6654E-01  1.2185E+00  5.9871E-01  7.3075E-01  2.1153E+00
             2.9102E+00
 PARAMETER:  1.5573E-01 -3.2519E-01  1.1027E+01  4.3926E-01  7.5047E-01 -4.3242E-02  2.9763E-01 -4.1297E-01 -2.1369E-01  8.4920E-01
             1.1682E+00
 GRADIENT:  -3.8037E+00  2.9062E-01 -3.3599E-03  9.8178E-01 -8.0262E-02  1.8200E-01 -2.2228E-02 -5.6748E-06  9.6075E-02  3.0262E-02
             1.8644E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2677.91992948905        NO. OF FUNC. EVALS.: 153
 CUMULATIVE NO. OF FUNC. EVALS.:     1567             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0571E+00  6.5344E-01  2.1284E+05  1.4042E+00  1.9164E+00  8.6549E-01  1.2234E+00  5.9395E-01  7.2821E-01  2.1163E+00
             2.9111E+00
 PARAMETER:  1.5552E-01 -3.2551E-01  1.2368E+01  4.3944E-01  7.5043E-01 -4.4456E-02  3.0167E-01 -4.2096E-01 -2.1717E-01  8.4968E-01
             1.1685E+00
 GRADIENT:  -4.4314E+02  4.4708E+01 -8.4519E-04  1.4407E+02 -1.0476E+01  5.9464E+00 -1.9001E+00  4.4106E-04 -6.2987E+00  6.2683E+00
             1.5794E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2677.92033346322        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1743             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0574E+00  6.5332E-01  2.1768E+05  1.4041E+00  1.9160E+00  8.6588E-01  1.2234E+00  5.9395E-01  7.2980E-01  2.1153E+00
             2.9101E+00
 PARAMETER:  1.5580E-01 -3.2569E-01  1.2391E+01  4.3938E-01  7.5024E-01 -4.4010E-02  3.0167E-01 -4.2095E-01 -2.1498E-01  8.4920E-01
             1.1682E+00
 GRADIENT:  -3.6039E+02  3.8610E+01 -8.3823E-04  1.2622E+02 -6.5065E+00 -1.4246E+01 -4.6674E-01  8.9964E-04 -1.0880E+00  2.9707E+00
             1.1745E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -2677.92050264409        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     1934
 NPARAMETR:  1.0571E+00  6.5337E-01  2.3415E+05  1.4042E+00  1.9163E+00  8.6585E-01  1.2231E+00  5.9406E-01  7.2953E-01  2.1149E+00
             2.9101E+00
 PARAMETER:  1.5554E-01 -3.2561E-01  1.2464E+01  4.3946E-01  7.5042E-01 -4.4041E-02  3.0135E-01 -4.2078E-01 -2.1536E-01  8.4900E-01
             1.1682E+00
 GRADIENT:  -2.2420E+02  2.3513E+01 -8.0006E-04  7.1921E+01 -5.2704E+00 -5.9345E+00 -1.3752E+00 -2.8712E-05 -4.2867E+00  9.4716E-01
             5.7249E+00

0ITERATION NO.:   78    OBJECTIVE VALUE:  -2677.92052336385        NO. OF FUNC. EVALS.: 109
 CUMULATIVE NO. OF FUNC. EVALS.:     2043
 NPARAMETR:  1.0571E+00  6.5337E-01  2.2839E+05  1.4042E+00  1.9161E+00  8.6596E-01  1.2231E+00  5.9409E-01  7.2968E-01  2.1152E+00
             2.9101E+00
 PARAMETER:  1.5556E-01 -3.2563E-01  1.2451E+01  4.3944E-01  7.5031E-01 -4.3917E-02  3.0138E-01 -4.2072E-01 -2.1514E-01  8.4914E-01
             1.1682E+00
 GRADIENT:  -2.6982E-01 -2.9383E-02  1.6025E-03 -1.2734E-01  8.4276E-03  6.2967E-02 -8.0289E-03  8.6904E-03  7.4169E-02 -1.8079E-02
            -1.2128E-02
 NUMSIGDIG:         3.4         3.2         1.9         4.1         4.7         3.6         4.4         4.2         3.5         4.2
                    5.4

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2043
 NO. OF SIG. DIGITS IN FINAL EST.:  1.9
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.4815E-04 -6.4760E-03  2.2831E-07 -8.0093E-03 -1.3859E-02
 SE:             2.9163E-02  1.3278E-02  1.4652E-07  2.4434E-02  2.7285E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8500E-01  6.2573E-01  1.1917E-01  7.4307E-01  6.1150E-01

 ETASHRINKSD(%)  2.3011E+00  5.5519E+01  1.0000E+02  1.8143E+01  8.5930E+00
 ETASHRINKVR(%)  4.5493E+00  8.0214E+01  1.0000E+02  3.2994E+01  1.6448E+01
 EBVSHRINKSD(%)  2.3089E+00  5.5490E+01  1.0000E+02  1.8321E+01  6.5127E+00
 EBVSHRINKVR(%)  4.5645E+00  8.0188E+01  1.0000E+02  3.3285E+01  1.2601E+01
 RELATIVEINF(%)  9.5236E+01  5.3239E-03  8.3148E-10  1.7926E-02  5.9692E+01
 EPSSHRINKSD(%)  1.5351E+01
 EPSSHRINKVR(%)  2.8346E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          889
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1633.8727120379081     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2677.9205233638495     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1044.0478113259414     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    51.69
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.46
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2677.921       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.06E+00  6.53E-01  2.31E+05  1.40E+00  1.92E+00  8.66E-01  1.22E+00  5.94E-01  7.30E-01  2.12E+00  2.91E+00
 


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
+        1.57E+03
 
 TH 2
+        5.39E+02  4.72E+02
 
 TH 3
+        2.93E-05 -1.49E-05  1.69E-12
 
 TH 4
+        5.30E+00  4.52E+02  5.63E-06  7.22E+02
 
 TH 5
+       -1.63E+00 -2.45E+01 -7.69E-07 -2.23E+01  6.57E+01
 
 TH 6
+        3.26E+02 -4.51E+02 -1.58E-05 -3.19E+01  7.44E+01  8.68E+02
 
 TH 7
+        1.76E+02 -5.53E+00 -3.31E-06  3.36E+01 -1.50E+01  2.20E+02  4.56E+01
 
 TH 8
+       -6.57E+01 -1.96E+01 -5.09E-06  5.01E+01 -1.54E+01 -3.12E+02 -4.43E+01  1.37E+02
 
 TH 9
+        3.87E+02  2.36E+02  2.86E-05  6.60E+01  6.90E+00  7.15E+02 -1.19E+02  1.32E+02  2.05E+02
 
 TH10
+        1.16E+01  1.98E+01 -2.94E-06  2.75E+00 -8.36E+00  3.29E+01  6.25E+00 -1.93E+01 -1.77E+01  3.30E+01
 
 TH11
+       -2.07E+01 -4.01E+01  1.18E-08 -1.70E+01  4.91E+00  6.56E+01 -5.39E+00 -4.03E+00 -2.14E+01  4.06E+00  1.41E+02
 
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
 #CPUT: Total CPU Time in Seconds,       65.243
Stop Time:
Sat Sep 25 01:59:16 CDT 2021
