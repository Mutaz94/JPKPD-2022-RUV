Wed Sep 29 10:03:24 CDT 2021
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
$DATA ../../../../data/int/D/dat92.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m92.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   57846.0085199435        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.1723E+03  9.2319E+02  3.5720E+00  7.0625E+02 -2.3315E+01 -3.1428E+03 -2.0043E+03 -7.3719E+01 -2.9709E+03 -6.1310E+02
            -1.1303E+05

0ITERATION NO.:    5    OBJECTIVE VALUE:  -441.665145850291        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1109E+00  2.1964E+00  8.8424E-01  1.9348E+00  9.5007E-01  4.1721E+00  5.1408E+00  9.8383E-01  2.2106E+00  1.2454E+00
             1.3219E+01
 PARAMETER:  2.0518E-01  8.8681E-01 -2.3029E-02  7.5998E-01  4.8784E-02  1.5284E+00  1.7372E+00  8.3700E-02  8.9327E-01  3.1946E-01
             2.6817E+00
 GRADIENT:  -1.2533E+01  5.0084E+01 -4.4399E+01  1.3403E+02 -1.8045E+01  1.5621E+02  5.9318E+01  4.5856E+00  1.2870E+01  2.3225E+01
             7.9976E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -503.832889283755        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0132E+00  3.6333E+00  1.7744E+01  2.3195E+00  1.8379E+00  2.5947E+00  1.6350E+01  7.2482E-01  1.4998E+00  8.5809E-01
             1.3572E+01
 PARAMETER:  1.1315E-01  1.3901E+00  2.9760E+00  9.4133E-01  7.0862E-01  1.0535E+00  2.8942E+00 -2.2183E-01  5.0535E-01 -5.3045E-02
             2.7080E+00
 GRADIENT:  -4.8871E+01  3.0431E+01  6.7895E+00  1.2584E+02 -7.4401E+01  6.6266E+01  5.8412E+01  5.6625E-02 -3.4898E+01  6.3632E+00
             6.8113E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -573.468227483764        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.1825E+00  1.9207E+00  4.4997E+00  1.3797E+00  1.7676E+00  2.3444E+00  3.6419E+00  2.1314E+00  1.7671E+00  8.5474E-01
             1.3673E+01
 PARAMETER:  2.6766E-01  7.5266E-01  1.6040E+00  4.2187E-01  6.6965E-01  9.5205E-01  1.3925E+00  8.5679E-01  6.6934E-01 -5.6961E-02
             2.7154E+00
 GRADIENT:  -2.0493E+01  3.0514E+01  2.8526E-01  6.7777E+01 -3.2356E+01  1.5218E+01 -3.8478E+01  3.3895E-01  9.9870E-01  1.2774E+01
             9.7048E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -598.767603633033        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  1.2232E+00  1.6368E+00  5.8799E+00  9.6655E-01  2.3084E+00  2.2716E+00  3.9912E+00  2.0111E+00  1.1090E+00  3.3833E-01
             1.3110E+01
 PARAMETER:  3.0146E-01  5.9273E-01  1.8715E+00  6.5975E-02  9.3654E-01  9.2050E-01  1.4841E+00  7.9870E-01  2.0347E-01 -9.8375E-01
             2.6734E+00
 GRADIENT:   5.8280E+00 -1.2055E+01 -3.7688E+00 -1.1667E+01  2.4920E+01  4.7262E+00  9.4837E+00  2.3434E-01  6.0450E+00  1.8853E+00
             9.7450E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -600.467375859374        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  1.1853E+00  1.6263E+00  5.4013E+00  9.7818E-01  2.1223E+00  2.2766E+00  3.9917E+00  1.5454E+00  1.1723E+00  2.6068E-01
             1.2672E+01
 PARAMETER:  2.7000E-01  5.8629E-01  1.7866E+00  7.7943E-02  8.5250E-01  9.2267E-01  1.4842E+00  5.3526E-01  2.5901E-01 -1.2445E+00
             2.6394E+00
 GRADIENT:  -6.0567E-01 -4.8384E+00 -6.3421E-01 -5.8191E+00  5.9665E+00  4.2766E+00  7.6352E+00 -6.6410E-03  6.3738E+00  1.1731E+00
             4.7212E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -600.470590228879        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      442
 NPARAMETR:  1.1822E+00  1.6555E+00  5.0827E+00  9.5786E-01  2.1038E+00  2.2730E+00  3.9448E+00  1.4404E+00  1.1299E+00  2.3986E-01
             1.2600E+01
 PARAMETER:  2.6735E-01  6.0412E-01  1.7258E+00  5.6945E-02  8.4375E-01  9.2111E-01  1.4724E+00  4.6489E-01  2.2209E-01 -1.3277E+00
             2.6337E+00
 GRADIENT:  -6.3842E-01 -3.8741E+00 -4.4457E-01 -5.2411E+00  4.7658E+00  3.5261E+00  6.3668E+00  3.8685E-02  5.5977E+00  1.0017E+00
             3.8139E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -603.127254401085        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      579
 NPARAMETR:  1.2182E+00  1.8096E+00  5.5676E+00  9.9566E-01  2.1427E+00  2.3883E+00  4.5221E+00  1.3334E+00  9.6129E-01  2.0315E-01
             1.2651E+01
 PARAMETER:  2.9741E-01  6.9312E-01  1.8170E+00  9.5652E-02  8.6205E-01  9.7058E-01  1.6090E+00  3.8777E-01  6.0516E-02 -1.4938E+00
             2.6378E+00
 GRADIENT:  -6.6515E-01  6.6603E+00 -1.0473E+00 -5.3757E+00  2.1173E+00  2.1067E+00  9.4705E+00  2.0993E-01  3.7472E+00  6.8094E-01
            -1.8240E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -604.523030291590        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      754
 NPARAMETR:  1.2268E+00  1.5083E+00  7.2645E+00  1.0908E+00  2.1751E+00  2.3624E+00  4.7572E+00  1.2433E+00  8.0112E-01  1.5313E-01
             1.2806E+01
 PARAMETER:  3.0444E-01  5.1095E-01  2.0830E+00  1.8687E-01  8.7705E-01  9.5968E-01  1.6597E+00  3.1781E-01 -1.2175E-01 -1.7765E+00
             2.6499E+00
 GRADIENT:   2.8043E-01  1.0697E-01 -5.2042E-02  6.8449E-01 -1.5355E-03  4.2588E-02 -1.2003E-01  2.2372E-01 -7.0210E-01  3.6229E-01
            -4.6486E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -604.755900172729        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      930
 NPARAMETR:  1.2281E+00  1.4722E+00  6.0649E+00  1.1045E+00  2.1074E+00  2.3606E+00  4.8031E+00  6.0175E-01  8.2409E-01  5.9483E-02
             1.2840E+01
 PARAMETER:  3.0548E-01  4.8678E-01  1.9025E+00  1.9937E-01  8.4546E-01  9.5891E-01  1.6693E+00 -4.0791E-01 -9.3474E-02 -2.7221E+00
             2.6526E+00
 GRADIENT:   1.0927E-01 -4.9174E-01 -7.0557E-01  7.6309E-01  1.0207E+00 -1.3296E-01  6.5654E-01  9.2814E-02 -9.2992E-02  5.5516E-02
             1.7522E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -604.847697956884        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1106
 NPARAMETR:  1.2275E+00  1.4934E+00  6.7328E+00  1.0997E+00  2.1461E+00  2.3615E+00  4.7747E+00  1.9558E-01  8.3257E-01  1.2988E-02
             1.2839E+01
 PARAMETER:  3.0495E-01  5.0108E-01  2.0070E+00  1.9502E-01  8.6367E-01  9.5928E-01  1.6633E+00 -1.5318E+00 -8.3237E-02 -4.2437E+00
             2.6525E+00
 GRADIENT:   3.1856E-02 -5.5726E-02 -2.0627E-02  2.2585E-02  4.9753E-03  1.6078E-02  1.3412E-01  7.9102E-03  7.9093E-02  2.5765E-03
            -1.8674E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -604.854090081086        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1293             RESET HESSIAN, TYPE I
 NPARAMETR:  1.2277E+00  1.4951E+00  6.6707E+00  1.0978E+00  2.1445E+00  2.3668E+00  4.7855E+00  1.1845E-02  8.2197E-01  1.0000E-02
             1.2835E+01
 PARAMETER:  3.0511E-01  5.0216E-01  1.9977E+00  1.9331E-01  8.6290E-01  9.6156E-01  1.6656E+00 -4.3358E+00 -9.6055E-02 -4.6108E+00
             2.6522E+00
 GRADIENT:   1.2300E+01  4.9021E+00  1.2230E-01  2.6066E+00  2.3817E+00  1.7891E+01  3.5688E+01  3.4297E-05  1.0221E-02  0.0000E+00
             4.7081E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -604.854794152119        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1477             RESET HESSIAN, TYPE I
 NPARAMETR:  1.2276E+00  1.4930E+00  6.6973E+00  1.0994E+00  2.1444E+00  2.3664E+00  4.7909E+00  1.0000E-02  8.2446E-01  1.0000E-02
             1.2836E+01
 PARAMETER:  3.0508E-01  5.0076E-01  2.0017E+00  1.9476E-01  8.6285E-01  9.6136E-01  1.6667E+00 -4.5489E+00 -9.3033E-02 -4.6108E+00
             2.6523E+00
 GRADIENT:   1.2277E+01  4.9207E+00  1.4612E-01  2.6660E+00  2.2450E+00  1.7833E+01  3.5781E+01  8.2330E-07  2.7181E-02  0.0000E+00
             4.7155E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -604.855046280351        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1660
 NPARAMETR:  1.2276E+00  1.4909E+00  6.7193E+00  1.1010E+00  2.1444E+00  2.3663E+00  4.7946E+00  1.0000E-02  8.2437E-01  1.0000E-02
             1.2835E+01
 PARAMETER:  3.0509E-01  4.9941E-01  2.0050E+00  1.9622E-01  8.6285E-01  9.6134E-01  1.6675E+00 -4.5717E+00 -9.3131E-02 -4.6108E+00
             2.6522E+00
 GRADIENT:   1.3637E-01  1.4461E-01 -5.9733E-03  2.2103E-01 -1.0992E-01  6.2522E-01  6.9928E-01  0.0000E+00 -9.7487E-02  0.0000E+00
            -8.8257E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -604.855468661717        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1847
 NPARAMETR:  1.2276E+00  1.4893E+00  6.7345E+00  1.1014E+00  2.1446E+00  2.3663E+00  4.7979E+00  1.0000E-02  8.2641E-01  1.0000E-02
             1.2836E+01
 PARAMETER:  3.0510E-01  4.9832E-01  2.0073E+00  1.9662E-01  8.6295E-01  9.6132E-01  1.6682E+00 -4.5717E+00 -9.0665E-02 -4.6108E+00
             2.6523E+00
 GRADIENT:   1.3748E-01  1.0213E-01 -3.3048E-04  2.2186E-02 -1.2690E-01  6.2316E-01  8.1524E-01  0.0000E+00 -4.9160E-02  0.0000E+00
            -6.7024E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -604.855811866012        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     2041             RESET HESSIAN, TYPE I
 NPARAMETR:  1.2276E+00  1.4867E+00  6.7489E+00  1.1023E+00  2.1453E+00  2.3662E+00  4.8027E+00  1.0000E-02  8.2923E-01  1.0000E-02
             1.2837E+01
 PARAMETER:  3.0507E-01  4.9655E-01  2.0094E+00  1.9737E-01  8.6330E-01  9.6130E-01  1.6692E+00 -4.5259E+00 -8.7258E-02 -4.6108E+00
             2.6523E+00
 GRADIENT:   1.2273E+01  4.8075E+00  1.5392E-01  2.5665E+00  2.2148E+00  1.7826E+01  3.6090E+01  1.2205E-05  8.3361E-02  0.0000E+00
             4.7360E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -604.855972732019        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2224
 NPARAMETR:  1.2276E+00  1.4855E+00  6.7603E+00  1.1029E+00  2.1456E+00  2.3662E+00  4.8048E+00  1.0000E-02  8.2951E-01  1.0000E-02
             1.2837E+01
 PARAMETER:  3.0507E-01  4.9577E-01  2.0111E+00  1.9796E-01  8.6343E-01  9.6129E-01  1.6696E+00 -4.5604E+00 -8.6915E-02 -4.6108E+00
             2.6523E+00
 GRADIENT:   1.2365E-01  2.9637E-02 -1.7057E-02 -1.5964E-01 -2.5845E-02  6.2288E-01  9.5245E-01  0.0000E+00  4.9563E-03  0.0000E+00
            -4.1294E-01

0ITERATION NO.:   81    OBJECTIVE VALUE:  -604.855972732019        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     2246
 NPARAMETR:  1.2276E+00  1.4855E+00  6.7603E+00  1.1029E+00  2.1456E+00  2.3662E+00  4.8048E+00  1.0000E-02  8.2951E-01  1.0000E-02
             1.2837E+01
 PARAMETER:  3.0507E-01  4.9577E-01  2.0111E+00  1.9796E-01  8.6343E-01  9.6129E-01  1.6696E+00 -4.5604E+00 -8.6915E-02 -4.6108E+00
             2.6523E+00
 GRADIENT:   1.2365E-01  2.9637E-02 -1.7057E-02 -1.5964E-01 -2.5845E-02  6.2288E-01  9.5245E-01  0.0000E+00  4.9563E-03  0.0000E+00
            -4.1294E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2246
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6476E-02  1.7345E-02 -1.1029E-05 -5.1307E-02  2.3616E-06
 SE:             2.8820E-02  2.5627E-02  1.6432E-05  1.0383E-02  9.0822E-05
 N:                     100         100         100         100         100

 P VAL.:         5.6754E-01  4.9851E-01  5.0212E-01  7.7604E-07  9.7926E-01

 ETASHRINKSD(%)  3.4497E+00  1.4146E+01  9.9945E+01  6.5216E+01  9.9696E+01
 ETASHRINKVR(%)  6.7805E+00  2.6290E+01  1.0000E+02  8.7901E+01  9.9999E+01
 EBVSHRINKSD(%)  4.4163E+00  9.0883E+00  9.9934E+01  7.2172E+01  9.9603E+01
 EBVSHRINKVR(%)  8.6376E+00  1.7351E+01  1.0000E+02  9.2256E+01  9.9998E+01
 RELATIVEINF(%)  9.1118E+01  3.9025E+01  8.6075E-06  3.7237E+00  3.1450E-04
 EPSSHRINKSD(%)  3.8802E+00
 EPSSHRINKVR(%)  7.6099E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -604.85597273201904     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       1049.2333870363918     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    72.75
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    16.34
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -604.856       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.23E+00  1.49E+00  6.76E+00  1.10E+00  2.15E+00  2.37E+00  4.80E+00  1.00E-02  8.30E-01  1.00E-02  1.28E+01
 


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
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.79E+02
 
 TH 2
+       -1.41E+01  7.11E+00
 
 TH 3
+       -3.40E-01  5.75E-01  1.95E-01
 
 TH 4
+       -1.04E+02  3.00E+01 -6.38E-01  1.96E+02
 
 TH 5
+        1.04E+01 -1.18E+01 -3.49E+00  2.56E+00  6.27E+01
 
 TH 6
+       -4.98E+01  5.58E+00  1.07E-02  3.90E+01 -1.32E+00  3.42E+01
 
 TH 7
+        1.19E+01 -3.21E+00  4.13E-03 -1.98E+01  8.44E-01 -3.66E+00  2.06E+00
 
 TH 8
+        6.50E-20  2.98E-20  1.30E-20 -1.05E-19 -2.31E-19 -5.66E-20  6.18E-21  1.48E-38
 
 TH 9
+        2.82E+01 -7.43E+00 -1.22E-01 -4.33E+01  4.11E+00 -1.13E+01  4.47E+00  1.20E-19  1.02E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -3.26E+00 -1.11E+00 -5.32E-03 -5.35E+00  3.71E-01  1.10E+00  5.16E-01 -3.08E-20  9.67E-01  0.00E+00  6.69E-01
 
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
+        1.14E+02
 
 TH 2
+       -6.97E-01  1.81E+01
 
 TH 3
+        9.02E-02  2.24E-01  1.41E-01
 
 TH 4
+       -3.94E+00  2.57E+01 -3.98E-01  1.44E+02
 
 TH 5
+       -6.02E-01 -5.54E+00 -1.93E+00 -9.58E-01  3.51E+01
 
 TH 6
+       -1.65E+00 -1.99E-02  2.59E-02  1.01E+00 -1.21E-01  2.66E+01
 
 TH 7
+        6.26E-01  2.23E+00 -1.06E-01 -1.43E+01  1.23E+00 -5.53E-02  5.32E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  2.17E-01
 
 TH 9
+       -1.39E-01 -2.10E+00 -2.10E-01 -2.58E+01  2.64E+00 -2.51E-01  2.54E+00  0.00E+00  1.49E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -5.18E+00 -2.43E+00 -3.10E-02 -8.62E+00  6.99E-01  8.83E-01  4.99E-01  0.00E+00  1.79E+00  0.00E+00  5.09E+00
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.21E+02
 
 TH 2
+        3.70E+01  1.82E+01
 
 TH 3
+        8.30E-02  2.09E-01  8.18E-02
 
 TH 4
+        5.73E+01  2.68E+01  1.04E-01  1.51E+02
 
 TH 5
+       -4.83E+00 -4.41E+00 -1.22E+00 -6.79E+00  2.00E+01
 
 TH 6
+        2.95E+01  8.03E+00  3.80E-02 -1.56E+01  2.44E-01  3.35E+01
 
 TH 7
+        2.89E+00  1.90E+00 -1.06E-01 -1.50E+01  2.01E+00  6.95E+00  5.07E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -7.19E+00 -2.93E+00 -9.01E-02 -3.84E+01  2.81E+00  4.71E+00  3.51E+00  0.00E+00  1.77E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.76E+01 -6.40E+00 -1.98E-01 -2.57E+01  5.26E+00 -5.16E+00  2.91E+00  0.00E+00  6.96E+00  0.00E+00  1.05E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,       89.196
Stop Time:
Wed Sep 29 10:04:55 CDT 2021
