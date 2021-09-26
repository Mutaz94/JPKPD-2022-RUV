Sat Sep 25 01:27:25 CDT 2021
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
$DATA ../../../../data/int/SL2/dat64.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      997
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

 TOT. NO. OF OBS RECS:      897
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
 RAW OUTPUT FILE (FILE): m64.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2406.36467916612        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.9667E+01 -5.8666E+01 -1.4077E+01  1.2679E+02  1.6930E+02 -2.8389E+01 -1.0015E+02 -1.2312E+02 -1.2248E+02 -1.7640E+00
            -2.6394E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3138.75389015840        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0164E+00  1.2933E+00  1.1033E+00  8.2229E-01  1.0804E+00  1.0996E+00  1.1074E+00  1.0410E+00  1.1623E+00  9.8456E-01
             1.8749E+00
 PARAMETER:  1.1626E-01  3.5719E-01  1.9835E-01 -9.5662E-02  1.7731E-01  1.9492E-01  2.0197E-01  1.4022E-01  2.5039E-01  8.4444E-02
             7.2857E-01
 GRADIENT:   1.7092E+01  2.5463E+01  1.6366E+00 -2.2795E+01 -3.7541E+01  1.1776E+01  1.1458E+01 -4.6655E+00 -6.6495E+00 -1.3888E+01
            -5.8843E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3147.26732248203        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0464E+00  1.6356E+00  1.3982E+00  6.8644E-01  1.4024E+00  1.1195E+00  7.0825E-01  1.9812E+00  1.5110E+00  1.3597E+00
             1.8722E+00
 PARAMETER:  1.4540E-01  5.9202E-01  4.3516E-01 -2.7623E-01  4.3817E-01  2.1285E-01 -2.4496E-01  7.8368E-01  5.1276E-01  4.0725E-01
             7.2712E-01
 GRADIENT:   7.0648E+01  9.7130E+01 -1.9789E+00  5.7176E+01  8.3169E+00  1.6594E+01 -2.2387E+00 -8.2401E+00  1.4031E+01 -6.1342E-01
            -4.1958E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3154.81504969502        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0020E+00  1.4318E+00  1.7345E+00  7.6741E-01  1.3457E+00  1.0744E+00  8.9201E-01  2.7117E+00  1.1150E+00  1.2134E+00
             1.8937E+00
 PARAMETER:  1.0198E-01  4.5896E-01  6.5072E-01 -1.6473E-01  3.9692E-01  1.7175E-01 -1.4281E-02  1.0976E+00  2.0883E-01  2.9344E-01
             7.3856E-01
 GRADIENT:  -1.3384E+01  3.1311E+00 -1.3299E+01 -5.2074E+00  1.0128E+01  1.9344E+00  2.7288E-01 -3.6607E+00  4.6860E+00 -8.7872E+00
             1.2886E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3162.70467189102        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:      327
 NPARAMETR:  1.0065E+00  1.0762E+00  2.9610E+00  1.0547E+00  1.2155E+00  1.0693E+00  1.3575E+00  3.0516E+00  6.2499E-01  1.0917E+00
             1.8849E+00
 PARAMETER:  1.0644E-01  1.7344E-01  1.1855E+00  1.5329E-01  2.9517E-01  1.6702E-01  4.0566E-01  1.2157E+00 -3.7002E-01  1.8769E-01
             7.3387E-01
 GRADIENT:  -1.8505E+01  5.4915E+01 -6.1845E+00  5.3564E+01 -5.1344E+01 -2.6408E+00  1.8631E-02 -3.1031E+01 -1.2294E+01  5.8610E+00
             3.4559E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3163.27464662941        NO. OF FUNC. EVALS.: 128
 CUMULATIVE NO. OF FUNC. EVALS.:      455
 NPARAMETR:  1.0107E+00  1.0405E+00  2.9606E+00  1.0418E+00  1.2455E+00  1.0721E+00  1.3518E+00  3.4680E+00  6.8448E-01  1.0755E+00
             1.8851E+00
 PARAMETER:  1.1063E-01  1.3972E-01  1.1854E+00  1.4098E-01  3.1954E-01  1.6960E-01  4.0146E-01  1.3436E+00 -2.7909E-01  1.7279E-01
             7.3396E-01
 GRADIENT:  -9.6623E+00  8.3408E+00 -1.9010E+01 -2.6539E+00 -2.8102E+01 -1.4878E+00  6.1885E+00 -1.0833E+01 -1.7739E+00 -1.1971E+00
             3.8667E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3163.56848212907        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:      619
 NPARAMETR:  1.0107E+00  1.0275E+00  2.9624E+00  1.0496E+00  1.2453E+00  1.0790E+00  1.2708E+00  3.4656E+00  7.3345E-01  1.0725E+00
             1.8857E+00
 PARAMETER:  1.1068E-01  1.2709E-01  1.1860E+00  1.4840E-01  3.1938E-01  1.7601E-01  3.3968E-01  1.3429E+00 -2.0999E-01  1.6997E-01
             7.3433E-01
 GRADIENT:  -9.4854E+00  5.1291E-01 -2.0249E+01 -1.3767E-01 -2.1992E+01  1.0024E+00  3.4302E-02 -1.0149E+01  6.6264E-02 -2.8314E-01
             4.0702E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3163.57099680853        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      795
 NPARAMETR:  1.0107E+00  1.0241E+00  2.9624E+00  1.0515E+00  1.2453E+00  1.0764E+00  1.2751E+00  3.4657E+00  7.3127E-01  1.0723E+00
             1.8857E+00
 PARAMETER:  1.1068E-01  1.2380E-01  1.1860E+00  1.5018E-01  3.1938E-01  1.7362E-01  3.4299E-01  1.3429E+00 -2.1298E-01  1.6984E-01
             7.3432E-01
 GRADIENT:  -9.5480E+00 -2.2170E-01 -2.0461E+01  1.6644E-01 -2.0575E+01  8.5898E-02  8.9253E-02 -1.0335E+01  2.7379E-02 -3.0107E-02
             4.0875E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3165.06780847333        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:      981
 NPARAMETR:  1.0110E+00  1.0048E+00  3.1499E+00  1.0650E+00  1.2482E+00  1.0651E+00  1.2873E+00  3.5356E+00  6.9283E-01  1.0229E+00
             1.8417E+00
 PARAMETER:  1.1093E-01  1.0476E-01  1.2474E+00  1.6297E-01  3.2172E-01  1.6311E-01  3.5256E-01  1.3629E+00 -2.6697E-01  1.2264E-01
             7.1071E-01
 GRADIENT:  -8.5686E+00  2.8049E+00 -1.5950E+01 -6.8344E-01 -1.3922E+01 -4.2267E+00 -4.1576E+00 -1.3190E+01 -4.8051E+00 -7.9764E+00
            -7.9588E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3165.16517639806        NO. OF FUNC. EVALS.: 215
 CUMULATIVE NO. OF FUNC. EVALS.:     1196             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0195E+00  1.0044E+00  3.1522E+00  1.0652E+00  1.2481E+00  1.0851E+00  1.2883E+00  3.5358E+00  7.5786E-01  1.0239E+00
             1.8414E+00
 PARAMETER:  1.1934E-01  1.0438E-01  1.2481E+00  1.6315E-01  3.2161E-01  1.8165E-01  3.5329E-01  1.3629E+00 -1.7725E-01  1.2364E-01
             7.1055E-01
 GRADIENT:   2.4987E+01  3.1578E+00 -1.4486E+01  3.3064E+00 -7.9365E+00  6.5335E+00  3.0972E+00 -1.0351E+01  3.9378E+00 -7.2997E+00
            -4.7190E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3165.25801642394        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:     1330
 NPARAMETR:  1.0155E+00  1.0044E+00  3.1522E+00  1.0652E+00  1.2481E+00  1.0765E+00  1.2883E+00  3.5357E+00  7.3051E-01  1.0239E+00
             1.8415E+00
 PARAMETER:  1.1537E-01  1.0438E-01  1.2481E+00  1.6314E-01  3.2160E-01  1.7373E-01  3.5329E-01  1.3629E+00 -2.1402E-01  1.2363E-01
             7.1057E-01
 GRADIENT:   1.8615E-01  1.7608E+00 -1.6035E+01 -2.2372E+00 -1.2623E+01 -5.2915E-02 -3.0586E-01 -1.2538E+01  1.4307E-02 -7.5215E+00
            -6.7901E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -3165.30390205402        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1515
 NPARAMETR:  1.0169E+00  1.0044E+00  3.1585E+00  1.0652E+00  1.2482E+00  1.0744E+00  1.2883E+00  3.5423E+00  7.1895E-01  1.0239E+00
             1.8424E+00
 PARAMETER:  1.1676E-01  1.0437E-01  1.2501E+00  1.6316E-01  3.2172E-01  1.7180E-01  3.5329E-01  1.3648E+00 -2.2997E-01  1.2365E-01
             7.1108E-01
 GRADIENT:   2.7556E+00  2.0215E+00 -1.5807E+01 -2.1445E+00 -1.3224E+01 -7.8482E-01 -1.3994E+00 -1.2529E+01 -1.4309E+00 -7.5715E+00
            -5.9869E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -3165.30939025128        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:     1680
 NPARAMETR:  1.0157E+00  1.0044E+00  3.1577E+00  1.0651E+00  1.2481E+00  1.0766E+00  1.2884E+00  3.5413E+00  7.2844E-01  1.0239E+00
             1.8428E+00
 PARAMETER:  1.1553E-01  1.0440E-01  1.2498E+00  1.6311E-01  3.2164E-01  1.7380E-01  3.5338E-01  1.3645E+00 -2.1685E-01  1.2362E-01
             7.1127E-01
 GRADIENT:   4.5789E-01  1.8015E+00 -1.5855E+01 -2.4458E+00 -1.2955E+01 -2.2042E-02 -4.7080E-01 -1.2363E+01 -2.1128E-01 -7.5121E+00
            -5.1884E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -3165.35447676894        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1866
 NPARAMETR:  1.0155E+00  1.0045E+00  3.1628E+00  1.0651E+00  1.2481E+00  1.0766E+00  1.2886E+00  3.5462E+00  7.2969E-01  1.0238E+00
             1.8443E+00
 PARAMETER:  1.1534E-01  1.0445E-01  1.2515E+00  1.6304E-01  3.2160E-01  1.7379E-01  3.5356E-01  1.3659E+00 -2.1513E-01  1.2357E-01
             7.1209E-01
 GRADIENT:   5.7000E-02  1.7695E+00 -1.5684E+01 -2.7954E+00 -1.3269E+01 -2.1991E-02 -3.0453E-01 -1.2129E+01 -3.8115E-03 -7.4851E+00
            -3.2074E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -3165.37476520733        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     2036
 NPARAMETR:  1.0156E+00  1.0042E+00  3.1654E+00  1.0650E+00  1.2481E+00  1.0801E+00  1.2886E+00  3.5488E+00  7.3224E-01  1.0238E+00
             1.8447E+00
 PARAMETER:  1.1549E-01  1.0423E-01  1.2523E+00  1.6302E-01  3.2161E-01  1.7706E-01  3.5354E-01  1.3666E+00 -2.1165E-01  1.2355E-01
             7.1232E-01
 GRADIENT:   1.7033E+01  3.5420E+00 -1.3916E+01  2.9944E+00 -9.2808E+00  4.7539E+00  8.6141E-01 -1.0335E+01  7.1604E-01 -7.3029E+00
            -1.2412E+00

0ITERATION NO.:   72    OBJECTIVE VALUE:  -3165.37476520733        NO. OF FUNC. EVALS.:  66
 CUMULATIVE NO. OF FUNC. EVALS.:     2102
 NPARAMETR:  1.0156E+00  1.0042E+00  3.1654E+00  1.0650E+00  1.2481E+00  1.0799E+00  1.2886E+00  3.5487E+00  7.3209E-01  1.0238E+00
             1.8447E+00
 PARAMETER:  1.1549E-01  1.0423E-01  1.2523E+00  1.6302E-01  3.2161E-01  1.7706E-01  3.5354E-01  1.3666E+00 -2.1165E-01  1.2355E-01
             7.1232E-01
 GRADIENT:   3.0395E-01 -1.0891E+04  8.8534E+02  6.9538E+03  3.5002E+03  1.2185E+00 -1.6061E+03  8.1945E+02  3.3340E-01  9.1694E+03
            -7.9745E+02
 NUMSIGDIG:         2.6         3.3         3.3         3.3         3.3         1.5         3.3         3.3         1.5         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2102
 NO. OF SIG. DIGITS IN FINAL EST.:  1.5
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.2452E-04 -9.9457E-03 -1.2632E-02  3.0103E-03 -3.4231E-02
 SE:             2.9618E-02  2.2811E-02  2.3114E-02  2.1511E-02  2.1957E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7779E-01  6.6283E-01  5.8472E-01  8.8870E-01  1.1898E-01

 ETASHRINKSD(%)  7.7450E-01  2.3580E+01  2.2566E+01  2.7936E+01  2.6443E+01
 ETASHRINKVR(%)  1.5430E+00  4.1600E+01  4.0040E+01  4.8068E+01  4.5894E+01
 EBVSHRINKSD(%)  7.3781E-01  2.4367E+01  2.5993E+01  3.0435E+01  2.5572E+01
 EBVSHRINKVR(%)  1.4702E+00  4.2796E+01  4.5229E+01  5.1607E+01  4.4605E+01
 RELATIVEINF(%)  9.8514E+01  1.5961E+01  3.6097E+01  1.3697E+01  3.6144E+01
 EPSSHRINKSD(%)  1.8797E+01
 EPSSHRINKVR(%)  3.4061E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          897
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1648.5757285691827     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3165.3747652073307     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1516.7990366381480     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    66.65
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    16.46
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3165.375       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.00E+00  3.17E+00  1.07E+00  1.25E+00  1.08E+00  1.29E+00  3.55E+00  7.32E-01  1.02E+00  1.84E+00
 


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
+        9.15E+02
 
 TH 2
+       -1.57E+03  2.59E+07
 
 TH 3
+        9.75E+00  2.67E+02  1.79E+04
 
 TH 4
+        4.61E+02  1.39E+04 -4.09E+02  9.40E+06
 
 TH 5
+        1.09E+02  3.00E+03 -1.20E+03 -4.82E+03  1.75E+06
 
 TH 6
+        2.79E+00  1.00E+03 -5.99E+00 -3.03E+02 -7.38E+01  1.64E+02
 
 TH 7
+       -4.51E+02  5.95E+06  4.33E+01  2.35E+03  5.25E+02  2.84E+02  1.37E+06
 
 TH 8
+        1.43E+01  4.08E+02  4.67E+01 -7.68E+02  4.18E+02 -9.12E+00  7.11E+01  1.21E+04
 
 TH 9
+        9.41E+00  1.92E+03 -2.97E+01 -8.36E+02 -2.81E+02 -1.26E+00  5.44E+02 -2.54E+01  8.39E+01
 
 TH10
+        6.18E+02  1.75E+04 -7.65E+02  1.29E+07 -9.11E+03 -3.73E+02  3.06E+03 -1.34E+03 -1.11E+03  1.77E+07
 
 TH11
+       -3.85E+01 -7.97E+02  5.80E+02  1.17E+03  3.51E+03  2.00E+01 -1.24E+02 -1.39E+02  9.55E+01  2.29E+03  1.64E+05
 
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
 #CPUT: Total CPU Time in Seconds,       83.211
Stop Time:
Sat Sep 25 01:28:50 CDT 2021
