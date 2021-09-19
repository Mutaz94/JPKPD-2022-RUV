Sat Sep 18 11:36:00 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat28.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m28.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1643.10926374170        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.6341E+01 -1.1212E+02 -5.2192E+01 -1.1287E+02  8.8245E+01  1.3154E+01  3.8081E-01  4.0112E+00 -1.6702E+01  2.1223E+01
            -9.3722E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1654.37788009615        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.7416E-01  1.0092E+00  1.1242E+00  1.0784E+00  9.4962E-01  9.3169E-01  8.9199E-01  1.0049E+00  1.1035E+00  7.1370E-01
             1.0704E+00
 PARAMETER:  7.3817E-02  1.0920E-01  2.1706E-01  1.7544E-01  4.8306E-02  2.9244E-02 -1.4298E-02  1.0491E-01  1.9850E-01 -2.3729E-01
             1.6806E-01
 GRADIENT:   1.0742E+01  1.0382E+01  2.4952E+01  8.6174E-02 -8.2024E+00 -1.4310E+01  3.1184E+00 -1.0067E+01  1.0751E+01 -1.6132E+01
             9.2909E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1656.59555905284        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.7202E-01  1.1092E+00  1.3733E+00  1.0189E+00  1.0866E+00  9.4310E-01  5.5654E-01  1.5387E+00  1.2000E+00  9.2198E-01
             1.0394E+00
 PARAMETER:  7.1620E-02  2.0365E-01  4.1719E-01  1.1877E-01  1.8304E-01  4.1420E-02 -4.8602E-01  5.3097E-01  2.8229E-01  1.8772E-02
             1.3862E-01
 GRADIENT:   9.5396E+00  6.6231E+00  3.2146E+00  2.5391E+00  2.7652E+00 -8.6881E+00  1.7875E+00  1.8773E+00  4.4551E+00  3.8145E-01
             3.4149E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1657.26599481894        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.6568E-01  9.3987E-01  1.2229E+00  1.1220E+00  9.7135E-01  9.6727E-01  5.7402E-01  1.2436E+00  1.0846E+00  8.4243E-01
             1.0255E+00
 PARAMETER:  6.5076E-02  3.7984E-02  3.0122E-01  2.1511E-01  7.0934E-02  6.6724E-02 -4.5509E-01  3.1799E-01  1.8118E-01 -7.1470E-02
             1.2523E-01
 GRADIENT:  -2.9049E+00  2.8306E+00  5.4596E-01  5.7762E+00 -4.5925E+00  1.2864E+00  4.1063E-02  1.0246E+00  2.5446E-01 -1.0947E+00
            -2.1630E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1657.88764591431        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      390
 NPARAMETR:  9.8551E-01  7.8631E-01  1.3352E+00  1.2369E+00  9.5734E-01  9.7263E-01  5.2942E-01  1.2467E+00  1.0134E+00  8.5092E-01
             1.0296E+00
 PARAMETER:  8.5406E-02 -1.4041E-01  3.8911E-01  3.1257E-01  5.6404E-02  7.2246E-02 -5.3597E-01  3.2054E-01  1.1328E-01 -6.1440E-02
             1.2921E-01
 GRADIENT:   1.2246E+01  8.8835E+00  3.1441E+00  1.4038E+01 -6.0077E+00  1.2961E+00 -1.0808E-01 -4.3821E-01  8.8340E-02 -2.9817E-01
            -9.6270E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1658.30387594790        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      565
 NPARAMETR:  9.7626E-01  5.5832E-01  1.3484E+00  1.3752E+00  8.9269E-01  9.6586E-01  3.9258E-01  1.1844E+00  9.1463E-01  8.2051E-01
             1.0324E+00
 PARAMETER:  7.5974E-02 -4.8283E-01  3.9895E-01  4.1860E-01 -1.3514E-02  6.5259E-02 -8.3501E-01  2.6921E-01  1.0759E-02 -9.7827E-02
             1.3185E-01
 GRADIENT:  -4.6474E+00  3.8529E+00  1.6382E+00  8.2821E+00 -3.6185E+00 -7.6884E-01 -1.1004E-01 -4.6911E-01 -1.2080E+00 -1.2012E-01
             1.6910E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1658.44985164989        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      740
 NPARAMETR:  9.7549E-01  3.9987E-01  1.3543E+00  1.4687E+00  8.5164E-01  9.6600E-01  2.7008E-01  1.1838E+00  8.5573E-01  7.9534E-01
             1.0321E+00
 PARAMETER:  7.5180E-02 -8.1662E-01  4.0328E-01  4.8440E-01 -6.0592E-02  6.5408E-02 -1.2090E+00  2.6872E-01 -5.5801E-02 -1.2898E-01
             1.3160E-01
 GRADIENT:  -1.6429E+00  1.5542E-01 -3.0807E-01  1.5998E-01  2.6287E-01 -2.4218E-01 -1.6829E-02  2.3478E-01  2.9108E-01  2.4168E-02
             1.5060E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1658.47789349623        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      916
 NPARAMETR:  9.7522E-01  3.0990E-01  1.3452E+00  1.5233E+00  8.2425E-01  9.6633E-01  2.0041E-01  1.1770E+00  8.2123E-01  7.7660E-01
             1.0320E+00
 PARAMETER:  7.4905E-02 -1.0715E+00  3.9653E-01  5.2087E-01 -9.3280E-02  6.5755E-02 -1.5074E+00  2.6295E-01 -9.6948E-02 -1.5283E-01
             1.3149E-01
 GRADIENT:   7.9827E-01  2.4941E-01  2.2062E-01  1.0261E+00 -3.1931E-01  1.5036E-01 -5.9465E-03 -6.2826E-02 -1.5517E-01 -1.6749E-03
            -7.4646E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1658.49082832208        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1091
 NPARAMETR:  9.7363E-01  2.1820E-01  1.3308E+00  1.5782E+00  7.9542E-01  9.6555E-01  1.2983E-01  1.1764E+00  7.8896E-01  7.5530E-01
             1.0325E+00
 PARAMETER:  7.3275E-02 -1.4223E+00  3.8575E-01  5.5628E-01 -1.2889E-01  6.4945E-02 -1.9416E+00  2.6243E-01 -1.3704E-01 -1.8064E-01
             1.3195E-01
 GRADIENT:   4.4815E-01  4.3292E-01  6.3536E-01  2.9071E+00 -1.7333E+00  1.0551E-01 -1.1091E-03  1.4392E-02 -3.1698E-01  6.8031E-02
            -7.4672E-03

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1658.50217710129        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1267
 NPARAMETR:  9.7217E-01  1.5737E-01  1.3430E+00  1.6156E+00  7.8511E-01  9.6463E-01  8.4956E-02  1.1965E+00  7.6839E-01  7.4544E-01
             1.0326E+00
 PARAMETER:  7.1777E-02 -1.7492E+00  3.9491E-01  5.7973E-01 -1.4194E-01  6.3985E-02 -2.3656E+00  2.7937E-01 -1.6346E-01 -1.9378E-01
             1.3204E-01
 GRADIENT:  -5.0061E-01  3.0366E-01  2.7845E-01  3.5355E+00 -6.1665E-01 -7.4549E-02 -1.4281E-04  3.8672E-02 -3.7726E-01 -3.5294E-02
            -3.5660E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1658.51518786635        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1445
 NPARAMETR:  9.7126E-01  9.9600E-02  1.3206E+00  1.6475E+00  7.6351E-01  9.6424E-01  4.6542E-02  1.1857E+00  7.5094E-01  7.3073E-01
             1.0328E+00
 PARAMETER:  7.0835E-02 -2.2066E+00  3.7805E-01  5.9925E-01 -1.6983E-01  6.3586E-02 -2.9674E+00  2.7034E-01 -1.8643E-01 -2.1371E-01
             1.3227E-01
 GRADIENT:  -1.6261E-01  1.8359E-01  6.1774E-01  2.6440E+00 -1.7720E+00 -4.4257E-02 -1.7426E-05  4.4348E-02 -2.7806E-01  1.1874E-01
            -3.2987E-03

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1658.52439020186        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1622
 NPARAMETR:  9.7063E-01  5.9228E-02  1.3327E+00  1.6716E+00  7.5882E-01  9.6411E-01  2.3062E-02  1.2025E+00  7.3870E-01  7.2474E-01
             1.0328E+00
 PARAMETER:  7.0186E-02 -2.7264E+00  3.8722E-01  6.1381E-01 -1.7599E-01  6.3453E-02 -3.6696E+00  2.8438E-01 -2.0286E-01 -2.2195E-01
             1.3226E-01
 GRADIENT:   2.1078E-01  4.7734E-02  1.6696E-01  1.2324E+00 -1.9570E-01  4.3170E-02 -6.2061E-07 -1.2394E-03 -1.4802E-01 -7.4077E-02
            -5.9044E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1658.52868192921        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1799
 NPARAMETR:  9.7005E-01  3.5841E-02  1.3186E+00  1.6840E+00  7.4864E-01  9.6381E-01  1.1649E-02  1.1940E+00  7.3221E-01  7.1864E-01
             1.0330E+00
 PARAMETER:  6.9591E-02 -3.2287E+00  3.7657E-01  6.2120E-01 -1.8950E-01  6.3143E-02 -4.3525E+00  2.7727E-01 -2.1169E-01 -2.3039E-01
             1.3251E-01
 GRADIENT:  -6.6837E-02  3.5169E-02  3.1529E-01  1.2269E+00 -9.2341E-01  5.1248E-03 -6.3471E-08  2.3241E-02 -1.3088E-01  7.5178E-02
             2.3704E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1658.53219954729        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1976
 NPARAMETR:  9.6957E-01  1.7378E-02  1.3266E+00  1.6950E+00  7.4759E-01  9.6343E-01  1.0000E-02  1.2036E+00  7.2701E-01  7.1669E-01
             1.0330E+00
 PARAMETER:  6.9094E-02 -3.9525E+00  3.8259E-01  6.2770E-01 -1.9090E-01  6.2747E-02 -5.3480E+00  2.8530E-01 -2.1882E-01 -2.3312E-01
             1.3246E-01
 GRADIENT:  -3.1912E-01  2.1424E-03 -7.7789E-02 -6.7670E-02  2.5319E-01 -8.7832E-02  0.0000E+00  7.2484E-04  3.9503E-02 -1.1257E-02
            -6.9344E-03

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1658.53334637914        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2152
 NPARAMETR:  9.6956E-01  1.0000E-02  1.3206E+00  1.6992E+00  7.4371E-01  9.6363E-01  1.0000E-02  1.1998E+00  7.2475E-01  7.1421E-01
             1.0331E+00
 PARAMETER:  6.9092E-02 -4.5219E+00  3.7808E-01  6.3014E-01 -1.9611E-01  6.2951E-02 -6.1313E+00  2.8212E-01 -2.2193E-01 -2.3658E-01
             1.3255E-01
 GRADIENT:   1.6722E-02  0.0000E+00  1.7903E-01  9.8259E-01 -5.0878E-01  1.9501E-02  0.0000E+00  8.3913E-03 -1.1447E-01  3.1012E-02
             1.5044E-03

0ITERATION NO.:   73    OBJECTIVE VALUE:  -1658.53358121567        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:     2253
 NPARAMETR:  9.6949E-01  1.0000E-02  1.3214E+00  1.6987E+00  7.4423E-01  9.6354E-01  1.0000E-02  1.2004E+00  7.2491E-01  7.1432E-01
             1.0331E+00
 PARAMETER:  6.9112E-02 -4.5169E+00  3.7878E-01  6.2995E-01 -1.9549E-01  6.2933E-02 -6.1248E+00  2.8276E-01 -2.2162E-01 -2.3639E-01
             1.3252E-01
 GRADIENT:   6.9989E-02  0.0000E+00  1.8759E-02  1.0844E-01 -4.5291E-02  1.0133E-02  0.0000E+00  2.4354E-03  9.2114E-03  1.0482E-03
             4.9814E-04

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2253
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.9171E-05 -2.7532E-06 -2.2958E-02 -4.8283E-03 -3.0342E-02
 SE:             2.9836E-02  1.8159E-06  1.9871E-02  2.9243E-02  1.8758E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9735E-01  1.2948E-01  2.4793E-01  8.6886E-01  1.0577E-01

 ETASHRINKSD(%)  4.6552E-02  9.9994E+01  3.3431E+01  2.0317E+00  3.7158E+01
 ETASHRINKVR(%)  9.3083E-02  1.0000E+02  5.5685E+01  4.0221E+00  6.0508E+01
 EBVSHRINKSD(%)  4.5238E-01  9.9994E+01  3.4586E+01  2.5299E+00  3.6110E+01
 EBVSHRINKVR(%)  9.0272E-01  1.0000E+02  5.7210E+01  4.9958E+00  5.9180E+01
 RELATIVEINF(%)  9.6243E+01  1.7225E-08  6.2528E+00  7.2219E+00  3.5121E+00
 EPSSHRINKSD(%)  4.4527E+01
 EPSSHRINKVR(%)  6.9227E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1658.5335812156698     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -923.38275465193158     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.28
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.53
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1658.534       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.70E-01  1.00E-02  1.32E+00  1.70E+00  7.44E-01  9.64E-01  1.00E-02  1.20E+00  7.25E-01  7.14E-01  1.03E+00
 


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
+        0.00E+00  0.00E+00
 
 TH 3
+       -9.48E-01  0.00E+00  1.94E+02
 
 TH 4
+       -1.10E+01  0.00E+00 -2.26E+01  6.96E+02
 
 TH 5
+        2.48E+00  0.00E+00 -4.90E+02 -1.14E+02  1.64E+03
 
 TH 6
+        3.89E+00  0.00E+00  3.12E-01 -3.74E+00 -5.77E+00  2.21E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        3.43E+00  0.00E+00 -3.05E+01 -2.74E+00 -2.11E+01 -5.98E+00  0.00E+00  4.08E+01
 
 TH 9
+        2.87E+00  0.00E+00  6.14E+00 -1.92E+00  5.39E+00 -1.68E+00  0.00E+00 -3.83E-02  3.47E+02
 
 TH10
+       -6.00E+00  0.00E+00  3.42E+00 -9.83E-01 -1.18E+02  8.79E+00  0.00E+00  2.20E+01  1.91E+00  9.40E+01
 
 TH11
+       -4.19E+00  0.00E+00 -5.35E+00 -7.19E+00 -1.20E+01  2.81E+00  0.00E+00  1.07E+01  1.32E+01  1.53E+01  1.96E+02
 
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
 #CPUT: Total CPU Time in Seconds,       30.855
Stop Time:
Sat Sep 18 11:36:33 CDT 2021
