Wed Sep 29 18:22:12 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat68.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m68.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1708.85724646088        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0883E+02  4.8752E+01 -2.2645E+01  1.3068E+02  4.8383E+01  4.9059E+01  2.1489E+01  3.9681E+00  4.5380E+01  2.2191E+00
             2.0807E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1714.32173605219        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      192
 NPARAMETR:  1.0118E+00  1.0133E+00  1.0028E+00  9.6739E-01  9.6121E-01  9.9103E-01  8.3948E-01  9.8222E-01  7.2355E-01  9.7426E-01
             9.1447E-01
 PARAMETER:  1.1177E-01  1.1318E-01  1.0283E-01  6.6843E-02  6.0440E-02  9.0989E-02 -7.4974E-02  8.2064E-02 -2.2359E-01  7.3923E-02
             1.0585E-02
 GRADIENT:  -1.2020E+01  2.8529E+01  1.9181E+01  2.3511E+00 -2.3281E+01 -1.5799E+00 -7.6690E+00 -2.8937E+00 -2.7047E+01 -2.4225E+00
            -2.4758E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1715.72976134454        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  1.0119E+00  9.0545E-01  8.9084E-01  1.0250E+00  8.7118E-01  9.9682E-01  7.2256E-01  7.0496E-01  7.8409E-01  9.5122E-01
             9.0354E-01
 PARAMETER:  1.1186E-01  6.7943E-04 -1.5590E-02  1.2468E-01 -3.7908E-02  9.6818E-02 -2.2496E-01 -2.4962E-01 -1.4324E-01  4.9985E-02
            -1.4354E-03
 GRADIENT:  -1.2722E+01  1.7119E+01  9.4075E+00  1.0965E+01 -1.5809E+01  2.7967E-01 -7.5871E+00 -2.1350E+00 -1.5139E+01 -5.4873E-01
            -2.8717E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1717.80190041516        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      549
 NPARAMETR:  1.0185E+00  8.8859E-01  8.6197E-01  1.0306E+00  8.5601E-01  9.9679E-01  9.1007E-01  7.0343E-01  7.7991E-01  9.0449E-01
             9.6454E-01
 PARAMETER:  1.1832E-01 -1.8117E-02 -4.8539E-02  1.3016E-01 -5.5469E-02  9.6788E-02  5.7675E-03 -2.5179E-01 -1.4858E-01 -3.8740E-04
             6.3897E-02
 GRADIENT:   5.9162E-01  2.5416E+00 -5.1951E-01  4.1443E+00 -1.9487E-01  5.3210E-01 -4.6948E-01  1.8383E-01  2.1802E-01  5.5277E-01
             1.5476E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1717.85760195200        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      727
 NPARAMETR:  1.0166E+00  7.6702E-01  9.0026E-01  1.1064E+00  8.2315E-01  9.9453E-01  1.0236E+00  6.9598E-01  7.3309E-01  8.8816E-01
             9.6672E-01
 PARAMETER:  1.1642E-01 -1.6525E-01 -5.0720E-03  2.0115E-01 -9.4615E-02  9.4518E-02  1.2335E-01 -2.6244E-01 -2.1049E-01 -1.8600E-02
             6.6149E-02
 GRADIENT:  -1.2886E+00  4.8747E+00  3.5372E+00  5.2449E+00 -6.1327E+00  1.1075E-01 -4.4077E-01 -7.3256E-01 -2.2912E-01 -5.4165E-01
             1.5215E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1717.92937538253        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      906
 NPARAMETR:  1.0151E+00  6.5923E-01  1.0474E+00  1.1845E+00  8.4959E-01  9.9248E-01  1.0529E+00  8.9017E-01  7.0926E-01  9.2272E-01
             9.7640E-01
 PARAMETER:  1.1495E-01 -3.1668E-01  1.4632E-01  2.6935E-01 -6.3003E-02  9.2447E-02  1.5159E-01 -1.6343E-02 -2.4353E-01  1.9574E-02
             7.6113E-02
 GRADIENT:   4.4953E-01  8.5346E+00  3.2689E+00  1.4033E+01 -1.0982E+01  3.4370E-01 -8.7241E-01  4.3935E-01 -1.0177E+00  9.1324E-01
             6.0460E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1718.36245824501        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1083
 NPARAMETR:  1.0120E+00  4.2556E-01  1.3558E+00  1.3460E+00  9.0125E-01  9.8678E-01  1.1737E+00  1.1769E+00  6.5890E-01  9.9390E-01
             9.6771E-01
 PARAMETER:  1.1195E-01 -7.5435E-01  4.0438E-01  3.9713E-01 -3.9701E-03  8.6696E-02  2.6017E-01  2.6289E-01 -3.1718E-01  9.3885E-02
             6.7173E-02
 GRADIENT:   4.3471E+00  7.8564E+00 -7.1806E-01  2.7038E+01 -2.1608E+00  7.0419E-02 -7.3941E-01  2.2740E-02 -2.1030E+00  1.5067E+00
             2.5279E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1718.63685575623        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1262
 NPARAMETR:  1.0068E+00  2.9487E-01  1.4990E+00  1.4328E+00  9.0733E-01  9.8329E-01  1.3438E+00  1.2820E+00  6.3182E-01  9.9320E-01
             9.4875E-01
 PARAMETER:  1.0683E-01 -1.1212E+00  5.0480E-01  4.5960E-01  2.7523E-03  8.3149E-02  3.9554E-01  3.4841E-01 -3.5915E-01  9.3175E-02
             4.7391E-02
 GRADIENT:  -1.9130E+00  6.9444E+00  3.6770E+00  3.2050E+01 -3.1606E+00 -4.2217E-01 -5.6931E-01 -3.3581E+00 -1.8214E+00 -1.3948E+00
            -6.2334E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1719.21081912134        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1440
 NPARAMETR:  1.0023E+00  1.4923E-01  1.8490E+00  1.5346E+00  9.5426E-01  9.8024E-01  1.7616E+00  1.6529E+00  6.0230E-01  1.0143E+00
             9.4790E-01
 PARAMETER:  1.0231E-01 -1.8023E+00  7.1463E-01  5.2825E-01  5.3184E-02  8.0047E-02  6.6620E-01  6.0251E-01 -4.0700E-01  1.1425E-01
             4.6494E-02
 GRADIENT:  -7.0816E+00  3.8785E+00  2.9658E+00  2.3918E+01 -6.6371E+00 -4.3367E-01 -1.7993E-01  1.7556E-01  2.5236E-01 -1.6226E+00
            -5.8535E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1719.70705015797        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1618
 NPARAMETR:  1.0040E+00  5.4589E-02  2.0896E+00  1.5987E+00  9.8084E-01  9.7905E-01  2.8798E+00  1.8313E+00  5.7612E-01  1.0387E+00
             9.6152E-01
 PARAMETER:  1.0397E-01 -2.8079E+00  8.3699E-01  5.6916E-01  8.0653E-02  7.8827E-02  1.1577E+00  7.0502E-01 -4.5143E-01  1.3801E-01
             6.0764E-02
 GRADIENT:  -6.3801E-01  1.3109E+00  3.7276E+00  1.1889E+01 -6.4696E+00 -1.4269E-02 -8.6342E-02 -4.7761E-01 -2.2064E+00 -7.8830E-01
             2.2228E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1720.08084624691        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1804
 NPARAMETR:  1.0043E+00  1.0000E-02  2.0446E+00  1.6132E+00  9.7554E-01  9.7864E-01  5.5406E+00  1.8390E+00  5.7296E-01  1.0456E+00
             9.6011E-01
 PARAMETER:  1.0433E-01 -4.8112E+00  8.1520E-01  5.7819E-01  7.5241E-02  7.8413E-02  1.8121E+00  7.0925E-01 -4.5693E-01  1.4462E-01
             5.9293E-02
 GRADIENT:   2.5564E+00  0.0000E+00 -3.4629E+00 -3.7306E+01  7.9015E+00  2.0276E-01 -2.5914E-03  2.3343E+00  5.9084E-01  4.6956E-01
             6.9995E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1720.14177300323        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1988
 NPARAMETR:  1.0037E+00  1.0000E-02  2.0519E+00  1.6178E+00  9.6878E-01  9.7847E-01  5.7914E+00  1.8090E+00  5.7222E-01  1.0383E+00
             9.5941E-01
 PARAMETER:  1.0371E-01 -4.8112E+00  8.1875E-01  5.8104E-01  6.8283E-02  7.8233E-02  1.8564E+00  6.9278E-01 -4.5823E-01  1.3761E-01
             5.8563E-02
 GRADIENT:   9.9499E-01  0.0000E+00  1.1796E+00 -1.7512E+01 -2.9458E-01  6.0688E-02 -5.0123E-03  3.5290E-03  1.6305E-01  7.5595E-02
            -1.0986E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1720.15740609001        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     2185             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0041E+00  1.0000E-02  2.0308E+00  1.6157E+00  9.6667E-01  9.7854E-01  6.3642E+00  1.8003E+00  5.7235E-01  1.0350E+00
             9.5944E-01
 PARAMETER:  1.0408E-01 -4.8112E+00  8.0845E-01  5.7974E-01  6.6097E-02  7.8306E-02  1.9507E+00  6.8795E-01 -4.5800E-01  1.3440E-01
             5.8594E-02
 GRADIENT:   5.0575E+02  0.0000E+00  9.7714E+00  1.2781E+03  9.2928E+00  4.7630E+01  2.5660E-01  3.9332E+00  3.0665E+01  1.1454E+00
             7.8492E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1720.16504337898        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     2376             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0041E+00  1.0000E-02  2.0196E+00  1.6151E+00  9.6454E-01  9.7857E-01  6.9117E+00  1.7912E+00  5.7242E-01  1.0326E+00
             9.5952E-01
 PARAMETER:  1.0407E-01 -4.8112E+00  8.0290E-01  5.7939E-01  6.3894E-02  7.8333E-02  2.0332E+00  6.8291E-01 -4.5789E-01  1.3205E-01
             5.8681E-02
 GRADIENT:   5.0554E+02  0.0000E+00  9.6806E+00  1.2759E+03  9.6400E+00  4.7618E+01  3.1055E-01  3.8332E+00  3.0662E+01  9.8854E-01
             8.1318E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1720.17210610671        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     2567             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0040E+00  1.0000E-02  2.0122E+00  1.6152E+00  9.6221E-01  9.7854E-01  7.3672E+00  1.7843E+00  5.7240E-01  1.0318E+00
             9.5945E-01
 PARAMETER:  1.0399E-01 -4.8112E+00  7.9923E-01  5.7945E-01  6.1480E-02  7.8307E-02  2.0970E+00  6.7900E-01 -4.5791E-01  1.3133E-01
             5.8605E-02
 GRADIENT:   5.0499E+02  0.0000E+00  1.0016E+01  1.2776E+03  8.6076E+00  4.7597E+01  3.6000E-01  3.7514E+00  3.0676E+01  1.1294E+00
             7.7865E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1720.17758373634        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:     2741
 NPARAMETR:  1.0040E+00  1.0000E-02  2.0102E+00  1.6145E+00  9.6006E-01  9.7855E-01  9.6227E+00  1.7800E+00  5.7244E-01  1.0333E+00
             9.5949E-01
 PARAMETER:  1.0399E-01 -4.8112E+00  7.9824E-01  5.7905E-01  5.9246E-02  7.8312E-02  2.3641E+00  6.7659E-01 -4.5785E-01  1.3273E-01
             5.8651E-02
 GRADIENT:   5.0487E+02  0.0000E+00  1.0964E+01  1.2752E+03  6.4599E+00  4.7590E+01  7.2739E-01  3.6805E+00  3.0927E+01  1.6236E+00
             8.6827E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1720.19552895387        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:     2812
 NPARAMETR:  1.0046E+00  1.0000E-02  1.9382E+00  1.6142E+00  9.4650E-01  9.7970E-01  5.0768E+00  1.7275E+00  5.7331E-01  1.0217E+00
             9.5950E-01
 PARAMETER:  1.0463E-01 -4.8112E+00  7.6178E-01  5.7883E-01  4.5017E-02  7.9487E-02  1.7247E+00  6.4669E-01 -4.5632E-01  1.2146E-01
             5.8658E-02
 GRADIENT:   5.0858E+02  0.0000E+00  9.4286E+00  1.2794E+03  7.8851E+00  4.7809E+01  1.4944E-01  3.5350E+00  3.0491E+01  1.1481E+00
             7.6866E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1720.20917272191        NO. OF FUNC. EVALS.: 145
 CUMULATIVE NO. OF FUNC. EVALS.:     2957
 NPARAMETR:  1.0037E+00  1.0000E-02  1.9269E+00  1.6109E+00  9.4267E-01  9.7788E-01  4.6498E+00  1.7124E+00  5.7401E-01  1.0211E+00
             9.5971E-01
 PARAMETER:  1.0365E-01 -4.8112E+00  7.5590E-01  5.7681E-01  4.0960E-02  7.7629E-02  1.6368E+00  6.3789E-01 -4.5510E-01  1.2093E-01
             5.8879E-02
 GRADIENT:   1.4274E+00  0.0000E+00  6.7665E-01 -2.3935E+01 -5.1817E-01 -1.9134E-01 -3.3456E-03  7.2463E-02  2.8305E-01  2.6308E-01
             1.5733E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1720.21194199352        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     3150             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0039E+00  1.0000E-02  1.9196E+00  1.6108E+00  9.4209E-01  9.7862E-01  5.0530E+00  1.7090E+00  5.7393E-01  1.0184E+00
             9.5937E-01
 PARAMETER:  1.0385E-01 -4.8112E+00  7.5214E-01  5.7671E-01  4.0342E-02  7.8388E-02  1.7200E+00  6.3593E-01 -4.5525E-01  1.1819E-01
             5.8523E-02
 GRADIENT:   5.0306E+02  0.0000E+00  9.9852E+00  1.2642E+03  8.0865E+00  4.7444E+01  1.4876E-01  3.2299E+00  3.0533E+01  1.0517E+00
             7.8120E-01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1720.21349271447        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     3343
 NPARAMETR:  1.0038E+00  1.0000E-02  1.9156E+00  1.6106E+00  9.4131E-01  9.7862E-01  5.4167E+00  1.7062E+00  5.7392E-01  1.0176E+00
             9.5938E-01
 PARAMETER:  1.0384E-01 -4.8112E+00  7.5002E-01  5.7661E-01  3.9515E-02  7.8391E-02  1.7895E+00  6.3425E-01 -4.5527E-01  1.1749E-01
             5.8530E-02
 GRADIENT:   1.9178E+00  0.0000E+00 -5.5295E-02 -2.3719E+01  1.0033E+00  1.0570E-01 -4.9230E-03  1.5169E-01  2.0768E-01 -1.0183E-01
            -1.5302E-02

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1720.21524787844        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     3540             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0038E+00  1.0000E-02  1.9123E+00  1.6104E+00  9.3898E-01  9.7863E-01  6.0437E+00  1.7014E+00  5.7392E-01  1.0187E+00
             9.5943E-01
 PARAMETER:  1.0383E-01 -4.8112E+00  7.4829E-01  5.7647E-01  3.7037E-02  7.8399E-02  1.8990E+00  6.3146E-01 -4.5527E-01  1.1857E-01
             5.8587E-02
 GRADIENT:   5.0280E+02  0.0000E+00  1.0750E+01  1.2633E+03  6.0747E+00  4.7430E+01  2.2639E-01  3.1799E+00  3.0560E+01  1.4694E+00
             8.4227E-01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1720.21793628234        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     3733
 NPARAMETR:  1.0038E+00  1.0000E-02  1.9078E+00  1.6102E+00  9.3846E-01  9.7864E-01  6.6507E+00  1.6983E+00  5.7387E-01  1.0180E+00
             9.5946E-01
 PARAMETER:  1.0382E-01 -4.8112E+00  7.4594E-01  5.7636E-01  3.6488E-02  7.8405E-02  1.9947E+00  6.2962E-01 -4.5534E-01  1.1782E-01
             5.8614E-02
 GRADIENT:   1.8232E+00  0.0000E+00  5.1722E-01 -2.3620E+01 -5.9264E-01  1.2235E-01 -6.9412E-03  1.2371E-01  2.5544E-01  2.6267E-01
             5.2592E-02

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1720.22062696678        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     3930             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0038E+00  1.0000E-02  1.9022E+00  1.6100E+00  9.3817E-01  9.7864E-01  7.3600E+00  1.6950E+00  5.7377E-01  1.0151E+00
             9.5938E-01
 PARAMETER:  1.0381E-01 -4.8112E+00  7.4303E-01  5.7624E-01  3.6175E-02  7.8407E-02  2.0961E+00  6.2766E-01 -4.5553E-01  1.1495E-01
             5.8530E-02
 GRADIENT:   5.0256E+02  0.0000E+00  9.9563E+00  1.2623E+03  8.0742E+00  4.7403E+01  3.5783E-01  3.1283E+00  3.0542E+01  9.6409E-01
             7.7537E-01

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1720.22248289384        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     4123
 NPARAMETR:  1.0038E+00  1.0000E-02  1.8990E+00  1.6099E+00  9.3760E-01  9.7864E-01  8.0104E+00  1.6926E+00  5.7367E-01  1.0143E+00
             9.5936E-01
 PARAMETER:  1.0381E-01 -4.8112E+00  7.4132E-01  5.7616E-01  3.5565E-02  7.8410E-02  2.1807E+00  6.2625E-01 -4.5569E-01  1.1419E-01
             5.8506E-02
 GRADIENT:   2.0114E+00  0.0000E+00 -1.0696E-01 -2.3806E+01  1.0452E+00  1.2115E-01 -6.4379E-03  9.4675E-02  2.3243E-01 -1.8595E-01
            -3.9725E-02

0ITERATION NO.:  120    OBJECTIVE VALUE:  -1720.22527213994        NO. OF FUNC. EVALS.: 147
 CUMULATIVE NO. OF FUNC. EVALS.:     4270             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0038E+00  1.0000E-02  1.8974E+00  1.6098E+00  9.3689E-01  9.7865E-01  9.6212E+00  1.6909E+00  5.7334E-01  1.0150E+00
             9.5937E-01
 PARAMETER:  1.0380E-01 -4.8112E+00  7.4050E-01  5.7611E-01  3.4812E-02  7.8417E-02  2.3640E+00  6.2526E-01 -4.5627E-01  1.1486E-01
             5.8518E-02
 GRADIENT:   5.0252E+02  0.0000E+00  1.0058E+01  1.2620E+03  7.7161E+00  4.7416E+01  7.2521E-01  3.1217E+00  3.0560E+01  1.0995E+00
             7.8576E-01

0ITERATION NO.:  125    OBJECTIVE VALUE:  -1720.22669768371        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     4459
 NPARAMETR:  1.0038E+00  1.0000E-02  1.8910E+00  1.6095E+00  9.3537E-01  9.7864E-01  9.6129E+00  1.6853E+00  5.7329E-01  1.0138E+00
             9.5948E-01
 PARAMETER:  1.0379E-01 -4.8112E+00  7.3713E-01  5.7592E-01  3.3187E-02  7.8410E-02  2.3631E+00  6.2194E-01 -4.5636E-01  1.1375E-01
             5.8634E-02
 GRADIENT:   1.8596E+00  0.0000E+00  1.1632E-01 -2.3721E+01  3.4702E-01  1.0768E-01  1.8229E-02  6.9702E-02  1.5797E-01 -8.9316E-03
             2.3020E-02

0ITERATION NO.:  130    OBJECTIVE VALUE:  -1720.22770431730        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     4653             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0038E+00  1.0000E-02  1.8826E+00  1.6091E+00  9.3281E-01  9.7866E-01  9.5032E+00  1.6791E+00  5.7345E-01  1.0145E+00
             9.5950E-01
 PARAMETER:  1.0377E-01 -4.8112E+00  7.3266E-01  5.7567E-01  3.0450E-02  7.8430E-02  2.3516E+00  6.1823E-01 -4.5608E-01  1.1441E-01
             5.8652E-02
 GRADIENT:   5.0200E+02  0.0000E+00  1.0348E+01  1.2598E+03  6.4786E+00  4.7372E+01  6.9619E-01  3.1541E+00  3.0490E+01  1.4641E+00
             8.8208E-01

0ITERATION NO.:  135    OBJECTIVE VALUE:  -1720.22857257235        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     4844
 NPARAMETR:  1.0038E+00  1.0000E-02  1.8800E+00  1.6089E+00  9.3235E-01  9.7866E-01  9.6135E+00  1.6759E+00  5.7361E-01  1.0131E+00
             9.5949E-01
 PARAMETER:  1.0376E-01 -4.8112E+00  7.3125E-01  5.7558E-01  2.9956E-02  7.8427E-02  2.3632E+00  6.1635E-01 -4.5581E-01  1.1302E-01
             5.8643E-02
 GRADIENT:   1.8402E+00  0.0000E+00  3.4098E-01 -2.3714E+01 -4.8496E-01  1.2108E-01  1.8552E-02  8.1137E-02  2.7270E-01  2.0179E-01
             5.4676E-02

0ITERATION NO.:  140    OBJECTIVE VALUE:  -1720.22904691610        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     5038             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0038E+00  1.0000E-02  1.8742E+00  1.6087E+00  9.3218E-01  9.7866E-01  9.5603E+00  1.6723E+00  5.7359E-01  1.0108E+00
             9.5938E-01
 PARAMETER:  1.0376E-01 -4.8112E+00  7.2820E-01  5.7545E-01  2.9774E-02  7.8424E-02  2.3576E+00  6.1422E-01 -4.5584E-01  1.1078E-01
             5.8534E-02
 GRADIENT:   5.0177E+02  0.0000E+00  9.6978E+00  1.2589E+03  8.3862E+00  4.7346E+01  7.1020E-01  2.9951E+00  3.0487E+01  9.0952E-01
             7.7424E-01

0ITERATION NO.:  145    OBJECTIVE VALUE:  -1720.22948426552        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     5227
 NPARAMETR:  1.0038E+00  1.0000E-02  1.8731E+00  1.6087E+00  9.3162E-01  9.7866E-01  9.5556E+00  1.6715E+00  5.7359E-01  1.0105E+00
             9.5939E-01
 PARAMETER:  1.0375E-01 -4.8112E+00  7.2760E-01  5.7542E-01  2.9165E-02  7.8430E-02  2.3571E+00  6.1374E-01 -4.5584E-01  1.1047E-01
             5.8545E-02
 GRADIENT:   1.8369E+00  0.0000E+00 -1.5640E-01 -2.3721E+01  7.0911E-01  1.0597E-01  1.5809E-02  7.7118E-02  1.9230E-01 -1.1104E-01
            -2.0812E-02

0ITERATION NO.:  150    OBJECTIVE VALUE:  -1720.22969967370        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     5421             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0038E+00  1.0000E-02  1.8723E+00  1.6086E+00  9.3055E-01  9.7867E-01  9.5662E+00  1.6695E+00  5.7364E-01  1.0115E+00
             9.5944E-01
 PARAMETER:  1.0375E-01 -4.8112E+00  7.2718E-01  5.7536E-01  2.8025E-02  7.8435E-02  2.3582E+00  6.1253E-01 -4.5575E-01  1.1141E-01
             5.8595E-02
 GRADIENT:   5.0173E+02  0.0000E+00  1.0349E+01  1.2585E+03  6.7806E+00  4.7373E+01  7.1112E-01  2.9712E+00  3.0510E+01  1.2110E+00
             8.2052E-01

0ITERATION NO.:  155    OBJECTIVE VALUE:  -1720.23002874341        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     5612
 NPARAMETR:  1.0038E+00  1.0000E-02  1.8706E+00  1.6085E+00  9.3041E-01  9.7867E-01  9.5637E+00  1.6686E+00  5.7365E-01  1.0109E+00
             9.5943E-01
 PARAMETER:  1.0375E-01 -4.8112E+00  7.2625E-01  5.7533E-01  2.7874E-02  7.8434E-02  2.3580E+00  6.1201E-01 -4.5574E-01  1.1084E-01
             5.8583E-02
 GRADIENT:   1.9632E+00  0.0000E+00  1.8361E-01 -2.3745E+01 -1.8103E-01  1.2647E-01  1.6190E-02  6.4968E-02  2.1548E-01  7.8581E-02
             1.0479E-02

0ITERATION NO.:  156    OBJECTIVE VALUE:  -1720.23002874341        NO. OF FUNC. EVALS.:  28
 CUMULATIVE NO. OF FUNC. EVALS.:     5640
 NPARAMETR:  1.0038E+00  1.0000E-02  1.8679E+00  1.6085E+00  9.3067E-01  9.7867E-01  9.5571E+00  1.6675E+00  5.7365E-01  1.0098E+00
             9.5938E-01
 PARAMETER:  1.0375E-01 -4.8112E+00  7.2625E-01  5.7533E-01  2.7874E-02  7.8434E-02  2.3580E+00  6.1201E-01 -4.5574E-01  1.1084E-01
             5.8583E-02
 GRADIENT:   3.0328E-03  0.0000E+00  2.3372E-01  1.5209E-01 -1.7368E-01 -4.4491E-04  1.0984E-04  4.0487E-02  7.0466E-04  3.9632E-02
             1.1067E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     5640
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.3013E-04  1.8980E-05 -3.6008E-02 -9.3084E-03 -4.3880E-02
 SE:             2.9869E-02  1.9782E-03  1.8906E-02  2.8840E-02  2.0304E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9652E-01  9.9234E-01  5.6834E-02  7.4688E-01  3.0686E-02

 ETASHRINKSD(%)  1.0000E-10  9.3373E+01  3.6663E+01  3.3817E+00  3.1979E+01
 ETASHRINKVR(%)  1.0000E-10  9.9561E+01  5.9884E+01  6.6490E+00  5.3731E+01
 EBVSHRINKSD(%)  3.8622E-01  9.3555E+01  4.0126E+01  3.8186E+00  2.8103E+01
 EBVSHRINKVR(%)  7.7094E-01  9.9585E+01  6.4151E+01  7.4914E+00  4.8309E+01
 RELATIVEINF(%)  9.5570E+01  9.0958E-03  8.8390E+00  2.2472E+00  8.6247E+00
 EPSSHRINKSD(%)  4.5173E+01
 EPSSHRINKVR(%)  6.9940E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1720.2300287434141     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -985.07920217967592     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    81.22
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     6.25
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1720.230       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.00E-02  1.87E+00  1.61E+00  9.30E-01  9.79E-01  9.56E+00  1.67E+00  5.74E-01  1.01E+00  9.59E-01
 


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
+        1.19E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -5.99E+00  0.00E+00  2.75E+01
 
 TH 4
+       -1.07E+02  0.00E+00 -1.16E+01  1.27E+03
 
 TH 5
+        4.43E+01  0.00E+00 -1.33E+02 -1.72E+02  7.02E+02
 
 TH 6
+       -1.88E+00  0.00E+00 -8.77E-01 -5.49E+01  7.80E+00  2.17E+02
 
 TH 7
+        3.80E-02  0.00E+00  1.66E-03  3.52E-02 -4.20E-03 -1.80E-02  3.64E-05
 
 TH 8
+        1.05E+00  0.00E+00  8.73E-01 -6.39E+00 -5.19E+00  7.12E-01 -4.45E-04  4.07E-01
 
 TH 9
+        1.74E+02  0.00E+00  4.64E+00  1.50E+02 -1.15E+00 -8.26E+01  1.65E-01 -2.40E+00  7.50E+02
 
 TH10
+       -3.30E-01  0.00E+00  1.49E+01 -1.74E+00 -7.84E+01  1.02E+00  2.61E-04  1.41E+00 -1.33E+00  1.06E+01
 
 TH11
+        3.55E-02  0.00E+00 -2.38E+01 -8.62E+01  7.73E+01  1.31E+01 -7.00E-03  9.60E+00 -3.79E+01  1.38E+01  3.10E+02
 
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
+        1.14E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.19E+00  0.00E+00  4.40E+01
 
 TH 4
+       -1.68E+01  0.00E+00 -2.75E+01  1.21E+03
 
 TH 5
+       -1.71E+00  0.00E+00 -1.27E+02 -1.02E+02  6.79E+02
 
 TH 6
+       -1.63E-01  0.00E+00  1.14E-01 -3.20E+00  2.57E-01  2.05E+02
 
 TH 7
+        4.24E-04  0.00E+00  2.37E-03  1.77E-03 -3.51E-03 -3.53E-04  1.69E-03
 
 TH 8
+        1.33E-01  0.00E+00 -1.55E+01 -4.65E+00 -1.02E+01  5.01E-02  1.78E-04  2.06E+01
 
 TH 9
+        2.85E+00  0.00E+00  5.27E+00 -1.62E+00  1.00E+00 -2.56E+00  1.15E-01  7.36E-01  5.23E+02
 
 TH10
+        6.64E-01  0.00E+00  1.23E-01 -3.41E+00 -7.98E+01  3.71E-01  5.20E-03  9.65E+00  3.26E+00  6.56E+01
 
 TH11
+       -8.65E+00  0.00E+00 -3.36E+00 -1.64E+01 -5.13E+00  3.20E+00  5.47E-03  4.92E+00  2.01E+01  1.03E+01  2.25E+02
 
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
+        1.14E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -2.02E+01  0.00E+00  4.27E+01
 
 TH 4
+        1.00E+02  0.00E+00 -3.80E+01  1.21E+03
 
 TH 5
+       -4.04E+01  0.00E+00 -1.29E+02 -5.59E+01  6.83E+02
 
 TH 6
+       -1.32E+01  0.00E+00  1.02E-01  2.90E+01 -3.75E+00  2.02E+02
 
 TH 7
+       -3.45E-02  0.00E+00  2.96E-03 -2.19E-02  3.93E-03  1.29E-02  2.70E-05
 
 TH 8
+        1.79E+01  0.00E+00 -1.21E+01 -2.35E+00 -1.28E+01 -1.21E+00 -6.91E-04  1.53E+01
 
 TH 9
+       -1.31E+02  0.00E+00  1.45E+01 -1.04E+02 -8.84E+00  5.06E+01  9.56E-02 -4.64E+00  4.08E+02
 
 TH10
+        2.76E+01  0.00E+00  1.97E+00 -2.92E+01 -7.53E+01  9.70E-01  2.08E-03  1.07E+01  1.09E+01  5.73E+01
 
 TH11
+       -2.48E+01  0.00E+00  6.98E+00  2.11E+01 -6.46E+01  2.31E+00  8.65E-03  3.76E+00  4.89E+01  1.79E+01  1.72E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.02
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       86.035
Stop Time:
Wed Sep 29 18:23:51 CDT 2021
